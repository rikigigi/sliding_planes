use superpoly::poly::Poly;
use rand::Rng;
use rand::distributions::{Distribution, Uniform};
use std::io;

fn read_number<T : std::str::FromStr >() -> T {
    let mut input = String::new();
    loop {
      io::stdin()
         .read_line(&mut input)
         .expect("Reading of the input failed");
      let n : T = match input.trim().parse() {
             Ok(num) => num,
             Err(_) => { println!("You must enter a number with correct type. You entered'{input}'");  continue },
          };
          return n;
    }
}



enum McMove {
    SingleShift{dx: f64, idx : usize},
    PairShift{dx: f64, idx : usize},
    TripletShift{dx : f64, idx : usize},
}

struct PlaneSystem {
    x : [f64;5] ,
    t : f64,
    p : Poly,
    beta : f64,
    p2 : f64
}


fn vp(x : f64 ) -> f64 {
   x - ((x+3.0)/6.0).floor()*6.0
}

impl PlaneSystem {

   fn new(p: Poly, t: f64, beta:f64, p2 : f64) -> PlaneSystem {
       PlaneSystem {
          x: [0.0,0.0,0.0,0.0,0.0],
          t,
          p,
          beta,
          p2
       }
   }

   fn abc(&self) -> [f64;5] {
      self.x.map(|a| a-((a+0.5)/3.0).floor()*3.0)
   }

   fn v(&self) -> (f64,f64,f64) {
       let (mut v0,mut v1,mut v2) = (0.0,0.0,0.0);
       for i in 1..5{
           let vs = self.p.calc(vp(self.x[i]-self.x[i-1]));
           v0+=vs.0;
           v1+=vs.1;
           v2+=vs.2;
       }
       for x in [self.x[0], self.t - self.x[4]] {
          let vs = self.p.calc(vp(x));
          v0+=vs.0;
          v1+=vs.1;
          v2+=vs.2;
       }
       // second neighbor interaction
       for i in 2..5{
           let vs = self.p.calc(vp(self.x[i]-self.x[i-2]));
           v0+=vs.0*self.p2;
           v1+=vs.1*self.p2;
           v2+=vs.2*self.p2;
       }
       for x in [self.x[1], self.x[0]+self.t - self.x[4], self.t - self.x[3]] {
          let vs = self.p.calc(vp(x));
          v0+=vs.0*self.p2;
          v1+=vs.1*self.p2;
          v2+=vs.2*self.p2;
       }
       (v0,v1,v2)
   }
   
   fn dv_dt(&self) -> f64 {
      self.p.calc(vp(self.t-self.x[4])).1+
      self.p2*self.p.calc(vp(self.t+self.x[0]-self.x[4])).1
   }

   fn dp(&self, idx1 : usize, idx2 : usize, dx : f64) -> f64 {
       //idx1 is the index of the plane in the real system, starting from 0 (that is the fixed plane)
       //idx2 is the index of the plane that is found going up. If it is less than idx1, the code
       //will go back throw the pbc, thus adding each time the tilt
       //it does something similar to p(x2-x1+dx)-p(x2-x1) , with pbc
       //put dx with the opposite sign if you want to translate plane idx1!
       let mut x2 = if idx2%6 == 0 {0.0} else {self.x[idx2%6-1]} + self.t*(idx2/6) as f64;
       let x1 = if idx1%6 == 0 {0.0} else {self.x[idx1%6-1]} + self.t*(idx1/6) as f64;

       //move idx2_ after idx1
       let mut idx2_ = idx2;
       while idx2_ < idx1 {
           idx2_ += 6;
           x2 += self.t;
       }
       self.p.calc(vp(x2-x1+dx)).0-self.p.calc(vp(x2-x1)).0
   }

   fn move_dv_new(&self, m : &Option<McMove>) -> Option<f64> {
       match m {
           Some(McMove::SingleShift{dx,idx}) => {
              let idx = *idx;
              let dx = *dx;
              Some(
                    self.dp(idx,idx+1,dx)+self.dp(idx+1,idx+2,-dx)+ //first neigh
                    self.p2*(self.dp(idx+6-1,idx+1+6,dx)+self.dp(idx+1,idx+3,-dx)) //second neigh
                 )
           },
           Some(McMove::PairShift{dx,idx}) => {
              let idx = *idx;
              let dx = *dx;
              Some(
                    self.dp(idx,idx+1,dx)+self.dp(idx+2,idx+3,-dx)+ //first neigh
                    self.p2*( // second neigh
                       self.dp(idx+6-1,idx+1+6,dx)+self.dp(idx+2,idx+4,-dx)+
                       self.dp(idx,idx+2,dx)+self.dp(idx+1,idx+3,-dx)
                    )
                  )
           },
           Some(McMove::TripletShift{dx,idx}) => {
              let idx = *idx;
              let dx = *dx;
              Some(
                    self.dp(idx,idx+1,dx)+self.dp(idx+3,idx+4,-dx)+ //first neigh
                    self.p2*(
                       self.dp(idx+2,idx+4,-dx)+self.dp(idx,idx+2,dx)+ //middle plane
                       self.dp(idx+6-1,idx+6+1,dx)+self.dp(idx+3,idx+5,-dx) //1st and 3rd
                    )
                  )
           },
           None => None
       }
   }

   fn move_dv(&self, m : &Option<McMove>) -> Option<f64> {
       match m {
          Some(McMove::SingleShift{dx,idx}) => {
              let idx = *idx;
              match idx {
                  0 => Some(  //moving the plane x_1, whose position is x[0]. The difference is due to p(x[0]) and p(x[1]-x[0])
                      self.p.calc(vp(self.x[idx]+dx)).0-self.p.calc(vp(self.x[idx])).0 + self.p.calc(vp(self.x[idx+1]-(self.x[idx]+dx) )).0-self.p.calc(vp(self.x[idx+1]-self.x[idx])).0
                      ),
                  1..=3 => Some( // difference due to p(x[i]-x[i-1]) and p(x[i+1]-x[i])
                      self.p.calc(vp(self.x[idx]+dx-(self.x[idx-1]) )).0 - self.p.calc(vp(self.x[idx]-(self.x[idx-1]) )).0 + self.p.calc(vp(self.x[idx+1]-(self.x[idx]+dx) )).0 - self.p.calc(vp(self.x[idx+1]-(self.x[idx]) )).0
                      ),
                  4 => Some (// moving last plane, with pbc, interacting with the first fixed plane
                      self.p.calc(vp(self.x[idx]+dx-self.x[idx-1])).0- self.p.calc(vp(self.x[idx]-self.x[idx-1])).0 + self.p.calc(vp(self.t-(self.x[idx]+dx))).0 - self.p.calc(vp(self.t-self.x[idx])).0 
                      ),
                  _ => None
              }
          },
          Some(McMove::PairShift{dx,idx}) => {
              let idx = *idx;
              match idx {
                  0 => Some ( //moving planes x_1 and x_2.
                        self.p.calc(vp(self.x[idx]+dx)).0-self.p.calc(vp(self.x[idx])).0 + self.p.calc(vp(self.x[idx+2]-(self.x[idx+1]+dx) )).0 - self.p.calc(vp(self.x[idx+2]-(self.x[idx+1]) )).0
                       ),
                 1..=2 => Some (
                       self.p.calc(vp(self.x[idx]+dx-(self.x[idx-1]) )).0 - self.p.calc(vp(self.x[idx]-(self.x[idx-1]) )).0 + self.p.calc(vp(self.x[idx+2]-(self.x[idx+1]+dx) )).0 - self.p.calc(vp(self.x[idx+2]-(self.x[idx+1]) )).0 
                       ),
                 3 => Some ( //the last one interact with the first fixed one due to pbc
                       self.p.calc(vp(self.x[idx]+dx-(self.x[idx-1]) )).0 - self.p.calc(vp(self.x[idx]-(self.x[idx-1]) )).0 +  self.p.calc(vp(self.t-(self.x[idx+1]+dx))).0 - self.p.calc(vp(self.t-self.x[idx+1])).0
                     ),
                 _ => None
              }
          },
          Some(McMove::TripletShift{dx,idx}) => {
              let idx = *idx;
              match idx {
                  0 => Some ( //moving planes x_1 and x_2 and x_3.
                        self.p.calc(vp(self.x[idx]+dx)).0-self.p.calc(vp(self.x[idx])).0 + self.p.calc(vp(self.x[idx+3]-(self.x[idx+2]+dx) )).0 - self.p.calc(vp(self.x[idx+3]-(self.x[idx+2]) )).0
                       ),
                 1 => Some ( // x_2, x_3 and x_4
                       self.p.calc(vp(self.x[idx]+dx-(self.x[idx-1]) )).0 - self.p.calc(vp(self.x[idx]-(self.x[idx-1]) )).0 + self.p.calc(vp(self.x[idx+2]-(self.x[idx+1]+dx) )).0 - self.p.calc(vp(self.x[idx+2]-(self.x[idx+1]) )).0 
                       ),
                 2 => Some ( // x_3, x_4, x_5, but the last one interact with the first fixed one due to pbc
                       self.p.calc(vp(self.x[idx]+dx-(self.x[idx-1]) )).0 - self.p.calc(vp(self.x[idx]-(self.x[idx-1]) )).0 +  self.p.calc(vp(self.t-(self.x[idx+2]+dx))).0 - self.p.calc(vp(self.t-self.x[idx+2])).0
                     ),
                 _ => None
              }
          }
          None => None
       }
   }

   fn apply_move(&mut self, m : Option<McMove> ) -> Option<McMove> {
       match m {
           Some(McMove::SingleShift{dx,idx}) => 
                match idx {
                   0..=4 => { 
                           self.x[idx]+=dx;
                           m
                        },
                   _ => None
                  },
           Some(McMove::PairShift{dx,idx}) =>
                match idx {
                    0..=3 =>  { 
                           self.x[idx]+=dx; self.x[idx+1]+=dx;
                           m
                        },
                    _ => None
               },
           Some(McMove::TripletShift{dx,idx}) =>
                match idx {
                    0..=2 =>  { 
                           self.x[idx]+=dx; self.x[idx+1]+=dx; self.x[idx+2]+=dx;
                           m
                        },
                    _ => None
               },
           None => None
      }
   }
}

fn generate_random_move(rng : &mut rand::rngs::ThreadRng, dxmax : f64) -> Option<McMove> {
    let moveChoice = Uniform::from(0..3);
    let planeChoice = Uniform::from(0..5);
    let dxChoice = Uniform::from(-dxmax..dxmax);
    let idx = moveChoice.sample(rng);
    let add = (rng.gen_range(0..=2) as f64)*2.0;
    match idx {
        0 => Some(McMove::SingleShift{
                dx : dxChoice.sample(rng)+add,
                idx : planeChoice.sample(rng)
             }),
        1 => Some(McMove::PairShift{
                dx : dxChoice.sample(rng)+add,
                idx : planeChoice.sample(rng)%4
             }),
        2 => Some(McMove::PairShift{
                dx : dxChoice.sample(rng)+add,
                idx : planeChoice.sample(rng)%3
             }),
        _ => None
    }
}

struct OnlineAverage {
    mean : f64,
    M2 : f64,
    count : u64
}

impl OnlineAverage {
   fn new() -> OnlineAverage {
     OnlineAverage{
       mean: 0.0,
       M2: 0.0,
       count: 0
     }
   }
   fn update(&mut self,newValue : f64) {
      self.count+=1;
      let delta = newValue - self.mean;
      self.mean += delta/self.count as f64;
      let delta2 = newValue - self.mean;
      self.M2 += delta*delta2
   }
   fn get_mean_variance_count(&self) -> (f64,f64,u64) {
      if self.count >= 2 {
         (self.mean,self.M2/self.count as f64, self.count)
      } else {
         (self.mean,self.M2,self.count)
      }
   }

}

fn main() {

    println!("#Please input h0, h1, second_neigh factor, tilt, beta and number of trials of the inner loop and the outer loop");
    let h0 = read_number::<f64>();
    let h1 = read_number::<f64>();
    let p2 = read_number::<f64>();
    let t = read_number::<f64>(); 
    let beta = read_number::<f64>(); 
    let nsteps = read_number::<u64>();
    let nsteps_out = read_number::<u64>();

    let p = Poly::new(vec![3.0,4.55,4.3,5.02,7.1],vec![0.0,0.0,0.0,9.0,90.0], h0, h1);

    let mut system = PlaneSystem::new(p,t, beta,p2);
    let mut rng = rand::thread_rng();
    let rnd_01 = Uniform::from(0.0..1.0);

    let mut v_mean = OnlineAverage::new();
    let mut de_mean = OnlineAverage::new();
    let mut acc_mean = OnlineAverage::new();
    let mut te_mean = OnlineAverage::new();

    for jstep in 0..nsteps_out { 
       let mut accepted = 0;
       let mut rejected = 0;
       for istep in 0..nsteps {
          let m = generate_random_move(&mut rng, 0.1);
          let de_ = system.move_dv(&m);
          let de_new = system.move_dv_new(&m);
          //assert!(
          //   match de_new {
          //       Some(xnew) => match de_ {
          //                Some(xold) => (xnew-xold).abs() < 1e-12,
          //                None => false
          //            },
          //       None => match de_ {
          //                None => true,
          //                Some(_) => false
          //           }
          //   }
          //);
          match de_ {
             Some(de) => {
                       te_mean.update(de);
                       let p = (-system.beta*de).exp();
                       if p >= 1.0 {
                          let m = system.apply_move(m);
                          match m { None => {rejected +=1;}, McMove => {accepted+=1;}}
                       } else {
                          let r = rnd_01.sample(&mut rng);
                          if r < p {
                              let m = system.apply_move(m);
                              match m { None => {rejected +=1;}, McMove => {accepted+=1;}}
                          } else {rejected += 1;}
                       }
                  },
             None => ()
          }
       }
       let v = system.v();
       let fupper = system.dv_dt();
       let abc = system.abc();
       let accr=accepted as f64 / (accepted + rejected) as f64;
       v_mean.update(v.0);
       de_mean.update(fupper);
       acc_mean.update(accr);
       println!("{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
            accr,
            fupper, v.0, v.1, v.2,
            system.x[0],system.x[1],system.x[2],system.x[3], system.x[4],
            abc[0],abc[1],abc[2],abc[3],abc[4] );
    }
    println!("#!delta_v {} {} {}",te_mean.get_mean_variance_count().0, te_mean.get_mean_variance_count().1, te_mean.get_mean_variance_count().2);
    println!("#!v {} {} {}",v_mean.get_mean_variance_count().0, v_mean.get_mean_variance_count().1, v_mean.get_mean_variance_count().2);
    println!("#!dv_dt {} {} {}",de_mean.get_mean_variance_count().0, de_mean.get_mean_variance_count().1, de_mean.get_mean_variance_count().2);
    println!("#!accr {} {} {}",acc_mean.get_mean_variance_count().0, acc_mean.get_mean_variance_count().1, acc_mean.get_mean_variance_count().2);

}
