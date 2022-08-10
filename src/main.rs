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
    PairShift{dx: f64, idx : usize}
}

struct PlaneSystem {
    x : [f64;5] ,
    t : f64,
    p : Poly,
    beta : f64
}


fn vp(x : f64 ) -> f64 {
   x - ((x+3.0)/6.0).floor()*6.0
}

impl PlaneSystem {

   fn new(p: Poly, t: f64, beta:f64) -> PlaneSystem {
       PlaneSystem {
          x: [0.0,0.0,0.0,0.0,0.0],
          t,
          p,
          beta
       }
   }


   fn v(&self) -> (f64,f64,f64) {
       let (mut v0,mut v1,mut v2) = (0.0,0.0,0.0);
       for i in 1..5{
           let vs = self.p.calc(vp(self.x[i]-self.x[i-1]));
           v0+=vs.0;
           v1+=vs.1;
           v2+=vs.2;
       }
       let vs = self.p.calc(vp(self.x[0]));
       v0+=vs.0;
       v1+=vs.1;
       v2+=vs.2;  
       let vs = self.p.calc(vp(self.t - self.x[4]));
       v0+=vs.0;
       v1+=vs.1;
       v2+=vs.2;  
       (v0,v1,v2)
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
           None => None
      }
   }
}

fn generate_random_move(rng : &mut rand::rngs::ThreadRng, dxmax : f64) -> Option<McMove> {
    let moveChoice = Uniform::from(0..2);
    let planeChoice = Uniform::from(0..5);
    let dxChoice = Uniform::from(-dxmax..dxmax);
    let idx = moveChoice.sample(rng);
    match idx {
        0 => Some(McMove::SingleShift{
                dx : dxChoice.sample(rng),
                idx : planeChoice.sample(rng)
             }),
        1 => Some(McMove::PairShift{
                dx : dxChoice.sample(rng),
                idx : planeChoice.sample(rng)
             }),
        _ => None
    }
}



fn main() {

    println!("#Please input h0, h1, tilt, beta and number of trials of the inner loop and the outer loop");
    let h0 = read_number::<f64>();
    let h1 = read_number::<f64>();
    let t = read_number::<f64>(); 
    let beta = read_number::<f64>(); 
    let nsteps = read_number::<u64>();
    let nsteps_out = read_number::<u64>();

    let p = Poly::new(vec![3.0,4.55,4.3,5.02,7.1],vec![0.0,0.0,0.0,9.0,90.0], h0, h1);

    let mut system = PlaneSystem::new(p,t, beta);
    let mut rng = rand::thread_rng();
    let rnd_01 = Uniform::from(0.0..1.0);

    for jstep in 0..nsteps_out { 
       let mut accepted = 0;
       let mut rejected = 0;
       for istep in 0..nsteps {
          let m = generate_random_move(&mut rng, 0.1);
          let de_ = system.move_dv(&m);
          match de_ {
             Some(de) => {
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
       println!("{} {} {} {} {} {} {} {} {}",
            accepted as f64 / (accepted + rejected) as f64,
            v.0, v.1, v.2,
            system.x[0],system.x[1],system.x[2],system.x[3], system.x[4] )
    }

}
