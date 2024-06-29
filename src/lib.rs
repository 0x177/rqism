use ndarray::{array,Array, Array1,Array2,linalg::kron};
use num::{complex::{Complex}};
use rand::prelude::*;

mod instruction;
use crate::instruction::Instruction;

const ONE_SQR_TWO: f32 = 0.7071067811865475;

#[derive(Clone,Debug)]
pub struct QuantumState {
    pub n: usize,
    pub state: Array1<Complex<f32>>,
    pub reg: usize,
    pub rng: ThreadRng
}

fn kron_expt(u: &Array2::<Complex<f32>>, n: usize) -> Array2::<Complex<f32>> {
    match n {
        0 => array![[Complex::<f32>::from(1.0)]],
        1 => u.clone(),
        _ => kron(&kron_expt(u, n - 1), u),
    }
}

fn lift(t: &Array2::<Complex<f32>>, i: usize, n: usize) -> Array2::<Complex<f32>> {
    let left = kron_expt(&Array2::<Complex<f32>>::eye(2), n - i - t.nrows().ilog2() as usize);
    let right = kron_expt(&Array2::<Complex<f32>>::eye(2), i);
    kron(&left, &kron(t, &right))
}

fn ptt(p: &Vec<usize>) -> Vec<(usize, usize)> {
    p.iter()
        .enumerate()
        .filter_map(|(mut src, &dest)| {
            while src < dest {
                src = p[src];
            }
            (src > dest).then_some((dest, src))
        })
        .collect()
}

fn tta(t: &Vec<(usize, usize)>) -> Vec<usize> {
    t.iter()
        .flat_map(|&(a, b)| {
            if b.saturating_sub(a) <= 1 {
                vec![a]
            } else {
                (a..b).chain((a..b - 1).rev()).collect()
            }
        })
        .collect()
}

impl QuantumState {
    pub fn new(n: usize) -> Self {
	let rng = rand::thread_rng();
        let mut state = Array::zeros(2usize.pow(n as u32));
        state[[0]] = Complex::new(1.0, 0.0);

        Self { n, state, reg: 0,rng}
    }
    
    pub fn gate_apply_sq(&self,t: &Array2<Complex<f32>>,i: usize) -> Self {
	Self {
	    n: self.n,
	    state: lift(t,i,self.n).dot(&self.state),
	    reg: 0,
	    rng: self.rng.clone()
	}
    }

    pub fn gate_apply_nq(&self, t: &Array2::<Complex<f32>>, i: &[usize]) -> Self {

	let swap: Array2::<Complex<f32>> = array![
            [Complex::<f32>::from(1.0), Complex::<f32>::from(0.0), Complex::<f32>::from(0.0), Complex::<f32>::from(0.0)],
            [Complex::<f32>::from(0.0), Complex::<f32>::from(0.0), Complex::<f32>::from(1.0), Complex::<f32>::from(0.0)],
            [Complex::<f32>::from(0.0), Complex::<f32>::from(1.0), Complex::<f32>::from(0.0), Complex::<f32>::from(0.0)],
            [Complex::<f32>::from(0.0), Complex::<f32>::from(0.0), Complex::<f32>::from(0.0), Complex::<f32>::from(1.0)]
	];
	
	let i_to_op = |i: &[usize]| {
	    i.iter()
		.fold(kron_expt(&Array2::<Complex<f32>>::eye(2),self.n as usize), |a,i| {
		    a.dot(&lift(&swap,*i,self.n))
		})
	};

	let idk = i
	    .iter()
	    .copied()
	    .rev()
	    .chain(0..self.n).filter(|ind| !i.contains(ind))
	    .collect::<Vec<_>>();

        let tr = tta(&ptt(&idk));
        let tf = i_to_op(&tr);
        let ft = i_to_op(&tr.into_iter().rev().collect::<Vec<_>>());
        let all = tf.dot(&lift(t,0,self.n).dot(&ft));
	
	Self {
            n: self.n,
            state: all.dot(&self.state),
	    reg: 0,
	    rng: self.rng.clone()
	}
    }

    pub fn execute_instruction(&self,ins: &Instruction) -> Self {
	match ins {
	    Instruction::Gate {matrix,indices} => {
		if indices.len() == 1 {
		    self.gate_apply_sq(&matrix,indices[0])
		} else {
		    self.gate_apply_nq(&matrix,indices)
		}
	    },
	    Instruction::Measure {indices} => {
		self.measure_qubits(indices)
	    }
	}
    }

    pub fn execute_circuit(&self,circuit: Vec<Instruction>,index: usize) -> Self {
	if index == circuit.len() {return self.clone()}

	(self.execute_instruction(&circuit[index])).execute_circuit(circuit,index+1)
    }

    pub fn measure(&self) -> Self {
	let mut s = self.clone();
	let mut rand = s.rng.gen_range(0.0..1.0);
	let mut k: i32 = -1;
	
	for (i,theta) in self.state.iter().enumerate() {
	    rand -= theta.norm_sqr();

	    if rand < 0.0 {
		k = i as i32;
	    }
	} 

	if k == -1 {
	    k = self.n as i32 - 1;
	}

	s.state.fill(0.0.into());
	s.state[0] = 1.0.into();
	s.reg = k as usize;

	return s;
    }

    pub fn measure_qubits(&self,indices: &[usize]) -> Self {
	let mut s = self.clone();

	for i in indices {
	    let rand = s.rng.gen_range(0.0..1.0);
	    s.state[*i] = if rand >= s.state[*i].re {Complex::new(1.0,0.0)} else {Complex::new(0.0,0.0)}
	}
	
	return s;
    }

    pub fn print(&self) -> Self {
	self.state.iter().enumerate().for_each(|(x,y)| {
	    println!("quibit {x}: {y}");
	}); 
	self.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn bell() {
	// bell
	let mut counts = vec![0; 4];
	let machine = QuantumState::new(2);

	let circuit = vec![
	    Instruction::hadamard(0),
	    Instruction::cnot(vec![0,1]),
	    Instruction::Measure {indices: vec![0,1]}
	];
	
	for _ in 0..1000 {
	    counts[machine
		   .execute_circuit(circuit.clone(),0)
		   .measure()
		   .reg
	    ] += 1;
	}

	println!("{:?}",counts);
    }
}

