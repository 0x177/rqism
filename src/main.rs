use ndarray::{Array, Array1,Array2};
use num::{integer::Roots,complex::Complex};

const one_sqr_two: f32 = 0.7071067811865475;

#[derive(Clone,Debug)]
struct QuantumState {
    n: u32,
    state: Array1<Complex<f32>>,
}

impl QuantumState {
    fn new(n: u32) -> Self {
        let mut state = Array::zeros(2usize.pow(n as u32));
        state[[0]] = Complex::new(1.0, 0.0);

        Self { n, state }
    }

    fn gate_apply(&self, t: Array2<Complex<f32>>, i: u32) -> Self {
        let left = Array::eye(2usize.pow(i));
        let right = Array::eye(2usize.pow(self.n - i - (t.shape()[0] as usize).sqrt() as u32));

        let all = ndarray::linalg::kron(&ndarray::linalg::kron(&left, &t.view()), &right);

	Self {
	    n: self.n,
	    state: all.dot(&self.state)
	}
    }

    fn hadamard(&self,i: u32) -> Self {
	let gate = Complex::new(one_sqr_two,0.0) * Array::from_shape_vec((2, 2), vec![
	    Complex::new(1.0,0.0),Complex::new(1.0,0.0),Complex::new(1.0,0.0),Complex::new(-1.0,0.0)]
	).unwrap();

	self.gate_apply(gate, i)
    }

    fn cnot(&self,i: u32) -> Self {
	let gate = Array::from_shape_vec((4, 4), vec![
	    Complex::new(1.0,0.0),Complex::new(0.0,0.0),Complex::new(0.0,0.0),Complex::new(0.0,0.0),
	    Complex::new(0.0,0.0),Complex::new(1.0,0.0),Complex::new(0.0,0.0),Complex::new(0.0,0.0),
	    Complex::new(0.0,0.0),Complex::new(0.0,0.0),Complex::new(0.0,0.0),Complex::new(1.0,0.0),
	    Complex::new(0.0,0.0),Complex::new(0.0,0.0),Complex::new(1.0,0.0),Complex::new(0.0,0.0),
	]
	).unwrap();

	self.gate_apply(gate,i)
    }

    fn t_gate(&self,i: u32) -> Self {
	let gate = Array::from_shape_vec((2, 2), vec![
	    Complex::new(1.0,0.0),Complex::new(0.0,0.0),Complex::new(0.0,0.0),Complex::new(one_sqr_two + one_sqr_two,1.0)]
	).unwrap();

	self.gate_apply(gate,i)
    }

    fn print(&self) -> Self {
	println!("{:?}",self);

	self.clone()
    }
}

fn main() {
    //quantum hello world
    QuantumState::new(2)
	.hadamard(0)
	.cnot(0)
	.print();
}
