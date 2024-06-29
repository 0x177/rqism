use ndarray::{Array,Array2};
use num::{complex::{Complex}};

const ONE_SQR_TWO: f32 = 0.7071067811865475;

#[derive(Clone)]
pub enum Instruction {
    Gate {
	matrix: Array2::<Complex<f32>>,
	indices: Vec<usize>
    },
    Measure {
	indices: Vec<usize>
    }
}

impl Instruction {
    pub fn hadamard(i: usize) -> Self {
	let gate = Complex::new(ONE_SQR_TWO,0.0) * Array::from_shape_vec((2, 2), vec![
	    Complex::new(1.0,0.0),Complex::new(1.0,0.0),Complex::new(1.0,0.0),Complex::new(-1.0,0.0)]
	).unwrap();

	Self::Gate {matrix: gate,indices: vec![i]}
    }

    pub fn not(i: usize) -> Self {
	let gate = Array::from_shape_vec((2,2), vec![
	    Complex::new(0.0,0.0),Complex::new(1.0,0.0),Complex::new(1.0,0.0),Complex::new(0.0,0.0)
	]).unwrap();

	Self::Gate {matrix: gate,indices: vec![i]}
    }

    pub fn cnot(i: Vec<usize>) -> Self {
	let gate = Array::from_shape_vec((4, 4), vec![
	    Complex::new(1.0,0.0),Complex::new(0.0,0.0),Complex::new(0.0,0.0),Complex::new(0.0,0.0),
	    Complex::new(0.0,0.0),Complex::new(1.0,0.0),Complex::new(0.0,0.0),Complex::new(0.0,0.0),
	    Complex::new(0.0,0.0),Complex::new(0.0,0.0),Complex::new(0.0,0.0),Complex::new(1.0,0.0),
	    Complex::new(0.0,0.0),Complex::new(0.0,0.0),Complex::new(1.0,0.0),Complex::new(0.0,0.0),
	]
	).unwrap();

	Self::Gate {matrix: gate,indices: i}
    }

    pub fn t_gate(i: usize) -> Self {
	let gate = Array::from_shape_vec((2, 2), vec![
	    Complex::new(1.0,0.0),Complex::new(0.0,0.0),Complex::new(0.0,0.0),Complex::new(ONE_SQR_TWO + ONE_SQR_TWO,1.0)]
	).unwrap();

	Self::Gate {matrix: gate,indices: vec![i]}
    }
}
