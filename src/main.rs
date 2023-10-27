use num::complex::ComplexFloat;
use rand::Rng;

const IDENT: [[f32; 2]; 2] = [
    [1.0,0.0],
    [0.0,1.0],
];

fn trv<const N:usize>(s:[[f32; N]; N]) -> Vec<Vec<f32>> {
    s.iter().map(|&e| e.to_vec()).collect()
}

fn operator_apply(matrix: &Vec<Vec<f32>>, vector: &Vec<num::complex::Complex<f32>>) -> Vec<num::complex::Complex<f32>> {
    let mut result = vec![num::complex::Complex::new(0.0,0.0); matrix.len()];
    
    for i in 0..matrix.len() {
        for j in 0..matrix[0].len() {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    
    result
}

fn operator_compose(mat1: &Vec<Vec<f32>>, mat2: &Vec<Vec<f32>>) -> Vec<Vec<f32>> {
    let rows1 = mat1.len();
    let cols1 = mat1[0].len();
    let cols2 = mat2[0].len();

    let mut result = vec![vec![0.0; cols2]; rows1];

    for i in 0..rows1 {
        for j in 0..cols2 {
            for k in 0..cols1 {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    result
}

fn kron_mult(mat1: &Vec<Vec<f32>>, mat2: &Vec<Vec<f32>>) -> Vec<Vec<f32>> {
    let m = mat1.len();
    let n = mat1[0].len();

    let p = mat2.len();
    let q = mat2[0].len();

    let mut temp = vec![vec![0.0; n*q]; m*p];

    for i in 0..m {
        for j in 0..n {
            let mat_one_i_j = mat1[i][j];
            let y = i * p;
            let x = j * q;

            for u in 0..p {
                for v in 0..q {
                    temp[y+u][x+v] = mat_one_i_j * mat2[u][v];
                }
            }
        }
    }
    
    temp
}

fn kron_exp(m:&Vec<Vec<f32>>,n:i32) -> Vec<Vec<f32>> {
    if n < 1 {return vec![vec![1.0,0.0],vec![0.0,0.0]]}
    if n == 1 {return m.to_vec()}
    return kron_mult(&kron_exp(&m,n-1),&m);
}

fn make_q_machine(n:u32) -> QuantumMachine {
    let mut qb = vec![num::complex::Complex::new(0.0,0.0); 2_u32.pow(n) as usize];
    qb[0] = num::complex::Complex::new(1.0,0.0);
    QuantumMachine {
	state: qb,
	measurement: 0.0,
    }
}

#[derive(Clone)]
struct QuantumMachine {
    state: Vec<num::complex::Complex<f32>>,
    measurement: f32, 
}

fn qubit_cnt(d: u32) -> u32 {
    return d.count_ones()+d.count_zeros();
}

fn lift(U: &Vec<Vec<f32>>, i: usize, n: usize) -> Vec<Vec<f32>> {
    kron_mult(
	&kron_exp(&trv(IDENT),(n-i-qubit_cnt(U[0].len() as u32) as usize) as i32),
	&kron_mult(&U,
			   &kron_exp(&trv(IDENT), i as i32)
	)
    )
}

fn apply_sq_gate(state: &Vec<num::complex::Complex<f32>>, U: &Vec<Vec<f32>>,q: usize) -> Vec<num::complex::Complex<f32>> {
    operator_apply(&lift(U,q,state.len()),state)
}

fn premut_transp(prem:&Vec<f32>) -> Vec<(f32,f32)> {
    let mut swaps = vec![];
    for dest in 0..prem.len() {
	let mut src = prem[dest];
	while src < dest as f32 {
	    src = prem[src as usize];
	}

	if src < dest as f32 {
	    swaps.push((src as f32,dest as f32));
	} else {
	    swaps.push((dest as f32,src as f32));
	}
    }
 
    swaps
}

fn transpositions_to_adjacent_transpositions(transpositions: &Vec<(f32, f32)>) -> Vec<f32> {
    let expand_cons = |c: &(f32, f32)| -> Vec<f32> {
        if c.1 - c.0 == 1.0 {
            vec![c.0]
        } else {
            let trans: Vec<f32> = ((c.0 as i32)..(c.1 as i32)).collect::<Vec<i32>>().iter().map(|x| {*x as f32}).collect();
            trans.into_iter().chain(trans.iter().cloned().rev().collect::<Vec<f32>>().pop()).collect()
        }
    };

    transpositions.iter().map(|c| expand_cons(c)).flatten().collect()
}

fn apply_nq_gate() {todo!()}

fn apply_gate(state:&Vec<num::complex::Complex<f32>>,u:&Vec<Vec<f32>>,qb:&Vec<usize>) -> Vec<num::complex::Complex<f32>> {
    assert!(qb.len() == qubit_cnt(u[0].len() as u32).try_into().unwrap(), "invalid attempt to apply a gate");

    if qb.len() == 1 {
	return apply_sq_gate(state, u, qb[0]);
    }

    apply_nq_gate(state, u,qb)
}

fn sample(state:&Vec<num::complex::Complex<f32>>) -> f32 {
    let mut rng = rand::thread_rng();
    let mut r = num::complex::Complex::new(rng.gen_range(0.0_f32..1.0_f32),0.0);

    for i in 0..state.len() {
	r -= (state[i as usize].abs()).powf(2.0);

	if r.re < 0.0 {return i as f32};
    }

    return 0.0;
}

fn collapse(state:&Vec<num::complex::Complex<f32>>,basis: i32) -> Vec<num::complex::Complex<f32>> {
    state
	.iter()
	.enumerate()
	.map(|(i,_)| {
	    if i == basis as usize {
		return num::complex::Complex::new(1.0,0.0);
	    }
	    num::complex::Complex::new(0.0,0.0)
	})
	.collect()
}

fn observe(machine:&QuantumMachine) -> QuantumMachine {
    //let mut m = machine.clone();
    //let b = sample(&m.state);
    //m.state = collapse(&m.state,b as i32);
    //m.measurement = b;

    //return m;

    let r = sample(&machine.state);
    QuantumMachine {
	state: collapse(&machine.state,r as i32),
	measurement: r
    }
}

fn run_program(program:
	       &Vec<(&str,(Vec<Vec<f32>>,Vec<usize>))>,
	       machine:&QuantumMachine) -> QuantumMachine
{
    let mut m = machine.clone();
    let i = program
	.iter()
	//.enumerate()
	.for_each(|(inst,data)| {
	    //TODO make GATE and MEASURE an enum
	    match inst {
		&"GATE" => {
		    let (gate,qubits) = data;
		    //for i in 0..qubits.len()-1 {
			//m.state[qubits[i as usize]] = apply_gate(&m.state,&gate,&qubits);
		    //}
		    m.state = apply_gate(&m.state,&gate,&qubits);
		},
		&"MEASURE" => {
		    observe(&m);
		}
	    }
	});
    m
}

fn main() {
    println!("Hello, world!");
}
