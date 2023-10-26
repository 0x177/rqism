
fn trv(s:[[f32]]) -> Vec<Vec<f32>> {
    s.iter().map(|&e| e.to_vec()).collect()
}

fn operator_apply(matrix: &Vec<Vec<f32>>, vector: &Vec<num::complex::Complex<f32>>) -> Vec<num::complex::Complex<f32>> {
    let mut result = vec![0.0; matrix.len()];
    
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
	measureument: 0,
    }
}

struct QuantumMachine {
    state: Vec<num::complex::Complex<f32>>,
    measurement: u32, 
}

fn qubit_cnt(d: u32) -> u32 {
    return d.count_ones()+d.count_zeros();
}

fn lift(U: &Vec<Vec<f32>>, i: usize, n: usize) -> Vec<Vec<f32>> {
    kronecker_multiply(
	kron_exp(trv(IDENT),n-i-qubit_cnt(U[0].len())),
	kronecker_multiply(U,
			   kron_exp(trv(IDENT), i)
	)
    )
}

fn apply_sq_gate(state: &num::complex::Complex<f32>, U: &Vec<Vec<f32>>,q: usize) -> num::complex::Complex<f32> {
    operator_apply(lift(U,q,state.length),state)
}

fn premut_transp(prem:&Vec<f32>) -> Vec<(f32,f32)> {
    let mut swaps = vec![];
    for dest in 0..prem.len {
	let mut src = prem[dest];
	while src < dest {
	    src = prem[src];
	}
    }
    if src < dest {
	swaps.push((src,dest));
    } else {
	swaps.push((dest,src));
    }
}

fn transpositions_to_adjacent_transpositions(transpositions: &Vec<(f32, f32)>) -> Vec<f32> {
    let expand_cons = |c: &(f32, i32)| -> Vec<f32> {
        if c.1 - c.0 == 1.0 {
            vec![c.0]
        } else {
            let trans: Vec<f32> = (c.0..c.1).collect().map(|x| {x as f32});
            trans.into_iter().chain(trans.iter().cloned().rev().collect().pop()).collect()
        }
    }

    transpositions.iter().map(|c| expand_cons(c)).flatten().collect()
}

fn apply_nq_gate() {todo!()}

fn apply_gate(state:&num::complex::Complex<f32>,u:&Vec<Vec<f32>>,qb:&Vec<num::complex::Complex<f32>>) -> num::complex::Complex<f32> {
    assert_eq!(qubits.len() == qubit_cnt(U[0].len()), "invalid attempt to apply a gate");

    if qubits.len() == 1 {
	apply_sq_gate(state, u, qb[0])
    }

    apply_nq_gate(state, u,qb);
}

fn sample(state:&Vec<num::complex::Complex<f32>>) -> f32 {
    let mut r = num::complex::Complex::new(rng.gen_range(0.0,1.0),0.0);

    for i in 0..state.len() {
	r -= (state[i as usize].abs()).powf(2.0);

	if r < 0.0 {return i as f32};
    }

    return 0.0;
}

fn collapse(state:&Vec<num::complex::Complex<f32>>,basis: i32) -> Vec<num::complex::Complex<f32>> {
    state
	.iter()
	.enumerate()
	.for_each(|i,_| {
	    if i == basis {
		num::complex::Complex::new(1.0,0.0)
	    }
	    num::complex::Complex::new(0.0,0.0)
	})
	.collect()
}

fn main() {
    println!("Hello, world!");
}
