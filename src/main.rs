//! gfe_bm_rust command line.
//!
//!   gfe_bm_rust table_s8 <dir> <first_sample> <last_sample> <alg1_starts> <vns_starts> <neighmax> <itermax> <seed>
//!
//! Reads inputMCXG10_BIC.txt, inputMCYG10_BIC.txt and assignmentG010G3.txt ... assignmentG010G15.txt from <dir>
//! (the authors' Monte Carlo files for Table S.VIII, G0 = 10, from their replication package; N = 90, T = 7, K = 2),
//! and for samples first..last (1-based, inclusive) and G = 3..15 prints, as CSV rows when each finishes:
//!   Q_authors       the objective of the authors' grouping (theta and alpha concentrated out)
//!   Q_local         a local search (single-unit moves) from the authors' grouping
//!   Q_vns_authors   the neighbourhood search of Algorithm 2 from the authors' grouping
//!   Q_alg1          Algorithm 1 with <alg1_starts> random starts
//!   Q_vns           Algorithm 2: <vns_starts> Algorithm 1 starts, then the neighbourhood search
//!   Q_best          the smallest of these
//! The replication files are the authors'; they are not part of this repository.
//!
//! Eliot R. Weinstein
//! eliotw@uchicago.edu
//! erweinstein.com

use gfe_bm_rust::Panel;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::io::Write;

fn read_matrix(path: &str) -> Vec<Vec<f64>> {
    std::fs::read_to_string(path).expect(path).lines().filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(|v| v.parse::<f64>().unwrap()).collect()).collect()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 10 || a[1] != "table_s8" {
        eprintln!("usage: gfe_bm_rust table_s8 <dir> <first_sample> <last_sample> <alg1_starts> <vns_starts> <neighmax> <itermax> <seed>");
        std::process::exit(2);
    }
    let dir = &a[2]; let first: usize = a[3].parse().unwrap(); let last: usize = a[4].parse().unwrap();
    let s1: usize = a[5].parse().unwrap(); let sv: usize = a[6].parse().unwrap();
    let nm: usize = a[7].parse().unwrap(); let im: usize = a[8].parse().unwrap(); let seed: u64 = a[9].parse().unwrap();
    let (n, t, k) = (90usize, 7usize, 2usize);
    let x = read_matrix(&format!("{}/inputMCXG10_BIC.txt", dir));
    let y = read_matrix(&format!("{}/inputMCYG10_BIC.txt", dir));
    let gs: Vec<usize> = (3..=15).collect();
    let assign: Vec<Vec<Vec<usize>>> = gs.iter().map(|g| read_matrix(&format!("{}/assignmentG010G{}.txt", dir, g)).into_iter()
        .map(|r| r.into_iter().map(|v| v as usize - 1).collect()).collect()).collect();
    if first == 1 { println!("sample,G,Q_authors,Q_local,Q_vns_authors,Q_alg1,Q_vns,Q_best,secs"); }
    let jobs: Vec<(usize, usize)> = (first - 1..last).flat_map(|r| (0..gs.len()).map(move |j| (r, j))).collect();
    let out = std::sync::Mutex::new(std::io::stdout());
    jobs.par_iter().for_each(|&(r, j)| {
        let t0 = std::time::Instant::now();
        let off = r * n * t;
        let p = Panel { n, t, k, y: (0..n * t).map(|i| y[off + i][0]).collect(), x: (0..n * t).flat_map(|i| x[off + i].clone()).collect() };
        let g = gs[j];
        let ga = &assign[j][r];
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed ^ ((r as u64) << 16) ^ g as u64);
        let q_auth = p.objective(ga, g);
        let q_loc = p.local_search(ga.clone(), g).objective;
        let q_va = p.vns_from(ga.clone(), g, nm, im, &mut rng).objective;
        let q_a1 = p.alg1(g, s1, &mut rng).objective;
        let q_vns = p.vns(g, sv, nm, im, &mut rng).objective;
        let best = [q_auth, q_loc, q_va, q_a1, q_vns].iter().cloned().fold(f64::INFINITY, f64::min);
        let line = format!("{},{},{:.10},{:.10},{:.10},{:.10},{:.10},{:.10},{:.2}", r + 1, g, q_auth, q_loc, q_va, q_a1, q_vns, best, t0.elapsed().as_secs_f64());
        let mut o = out.lock().unwrap(); writeln!(o, "{}", line).unwrap(); o.flush().unwrap();
    });
}
