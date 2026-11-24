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
//! Options (after the nine arguments; added 24 November 2026):
//!   --out <file>     append the rows to <file> (header first if it is new), skipping (sample, G) pairs already in it,
//!                    so the whole table can be run in one call and the call stopped and restarted
//!   --incremental    use the local search through sufficient statistics (`local_search_incremental`) in place of
//!                    the full recomputation, in Q_local and inside both neighbourhood searches
//! The jobs (sample, G) go to the threads from one queue, largest G first, so a call over many samples keeps every
//! thread busy until the last few jobs; RAYON_NUM_THREADS sets the number of threads. Each job's random stream
//! depends only on its seed, sample and G, so the results do not depend on the threads or the order.
//! The replication files are the authors'; they are not part of this repository.
//!
//! Eliot R. Weinstein
//! eliotw@uchicago.edu
//! erweinstein.com

use gfe_bm_rust::Panel;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};

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
    let incremental = a.iter().any(|v| v == "--incremental");
    let out_file = a.iter().position(|v| v == "--out").map(|p| a[p + 1].clone());
    let header = "sample,G,Q_authors,Q_local,Q_vns_authors,Q_alg1,Q_vns,Q_best,secs";
    let mut done = std::collections::HashSet::new();
    let sink: Box<dyn Write + Send> = match &out_file {
        Some(f) => {
            if let Ok(txt) = std::fs::read_to_string(f) {
                for l in txt.lines().skip(1) { let v: Vec<&str> = l.split(',').collect();
                    if v.len() == 9 { done.insert((v[0].parse::<usize>().unwrap(), v[1].parse::<usize>().unwrap())); } }
            }
            let new = !std::path::Path::new(f).exists();
            let mut h = std::fs::OpenOptions::new().create(true).append(true).open(f).expect(f);
            if new { writeln!(h, "{}", header).unwrap(); }
            Box::new(h)
        }
        None => { if first == 1 { println!("{}", header); } Box::new(std::io::stdout()) }
    };
    let mut jobs: Vec<(usize, usize)> = (first - 1..last).flat_map(|r| (0..gs.len()).map(move |j| (r, j)))
        .filter(|&(r, j)| !done.contains(&(r + 1, gs[j]))).collect();
    jobs.sort_by(|x, y| y.1.cmp(&x.1).then(x.0.cmp(&y.0)));           // largest G first
    if out_file.is_some() { eprintln!("{} jobs to run ({} already in the file)", jobs.len(), done.len()); }
    let out = std::sync::Mutex::new(sink);
    let next = AtomicUsize::new(0);
    let run = |&(r, j): &(usize, usize)| {
        let t0 = std::time::Instant::now();
        let off = r * n * t;
        let p = Panel { n, t, k, y: (0..n * t).map(|i| y[off + i][0]).collect(), x: (0..n * t).flat_map(|i| x[off + i].clone()).collect() };
        let g = gs[j];
        let ga = &assign[j][r];
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed ^ ((r as u64) << 16) ^ g as u64);
        let q_auth = p.objective(ga, g);
        let q_loc = if incremental { p.local_search_incremental(ga.clone(), g) } else { p.local_search(ga.clone(), g) }.objective;
        let q_va = p.vns_from_with(ga.clone(), g, nm, im, &mut rng, incremental).objective;
        let q_a1 = p.alg1(g, s1, &mut rng).objective;
        let q_vns = p.vns_with(g, sv, nm, im, &mut rng, incremental).objective;
        let best = [q_auth, q_loc, q_va, q_a1, q_vns].iter().cloned().fold(f64::INFINITY, f64::min);
        let line = format!("{},{},{:.10},{:.10},{:.10},{:.10},{:.10},{:.10},{:.2}", r + 1, g, q_auth, q_loc, q_va, q_a1, q_vns, best, t0.elapsed().as_secs_f64());
        let mut o = out.lock().unwrap(); writeln!(o, "{}", line).unwrap(); o.flush().unwrap();
    };
    rayon::broadcast(|_| loop {
        let i = next.fetch_add(1, Ordering::SeqCst); if i >= jobs.len() { break; }
        run(&jobs[i]);
    });
}
