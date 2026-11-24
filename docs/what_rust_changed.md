# What the Rust version changed

Eliot R. Weinstein · eliotw@uchicago.edu · erweinstein.com · 18 November 2026; updated 24 November 2026

This note records what writing grouped fixed-effects (GFE) estimation in Rust bought, measured against the
Python translation ([gfe_bm_python](https://github.com/erweinstein/gfe_bm_python)) and against a C++ port of
the same code through Rcpp. The estimator and its algorithms are Bonhomme and Manresa's (2015); see the
README for credit.

## The result it made possible

The Monte Carlo of the paper's Table S.VIII chooses the number of groups by the BIC of the supplement's (S20):

BIC(G) = Q(G)/NT + σ̂² (GT + N + K) ln(NT)/NT,  with σ̂² = Q(G_max)/(NT − G_max T − N − K) from the G = 15 fit.

So a search that stops short at G = 15 raises σ̂², inflates the penalty and pushes the chosen G down. The
objective Q(θ, α, g) = Σᵢ Σₜ (yᵢₜ − xᵢₜ′θ − α_{gᵢ t})² is minimized over all G-partitions of the N = 90 units;
Algorithm 1 is a k-means-type descent, so it needs many random starts, and the variable neighbourhood search adds
jumps and single-unit moves. The Rust code ran 1,300 independent (sample, G) jobs, each with 1,000 Algorithm 1
starts and two neighbourhood searches, in parallel (about 3.3 minutes per sample on 2 cores; with the incremental
local search added on 24 November, about 5 seconds per sample on one thread). That depth showed:

- some search beat the authors' grouping in 6% of samples at G = 3, 26–37% at G = 6–14 and 1% at G = 15 (mean
  gain 0.004–0.044; largest 0.46);
- because the G = 15 fits barely move, σ̂² and the published BIC choices survive (G = 9 in 46 rather than 42
  samples, G = 10 in 40 rather than 42);
- Algorithm 1 alone never reaches the best objective from G = 13 on, so it overstates σ̂² by 3.5% and shifts the
  choice down (G = 10 in 32 samples).

The README's results section has the full table.

![Share of samples in which a search found a lower objective than the authors' grouping, by G](fig_gfe.png)

## Measured speed

One thread, on an idle 2-core Linux container, all implementations in one session; medians over 5 samples × 3
runs of the compute time each program reports (reading the input files dominates the whole-process time).
G = 10. Algorithm 1 was timed on 18 November and the local search on 24 November.

| Task (unit of work) | Python (NumPy) | Rcpp | Rust | Rust, incremental |
|---|---|---|---|---|
| Algorithm 1 (one random start; 1,000 starts per run) | 1.81 ms | 0.165 ms | 0.209 ms | |
| Local search from the authors' grouping (one search) | 12.0 ms (incremental) | 6.23 ms | 12.2 ms | 0.148 ms |
| The Table S.VIII search for one sample (13 jobs, G = 3–15; samples 1 and 2, one run each) | | | 63.6 s | 5.0–5.4 s |

Versions: Python 3.13 with NumPy 2.5; Rcpp 1.1.22 with g++ 13.3 at R's default -O2; Rust 1.97 with opt-level 3
and link-time optimization.

- **The local search, like for like.** Python's local search updates sufficient statistics for each trial move;
  the Rcpp port and the original Rust search recompute the objective. On 18 November that difference was enough
  for Python to keep pace with Rust (10.1 against 11.3 ms then; 12.0 against 12.2 ms on 24 November). The
  incremental search in Rust (`local_search_incremental`, added 24 November) runs the same algorithm as Python's
  and takes 0.148 ms: 81 times faster than Python's, 83 times faster than Rust's full recomputation and 42 times
  faster than the Rcpp port. All four reach the same objectives on the 5 samples (8.8795721449 on sample 1, for
  example). A trial move costs O(T(K + 1)²) arithmetic, a few hundred operations, so Python spends nearly all its
  time on interpreter and NumPy call overhead; with the algorithm equal, that overhead is the 81×.
- **The whole search.** The neighbourhood search calls the local search after every jump, so the incremental
  search cuts a sample's 13 jobs from 63.6 s to about 5 s on one thread (12×), and Algorithm 1's 1,000 starts now
  take most of the time. Over all 1,300 jobs the six objectives agree with `output/table_s8_search.csv` in every
  printed digit.
- **Algorithm 1.** Rust is 8.7× faster than Python, and the Rcpp version 1.3× faster than Rust. The Rcpp
  version is a line-for-line port of `src/lib.rs`, written for this comparison and not in a repository. From
  1,000 starts its Algorithm 1 finds the same best objective as Rust on 10 of 20 (sample, G) cases, a lower one on
  4 and a higher one on 6, because its random starts come from R's generator. I have not diagnosed why the Rust
  version is slower; bounds checks on the index arithmetic in the inner loops are one candidate. R's -O3 was
  slower than its default -O2 here (0.204 against 0.165 ms).

So, for the full-recomputation kernels, the speed is that of compiled code, not of Rust in particular. With a
cheap per-step update, the gap between a compiled language and an interpreted one is far larger than for
vectorized work, and that is where the 81× comes from.

What Rust added beyond speed was the parallel driver: rayon runs work only on types the compiler can prove safe to
share (Send and Sync), so the 1,300 jobs run in parallel without data races (Jung et al. 2018), and each
(sample, G) job owns a Xoshiro256++ stream (Blackman and Vigna 2021) seeded from its indices, so results do not
depend on the number of threads. The rows are written as jobs finish, so only their order varies between runs.

**More cores.** `run_table_s8.sh` calls the program once per sample, so a call has only 13 jobs and waits for
its longest (G = 15). Simulating that schedule from the measured job times gives at most about 5.9× on any
number of cores. Since 24 November, `run_table_s8_all.sh` puts all 1,300 jobs in one queue, largest G first;
the same simulation then gives 2.0× on 2 cores, 4.0× on 4, 15.9× on 16 and 63× on 64, because the longest job
(under a second with the incremental search) is small next to the total (about 8.5 minutes on one thread). The
simulation leaves out memory bandwidth, caches and clock speed, so it is an upper bound; beyond 2 cores nothing
was measured.

## Possible problems

- **The two Algorithm 1 versions are not identical.** When a group would be empty, Rust moves the unit farthest
  from its group path into it, while the Python translation abandons the start. From 1,000 starts at G = 10 the two
  reach the same objective on 2 of 5 samples and Rust a lower one on the other 3 (for example 8.5261 against
  8.5334). That could be the repair or the different random streams; I have not separated them.
- **The benchmarks are small**: one machine, 2 cores, short runs. The job times behind the 1,300-job simulation
  come from the check run of the incremental search, which ran on one thread while another program used both
  cores, so they are slower than idle times (8.5 minutes in all, against about 5.2 s × 100 samples idle); the
  simulation needs only their relative sizes.
- **The 81× is against my Python code as written.** Python's local search loops over units and groups in Python
  and calls NumPy for each 3 × 3 update. Evaluating all G − 1 candidate groups for a unit in one vectorized call,
  or compiling the loop with Numba, would narrow the gap; I have not tried either.
- **A metadata inconsistency.** `Cargo.toml` declares `license = "MIT"`, while `LICENSE` and `NOTICE.md` say BSD
  3-Clause. One of them should be fixed before anyone depends on the crate.

## References

- Blackman, D. and S. Vigna (2021). Scrambled linear pseudorandom number generators. *ACM Transactions on
  Mathematical Software* 47(4), 1–32. doi:10.1145/3460772
- Bonhomme, S. and E. Manresa (2015). Grouped patterns of heterogeneity in panel data. *Econometrica* 83(3),
  1147–1184. doi:10.3982/ECTA11319
- Jung, R., J.-H. Jourdan, R. Krebbers and D. Dreyer (2018). RustBelt: securing the foundations of the Rust
  programming language. *Proceedings of the ACM on Programming Languages* 2(POPL), article 66.
  doi:10.1145/3158154
