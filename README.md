# gfe_bm_rust: grouped fixed-effects estimation in Rust

Eliot R. Weinstein · eliotw@uchicago.edu · erweinstein.com · November 2026

> Bonhomme, Stéphane and Elena Manresa (2015), "Grouped Patterns of Heterogeneity in Panel Data,"
> *Econometrica* 83(3), 1147–1184. doi:10.3982/ECTA11319

**Credit.** The grouped fixed-effects (GFE) estimator, its theory and its algorithms (Algorithm 1 and the variable
neighborhood search of Algorithm 2) are the work of **Stéphane Bonhomme and Elena Manresa**. This repository is a
Rust implementation of those algorithms, written after the R and Python translations of their replication code
([gfe_bm_R](https://github.com/erweinstein/gfe_bm_R), [gfe_bm_python](https://github.com/erweinstein/gfe_bm_python)),
which follow the original files line by line. It was not reviewed by the original authors; errors are mine. Please
cite their paper when you use it.

## Why a Rust version

For the Monte Carlo of the BIC choice of the number of groups (the paper's Table S.VIII), how good the search is
matters: a search that stops short at large G fits worse there, which raises the error variance estimated from
the G = 15 fit and with it the BIC's penalty, so the BIC chooses too few groups. The Python translation with
Algorithm 1 and 100 random starts chose about half a group fewer than published; on the authors' own samples and
groupings (`scripts/04d_authors_groupings.py` in gfe_bm_python), their groupings reproduce the published
distribution exactly and ours fall short at large G. Settling how good a fit is takes many more starts than an
interpreted language allows for 100 samples × 13 values of G; this code runs Algorithm 1 with 1,000 starts and
two variable neighborhood searches per sample and G in seconds to a minute.

## The model

For a balanced panel of N units over T periods with K covariates,

  y_it = x_it′θ + α_{g_i t} + v_it,

and the GFE estimate minimizes Q = Σ_i Σ_t (y_it − x_it′θ − α_{g_i t})² over θ, the group-time effects α and
the grouping g. Given g, θ comes from the within-(group, time) regression and α from the group-time means; given
(θ, α), each unit goes to its nearest group path.

## What is here

- `src/lib.rs`: the panel type; least squares given a grouping (`update`, `objective`); Algorithm 1 from a start
  (`alg1_run`, stopping when the assignment no longer changes) and with random starts (`alg1`; θ ~ N(0, I), α from
  G random units, as in the authors' code); a local search of single-unit moves (`local_search`); and the
  neighborhood search of Algorithm 2 (`vns`, or `vns_from` a given grouping).
- `src/main.rs`: `gfe_bm_rust table_s8`, which reads the authors' Monte Carlo files for G0 = 10 and reports, for
  each sample and G = 3 to 15, the objective of their grouping, of a local search and a neighborhood search from
  it, of Algorithm 1, and of Algorithm 2.
- `run_table_s8.sh`: runs it sample by sample, resumably, into `output/table_s8_search.csv`.
- `docs/what_rust_changed.md`: what the Rust version bought, with timings against the Python translation and a
  C++ (Rcpp) port of the same code.

Checks: on the authors' first sample, the least-squares step gives their grouping's objective to 13 digits
(18.7959402077 at G = 3, as in the Python translation), and from the Python translation's random starts Algorithm 1
converges to the same objectives (18.9191, 20.0658, 18.8115).

## Running it

Rust 1.97 (cargo); crates rand 0.8, rand_xoshiro 0.6 and rayon 1.10, fetched by cargo.

```
cargo build --release
./run_table_s8.sh <folder with inputMCXG10_BIC.txt, inputMCYG10_BIC.txt and assignmentG010G3.txt ... G15.txt>
```

The replication files are the authors' and are not included (see `NOTICE.md`); they are in the `MonteCarlo`
folder of `Bonhomme_Manresa_codes`, on [Stéphane Bonhomme's website](https://sites.google.com/site/stephanebonhommeresearch/).

## Results (Table S.VIII, G0 = 10; 100 samples; `scripts/summarize_table_s8.py`)

For each of the authors' 100 simulated samples and G = 3 to 15, the run computed the objective of the authors'
grouping, of a local search and a neighbourhood search started from it, of Algorithm 1 with 1,000 random starts,
and of Algorithm 2 from random starts (about 3.3 minutes per sample on 2 cores). The BIC choice of G, computed as
`BIC_G10_MC.m` does (σ̂² from each sample's G = 15 fit, averaged over the samples), with each set of objectives
(`output/table_s8_bic_choices.csv`):

| Objective used | σ̂² | G = 7 | G = 8 | G = 9 | G = 10 | G = 11 | G = 12 |
|---|---|---|---|---|---|---|---|
| The authors' groupings | 0.01556 | 0 | 10 | 42 | 42 | 6 | 0 |
| Algorithm 1, 1,000 starts | 0.01611 | 1 | 19 | 45 | 32 | 3 | 0 |
| Algorithm 2 from random starts | 0.01567 | 0 | 13 | 47 | 33 | 6 | 1 |
| The best of all searches | 0.01556 | 0 | 9 | 46 | 40 | 5 | 0 |
| Published (%) | | 0 | 10 | 42 | 42 | 6 | 0 |

- **The authors' groupings are good but not global optima.** Some search found a lower objective than theirs in
  6% of samples at G = 3, 26–37% at G = 6 to 14, and 1% at G = 15 (`output/table_s8_improvement.csv`; mean
  improvement 0.004–0.044, largest 0.46).
- **The published distribution survives better fits.** Because the improvements at G = 15 are negligible, σ̂²
  and the penalty are unchanged, and the best fits move only a few samples (G = 9 in 46 instead of 42, G = 10 in
  40 instead of 42).
- **The search matters through σ̂².** Algorithm 1 even with 1,000 starts never reaches the best objective from
  G = 13 on, so its σ̂² is 3.5% too high and its choices shift down (G = 10 in 32 samples). The neighbourhood
  search from random starts finds the best at G = 15 in only a third of samples. A strong search at the largest
  G, where σ̂² comes from, matters more for the BIC than a strong search at the candidate values.

The CSV holds objective values only; the authors' samples and groupings are not redistributed (`NOTICE.md`).
