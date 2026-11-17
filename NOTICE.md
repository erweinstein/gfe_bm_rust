# Notice: scope of the license and attribution

gfe_bm_rust: grouped fixed-effects estimation in Rust, after the replication code for Bonhomme and Manresa (2015),
"Grouped Patterns of Heterogeneity in Panel Data"

Copyright (c) 2026 Eliot R. Weinstein
eliotw@uchicago.edu
https://erweinstein.com

The Rust code is licensed under the BSD 3-Clause License (see `LICENSE`), as is the Python translation,
gfe_bm_python (https://github.com/erweinstein/gfe_bm_python).

That license covers only the material written by Eliot R. Weinstein in this repository: the Rust code and the
documentation. It does not cover, and grants no rights in, the following, which remain the work and property of
their authors and are subject to their own terms:

- The grouped fixed-effects estimator, its theory, the computational algorithms (Algorithm 1 and the variable
  neighborhood search), and the original MATLAB, Fortran, and Stata replication code, by Stéphane Bonhomme and
  Elena Manresa: Bonhomme, S. and E. Manresa (2015), "Grouped Patterns of Heterogeneity in Panel Data,"
  *Econometrica* 83(3), 1147–1184. doi:10.3982/ECTA11319
- The simulated samples and groupings of their replication package, which the `table_s8` command reads from a
  local copy. They are not included in this repository.

If you use this code, please cite Bonhomme and Manresa (2015).
