//! gfe_bm_rust -- grouped fixed-effects estimation (Bonhomme and Manresa 2015, "Grouped patterns of heterogeneity in
//! panel data", Econometrica 83(3), 1147-1184) for balanced panels with covariates:
//!   y_it = x_it' theta + alpha_{g_i t} + v_it,   minimize Q = sum_i sum_t (y_it - x_it' theta - alpha_{g_i t})^2
//! over theta, the group-time effects alpha and the grouping g. Algorithm 1 (alternate assignment and least
//! squares, many random starts), a local search (single-unit moves to a local minimum), and Algorithm 2 (variable
//! neighborhood search). A Rust companion to the R and Python translations (gfe_bm_R, gfe_bm_python), written for
//! speed: thousands of starts per sample for the Monte Carlo of the BIC choice of G (Table S.VIII).
//!
//! Eliot R. Weinstein
//! eliotw@uchicago.edu
//! erweinstein.com

use rand::prelude::*;

/// A balanced panel: y[i * t + s], x[(i * t + s) * k + j]
#[derive(Clone)]
pub struct Panel { pub n: usize, pub t: usize, pub k: usize, pub y: Vec<f64>, pub x: Vec<f64> }

#[derive(Clone, Debug)]
pub struct Fit { pub groups: Vec<usize>, pub theta: Vec<f64>, pub alpha: Vec<f64>, pub objective: f64 }

fn solve(a: &mut [f64], b: &mut [f64], k: usize) -> Option<Vec<f64>> {
    // Gaussian elimination with partial pivoting on the k x k system a x = b
    for c in 0..k {
        let mut p = c; for r in c..k { if a[r * k + c].abs() > a[p * k + c].abs() { p = r; } }
        if a[p * k + c].abs() < 1e-12 { return None; }
        if p != c { for j in 0..k { a.swap(c * k + j, p * k + j); } b.swap(c, p); }
        for r in (c + 1)..k { let f = a[r * k + c] / a[c * k + c];
            for j in c..k { a[r * k + j] -= f * a[c * k + j]; } b[r] -= f * b[c]; }
    }
    let mut x = vec![0.0; k];
    for c in (0..k).rev() { let mut s = b[c]; for j in (c + 1)..k { s -= a[c * k + j] * x[j]; } x[c] = s / a[c * k + c]; }
    Some(x)
}

impl Panel {
    /// least squares given a grouping: theta from the within-group-time regression, alpha = group-time means of
    /// y - x' theta; returns the fit (objective included). Empty groups get alpha = 0.
    pub fn update(&self, groups: &[usize], g: usize) -> Option<Fit> {
        let (n, t, k) = (self.n, self.t, self.k);
        let mut cnt = vec![0.0; g]; for &gi in groups { cnt[gi] += 1.0; }
        let mut ym = vec![0.0; g * t]; let mut xm = vec![0.0; g * t * k];
        for i in 0..n { let gi = groups[i]; for s in 0..t {
            ym[gi * t + s] += self.y[i * t + s];
            for j in 0..k { xm[(gi * t + s) * k + j] += self.x[(i * t + s) * k + j]; } } }
        for gg in 0..g { if cnt[gg] > 0.0 { for s in 0..t { ym[gg * t + s] /= cnt[gg]; for j in 0..k { xm[(gg * t + s) * k + j] /= cnt[gg]; } } } }
        let mut a = vec![0.0; k * k]; let mut b = vec![0.0; k];
        let mut dx = vec![0.0; k];
        for i in 0..n { let gi = groups[i]; for s in 0..t {
            for j in 0..k { dx[j] = self.x[(i * t + s) * k + j] - xm[(gi * t + s) * k + j]; }
            let dy = self.y[i * t + s] - ym[gi * t + s];
            for j in 0..k { b[j] += dx[j] * dy; for l in 0..k { a[j * k + l] += dx[j] * dx[l]; } } } }
        let theta = if k > 0 { solve(&mut a, &mut b, k)? } else { Vec::new() };
        let mut alpha = vec![0.0; g * t];
        for gg in 0..g { for s in 0..t { let mut v = ym[gg * t + s]; for j in 0..k { v -= xm[(gg * t + s) * k + j] * theta[j]; } alpha[gg * t + s] = v; } }
        let objective = self.objective_given(groups, &theta, &alpha, t);
        Some(Fit { groups: groups.to_vec(), theta, alpha, objective })
    }
    fn resid(&self, i: usize, s: usize, theta: &[f64]) -> f64 {
        let mut r = self.y[i * self.t + s]; for j in 0..self.k { r -= self.x[(i * self.t + s) * self.k + j] * theta[j]; } r
    }
    fn objective_given(&self, groups: &[usize], theta: &[f64], alpha: &[f64], t: usize) -> f64 {
        let mut q = 0.0;
        for i in 0..self.n { for s in 0..t { let e = self.resid(i, s, theta) - alpha[groups[i] * t + s]; q += e * e; } }
        q
    }
    /// the objective of a grouping with theta and alpha concentrated out
    pub fn objective(&self, groups: &[usize], g: usize) -> f64 { self.update(groups, g).map_or(f64::INFINITY, |f| f.objective) }

    /// assign every unit to its nearest group-time path given theta and alpha; never leave a group empty (a group
    /// that would be empty takes the unit farthest from its assigned path)
    fn assign(&self, theta: &[f64], alpha: &[f64], g: usize) -> Vec<usize> {
        let (n, t) = (self.n, self.t);
        let mut groups = vec![0; n]; let mut dist = vec![0.0; n];
        for i in 0..n {
            let mut best = (f64::INFINITY, 0);
            for gg in 0..g { let mut d = 0.0; for s in 0..t { let e = self.resid(i, s, theta) - alpha[gg * t + s]; d += e * e; }
                if d < best.0 { best = (d, gg); } }
            groups[i] = best.1; dist[i] = best.0;
        }
        let mut cnt = vec![0usize; g]; for &gi in &groups { cnt[gi] += 1; }
        for gg in 0..g { if cnt[gg] == 0 {
            let far = (0..n).filter(|&i| cnt[groups[i]] > 1).max_by(|&a, &b| dist[a].partial_cmp(&dist[b]).unwrap()).unwrap();
            cnt[groups[far]] -= 1; groups[far] = gg; cnt[gg] = 1; dist[far] = 0.0; } }
        groups
    }

    /// Algorithm 1 from given theta and alpha, to convergence
    pub fn alg1_run(&self, mut theta: Vec<f64>, mut alpha: Vec<f64>, g: usize, max_iter: usize) -> Option<Fit> {
        // stop when the assignment no longer changes (the original code's deltapar == 0)
        let mut prev: Option<Vec<usize>> = None; let mut fit = None;
        for _ in 0..max_iter {
            let groups = self.assign(&theta, &alpha, g);
            if prev.as_ref() == Some(&groups) { break; }
            let f = self.update(&groups, g)?;
            theta = f.theta.clone(); alpha = f.alpha.clone();
            prev = Some(groups); fit = Some(f);
        }
        fit
    }

    /// a random start as in the authors' code: theta ~ N(0, sd^2) per component, alpha = y - x' theta of G random units
    pub fn random_start<R: Rng>(&self, g: usize, rng: &mut R, theta_sd: f64) -> (Vec<f64>, Vec<f64>) {
        let theta: Vec<f64> = (0..self.k).map(|_| theta_sd * rng.sample::<f64, _>(rand_distr_normal())).collect();
        let units: Vec<usize> = rand::seq::index::sample(rng, self.n, g).into_vec();
        let mut alpha = vec![0.0; g * self.t];
        for (gg, &i) in units.iter().enumerate() { for s in 0..self.t { alpha[gg * self.t + s] = self.resid(i, s, &theta); } }
        (theta, alpha)
    }

    /// Algorithm 1 with many random starts; the best fit
    pub fn alg1<R: Rng>(&self, g: usize, starts: usize, rng: &mut R) -> Fit {
        let mut best: Option<Fit> = None;
        for _ in 0..starts {
            let (th, al) = self.random_start(g, rng, 1.0);
            if let Some(f) = self.alg1_run(th, al, g, 1000) {
                if best.as_ref().map_or(true, |b| f.objective < b.objective) { best = Some(f); }
            }
        }
        best.expect("no start converged")
    }

    /// local search: move single units between groups while the concentrated objective falls (no empty groups)
    pub fn local_search(&self, mut groups: Vec<usize>, g: usize) -> Fit {
        let mut obj = self.objective(&groups, g);
        let mut cnt = vec![0usize; g]; for &gi in &groups { cnt[gi] += 1; }
        loop {
            let mut improved = false;
            for i in 0..self.n {
                let g0 = groups[i]; if cnt[g0] <= 1 { continue; }
                let mut best = (obj, g0);
                for gg in 0..g { if gg == g0 { continue; }
                    groups[i] = gg; let o = self.objective(&groups, g);
                    if o < best.0 - 1e-10 * best.0.abs().max(1.0) { best = (o, gg); } }
                groups[i] = best.1;
                if best.1 != g0 { cnt[g0] -= 1; cnt[best.1] += 1; obj = best.0; improved = true; }
            }
            if !improved { break; }
        }
        self.update(&groups, g).unwrap()
    }

    /// the same local search with the objective computed from sufficient statistics, as in gfe_bm_python's
    /// gfe_local_search (and the authors' Brusco_cov.f90): with z_it = (y_it, x_it')', S_gs = sum_{i in g} z_is and
    /// n_gs the group sizes, the within-(group, time) cross products are M = sum_it z_it z_it' - sum_gs S_gs S_gs' / n_gs
    /// and Q = M_yy - M_yx M_xx^{-1} M_xy. A move changes two groups' S, so it costs O(T (K+1)^2) operations
    /// instead of a pass over the whole panel. Same moves, same order and same tolerance as `local_search`; the
    /// objective it returns is recomputed directly (by `update`) at the end. Added 24 November 2026.
    pub fn local_search_incremental(&self, mut groups: Vec<usize>, g: usize) -> Fit {
        let (n, t, k) = (self.n, self.t, self.k); let k1 = k + 1; let kk = k1 * k1;
        let z = |i: usize, s: usize, c: usize| if c == 0 { self.y[i * t + s] } else { self.x[(i * t + s) * k + c - 1] };
        let mut w = vec![0.0; kk];
        for i in 0..n { for s in 0..t { for a in 0..k1 { for b in 0..k1 { w[a * k1 + b] += z(i, s, a) * z(i, s, b); } } } }
        let mut sg = vec![0.0; g * t * k1]; let mut cnt = vec![0usize; g];
        for i in 0..n { cnt[groups[i]] += 1; for s in 0..t { for c in 0..k1 { sg[(groups[i] * t + s) * k1 + c] += z(i, s, c); } } }
        // between cross products of one group from its sums (t x k1) and size
        let between = |sums: &[f64], m: usize, out: &mut [f64]| {
            for v in out.iter_mut() { *v = 0.0; }
            if m == 0 { return; }
            let inv = 1.0 / m as f64;
            for s in 0..t { let r = &sums[s * k1..(s + 1) * k1];
                for a in 0..k1 { let ra = r[a] * inv; for b in a..k1 { out[a * k1 + b] += ra * r[b]; } } }
            for a in 0..k1 { for b in 0..a { out[a * k1 + b] = out[b * k1 + a]; } }
        };
        let mut bg = vec![0.0; g * kk];
        for gg in 0..g { let (lo, hi) = (gg * t * k1, (gg + 1) * t * k1); let sums = sg[lo..hi].to_vec();
            between(&sums, cnt[gg], &mut bg[gg * kk..(gg + 1) * kk]); }
        let mut m = w.clone(); for gg in 0..g { for e in 0..kk { m[e] -= bg[gg * kk + e]; } }
        let obj_of = |m: &[f64]| -> f64 {
            match k {
                0 => m[0],
                1 => if m[k1 + 1] <= 1e-12 * m[0].abs() { f64::INFINITY } else { m[0] - m[k1] * m[k1] / m[k1 + 1] },
                2 => { let (a, b, c) = (m[k1 + 1], m[k1 + 2], m[2 * k1 + 2]); let det = a * c - b * b;
                       if det <= 1e-12 * a * c { return f64::INFINITY; }
                       let (u, v) = (m[k1], m[2 * k1]); m[0] - (c * u * u - 2.0 * b * u * v + a * v * v) / det }
                _ => { let mut a = vec![0.0; k * k]; let mut b = vec![0.0; k];
                       for r in 0..k { b[r] = m[(r + 1) * k1]; for c in 0..k { a[r * k + c] = m[(r + 1) * k1 + c + 1]; } }
                       let bb = b.clone();
                       match solve(&mut a, &mut b, k) { Some(x) => m[0] - (0..k).map(|r| bb[r] * x[r]).sum::<f64>(), None => f64::INFINITY } }
            }
        };
        let mut obj = obj_of(&m);
        let mut zi = vec![0.0; t * k1]; let mut s0 = vec![0.0; t * k1]; let mut s1 = vec![0.0; t * k1];
        let mut b0 = vec![0.0; kk]; let mut b1 = vec![0.0; kk]; let mut mn = vec![0.0; kk];
        let mut keep_s = vec![0.0; t * k1]; let mut keep_b = vec![0.0; kk]; let mut keep_m = vec![0.0; kk];
        loop {
            let mut improved = false;
            for i in 0..n {
                let g0 = groups[i]; if cnt[g0] <= 1 { continue; }
                for s in 0..t { for c in 0..k1 { zi[s * k1 + c] = z(i, s, c); } }
                for e in 0..t * k1 { s0[e] = sg[g0 * t * k1 + e] - zi[e]; }
                between(&s0, cnt[g0] - 1, &mut b0);
                let mut best = (obj, g0);
                for gg in 0..g { if gg == g0 { continue; }
                    for e in 0..t * k1 { s1[e] = sg[gg * t * k1 + e] + zi[e]; }
                    between(&s1, cnt[gg] + 1, &mut b1);
                    for e in 0..kk { mn[e] = m[e] + bg[g0 * kk + e] + bg[gg * kk + e] - b0[e] - b1[e]; }
                    let o = obj_of(&mn);
                    if o < best.0 - 1e-10 * best.0.abs().max(1.0) { best = (o, gg);
                        keep_s.copy_from_slice(&s1); keep_b.copy_from_slice(&b1); keep_m.copy_from_slice(&mn); }
                }
                if best.1 != g0 {
                    let g1 = best.1;
                    sg[g0 * t * k1..(g0 + 1) * t * k1].copy_from_slice(&s0); bg[g0 * kk..(g0 + 1) * kk].copy_from_slice(&b0);
                    sg[g1 * t * k1..(g1 + 1) * t * k1].copy_from_slice(&keep_s); bg[g1 * kk..(g1 + 1) * kk].copy_from_slice(&keep_b);
                    m.copy_from_slice(&keep_m); obj = best.0;
                    groups[i] = g1; cnt[g0] -= 1; cnt[g1] += 1; improved = true;
                }
            }
            if !improved { break; }
        }
        self.update(&groups, g).unwrap()
    }

    /// `vns_from` with a choice of local search (`incremental` = true uses `local_search_incremental`)
    pub fn vns_from_with<R: Rng>(&self, groups: Vec<usize>, g: usize, neighmax: usize, itermax: usize, rng: &mut R, incremental: bool) -> Fit {
        let ls = |gr: Vec<usize>| if incremental { self.local_search_incremental(gr, g) } else { self.local_search(gr, g) };
        let mut best = ls(groups);
        let mut stale = 0;
        while stale < itermax {
            let mut improved = false;
            let mut nb = 1;
            while nb <= neighmax {
                let mut cand = best.groups.clone();
                let mut cnt = vec![0usize; g]; for &gi in &cand { cnt[gi] += 1; }
                for _ in 0..nb {
                    let i = rng.gen_range(0..self.n); if cnt[cand[i]] <= 1 { continue; }
                    let mut gg = rng.gen_range(0..g); while gg == cand[i] { gg = rng.gen_range(0..g); }
                    cnt[cand[i]] -= 1; cand[i] = gg; cnt[gg] += 1;
                }
                let f0 = self.update(&cand, g).unwrap();
                let f1 = self.alg1_run(f0.theta, f0.alpha, g, 1000).unwrap_or(best.clone());
                let f2 = ls(f1.groups);
                if f2.objective < best.objective - 1e-10 * best.objective.abs().max(1.0) { best = f2; improved = true; nb = 1; }
                else { nb += 1; }
            }
            if improved { stale = 0; } else { stale += 1; }
        }
        best
    }

    /// `vns` with a choice of local search
    pub fn vns_with<R: Rng>(&self, g: usize, starts: usize, neighmax: usize, itermax: usize, rng: &mut R, incremental: bool) -> Fit {
        let a1 = self.alg1(g, starts, rng);
        self.vns_from_with(a1.groups, g, neighmax, itermax, rng, incremental)
    }

    /// Algorithm 2 (variable neighborhood search), as in the authors' Fortran: from the best of `starts` Algorithm 1
    /// runs, repeat: jump to a random neighbour (relocate n units, n = 1..neighmax), run Algorithm 1 and the local
    /// search from it, keep it if better, else widen the neighbourhood; stop after itermax rounds without improvement
    pub fn vns<R: Rng>(&self, g: usize, starts: usize, neighmax: usize, itermax: usize, rng: &mut R) -> Fit {
        let a1 = self.alg1(g, starts, rng);
        self.vns_from(a1.groups, g, neighmax, itermax, rng)
    }

    /// the neighbourhood search of Algorithm 2 started from a given grouping (e.g. someone else's solution)
    pub fn vns_from<R: Rng>(&self, groups: Vec<usize>, g: usize, neighmax: usize, itermax: usize, rng: &mut R) -> Fit {
        let mut best = self.local_search(groups, g);
        let mut stale = 0;
        while stale < itermax {
            let mut improved = false;
            let mut nb = 1;
            while nb <= neighmax {
                let mut cand = best.groups.clone();
                let mut cnt = vec![0usize; g]; for &gi in &cand { cnt[gi] += 1; }
                for _ in 0..nb {
                    let i = rng.gen_range(0..self.n); if cnt[cand[i]] <= 1 { continue; }
                    let mut gg = rng.gen_range(0..g); while gg == cand[i] { gg = rng.gen_range(0..g); }
                    cnt[cand[i]] -= 1; cand[i] = gg; cnt[gg] += 1;
                }
                let f0 = self.update(&cand, g).unwrap();
                let f1 = self.alg1_run(f0.theta, f0.alpha, g, 1000).unwrap_or(best.clone());
                let f2 = self.local_search(f1.groups, g);
                if f2.objective < best.objective - 1e-10 * best.objective.abs().max(1.0) { best = f2; improved = true; nb = 1; }
                else { nb += 1; }
            }
            if improved { stale = 0; } else { stale += 1; }
        }
        best
    }
}

fn rand_distr_normal() -> StandardNormalBoxMuller { StandardNormalBoxMuller }
/// a standard normal draw (Box-Muller), to avoid a dependency for one distribution
pub struct StandardNormalBoxMuller;
impl Distribution<f64> for StandardNormalBoxMuller {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        let u1: f64 = rng.gen_range(f64::MIN_POSITIVE..1.0); let u2: f64 = rng.gen();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}
