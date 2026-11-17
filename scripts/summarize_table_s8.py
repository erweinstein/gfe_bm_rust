"""summarize_table_s8.py -- the BIC choice of G (Bonhomme and Manresa 2015, Supplementary Appendix eq. S20) from the
objectives in output/table_s8_search.csv, as BIC_G10_MC.m computes it: sigma2_hat = Q(15)/(NT - 15T - N - K) in
each sample, averaged over the samples; BIC(G) = Q(G)/NT + sigma2_hat (GT + N + K)/NT ln(NT); G = 3..15.
Writes output/table_s8_bic_choices.csv and output/table_s8_improvement.csv.

Eliot R. Weinstein · eliotw@uchicago.edu · erweinstein.com"""
import numpy as np
import pandas as pd

N, T, K = 90, 7, 2
NT = N * T
d = pd.read_csv("output/table_s8_search.csv")
rows = []
labels = {"Q_authors": "The authors' groupings", "Q_alg1": "Algorithm 1, 1,000 starts",
          "Q_vns": "Algorithm 2 (VNS) from random starts", "Q_best": "The best of all searches"}
for col, lab in labels.items():
    s2 = np.mean([g.set_index("G").loc[15, col] / (NT - 15 * T - N - K) for _, g in d.groupby("sample")])
    ch = []
    for _, g in d.groupby("sample"):
        g = g.set_index("G")
        bic = {G: g.loc[G, col] / NT + s2 * (G * T + N + K) / NT * np.log(NT) for G in g.index}
        ch.append(min(bic, key=bic.get))
    vc = pd.Series(ch).value_counts()
    rows.append({"objective": lab, "sigma2_hat": round(s2, 5), **{f"G{G}": int(vc.get(G, 0)) for G in range(6, 13)}})
rows.append({"objective": "Published (Table S.VIII, %)", "sigma2_hat": np.nan, "G6": 0, "G7": 0, "G8": 10, "G9": 42, "G10": 42,
             "G11": 6, "G12": 0})
out = pd.DataFrame(rows)
out.to_csv("output/table_s8_bic_choices.csv", index=False)
print(out.to_string(index=False))
d["gain"] = d.Q_authors - d.Q_best
imp = d.groupby("G").agg(share_improved=("gain", lambda x: (x > 1e-8).mean()), mean_gain=("gain", "mean"),
                         max_gain=("gain", "max"))
for c in ["Q_local", "Q_vns_authors", "Q_alg1", "Q_vns"]:
    imp["best_" + c] = ((d[c] - d.Q_best).abs() < 1e-8).groupby(d["G"]).mean()
imp = imp.round(4)
imp.to_csv("output/table_s8_improvement.csv")
print(imp.to_string())
