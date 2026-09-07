"""
Extends the Dyson-repulsion eigenvalue-trajectory diagnostic (scripts
05 and 17) from TNBC-only to LumA -- the "natural next step" flagged
explicitly in the Discussion and Conclusions: the reproducibility check
closed the single-cohort/single-order gap for TNBC, but the diagnostic
itself had never been run on the other real cohort in this paper.

Same method exactly, no new machinery: LumA's own 60 most-variable
genes, one fixed random patient-inclusion order first (comparable to
Fig. dyson/avoided), then N_ORDERS=20 independent random orders for a
reproducibility check (comparable to Fig. repro), using the identical
closest-approach/reopening criterion throughout. LumA has n=430
patients (vs. TNBC's 119), so this also checks whether the repulsion
pattern is a small-n-cohort artifact or holds at a much larger, more
Marchenko-Pastur-typical sample size (LumA's q=3.49 vs. TNBC's q=12.61).
"""
import numpy as np
import pandas as pd
import os

TOP_GENES = 60
K_EIGS = 8
N_START = 15
N_ORDERS = 20


def gap_pattern(traj, pair_idx, window=5):
    gaps = np.abs(traj[:, pair_idx] - traj[:, pair_idx + 1])
    min_idx = int(np.argmin(gaps))
    min_gap = gaps[min_idx]
    lo = max(0, min_idx - window)
    hi = min(len(gaps) - 1, min_idx + window)
    reopens_after = (hi > min_idx) and (gaps[hi] > 1.5 * min_gap)
    reopens_before = (lo < min_idx) and (gaps[lo] > 1.5 * min_gap)
    at_boundary = (min_idx == len(gaps) - 1)
    if at_boundary:
        return "inconclusive", min_idx, min_gap
    elif reopens_after and reopens_before:
        return "repulsion", min_idx, min_gap
    elif reopens_after or reopens_before:
        return "partial", min_idx, min_gap
    else:
        return "none", min_idx, min_gap


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    expr = pd.read_csv("../data/TCGA-BRCA/luma_expr.csv", index_col=0)
    var = expr.var(axis=1)
    sub_genes = var.sort_values(ascending=False).index[:TOP_GENES]
    X = expr.loc[sub_genes].to_numpy()
    p, n_total = X.shape
    print(f"LumA, top {p} most-variable genes, n_total={n_total} patients "
          f"(vs. TNBC's n=119 in the original Dyson-trajectory result).")

    # --- Single fixed order (seed 1, same convention as script 05) ---
    rng = np.random.default_rng(1)
    order = rng.permutation(n_total)
    Xp = X[:, order]
    ns = list(range(N_START, n_total + 1))
    traj = np.full((len(ns), K_EIGS), np.nan)
    for i, n in enumerate(ns):
        C = np.corrcoef(Xp[:, :n])
        traj[i] = np.sort(np.linalg.eigvalsh(C))[::-1][:K_EIGS]
    pd.DataFrame(traj, columns=[f"eig{k+1}" for k in range(K_EIGS)]).to_csv(
        "../results/luma_dyson_trajectories.csv", index=False)

    print(f"\nSingle-order scan (seed 1), {len(ns)} steps (n={N_START}..{n_total}):")
    for pair_idx in range(K_EIGS - 1):
        pattern, min_idx, min_gap = gap_pattern(traj, pair_idx)
        print(f"  lambda{pair_idx+1}-lambda{pair_idx+2}: {pattern} "
              f"(closest approach n={ns[min_idx]}, gap={min_gap:.4f})")

    # --- 20-order reproducibility check (same as script 17) ---
    print(f"\nReproducibility check: {N_ORDERS} independent random patient orders...")
    rows = []
    for seed in range(N_ORDERS):
        rng = np.random.default_rng(seed)
        perm = rng.permutation(n_total)
        Xr = X[:, perm]
        traj_r = np.full((len(ns), K_EIGS), np.nan)
        for i, n in enumerate(ns):
            C = np.corrcoef(Xr[:, :n])
            traj_r[i] = np.sort(np.linalg.eigvalsh(C))[::-1][:K_EIGS]
        for pair_idx in range(K_EIGS - 1):
            pattern, _, _ = gap_pattern(traj_r, pair_idx)
            rows.append(dict(order_seed=seed, pair=f"lambda{pair_idx+1}-lambda{pair_idx+2}",
                              pattern=pattern))
    df = pd.DataFrame(rows)
    df["pattern"] = df["pattern"].replace({
        "repulsion": "repulsion (both sides reopen)",
        "partial": "partial repulsion (one side reopens)",
        "inconclusive": "inconclusive (still narrowing at n_max)",
        "none": "no clear repulsion",
    })
    df.to_csv("../results/luma_dyson_repulsion_multi_order.csv", index=False)

    print("\nReproducibility across 20 random orders, by pair (LumA):")
    summary_rows = []
    for pair in df["pair"].unique():
        sub = df[df["pair"] == pair]
        counts = sub["pattern"].value_counts()
        n_rep = counts.get("repulsion (both sides reopen)", 0)
        n_part = counts.get("partial repulsion (one side reopens)", 0)
        n_inc = counts.get("inconclusive (still narrowing at n_max)", 0)
        n_none = counts.get("no clear repulsion", 0)
        print(f"  {pair}: repulsion={n_rep}/20, partial={n_part}/20, "
              f"inconclusive={n_inc}/20, none={n_none}/20")
        summary_rows.append(dict(pair=pair, repulsion=n_rep, partial=n_part,
                                  inconclusive=n_inc, none=n_none))
    pd.DataFrame(summary_rows).to_csv(
        "../results/luma_dyson_repulsion_multi_order_summary.csv", index=False)
    print("\nSaved luma_dyson_trajectories.csv, "
          "luma_dyson_repulsion_multi_order.csv and _summary.csv")
