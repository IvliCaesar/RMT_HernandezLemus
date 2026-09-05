"""
Closes a gap explicitly flagged but not yet attempted in the article's
Discussion/Conclusions: script 05 tracked eigenvalue trajectories under
ONE fixed random patient-inclusion order, finding several pairs with a
clear approach-then-reopen (repulsion) pattern and one pair
(lambda_2-lambda_3) still narrowing at the end of the real cohort
(n=119), reported as inconclusive rather than forced into a
conclusion. This script reruns the identical trajectory scan under
N_ORDERS=20 independent random patient orders (same 60-gene TNBC
subset, same method) and asks, for each adjacent pair: in what
fraction of random orders does it show a clear local-minimum-then-
reopen pattern, versus stay inconclusive (still narrowing at n=119) or
show no clear minimum at all in this range?

This is a real, direct test of reproducibility, not a re-assertion:
whichever pairs repel reliably across orders are the ones the "Dyson
repulsion" claim can actually stand on; pairs that only repel under
the original seed would mean that finding was an artifact of one
arbitrary ordering choice.
"""
import numpy as np
import pandas as pd
import os

N_SAMPLES_MIN = 15
K_EIGS = 8
N_ORDERS = 20


def gap_pattern(traj, ns, pair_idx):
    """For eigenvalue pair (pair_idx, pair_idx+1), find the closest
    approach and check whether the gap reopens by >=1.5x within 5 steps
    on both sides where available -- same criterion as script 05."""
    gaps = np.abs(traj[:, pair_idx] - traj[:, pair_idx + 1])
    min_idx = int(np.argmin(gaps))
    min_gap = gaps[min_idx]
    lo = max(0, min_idx - 5)
    hi = min(len(gaps) - 1, min_idx + 5)
    reopens_after = (hi > min_idx) and (gaps[hi] > 1.5 * min_gap)
    reopens_before = (lo < min_idx) and (gaps[lo] > 1.5 * min_gap)
    at_boundary = (min_idx == len(gaps) - 1)  # still narrowing at n=119
    if at_boundary:
        return "inconclusive (still narrowing at n_max)"
    elif reopens_after and reopens_before:
        return "repulsion (both sides reopen)"
    elif reopens_after or reopens_before:
        return "partial repulsion (one side reopens)"
    else:
        return "no clear repulsion"


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    expr = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
    var = expr.var(axis=1)
    sub_genes = var.sort_values(ascending=False).index[:60]
    X = expr.loc[sub_genes].to_numpy()
    p, n_total = X.shape
    ns = list(range(N_SAMPLES_MIN, n_total + 1))
    print(f"TNBC, top {p} genes, n_total={n_total}, {N_ORDERS} random orders")

    rows = []
    for order_seed in range(N_ORDERS):
        rng = np.random.default_rng(order_seed)
        perm = rng.permutation(n_total)
        Xp = X[:, perm]
        traj = np.full((len(ns), K_EIGS), np.nan)
        for i, n in enumerate(ns):
            C = np.corrcoef(Xp[:, :n])
            traj[i] = np.sort(np.linalg.eigvalsh(C))[::-1][:K_EIGS]
        for pair_idx in range(K_EIGS - 1):
            pattern = gap_pattern(traj, ns, pair_idx)
            rows.append(dict(order_seed=order_seed,
                              pair=f"lambda{pair_idx+1}-lambda{pair_idx+2}",
                              pattern=pattern))

    df = pd.DataFrame(rows)
    df.to_csv("../results/dyson_repulsion_multi_order.csv", index=False)

    print("\nReproducibility across 20 random orders, by pair:")
    summary_rows = []
    for pair in df["pair"].unique():
        sub = df[df["pair"] == pair]
        counts = sub["pattern"].value_counts()
        n_repulsion = counts.get("repulsion (both sides reopen)", 0)
        n_partial = counts.get("partial repulsion (one side reopens)", 0)
        n_inconclusive = counts.get("inconclusive (still narrowing at n_max)", 0)
        n_none = counts.get("no clear repulsion", 0)
        print(f"  {pair}: repulsion={n_repulsion}/20, partial={n_partial}/20, "
              f"inconclusive={n_inconclusive}/20, none={n_none}/20")
        summary_rows.append(dict(pair=pair, repulsion=n_repulsion,
                                  partial=n_partial,
                                  inconclusive=n_inconclusive, none=n_none))
    pd.DataFrame(summary_rows).to_csv(
        "../results/dyson_repulsion_multi_order_summary.csv", index=False)
    print("\nSaved results/dyson_repulsion_multi_order.csv and _summary.csv")
