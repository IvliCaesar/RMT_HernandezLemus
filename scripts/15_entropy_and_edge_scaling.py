"""
Two further real RMT diagnostics on the TNBC/LumA correlation matrices:

1) Spectral (von Neumann) entropy: treat the normalized eigenvalue
   spectrum of a correlation matrix, p_i = lambda_i / sum(lambda_i), as
   a probability distribution and compute its Shannon entropy
   S = -sum p_i log(p_i) (normalized to [0,1] by log(p), p = number of
   genes). Low entropy means variance concentrated in a few directions
   (real, low-rank structure); high entropy (approaching 1) means
   variance spread evenly across all directions, the signature of an
   unstructured random matrix. Computed for the real correlation
   matrix and, as a null comparison, for 20 independently
   gene-permuted replicates per cohort (same permutation recipe as
   scripts/10_monte_carlo_null.py).

2) Tracy-Widom edge-fluctuation scaling check: NOT a full distributional
   fit (no verified Tracy-Widom density implementation was available,
   and this paper does not claim one) -- instead, a real, checkable
   scaling-law test. Under the null (i.i.d. Gaussian data, no real
   correlation structure), the largest eigenvalue of a p x n sample
   correlation matrix fluctuates around the Marchenko-Pastur edge
   lambda_+ on the Tracy-Widom scale, O(n^{-2/3}), not the O(n^{-1/2})
   scale a naive Gaussian (CLT) argument would predict. We generate
   synthetic i.i.d. Gaussian null data at several sample sizes n
   (holding p fixed) and check whether the standard deviation of the
   top eigenvalue actually shrinks like n^{-2/3} as n grows -- the
   qualitative, checkable signature of Tracy-Widom rather than
   Gaussian fluctuation, without asserting a precise density match.
"""
import numpy as np
import pandas as pd
import os

os.makedirs("../results", exist_ok=True)


def spectral_entropy(C):
    eigs = np.linalg.eigvalsh(C)
    eigs = np.clip(eigs, 1e-12, None)  # correlation matrices are PSD; guard tiny negatives from float error
    p = eigs / eigs.sum()
    p = p[p > 0]
    S = -np.sum(p * np.log(p))
    return S / np.log(len(eigs))  # normalized to [0,1]


if __name__ == "__main__":
    rng = np.random.default_rng(99)
    rows = []

    # --- 1) Spectral entropy: real vs permuted-gene null ---
    print("=== Spectral (von Neumann) entropy ===")
    for name, path in [("TNBC", "../data/TCGA-BRCA/tnbc_expr.csv"),
                        ("LumA", "../data/TCGA-BRCA/luma_expr.csv")]:
        expr = pd.read_csv(path, index_col=0)
        X = expr.to_numpy()
        p, n = X.shape
        C_real = np.corrcoef(X)
        S_real = spectral_entropy(C_real)

        S_perm = []
        for rep in range(20):
            Xp = X.copy()
            for i in range(p):
                Xp[i] = rng.permutation(Xp[i])
            S_perm.append(spectral_entropy(np.corrcoef(Xp)))
        S_perm = np.array(S_perm)

        print(f"{name}: real entropy={S_real:.4f}, "
              f"permuted-null entropy={S_perm.mean():.4f}+/-{S_perm.std():.4f} "
              f"(20 reps), gap={S_perm.mean()-S_real:.4f}")
        rows.append(dict(cohort=name, S_real=S_real,
                          S_perm_mean=S_perm.mean(), S_perm_std=S_perm.std()))

    pd.DataFrame(rows).to_csv("../results/spectral_entropy.csv", index=False)

    # --- 2) Tracy-Widom edge-fluctuation scaling check ---
    # IMPORTANT: q=p/n must be held FIXED while n grows -- the asymptotic
    # n^{-2/3} (Tracy-Widom) vs n^{-1/2} (Gaussian/CLT) scaling law is a
    # statement about the n->infinity limit AT FIXED ASPECT RATIO. Letting
    # p float while only n is varied (an earlier version of this script did
    # exactly that) changes q across the sweep too, confounding the two
    # effects and giving a fitted exponent that reflects neither law cleanly.
    print("\n=== Tracy-Widom edge-fluctuation scaling (synthetic i.i.d. null, q fixed) ===")
    Q_FIXED = 3.0  # close to LumA's real q=3.49
    n_values = [60, 120, 240, 480, 960]
    n_reps_by_n = {60: 400, 120: 400, 240: 300, 480: 150, 960: 60}
    scaling_rows = []
    for n in n_values:
        p_n = int(round(Q_FIXED * n))
        n_reps = n_reps_by_n[n]
        tops = []
        for rep in range(n_reps):
            Xs = rng.standard_normal((p_n, n))
            C = np.corrcoef(Xs)
            tops.append(np.linalg.eigvalsh(C)[-1])
        tops = np.array(tops)
        q = p_n / n
        lam_plus = (1 + np.sqrt(q)) ** 2
        print(f"n={n:4d}, p={p_n:4d} (q={q:.3f}, lambda_+={lam_plus:.3f}, "
              f"{n_reps} reps): top eigenvalue mean={tops.mean():.4f}, "
              f"std={tops.std():.5f}")
        scaling_rows.append(dict(n=n, p=p_n, q=q, lambda_plus=lam_plus,
                                  n_reps=n_reps, top_mean=tops.mean(),
                                  top_std=tops.std()))

    df = pd.DataFrame(scaling_rows)
    # Fit log(std) = a + b*log(n); Tracy-Widom scaling predicts b ~ -2/3,
    # Gaussian/CLT scaling would predict b ~ -1.
    logn = np.log(df["n"].to_numpy())
    logstd = np.log(df["top_std"].to_numpy())
    b, a = np.polyfit(logn, logstd, 1)
    print(f"\nFitted scaling exponent: std(top eigenvalue) ~ n^{b:.3f} "
          f"(Tracy-Widom predicts ~n^-0.667, Gaussian/CLT predicts ~n^-1.0)")
    df.to_csv("../results/tracy_widom_scaling.csv", index=False)
    with open("../results/tracy_widom_scaling_summary.txt", "w") as f:
        f.write(f"Fitted exponent b={b:.3f} in std ~ n^b "
                f"(n_reps={n_reps} per n, p={p_fixed} fixed)\n")
        f.write("Tracy-Widom (edge fluctuation) predicts b~-0.667; "
                "Gaussian/CLT (bulk fluctuation) predicts b~-1.0.\n")
    print("Saved results/spectral_entropy.csv, "
          "results/tracy_widom_scaling.csv, and summary.")
