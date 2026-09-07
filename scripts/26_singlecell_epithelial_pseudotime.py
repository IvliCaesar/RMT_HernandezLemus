"""
Closes the corrected Line 5/Line 4(ii) retest the Roadmap and Discussion
both point to explicitly: the global single-cell PC1 (Sec. wavelets,
script 23) carries only 9.7% of variance, no better than bulk PC1's
10.2%, because the full 9-cell-type GSE176078 sample is a categorical
mixture (T-cells vs. everything else), not a continuous trajectory.
The corrected next step, stated but not yet attempted, is to restrict
pseudo-time inference to a single, continuously-varying lineage -- the
894 Cancer Epithelial cells already present in this same sample -- and
retest both open items on that restricted population:

  (a) Line 4(ii): does the Dyson-repulsion eigenvalue-trajectory
      diagnostic (script 05, applied there to TNBC patients) show the
      same approach-then-reopen pattern when "time" is real single-cell
      pseudo-time within one lineage, instead of patient sample size?

  (b) Line 5 retest: does the wavelet coarse/fine decomposition (script
      04, whose real-data application failed a random-ordering control
      on bulk PC1) recover different coarse vs. fine hub genes along
      this restricted, higher-variance ordering axis, and does a
      random-order control still destroy that result the way it did on
      bulk data (script 04b)?

No new methodology anywhere in this script: the same PC1-ordering,
same Dyson-trajectory scan (script 05/17), and same wavelet pipeline
(script 04/04b) already used elsewhere in this paper, applied to a
population subset chosen specifically to fix the ordering-axis
weakness this Roadmap item already diagnosed.
"""
import sys
import os
import numpy as np
import pandas as pd
import pywt
from scipy.io import mmread

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score

DATA_DIR = "../data/GSE176078_singlecell"
N_GENES = 60          # match the TOP_GENES subset used for the bulk Dyson
                      # trajectory scan (script 05), for direct comparability
K_EIGS = 8
N_START = 30
N_RANDOM_ORDERS = 20
WAVELET_THRESHOLDS = np.round(np.arange(0.10, 0.80, 0.05), 2)


def load_epithelial_subset():
    mat = mmread(f"{DATA_DIR}/count_matrix_sparse.mtx").tocsr()  # genes x cells
    genes = pd.read_csv(f"{DATA_DIR}/count_matrix_genes.tsv", header=None)[0].to_numpy()
    barcodes = pd.read_csv(f"{DATA_DIR}/count_matrix_barcodes.tsv", header=None)[0].to_numpy()
    meta = pd.read_csv(f"{DATA_DIR}/metadata.csv", index_col=0).reindex(barcodes)
    is_epi = (meta["celltype_major"] == "Cancer Epithelial").to_numpy()
    print(f"Cancer Epithelial cells: {is_epi.sum()} / {len(barcodes)} "
          f"(real annotation from Wu et al. 2021, not used below).")
    mat_epi = mat[:, is_epi]
    return mat_epi, genes


def normalize_cp10k_log1p(mat):
    cell_totals = np.asarray(mat.sum(axis=0)).ravel()
    cell_totals[cell_totals == 0] = 1.0
    scale = 1e4 / cell_totals
    mat_norm = mat.multiply(scale[np.newaxis, :]).tocsr()
    mat_norm.data = np.log1p(mat_norm.data)
    return mat_norm


def gap_pattern(traj, pair_idx, window=5):
    """Same criterion as scripts 05/17: closest approach, then check
    whether the gap reopens >=1.5x within `window` steps on each side."""
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


def find_tau_star_and_hubs(C, genes, thresholds):
    results = []
    for tau in thresholds:
        A = C.copy()
        A[np.abs(A) < tau] = 0.0
        np.fill_diagonal(A, 0.0)
        eigs = np.linalg.eigvalsh(A)
        if len(np.unique(np.round(eigs, 8))) < 10:
            continue
        spacings = unfold_spectrum(eigs)
        ks_p, ks_g = poisson_goe_fit_score(spacings)
        regime = "GOE" if ks_g < ks_p else "Poisson"
        results.append((tau, regime))
    tau_star = None
    for i in range(len(results) - 1, 0, -1):
        if results[i][1] == "Poisson" and results[i - 1][1] == "GOE":
            tau_star = results[i - 1][0]
            break
    if tau_star is None:
        return None, []
    A = np.abs(C.copy())
    A[A < tau_star] = 0.0
    np.fill_diagonal(A, 0.0)
    degree = (A > 0).sum(axis=1)
    hubs = list(genes[np.argsort(degree)[::-1][:10]])
    return tau_star, hubs


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    mat_epi, genes = load_epithelial_subset()
    norm = normalize_cp10k_log1p(mat_epi)
    n_cells = norm.shape[1]

    detected_frac = np.asarray((norm > 0).sum(axis=1)).ravel() / n_cells
    keep = detected_frac >= 0.03
    dense = np.asarray(norm[keep].todense())
    genes_kept = genes[keep]
    print(f"Genes detected in >=3% of epithelial cells: {keep.sum()} / {len(genes)}")

    # --- Pseudo-time: PC1 restricted to this single lineage ---
    var_full = dense.var(axis=1)
    pt_gene_idx = np.argsort(var_full)[::-1][:1500]  # same p as script 23
    Xpt = dense[pt_gene_idx]
    Xc = Xpt - Xpt.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    pc1 = Vt[0]
    pc1_var_frac = 100 * S[0] ** 2 / np.sum(S ** 2)
    print(f"\nCancer-Epithelial-only PC1 explains {pc1_var_frac:.1f}% of "
          f"variance (n={n_cells} cells, p=1500 genes) -- compare against "
          f"global single-cell PC1 (9.7%, script 23/wavelets discussion) "
          f"and bulk-patient PC1 (10.2%, script 04).")
    order = np.argsort(pc1)

    # --- Dyson eigenvalue trajectory along real single-lineage pseudo-time ---
    dyson_gene_idx = np.argsort(var_full)[::-1][:N_GENES]
    Xd = dense[dyson_gene_idx][:, order]
    ns = list(range(N_START, n_cells + 1))
    traj = np.full((len(ns), K_EIGS), np.nan)
    for i, n in enumerate(ns):
        C = np.corrcoef(Xd[:, :n])
        traj[i] = np.sort(np.linalg.eigvalsh(C))[::-1][:K_EIGS]
    pd.DataFrame(traj, columns=[f"eig{k+1}" for k in range(K_EIGS)]).to_csv(
        "../results/singlecell_epithelial_dyson_trajectories.csv", index=False)

    print(f"\nDyson trajectory scan: top-{K_EIGS} eigenvalues of a "
          f"{N_GENES}-gene correlation matrix, {n_cells - N_START + 1} "
          f"steps (n={N_START}..{n_cells} cells) along real pseudo-time.")
    real_patterns = {}
    for pair_idx in range(K_EIGS - 1):
        pattern, min_idx, min_gap = gap_pattern(traj, pair_idx)
        real_patterns[pair_idx] = pattern
        print(f"  lambda{pair_idx+1}-lambda{pair_idx+2}: {pattern} "
              f"(closest approach at step n={ns[min_idx]}, gap={min_gap:.4f})")

    # Robustness check: is the real pseudo-time order doing anything a
    # random cell order wouldn't? Same reproducibility logic as script 17,
    # but here the comparison is real-order vs. random-order (not
    # real-order vs. real-order under different seeds), since pseudo-time
    # is the whole point being tested, not an arbitrary nuisance choice.
    print(f"\nRandom-order control: rerunning the identical scan under "
          f"{N_RANDOM_ORDERS} random (non-pseudo-time) cell orderings...")
    random_rows = []
    for seed in range(N_RANDOM_ORDERS):
        rng = np.random.default_rng(seed)
        perm = rng.permutation(n_cells)
        Xr = dense[dyson_gene_idx][:, perm]
        traj_r = np.full((len(ns), K_EIGS), np.nan)
        for i, n in enumerate(ns):
            C = np.corrcoef(Xr[:, :n])
            traj_r[i] = np.sort(np.linalg.eigvalsh(C))[::-1][:K_EIGS]
        for pair_idx in range(K_EIGS - 1):
            pattern, min_idx, min_gap = gap_pattern(traj_r, pair_idx)
            random_rows.append(dict(seed=seed, pair=f"lambda{pair_idx+1}-lambda{pair_idx+2}",
                                     pattern=pattern))
    df_rand = pd.DataFrame(random_rows)
    df_rand.to_csv("../results/singlecell_epithelial_dyson_random_control.csv", index=False)
    print("\nReal pseudo-time vs. random-order repulsion rate, by pair:")
    summary_rows = []
    for pair_idx in range(K_EIGS - 1):
        pair = f"lambda{pair_idx+1}-lambda{pair_idx+2}"
        counts = df_rand[df_rand["pair"] == pair]["pattern"].value_counts()
        n_rep = counts.get("repulsion", 0)
        n_part = counts.get("partial", 0)
        print(f"  {pair}: real={real_patterns[pair_idx]:>11s}  "
              f"random-order repulsion={n_rep}/{N_RANDOM_ORDERS}, "
              f"partial={n_part}/{N_RANDOM_ORDERS}")
        summary_rows.append(dict(pair=pair, real_pattern=real_patterns[pair_idx],
                                  random_repulsion=n_rep, random_partial=n_part))
    pd.DataFrame(summary_rows).to_csv(
        "../results/singlecell_epithelial_dyson_summary.csv", index=False)

    # --- Line 5 retest: wavelet coarse/fine decomposition on real pseudo-time ---
    print("\n--- Line 5 retest: wavelet decomposition along epithelial pseudo-time ---")
    wave_gene_idx = np.argsort(var_full)[::-1][:1500]
    Xw = dense[wave_gene_idx][:, order]
    p_w = Xw.shape[0]

    def wavelet_coarse_fine(Xo, n):
        coarse = np.zeros_like(Xo)
        fine = np.zeros_like(Xo)
        wavelet = "db4"
        level = pywt.dwt_max_level(n, pywt.Wavelet(wavelet).dec_len)
        level = max(1, min(level, 3))
        for i in range(Xo.shape[0]):
            coeffs = pywt.wavedec(Xo[i], wavelet, level=level)
            coarse_coeffs = [coeffs[0]] + [np.zeros_like(c) for c in coeffs[1:]]
            coarse[i] = pywt.waverec(coarse_coeffs, wavelet)[:n]
            fine_coeffs = [np.zeros_like(coeffs[0])] + \
                [np.zeros_like(c) for c in coeffs[1:-1]] + [coeffs[-1]]
            fine[i] = pywt.waverec(fine_coeffs, wavelet)[:n]
        return coarse, fine

    coarse, fine = wavelet_coarse_fine(Xw, n_cells)
    C_coarse = np.corrcoef(coarse)
    C_fine = np.corrcoef(fine)
    genes_w = genes_kept[wave_gene_idx]

    tau_coarse, hubs_coarse = find_tau_star_and_hubs(C_coarse, genes_w, WAVELET_THRESHOLDS)
    tau_fine, hubs_fine = find_tau_star_and_hubs(C_fine, genes_w, WAVELET_THRESHOLDS)
    print(f"Coarse-scale network: tau*={tau_coarse}, hubs={hubs_coarse}")
    print(f"Fine-scale network:   tau*={tau_fine}, hubs={hubs_fine}")
    real_overlap = len(set(hubs_coarse or []) & set(hubs_fine or []))
    print(f"Real hub-gene overlap (coarse vs fine, top-10 each): {real_overlap}/10")

    print(f"\nRandom-order control on the wavelet result "
          f"({N_RANDOM_ORDERS} replicates, same criterion as script 04b: "
          f"is the real overlap distinguishable from a meaningless ordering?)...")
    rand_overlaps = []
    for seed in range(N_RANDOM_ORDERS):
        rng = np.random.default_rng(1000 + seed)
        perm = rng.permutation(n_cells)
        Xr = Xw[:, perm]
        coarse_r, fine_r = wavelet_coarse_fine(Xr, n_cells)
        Cc_r = np.corrcoef(coarse_r)
        Cf_r = np.corrcoef(fine_r)
        _, hc_r = find_tau_star_and_hubs(Cc_r, genes_w, WAVELET_THRESHOLDS)
        _, hf_r = find_tau_star_and_hubs(Cf_r, genes_w, WAVELET_THRESHOLDS)
        rand_overlaps.append(len(set(hc_r or []) & set(hf_r or [])))
    rand_overlaps = np.array(rand_overlaps)
    print(f"Random-order overlaps (n={N_RANDOM_ORDERS}): "
          f"{rand_overlaps.tolist()}")
    n_ge_real = int((rand_overlaps >= real_overlap).sum())
    print(f"Real overlap ({real_overlap}) matched or exceeded by "
          f"{n_ge_real}/{N_RANDOM_ORDERS} random orderings.")

    pd.DataFrame({
        "seed": list(range(N_RANDOM_ORDERS)),
        "random_overlap": rand_overlaps,
    }).assign(real_overlap=real_overlap, pc1_var_frac=pc1_var_frac).to_csv(
        "../results/singlecell_epithelial_wavelet_control.csv", index=False)

    print("\nSaved: singlecell_epithelial_dyson_trajectories.csv, "
          "singlecell_epithelial_dyson_random_control.csv, "
          "singlecell_epithelial_dyson_summary.csv, "
          "singlecell_epithelial_wavelet_control.csv")
