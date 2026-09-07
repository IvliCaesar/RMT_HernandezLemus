"""
The corrected-but-still-linear pseudo-time tried in script 26 (PC1
restricted to the 894 Cancer Epithelial cells) explained only 8.6% of
variance, no better than the whole-sample figure -- and the Roadmap's
own next-step language after that result was explicit: a genuinely
different, non-PC1 pseudo-time estimator, not another restriction of
the same linear approach. This script builds one: diffusion pseudotime
(DPT, Haghverdi et al. 2016, Nature Methods), a standard single-cell
trajectory-inference method built from a cell-cell diffusion process
rather than a single global linear axis, implemented here directly in
numpy/scipy (kNN Gaussian-kernel affinity -> symmetric-normalized
diffusion operator -> eigendecomposition -> the closed-form DPT
distance from a root cell), since no external single-cell package
(scanpy, destiny) was available in this environment.

Once this genuinely non-linear ordering is in hand, it is used exactly
as PC1 was in script 26: as the "time" axis for the identical
Dyson-trajectory scan (with the same real-vs-random-order control) and
the identical wavelet coarse/fine retest, so the comparison to script
26's PC1-based results is apples-to-apples.
"""
import sys
import os
import numpy as np
import pandas as pd
import pywt
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score
from importlib import import_module
ep26 = import_module("26_singlecell_epithelial_pseudotime")

N_GENES_DPT = 1500     # genes used to build the diffusion graph itself
N_GENES_DYSON = 60     # genes used for the Dyson trajectory scan (as in script 26)
K_NEIGHBORS = 15
N_DIFFUSION_COMPONENTS = 15
K_EIGS = 8
N_START = 30
N_RANDOM_ORDERS = 20
WAVELET_THRESHOLDS = np.round(np.arange(0.10, 0.80, 0.05), 2)


def diffusion_pseudotime(X, k=K_NEIGHBORS, n_comp=N_DIFFUSION_COMPONENTS):
    """X: cells x genes (already normalized/log-transformed). Returns
    (pseudotime array, root cell index, eigenvalues used) via the
    standard symmetric-normalization trick: build the kNN Gaussian
    affinity W, form S = D^-1/2 W D^-1/2 (symmetric, real eigenvalues
    equal to the Markov transition operator's own), eigendecompose S
    directly with eigh, then apply Haghverdi et al. (2016)'s closed-form
    diffusion-pseudotime distance from a root cell using the transition
    operator's eigenvalues/eigenvectors (recovered from S's via the
    standard D^{-1/2} correction)."""
    n = X.shape[0]
    D2 = cdist(X, X, metric="sqeuclidean")
    # median-heuristic bandwidth over k-th nearest neighbor distances
    knn_dist = np.sort(D2, axis=1)[:, 1:k + 1]
    sigma2 = np.median(knn_dist)
    W = np.exp(-D2 / (2 * sigma2))
    # restrict to mutual kNN (sparsify, then symmetrize)
    idx_knn = np.argsort(D2, axis=1)[:, 1:k + 1]
    mask = np.zeros((n, n), dtype=bool)
    rows = np.repeat(np.arange(n), k)
    mask[rows, idx_knn.ravel()] = True
    mask = mask | mask.T
    W = W * mask
    np.fill_diagonal(W, 1.0)

    deg = W.sum(axis=1)
    d_inv_sqrt = 1.0 / np.sqrt(deg)
    S = W * d_inv_sqrt[:, None] * d_inv_sqrt[None, :]
    S = (S + S.T) / 2  # symmetrize away numerical asymmetry

    eigvals, eigvecs = np.linalg.eigh(S)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order][1:n_comp + 1]       # drop the trivial top eigenvalue (~1)
    eigvecs_S = eigvecs[:, order][:, 1:n_comp + 1]
    # recover the (non-symmetric) transition operator's right eigenvectors
    psi = eigvecs_S * d_inv_sqrt[:, None]
    psi = psi / np.linalg.norm(psi, axis=0, keepdims=True)

    # root cell: the one furthest, in diffusion-component space, from the
    # population centroid -- a standard, biology-free extremal-point
    # heuristic when no marker-based root is specified
    diff_coords = psi * eigvals[None, :]
    centroid = diff_coords.mean(axis=0)
    root = int(np.argmax(np.linalg.norm(diff_coords - centroid, axis=1)))

    # Haghverdi et al. (2016) closed-form DPT distance from the root
    weights = eigvals / (1.0 - eigvals)
    dpt = np.sqrt(np.sum((weights[None, :] * (psi - psi[root])) ** 2, axis=1))
    return dpt, root, eigvals


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    mat_epi, genes = ep26.load_epithelial_subset()
    norm = ep26.normalize_cp10k_log1p(mat_epi)
    n_cells = norm.shape[1]

    detected_frac = np.asarray((norm > 0).sum(axis=1)).ravel() / n_cells
    keep = detected_frac >= 0.03
    dense = np.asarray(norm[keep].todense())
    genes_kept = genes[keep]

    var_full = dense.var(axis=1)
    dpt_gene_idx = np.argsort(var_full)[::-1][:N_GENES_DPT]
    X = dense[dpt_gene_idx].T  # cells x genes

    print(f"Building diffusion map on {n_cells} Cancer Epithelial cells, "
          f"{N_GENES_DPT} genes, k={K_NEIGHBORS} nearest neighbors...")
    dpt, root, eigvals = diffusion_pseudotime(X)
    print(f"Root cell: index {root}. Top diffusion eigenvalues: "
          f"{np.round(eigvals[:5], 4)}")
    print(f"DPT range: [{dpt.min():.3f}, {dpt.max():.3f}], "
          f"mean={dpt.mean():.3f}, std={dpt.std():.3f}")

    # PC1 ordering from script 26's own already-computed result, for a
    # direct rank-correlation comparison (are these two axes actually
    # different, or coincidentally the same ordering under another name?)
    pc1_gene_idx = np.argsort(var_full)[::-1][:1500]
    Xpt = dense[pc1_gene_idx]
    Xc = Xpt - Xpt.mean(axis=1, keepdims=True)
    _, S_svd, Vt = np.linalg.svd(Xc, full_matrices=False)
    pc1 = Vt[0]
    rho, p_rho = spearmanr(dpt, pc1)
    print(f"\nSpearman rank correlation, DPT vs.\ PC1 ordering: "
          f"rho={rho:.3f} (p={p_rho:.2e}) -- "
          f"{'a genuinely different ordering' if abs(rho) < 0.5 else 'a similar ordering to PC1'}.")

    order = np.argsort(dpt)
    pd.DataFrame({"cell_rank": np.arange(n_cells), "dpt": dpt[order]}).to_csv(
        "../results/singlecell_epithelial_dpt.csv", index=False)

    # --- Dyson trajectory scan along DPT ordering (same method as script 26) ---
    dyson_gene_idx = np.argsort(var_full)[::-1][:N_GENES_DYSON]
    Xd = dense[dyson_gene_idx][:, order]
    ns = list(range(N_START, n_cells + 1))
    traj = np.full((len(ns), K_EIGS), np.nan)
    for i, n in enumerate(ns):
        C = np.corrcoef(Xd[:, :n])
        traj[i] = np.sort(np.linalg.eigvalsh(C))[::-1][:K_EIGS]

    print(f"\nDyson trajectory scan along DPT ordering "
          f"({n_cells - N_START + 1} steps):")
    real_patterns = {}
    for pair_idx in range(K_EIGS - 1):
        pattern, min_idx, min_gap = ep26.gap_pattern(traj, pair_idx)
        real_patterns[pair_idx] = pattern
        print(f"  lambda{pair_idx+1}-lambda{pair_idx+2}: {pattern} "
              f"(closest approach at step n={ns[min_idx]}, gap={min_gap:.4f})")

    print(f"\nRandom-order control ({N_RANDOM_ORDERS} replicates)...")
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
            pattern, _, _ = ep26.gap_pattern(traj_r, pair_idx)
            random_rows.append(dict(seed=seed, pair=f"lambda{pair_idx+1}-lambda{pair_idx+2}",
                                     pattern=pattern))
    df_rand = pd.DataFrame(random_rows)
    summary_rows = []
    print("\nDPT ordering vs. random-order repulsion rate, by pair:")
    for pair_idx in range(K_EIGS - 1):
        pair = f"lambda{pair_idx+1}-lambda{pair_idx+2}"
        counts = df_rand[df_rand["pair"] == pair]["pattern"].value_counts()
        n_rep = counts.get("repulsion", 0)
        print(f"  {pair}: DPT={real_patterns[pair_idx]:>11s}  "
              f"random-order repulsion={n_rep}/{N_RANDOM_ORDERS}")
        summary_rows.append(dict(pair=pair, dpt_pattern=real_patterns[pair_idx],
                                  random_repulsion=n_rep))
    pd.DataFrame(summary_rows).to_csv(
        "../results/singlecell_epithelial_dpt_dyson_summary.csv", index=False)

    # --- Wavelet retest along DPT ordering (same method as script 26) ---
    print("\n--- Wavelet retest along DPT ordering ---")
    wave_gene_idx = np.argsort(var_full)[::-1][:1500]
    Xw = dense[wave_gene_idx][:, order]
    genes_w = genes_kept[wave_gene_idx]

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

    tau_coarse, hubs_coarse = ep26.find_tau_star_and_hubs(C_coarse, genes_w, WAVELET_THRESHOLDS)
    tau_fine, hubs_fine = ep26.find_tau_star_and_hubs(C_fine, genes_w, WAVELET_THRESHOLDS)
    print(f"Coarse-scale network: tau*={tau_coarse}, hubs={hubs_coarse}")
    print(f"Fine-scale network:   tau*={tau_fine}, hubs={hubs_fine}")
    real_overlap = len(set(hubs_coarse or []) & set(hubs_fine or []))
    print(f"Real hub-gene overlap (coarse vs fine, top-10 each): {real_overlap}/10")

    rand_overlaps = []
    for seed in range(N_RANDOM_ORDERS):
        rng = np.random.default_rng(2000 + seed)
        perm = rng.permutation(n_cells)
        Xr = Xw[:, perm]
        coarse_r, fine_r = wavelet_coarse_fine(Xr, n_cells)
        Cc_r = np.corrcoef(coarse_r)
        Cf_r = np.corrcoef(fine_r)
        _, hc_r = ep26.find_tau_star_and_hubs(Cc_r, genes_w, WAVELET_THRESHOLDS)
        _, hf_r = ep26.find_tau_star_and_hubs(Cf_r, genes_w, WAVELET_THRESHOLDS)
        rand_overlaps.append(len(set(hc_r or []) & set(hf_r or [])))
    rand_overlaps = np.array(rand_overlaps)
    n_ge_real = int((rand_overlaps >= real_overlap).sum())
    print(f"Random-order overlaps: {rand_overlaps.tolist()}")
    print(f"Real overlap ({real_overlap}) matched or exceeded by "
          f"{n_ge_real}/{N_RANDOM_ORDERS} random orderings.")

    pd.DataFrame({"seed": list(range(N_RANDOM_ORDERS)), "random_overlap": rand_overlaps}).assign(
        real_overlap=real_overlap, dpt_pc1_spearman=rho).to_csv(
        "../results/singlecell_epithelial_dpt_wavelet_control.csv", index=False)

    print("\nSaved singlecell_epithelial_dpt.csv, "
          "singlecell_epithelial_dpt_dyson_summary.csv, "
          "singlecell_epithelial_dpt_wavelet_control.csv")
