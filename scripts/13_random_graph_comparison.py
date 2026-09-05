"""
Graph-theoretic characterization of the RMT-thresholded networks
(Line 2 of the program: "large-scale networks"), comparing each real
network against two standard random-graph nulls of the same size:
  - Erdos-Renyi G(n,p): same node/edge count, no other structure.
  - Configuration model: same node count AND same degree sequence
    (preserves the real hub-degree distribution, tests whether
    clustering/modularity beyond degree alone is real).
Reports clustering coefficient, largest connected component size, and
Louvain-community count/modularity for the real network and for the
mean over random-graph replicates of each null -- the standard
small-world/modular-network check (Luo et al. 2007 and the biological
network literature generally use exactly this kind of comparison).
Louvain (not greedy-modularity) is used for tractability on LumA's
~146K-edge network; both are standard modularity-maximization
heuristics.
"""
import numpy as np
import pandas as pd
import networkx as nx
import os

N_RANDOM_REPS = 5


def network_stats(G):
    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        return dict(clustering=0.0, giant_frac=0.0, n_communities=0,
                    modularity=0.0)
    clustering = nx.average_clustering(G)
    giant = max(nx.connected_components(G), key=len)
    giant_frac = len(giant) / G.number_of_nodes()
    communities = nx.community.louvain_communities(G, seed=0)
    modularity = nx.community.modularity(G, communities)
    return dict(clustering=clustering, giant_frac=giant_frac,
                n_communities=len(communities), modularity=modularity)


def build_real_network(name, path, tau_star):
    expr = pd.read_csv(path, index_col=0)
    genes = expr.index.to_numpy()
    C = np.corrcoef(expr.to_numpy())
    A = np.abs(C.copy())
    A[A < tau_star] = 0.0
    np.fill_diagonal(A, 0.0)
    G = nx.from_numpy_array(A > 0)
    G.remove_nodes_from(list(nx.isolates(G)))
    return G


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    rows = []

    for name, path, tau_star in [
        ("TNBC", "../data/TCGA-BRCA/tnbc_expr.csv", 0.45),
        ("LumA", "../data/TCGA-BRCA/luma_expr.csv", 0.25),
    ]:
        G = build_real_network(name, path, tau_star)
        n, m = G.number_of_nodes(), G.number_of_edges()
        degree_seq = [d for _, d in G.degree()]
        print(f"=== {name}: real network n={n} nodes, m={m} edges, "
              f"mean degree={2*m/n:.1f} ===")

        real_stats = network_stats(G)
        print(f"  REAL:  clustering={real_stats['clustering']:.4f}, "
              f"giant component={100*real_stats['giant_frac']:.1f}%, "
              f"{real_stats['n_communities']} communities, "
              f"modularity={real_stats['modularity']:.4f}")
        rows.append(dict(cohort=name, model="real", **real_stats))

        p = 2 * m / (n * (n - 1))
        er_stats_list = []
        for rep in range(N_RANDOM_REPS):
            G_er = nx.gnp_random_graph(n, p, seed=rep)
            er_stats_list.append(network_stats(G_er))
        er_mean = {k: np.mean([s[k] for s in er_stats_list])
                   for k in er_stats_list[0]}
        print(f"  ER({N_RANDOM_REPS} reps mean): clustering={er_mean['clustering']:.4f}, "
              f"giant component={100*er_mean['giant_frac']:.1f}%, "
              f"{er_mean['n_communities']:.1f} communities, "
              f"modularity={er_mean['modularity']:.4f}")
        rows.append(dict(cohort=name, model="erdos_renyi", **er_mean))

        cm_stats_list = []
        for rep in range(N_RANDOM_REPS):
            try:
                G_cm = nx.random_degree_sequence_graph(
                    degree_seq, seed=rep, tries=50)
                G_cm = nx.Graph(G_cm)
                G_cm.remove_edges_from(nx.selfloop_edges(G_cm))
            except nx.NetworkXError:
                continue
            cm_stats_list.append(network_stats(G_cm))
        if cm_stats_list:
            cm_mean = {k: np.mean([s[k] for s in cm_stats_list])
                       for k in cm_stats_list[0]}
            print(f"  ConfigModel({len(cm_stats_list)} reps mean): "
                  f"clustering={cm_mean['clustering']:.4f}, "
                  f"giant component={100*cm_mean['giant_frac']:.1f}%, "
                  f"{cm_mean['n_communities']:.1f} communities, "
                  f"modularity={cm_mean['modularity']:.4f}")
            rows.append(dict(cohort=name, model="config_model", **cm_mean))
        print()

    pd.DataFrame(rows).to_csv("../results/random_graph_comparison.csv",
                               index=False)
    print("Saved results/random_graph_comparison.csv")
