"""
Genetic-algorithm feature selection for TNBC-vs-LumA classification,
added as a fourth feature-selection arm alongside RANDOM/TOP-VAR/
RMT-HUB in 06_ml_recognition.py: does an evolutionary search over
20-gene subsets, directly optimizing cross-validated accuracy, do any
better than the three simpler selectors on the same task? Given that
task's already-documented ceiling effect (Table "tab:ml": all three
simpler selectors hit ~98% with 20 genes, even a RANDOM set), this is
a real test of whether a more expensive search finds real headroom or
just confirms the ceiling from yet another angle -- reported honestly
either way.

DEAP genetic algorithm: chromosome = a fixed-length list of 20 gene
indices (repeated indices collapse via set() before evaluation, so
effective panel size can be <=20); fitness = mean 5-fold stratified
CV accuracy of a standardized logistic regression on that gene subset;
tournament selection, two-point crossover, index-mutation.
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from deap import base, creator, tools, algorithms
import random

N_GENES_PANEL = 20
POP_SIZE = 40
N_GEN = 15
RNG_SEED = 0

if __name__ == "__main__":
    tnbc = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
    luma = pd.read_csv("../data/TCGA-BRCA/luma_expr.csv", index_col=0)
    combined = pd.concat([tnbc, luma], axis=1)
    y = np.array([1] * tnbc.shape[1] + [0] * luma.shape[1])
    X_all = combined.to_numpy()
    n_candidate_genes = X_all.shape[0]
    genes = combined.index.to_numpy()

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)

    def fitness(individual):
        idx = sorted(set(individual))
        Xf = X_all[idx].T
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
        acc = cross_val_score(clf, Xf, y, cv=cv, scoring="accuracy").mean()
        return (acc,)

    random.seed(RNG_SEED)
    np.random.seed(RNG_SEED)

    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()
    toolbox.register("gene_idx", random.randint, 0, n_candidate_genes - 1)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
                      toolbox.gene_idx, n=N_GENES_PANEL)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", fitness)
    toolbox.register("mate", tools.cxTwoPoint)

    def mutate_gene(individual, indpb=0.1):
        for i in range(len(individual)):
            if random.random() < indpb:
                individual[i] = random.randint(0, n_candidate_genes - 1)
        return (individual,)

    toolbox.register("mutate", mutate_gene, indpb=0.1)
    toolbox.register("select", tools.selTournament, tournsize=3)

    pop = toolbox.population(n=POP_SIZE)
    print(f"Running GA: pop={POP_SIZE}, generations={N_GEN}, "
          f"panel size<={N_GENES_PANEL}, candidate genes={n_candidate_genes}")

    best_per_gen = []
    for gen in range(N_GEN):
        fits = list(map(toolbox.evaluate, pop))
        for ind, fit in zip(pop, fits):
            ind.fitness.values = fit
        best = max(pop, key=lambda i: i.fitness.values[0])
        best_per_gen.append(best.fitness.values[0])
        print(f"  gen {gen}: best CV accuracy = {best.fitness.values[0]:.4f}")
        offspring = toolbox.select(pop, len(pop))
        offspring = list(map(toolbox.clone, offspring))
        offspring = algorithms.varAnd(offspring, toolbox, cxpb=0.5, mutpb=0.3)
        pop[:] = offspring

    fits = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fits):
        ind.fitness.values = fit
    best = max(pop, key=lambda i: i.fitness.values[0])
    best_idx = sorted(set(best))
    best_genes = list(genes[best_idx])

    print(f"\nFinal best: {len(best_idx)} unique genes, "
          f"CV accuracy = {best.fitness.values[0]:.4f}")
    print(f"Genes: {best_genes}")

    # ROC-AUC of the final best panel, for direct comparison with Table tab:ml
    Xf = X_all[best_idx].T
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
    auc = cross_val_score(clf, Xf, y, cv=cv, scoring="roc_auc").mean()
    print(f"ROC-AUC of final GA panel: {auc:.4f}")

    pd.DataFrame({"generation": list(range(N_GEN)),
                  "best_cv_accuracy": best_per_gen}).to_csv(
        "../results/ga_feature_selection_history.csv", index=False)
    with open("../results/ga_feature_selection_summary.txt", "w") as f:
        f.write(f"Final panel ({len(best_idx)} genes): {best_genes}\n")
        f.write(f"CV accuracy: {best.fitness.values[0]:.4f}\n")
        f.write(f"ROC-AUC: {auc:.4f}\n")
    print("Saved results/ga_feature_selection_history.csv and _summary.txt")
