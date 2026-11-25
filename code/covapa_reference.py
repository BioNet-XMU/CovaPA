"""Reference implementation of the pathway statistic and cohort inference.

Pass a participant × metabolite array already on the analysis scale, a
two-group label vector, and a dict mapping each tested pathway to its
measured column indices. Returns pathway, squared distance, Monte Carlo p,
and BH q. Input order and RNG affect the finite permutation sample.
Exchangeability of participant profiles is required for permutation inference.
"""
import numpy as np


def _sqrt(a):
    values, vectors = np.linalg.eigh(a)
    return (vectors * np.sqrt(np.maximum(values, 0))) @ vectors.T


def squared_bures(a, b):
    root = _sqrt(a)
    return float(max(0, np.trace(a + b) - 2 * np.trace(_sqrt(root @ b @ root))))


def covapa(x, labels, pathways, *, permutations=999, seed=20261119, ridge=1e-4):
    x = np.asarray(x, float)
    labels = np.asarray(labels)
    group = np.unique(labels)
    if x.ndim != 2 or not np.isfinite(x).all() or len(group) != 2 or len(labels) != len(x):
        raise ValueError("Expected finite participant × metabolite input and exactly two groups.")
    if any(sum(labels == g) < 2 for g in group):
        raise ValueError("Each group needs at least two participants.")
    sets = [np.array(sorted(set(ids)), dtype=int) for ids in pathways.values()]
    if not sets or any(len(ids) < 3 or ids.min() < 0 or ids.max() >= x.shape[1] for ids in sets):
        raise ValueError("Provide valid column indices for at least three metabolites per tested pathway.")
    degree = np.bincount(np.concatenate(sets), minlength=x.shape[1])
    features = [x[:, ids] / degree[ids] for ids in sets]
    identities = [ridge * np.eye(z.shape[1]) for z in features]

    def scores(case):
        return np.array([squared_bures(np.cov(z[case], rowvar=False) + i,
                                       np.cov(z[~case], rowvar=False) + i)
                         for z, i in zip(features, identities, strict=True)])

    case = labels == group[0]
    observed = scores(case)
    rng = np.random.default_rng(seed)
    exceed = np.zeros_like(observed, dtype=int)
    for _ in range(permutations):
        exceed += scores(rng.permutation(case)) >= observed
    p = (1 + exceed) / (permutations + 1)
    order = np.argsort(p)
    q = np.empty_like(p)
    q[order] = np.minimum(1, np.minimum.accumulate(
        (len(p) * p[order] / np.arange(1, len(p) + 1))[::-1])[::-1])
    return list(zip(pathways, observed, p, q, strict=True))
