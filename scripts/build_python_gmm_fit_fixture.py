#!/usr/bin/env python3
"""Build a deterministic full-fit GMM fixture independently with NumPy.

The starting centroids are fixed deliberately: this isolates the EM trajectory
from the documented sklearn-versus-R k-means initialization difference.
"""

import json
from pathlib import Path

import numpy as np


X = np.asarray([
    [-2.3, -1.7], [-2.0, -2.2], [-1.6, -1.8], [-2.2, -2.4],
    [1.7, 2.4], [2.1, 1.8], [2.5, 2.2], [1.8, 1.6],
], dtype=np.float64)
CENTERS = np.asarray([[-2.3, -1.7], [2.5, 2.2]], dtype=np.float64)
REG = 1e-6
TOL = 1e-8
MAX_ITER = 30


def log_scores(x, means, covariances, weights):
    scores = []
    for mean, covariance, weight in zip(means, covariances, weights):
        adjusted = covariance + REG * np.eye(x.shape[1])
        sign, logdet = np.linalg.slogdet(adjusted)
        assert sign > 0
        centered = x - mean
        mahalanobis = np.einsum(
            "ni,ij,nj->n", centered, np.linalg.inv(adjusted), centered
        )
        scores.append(
            np.log(weight) - 0.5 *
            (x.shape[1] * np.log(2 * np.pi) + logdet + mahalanobis)
        )
    return np.column_stack(scores)


def responsibilities(scores):
    row_max = np.max(scores, axis=1, keepdims=True)
    unscaled = np.exp(scores - row_max)
    return unscaled / np.sum(unscaled, axis=1, keepdims=True)


# Match PyCave's k-means initializer: hard nearest-centroid assignments, then
# empirical component probabilities, means, and maximum-likelihood covariance.
assignment = np.argmin(
    np.sum((X[:, None, :] - CENTERS[None, :, :]) ** 2, axis=2), axis=1
)
means = np.vstack([X[assignment == component].mean(axis=0) for component in range(2)])
weights = np.bincount(assignment, minlength=2).astype(float) / X.shape[0]
covariances = np.stack([
    ((X[assignment == component] - means[component]).T @
     (X[assignment == component] - means[component])) /
    np.sum(assignment == component)
    for component in range(2)
])

best_nll = np.inf
wait = 0
nll_history = []
for iteration in range(1, MAX_ITER + 1):
    scores = log_scores(X, means, covariances, weights)
    row_max = np.max(scores, axis=1)
    nll = -np.mean(row_max + np.log(np.exp(scores - row_max[:, None]).sum(axis=1)))
    nll_history.append(float(nll))
    resp = responsibilities(scores)
    nk = resp.sum(axis=0)
    means = (resp.T @ X) / nk[:, None]
    covariances = np.stack([
        ((resp[:, component, None] * (X - means[component])).T @
         (X - means[component])) / nk[component]
        for component in range(2)
    ])
    weights = nk / X.shape[0]
    if nll < best_nll - TOL:
        best_nll = nll
        wait = 0
    else:
        wait += 1
        if wait >= 2:
            break

final_scores = log_scores(X, means, covariances, weights)
row_max = np.max(final_scores, axis=1)
final_nll = -np.mean(
    row_max + np.log(np.exp(final_scores - row_max[:, None]).sum(axis=1))
)
payload = {
    "x": X.tolist(),
    "initial_centers": CENTERS.tolist(),
    "weights": weights.tolist(),
    "means": means.tolist(),
    "covariances": covariances.tolist(),
    "nll": float(final_nll),
    "labels": np.argmax(final_scores, axis=1).tolist(),
    "iterations": iteration,
    "nll_history": nll_history,
    "parameters": {
        "max_iter": MAX_ITER, "tolerance": TOL,
        "covariance_regularization": REG,
    },
    "provenance": {
        "reference": "independent NumPy implementation of PyCave 3.2.1 full-batch EM",
        "controlled_initialization": True,
        "purpose": "fit-to-fit R/Python EM trajectory parity",
    },
}
output = Path(__file__).resolve().parents[1] / "inst/extdata/python_gmm_fit_fixture.json"
output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(output)
