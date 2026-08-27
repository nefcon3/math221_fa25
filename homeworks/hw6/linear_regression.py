#!/usr/bin/env python3
# Linear regression fit to g(x1,x2)=x1*x2
# Using direct solve (OLS/Ridge) and also LM for iteration plots

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import os

#!/usr/bin/env python3
# Part (d): Fit the same dataset with a linear model f_lin(x; β, ν) = x^T β + ν
# and report the RMS fitting error. (No execution here.)

def make_data(N=200, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.uniform(-1.0, 1.0, size=(N, 2))
    y = X[:, 0] * X[:, 1]
    return X, y

def fit_linear_ols(X, y):
    Z = np.column_stack([X, np.ones(len(X))])  # design: [x1, x2, 1]
    w, *_ = np.linalg.lstsq(Z, y, rcond=None)  # w = [β1, β2, ν]
    yhat = Z @ w
    rms = np.sqrt(np.mean((yhat - y) ** 2))
    beta, nu = w[:2], w[2]
    return beta, nu, yhat, rms

def fit_linear_ridge(X, y, lam=1e-5):
    Z = np.column_stack([X, np.ones(len(X))])
    I = np.eye(Z.shape[1])
    w = np.linalg.solve(Z.T @ Z + lam * I, Z.T @ y)
    yhat = Z @ w
    rms = np.sqrt(np.mean((yhat - y) ** 2))
    beta, nu = w[:2], w[2]
    return beta, nu, yhat, rms

if __name__ == "__main__":
    X, y = make_data(N=200, seed=0)
    beta_ols, nu_ols, yhat_ols, rms_ols = fit_linear_ols(X, y)
    beta_ridge, nu_ridge, yhat_ridge, rms_ridge = fit_linear_ridge(X, y, lam=1e-5)

    print("OLS  : beta =", beta_ols, "nu =", nu_ols, "RMS =", rms_ols)
    print("Ridge: beta =", beta_ridge, "nu =", nu_ridge, "RMS =", rms_ridge)