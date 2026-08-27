#!/usr/bin/env python3
# Part (c): Fit the 2-2-3-1 network to g(x1,x2)=x1*x2 with Levenberg–Marquardt,
# plot f(θ)=‖r(θ)‖²+γ‖θ‖² and ‖∇f(θ)‖ versus iteration, and report RMS error.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import os

# ----- model pieces -----
def sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))

def sigmoid_prime(z):
    s = sigmoid(z)
    return s * (1.0 - s)

def fhat(theta, X):
    # theta = [θ1,...,θ13]; X shape (N,2) with columns x1,x2
    x1, x2 = X[:, 0], X[:, 1]
    s1 = theta[1]*x1 + theta[2]*x2 + theta[3]     # θ2,θ3,θ4
    s2 = theta[5]*x1 + theta[6]*x2 + theta[7]     # θ6,θ7,θ8
    s3 = theta[9]*x1 + theta[10]*x2 + theta[11]   # θ10,θ11,θ12
    a1, a2, a3 = sigmoid(s1), sigmoid(s2), sigmoid(s3)
    yhat = theta[0]*a1 + theta[4]*a2 + theta[8]*a3 + theta[12]  # θ1,θ5,θ9,θ13
    return yhat

def grad_fhat(theta, X):
    # returns Jacobian of fhat wrt theta: shape (N, 13)
    x1, x2 = X[:, 0], X[:, 1]
    s1 = theta[1]*x1 + theta[2]*x2 + theta[3]
    s2 = theta[5]*x1 + theta[6]*x2 + theta[7]
    s3 = theta[9]*x1 + theta[10]*x2 + theta[11]
    a1, a2, a3 = sigmoid(s1), sigmoid(s2), sigmoid(s3)
    g1, g2, g3 = sigmoid_prime(s1), sigmoid_prime(s2), sigmoid_prime(s3)

    N = X.shape[0]
    J = np.zeros((N, 13))
    J[:, 0]  = a1                    # d/dθ1
    J[:, 1]  = theta[0]*g1*x1        # d/dθ2
    J[:, 2]  = theta[0]*g1*x2        # d/dθ3
    J[:, 3]  = theta[0]*g1           # d/dθ4
    J[:, 4]  = a2                    # d/dθ5
    J[:, 5]  = theta[4]*g2*x1        # d/dθ6
    J[:, 6]  = theta[4]*g2*x2        # d/dθ7
    J[:, 7]  = theta[4]*g2           # d/dθ8
    J[:, 8]  = a3                    # d/dθ9
    J[:, 9]  = theta[8]*g3*x1        # d/dθ10
    J[:,10]  = theta[8]*g3*x2        # d/dθ11
    J[:,11]  = theta[8]*g3           # d/dθ12
    J[:,12]  = 1.0                   # d/dθ13
    return J

# ----- augmented residuals for LM with L2 regularization -----
def aug_residual(theta, X, y, gamma, trace=None):
    r = fhat(theta, X) - y                             # shape (N,)
    res = np.concatenate([r, np.sqrt(gamma)*theta])    # shape (N+13,)
    if trace is not None:
        J = np.vstack([grad_fhat(theta, X), np.sqrt(gamma)*np.eye(theta.size)])
        fval = np.dot(res, res)                        # = ‖r‖² + γ‖θ‖²
        grad = 2.0 * J.T @ res                         # ∇f
        trace["f"].append(fval)
        trace["grad_norm"].append(np.linalg.norm(grad))
        trace["theta"].append(theta.copy())
    return res

def aug_jacobian(theta, X, y, gamma):
    return np.vstack([grad_fhat(theta, X), np.sqrt(gamma)*np.eye(theta.size)])

# ----- one LM run -----
def run_lm(seed, gamma=1e-5, init_scale=0.5):
    rng = np.random.default_rng(seed)
    N = 200
    X = rng.uniform(-1.0, 1.0, size=(N, 2))
    y = X[:, 0] * X[:, 1]
    theta0 = rng.normal(scale=init_scale, size=13)
    trace = {"f": [], "grad_norm": [], "theta": []}
    fun = lambda th: aug_residual(th, X, y, gamma, trace)
    jac = lambda th: aug_jacobian(th, X, y, gamma)
    res = least_squares(fun, theta0, jac=jac, method="lm", max_nfev=500)
    theta_hat = res.x
    rms = np.sqrt(np.mean((fhat(theta_hat, X) - y) ** 2))
    return X, y, theta0, theta_hat, rms, trace

# ----- experiment with multiple starts -----
seeds = [0, 1, 2]
results = [run_lm(s) for s in seeds]

# Get script directory for saving outputs
script_dir = os.path.dirname(os.path.abspath(__file__))

# ----- plots -----
plt.figure(figsize=(7.2, 5.2))
for (X, y, th0, th, rms, tr), label in zip(results, [f"start {k}" for k in range(len(seeds))]):
    # Start from iteration 10 to avoid initial spike
    plt.plot(range(10, len(tr["f"])), tr["f"][10:], label=f"f(θ) — {label} (RMS={rms:.6g})")
plt.xlabel("iteration (residual evaluations)")
plt.ylabel(r"$f(\theta)=\|r(\theta)\|^2+\gamma\|\theta\|^2$")
plt.title("Objective value vs iteration (LM)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(script_dir, "nn_lm_f_vs_iter.png"), dpi=180)

plt.figure(figsize=(7.2, 5.2))
for (X, y, th0, th, rms, tr), label in zip(results, [f"start {k}" for k in range(len(seeds))]):
    plt.semilogy(tr["grad_norm"], label=f"‖∇f‖ — {label}")
plt.xlabel("iteration (residual evaluations)")
plt.ylabel(r"$\|\nabla f(\theta)\|_2$")
plt.title("Gradient norm vs iteration (LM)")
plt.legend()
plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(script_dir, "nn_lm_gradnorm_vs_iter.png"), dpi=180)

# RMS errors are displayed in the plot legend