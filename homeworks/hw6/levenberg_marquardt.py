import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import os

x = np.arange(0, 6, dtype=float)
y = np.array([5.2, 4.5, 2.7, 2.5, 2.1, 1.9], float)

def f_model(x, theta):
    return theta[0] * np.exp(theta[1] * x)

def residuals(theta, x, y):
    return f_model(x, theta) - y

logy = np.log(y)
A = np.vstack([np.ones_like(x), x]).T
beta, *_ = np.linalg.lstsq(A, logy, rcond=None)
theta0 = np.array([np.exp(beta[0]), beta[1]])

res = least_squares(residuals, x0=theta0, args=(x, y), method="lm", jac="2-point")
theta_hat = res.x

x_plot = np.linspace(x.min(), x.max(), 400)
y_plot = f_model(x_plot, theta_hat)

plt.figure(figsize=(6.0, 4.2))
plt.scatter(x, y, label="data", zorder=3)
plt.plot(
    x_plot,
    y_plot,
    label=rf"LM fit: $\hat f(x)=\hat\theta_1 e^{{\hat\theta_2 x}}$"
          rf"\n$\hat\theta_1={theta_hat[0]:.3g},\ \hat\theta_2={theta_hat[1]:.3g}$"
)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Exponential fit by Levenberg–Marquardt")
plt.grid(True, alpha=0.25)
plt.legend(frameon=True)
plt.tight_layout()

# Save plot in the same directory as this script
script_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(script_dir, "exp_fit_lm.png")
plt.savefig(output_path, dpi=180)
print(f"Plot saved to: {output_path}")