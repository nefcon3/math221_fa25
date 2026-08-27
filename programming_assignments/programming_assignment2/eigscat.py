import numpy as np
import matplotlib.pyplot as plt
import os

def eigscat(a, err, m, save_path=None, test_label=None):
    """
    Plot eigenvalues of a matrix and eigenvalues of random perturbations
    (real and complex), illustrating sensitivity of eigenvalues.

    Parameters
    ----------
    a : (n, n) array_like
        Input matrix.
    err : float
        Size (Frobenius norm) of the perturbations.
    m : int
        Number of perturbed matrices to compute for each case.
    save_path : str, optional
        Path to save the plot image. If None, displays the plot instead.
    test_label : str, optional
        Label to identify the test (e.g., "Test 1", "Test 6").
    """

    a = np.asarray(a, dtype=complex)

    # --- Compute eigenvalues of unperturbed and perturbed matrices ---
    # ea: eigenvalues of a (conjugated to mimic MATLAB code)
    w = np.linalg.eigvals(a)
    ea = np.conj(w)
    er = []  # real-perturbed eigenvalues
    ei = []  # complex-perturbed eigenvalues

    n = a.shape[0]

    # Real perturbations
    for _ in range(m):
        r = np.random.rand(*a.shape) - 0.5
        r = (err / np.linalg.norm(r)) * r  # scale to have norm 'err'
        w_r = np.linalg.eigvals(a + r)
        er.append(np.conj(w_r))

    # Complex perturbations
    j = 1j
    for _ in range(m):
        r_real = np.random.rand(*a.shape) - 0.5
        r_imag = np.random.rand(*a.shape) - 0.5
        r = r_real + j * r_imag
        r = (err / np.linalg.norm(r)) * r  # scale to have norm 'err'
        w_c = np.linalg.eigvals(a + r)
        ei.append(np.conj(w_c))

    # Convert lists to arrays
    er = np.concatenate(er) if er else np.array([])
    ei = np.concatenate(ei) if ei else np.array([])

    # --- Plot data ---

    # Create a single figure with 3 subplots side-by-side
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    title = f'Eigenvalue Perturbations (err={err:.2e})'
    if test_label:
        title = f'{test_label}: {title}'
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # Calculate combined range for all data to make square plots
    all_real = np.concatenate([ea.real, er.real if len(er) > 0 else [], ei.real if len(ei) > 0 else []])
    all_imag = np.concatenate([ea.imag, er.imag if len(er) > 0 else [], ei.imag if len(ei) > 0 else []])

    # Find the max range and center
    real_range = all_real.max() - all_real.min()
    imag_range = all_imag.max() - all_imag.min()
    max_range = max(real_range, imag_range) * 1.1  # Add 10% padding

    real_center = (all_real.max() + all_real.min()) / 2
    imag_center = (all_imag.max() + all_imag.min()) / 2

    # Set square limits for all plots
    xlim = [real_center - max_range/2, real_center + max_range/2]
    ylim = [imag_center - max_range/2, imag_center + max_range/2]

    # 1. Unperturbed eigenvalues only
    axes[0].plot(ea.real, ea.imag, 'ok', fillstyle='none', label='Unperturbed eigenvalues')
    axes[0].grid(True)
    axes[0].set_aspect('equal')
    axes[0].set_xlim(xlim)
    axes[0].set_ylim(ylim)
    axes[0].set_title('Unperturbed Eigenvalues')
    axes[0].set_xlabel('Real part')
    axes[0].set_ylabel('Imaginary part')
    axes[0].legend()

    # 2. Real perturbations only
    axes[1].plot(er.real, er.imag, 'xr', label='Perturbed eigenvalues (real)')
    axes[1].grid(True)
    axes[1].set_aspect('equal')
    axes[1].set_xlim(xlim)
    axes[1].set_ylim(ylim)
    axes[1].set_title('Real Perturbations')
    axes[1].set_xlabel('Real part')
    axes[1].set_ylabel('Imaginary part')
    axes[1].legend()

    # 3. Complex perturbations only
    axes[2].plot(ei.real, ei.imag, '.b', label='Perturbed eigenvalues (complex)')
    axes[2].grid(True)
    axes[2].set_aspect('equal')
    axes[2].set_xlim(xlim)
    axes[2].set_ylim(ylim)
    axes[2].set_title('Complex Perturbations')
    axes[2].set_xlabel('Real part')
    axes[2].set_ylabel('Imaginary part')
    axes[2].legend()

    # --- Compute "condition numbers" for eigenvalues (as in MATLAB) ---
    # v: eigenvectors, d: diag matrix of eigenvalues
    eigvals, v = np.linalg.eig(a)

    # Normalize each eigenvector column to have unit 2-norm
    for i in range(n):
        v[:, i] /= np.linalg.norm(v[:, i])

    # Try to compute condition numbers; if eigenvector matrix is singular, use inf
    try:
        vi = np.linalg.inv(v)
        cnd = np.zeros(n)
        for i in range(n):
            cnd[i] = np.linalg.norm(vi[i, :])
    except np.linalg.LinAlgError:
        # Eigenvector matrix is singular (defective matrix)
        cnd = np.full(n, np.inf)

    # Format the table as text
    table_text = "       Eigenvalues               Condition Numbers\n"
    table_text += "-" * 60 + "\n"
    for lam, cond in zip(eigvals, cnd):
        table_text += f"{lam: .6e}    {cond: .6e}\n"

    # Add table to the bottom of the figure
    fig.subplots_adjust(bottom=0.3)  # Make room for the table
    fig.text(0.5, 0.12, table_text, ha='center', va='top',
             fontfamily='monospace', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()


# Example usage:
if __name__ == "__main__":
    np.random.seed(0)  # for reproducibility
    m = 1000            # number of perturbed matrices

    # Create output directory
    output_dir = os.path.expanduser("~/Downloads/math221/programming_assignments/programming_assignment2/pseudospectra")
    os.makedirs(output_dir, exist_ok=True)
    print(f"Saving plots to: {output_dir}")

    # (1) a = randn(5)
    while True:
        A1 = np.random.randn(5, 5)
        eig1 = np.linalg.eigvals(A1)
        if np.any(np.abs(eig1.imag) > 1e-12):  # ensure at least one complex eigenvalue
            break
    for i, err in enumerate([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 2e-1]):
        print(f"\n=== Test (1): random 5x5, err={err} ===")
        err_str = f"{err:.2e}".replace(".", "_").replace("-", "m").replace("+", "p")
        save_file = os.path.join(output_dir, f"test1_randn5x5_err{err_str}.png")
        eigscat(A1, err=err, m=m, save_path=save_file, test_label="Test 1")
        print(f"Saved to: {save_file}")

    # (2) a = diag(ones(4,1),1)
    A2 = np.diag(np.ones(4), 1)
    for i, err in enumerate([1e-12, 1e-10, 1e-8]):
        print(f"\n=== Test (2): diag(ones(4),1), err={err} ===")
        err_str = f"{err:.2e}".replace(".", "_").replace("-", "m").replace("+", "p")
        save_file = os.path.join(output_dir, f"test2_diag_err{err_str}.png")
        eigscat(A2, err=err, m=m, save_path=save_file, test_label="Test 2")
        print(f"Saved to: {save_file}")

    # (3) specific 4x4 upper-triangular-ish matrix
    A3 = np.array([
        [1.0, 1e6, 0.0, 0.0],
        [0.0, 2.0, 1e-3, 0.0],
        [0.0, 0.0, 3.0, 10.0],
        [0.0, 0.0, -1.0, 4.0],
    ])
    for i, err in enumerate([1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]):
        print(f"\n=== Test (3): 4x4 with large/small entries, err={err} ===")
        err_str = f"{err:.2e}".replace(".", "_").replace("-", "m").replace("+", "p")
        save_file = os.path.join(output_dir, f"test3_4x4_err{err_str}.png")
        eigscat(A3, err=err, m=m, save_path=save_file, test_label="Test 3")
        print(f"Saved to: {save_file}")

    # (4) a = q * diag(ones(3,1),1) * q'
    Q, _ = np.linalg.qr(np.random.randn(4, 4))
    B4 = np.diag(np.ones(3), 1)   # 4x4 with ones on first superdiagonal
    A4 = Q @ B4 @ Q.T
    for i, err in enumerate([1e-16, 1e-14, 1e-12, 1e-10, 1e-8]):
        print(f"\n=== Test (4): q*diag(ones(3),1)*q^T, err={err} ===")
        err_str = f"{err:.2e}".replace(".", "_").replace("-", "m").replace("+", "p")
        save_file = os.path.join(output_dir, f"test4_qr_err{err_str}.png")
        eigscat(A4, err=err, m=m, save_path=save_file, test_label="Test 4")
        print(f"Saved to: {save_file}")

    # (5) 3x3 upper triangular with large off-diagonals
    A5 = np.array([
        [1.0, 1e3, 1e6],
        [0.0, 1.0, 1e3],
        [0.0, 0.0, 1.0],
    ])
    for i, err in enumerate([1e-7, 1e-6, 5e-6, 8e-6, 1e-5, 1.5e-5, 2e-5]):
        print(f"\n=== Test (5): 3x3 with large off-diagonals, err={err} ===")
        err_str = f"{err:.2e}".replace(".", "_").replace("-", "m").replace("+", "p")
        save_file = os.path.join(output_dir, f"test5_3x3_err{err_str}.png")
        eigscat(A5, err=err, m=m, save_path=save_file, test_label="Test 5")
        print(f"Saved to: {save_file}")

    # (6) 6x6 block upper-triangular matrix
    A6 = np.array([
        [1.0, 0.0, 0.0,   0.0,   0.0,   0.0],
        [0.0, 2.0, 1.0,   0.0,   0.0,   0.0],
        [0.0, 0.0, 2.0,   0.0,   0.0,   0.0],
        [0.0, 0.0, 0.0,   3.0, 1e2, 1e4],
        [0.0, 0.0, 0.0,   0.0,   3.0, 1e2],
        [0.0, 0.0, 0.0,   0.0,   0.0,   3.0],
    ])
    for i, err in enumerate([1e-10, 1e-8, 1e-6, 1e-4, 1e-3]):
        print(f"\n=== Test (6): 6x6 block upper-triangular, err={err} ===")
        err_str = f"{err:.2e}".replace(".", "_").replace("-", "m").replace("+", "p")
        save_file = os.path.join(output_dir, f"test6_6x6_err{err_str}.png")
        eigscat(A6, err=err, m=m, save_path=save_file, test_label="Test 6")
        print(f"Saved to: {save_file}")

    print(f"\nAll plots saved to: {output_dir}")