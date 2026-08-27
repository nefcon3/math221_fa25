import numpy as np
import time

# Global machine epsilon for single precision (same as LAPACK/Fortran)
macheps = np.finfo(np.float32).eps


def generate_system(dimension: int, well_conditioned: bool = True, exact: bool = False, seed: int | None = 0):
    """
    Returns (n, A, x, b, c) where:
      n = dimension
      A in R^{n x n}, x in R^n, b = A @ x, c = cond(A) using inf-norm via exact inverse.

    Parameters
    ----------
    dimension : int
        Size n of the square matrix.
    well_conditioned : bool, default True
        If exact is False, controls the random matrix family.
    exact : bool, default False
        If True, build the special deterministic matrix with ones on diagonal and last column,
        -1 below diagonal, and zeros elsewhere.
    seed : int | None, default 0
        Random seed to make x and (when exact is False) A reproducible.
    """
    n = int(dimension)
    if n <= 0:
        raise ValueError("dimension must be a positive integer")

    rng = np.random.default_rng(seed)

    # Random x (standard normal), single precision
    x = rng.standard_normal(n).astype(np.float32)

    if exact:
        # Build A with: ones on diagonal and last column; -1 below diagonal; zeros elsewhere.
        A = np.zeros((n, n), dtype=np.float32)
        np.fill_diagonal(A, np.float32(1.0))
        A[:, -1] = np.float32(1.0)
        lower_mask = np.tri(n, n, k=-1, dtype=bool)
        lower_mask[:, -1] = False  # exclude the last column from being set to -1
        A[lower_mask] = np.float32(-1.0)
    else:
        if well_conditioned:
            # Random permutation matrix + small Gaussian noise
            perm = rng.permutation(n)
            P = np.eye(n, dtype=np.float32)[perm]
            # Gaussian noise with mean roughly in [1e-5, 1e-2] range
            E = np.abs(rng.normal(0.005, 0.003, size=(n, n)).astype(np.float32))
            E = np.clip(E, 1e-5, 1e-2)
            A = P + E
        else:
            # Ill-conditioned: A = L @ U with tiny diagonals and larger off-diagonals
            # Lower triangular L with Gaussian entries
            L = np.zeros((n, n), dtype=np.float32)
            # Diagonal: Gaussian with mean in [1e-8, 1e-5] range
            diag_L = np.abs(rng.normal(5e-6, 3e-6, size=n).astype(np.float32))
            diag_L = np.clip(diag_L, 1e-8, 1e-5)
            L[np.diag_indices(n)] = diag_L
            # Off-diagonal: Gaussian with mean in [1e-4, 1e-1] range
            sub_L = np.abs(rng.normal(0.05, 0.03, size=(n, n)).astype(np.float32))
            sub_L = np.clip(sub_L, 1e-4, 1e-1)
            L += np.tril(sub_L, k=-1)

            # Upper triangular U with Gaussian entries
            U = np.zeros((n, n), dtype=np.float32)
            diag_U = np.abs(rng.normal(5e-6, 3e-6, size=n).astype(np.float32))
            diag_U = np.clip(diag_U, 1e-8, 1e-5)
            U[np.diag_indices(n)] = diag_U
            sup_U = np.abs(rng.normal(0.05, 0.03, size=(n, n)).astype(np.float32))
            sup_U = np.clip(sup_U, 1e-4, 1e-1)
            U += np.triu(sup_U, k=1)

            A = L @ U  # remains float32

    b = A @ x  # float32

    # Condition number using exact inverse and infinity norm (single precision)
    A_inv = np.linalg.inv(A).astype(np.float32)
    normA = np.float32(np.linalg.norm(A, np.inf))
    normAinv = np.float32(np.linalg.norm(A_inv, np.inf))
    c = np.float32(normA * normAinv)

    return n, A, x, b, c

def gepp_outputs(A, b):
    """
    Solve Ax=b using LAPACK SGESVX via f2py.

    Parameters
    ----------
    A : array_like
        Coefficient matrix (will be converted to float32)
    b : array_like
        Right-hand side vector (will be converted to float32)

    Returns
    -------
    xprime : ndarray
        Solution vector (float32)
    pivot_growth : float
        Pivot growth factor = 1/WORK(1)
    cond_estimate : float
        Condition number estimate = 1/RCOND
    ferr : float
        Forward error bound
    berr : float
        Backward error bound
    """
    import gepp_wrapper

    # Convert to single precision (float32) and ensure proper layout
    A32 = np.array(A, dtype=np.float32, order='F')
    b32 = np.array(b, dtype=np.float32).reshape(-1)
    n = A32.shape[0]

    # Call the Fortran subroutine
    x, rcond, ferr, berr, work1 = gepp_wrapper.gepp_wrapper(A32, b32)

    # Compute outputs as specified
    xprime = x.astype(np.float32)
    pivot_growth = np.float32(1.0 / work1) if work1 != 0 else np.float32(np.inf)
    cond_estimate = np.float32(1.0 / rcond) if rcond != 0 else np.float32(np.inf)
    ferr = np.float32(ferr)
    berr = np.float32(berr)

    return xprime, pivot_growth, cond_estimate, ferr, berr
def gecp_outputs(A, b):
    """
    Solve Ax=b using LAPACK complete pivoting (SGETC2/SGESC2) via f2py.

    Parameters
    ----------
    A : array_like
        Coefficient matrix (will be converted to float32)
    b : array_like
        Right-hand side vector (will be converted to float32)

    Returns
    -------
    xprime : ndarray
        Solution vector (float32)
    pivot_growth : float
        Pivot growth factor = 1/WORK(1)
    cond_estimate : float
        Condition number estimate = 1/RCOND
    ferr : float
        Forward error bound
    berr : float
        Backward error bound
    """
    import gecp_wrapper

    # Convert to single precision (float32) and ensure proper layout
    A32 = np.array(A, dtype=np.float32, order='F')
    b32 = np.array(b, dtype=np.float32).reshape(-1)
    n = A32.shape[0]

    # Call the Fortran subroutine
    x, rcond, ferr, berr, work1 = gecp_wrapper.gecp_wrapper(A32, b32)

    # Compute outputs as specified
    xprime = x.astype(np.float32)
    pivot_growth = np.float32(1.0 / work1) if work1 != 0 else np.float32(np.inf)
    cond_estimate = np.float32(1.0 / rcond) if rcond != 0 else np.float32(np.inf)
    ferr = np.float32(ferr)
    berr = np.float32(berr)

    return xprime, pivot_growth, cond_estimate, ferr, berr
def postprocessing(xprime, pivot_growth, cond_estimate, ferr, berr, A, x, b, c):
    """
    Compute full diagnostic report from solver outputs and true solution.

    Parameters
    ----------
    xprime : array_like
        Computed solution vector from solver
    pivot_growth : float
        Pivot growth factor from solver
    cond_estimate : float
        Estimated condition number (1/RCOND) from solver
    ferr : float
        Forward error bound from solver
    berr : float
        Backward error bound from solver
    A : array_like
        Coefficient matrix
    x : array_like
        True solution vector (for error comparison)
    b : array_like
        Right-hand side vector
    c : float
        True condition number of A (infinity norm)

    Returns
    -------
    dict
        Dictionary containing:
        - pivot_growth: Pivot growth factor
        - cond_estimate: Estimated condition number (1/RCOND)
        - cond_to_tcond: Ratio of estimated to true condition number
        - ferr: Forward error bound from solver
        - ferr_to_terr: Ratio of FERR to actual error
        - terr_to_cond_macheps: Scaled true error
        - scaled_berr: Backward error scaled by machine epsilon
        - berr_to_macheps: BERR divided by machine epsilon
    """
    # Ensure all inputs are single precision
    n = A.shape[0]
    A32 = np.array(A, dtype=np.float32, order='F')
    x32 = np.array(x, dtype=np.float32).reshape(-1)
    b32 = np.array(b, dtype=np.float32).reshape(-1)
    xprime32 = np.array(xprime, dtype=np.float32).reshape(-1)

    # Direct outputs
    pivot_growth = np.float32(pivot_growth)
    cond_estimate = np.float32(cond_estimate)
    ferr = np.float32(ferr)
    berr = np.float32(berr)
    c = np.float32(c)

    # cond_to_tcond = cond_estimate / c
    cond_to_tcond = (cond_estimate / c) if c != 0 else np.float32(np.inf)

    # terr = ||xprime - x||_inf / ||xprime||_inf
    num_terr = np.float32(np.linalg.norm(xprime32 - x32, ord=np.inf))
    den_terr = np.float32(np.linalg.norm(xprime32, ord=np.inf))
    terr = (num_terr / den_terr) if den_terr != 0.0 else (0.0 if num_terr == 0.0 else np.float32(np.inf))

    # ferr_to_terr = ferr / terr
    ferr_to_terr = (ferr / terr) if np.isfinite(terr) and terr != 0.0 else np.float32(np.inf)

    # terr_to_cond_macheps = terr / (macheps * cond_estimate)
    if np.isfinite(cond_estimate) and cond_estimate > 0.0:
        terr_to_cond_macheps = np.float32(terr / (macheps * cond_estimate))
    else:
        terr_to_cond_macheps = np.float32(np.inf)

    # scaled_berr = ||A*xprime - b||_inf / ((||A||_inf * ||xprime||_inf + ||b||_inf) * macheps)
    res = A32 @ xprime32 - b32
    num_scaled_berr = np.float32(np.linalg.norm(res, ord=np.inf))
    norm_A = np.float32(np.linalg.norm(A32, ord=np.inf))
    norm_xprime = np.float32(np.linalg.norm(xprime32, ord=np.inf))
    norm_b = np.float32(np.linalg.norm(b32, ord=np.inf))
    den_scaled_berr = (norm_A * norm_xprime + norm_b) * macheps
    scaled_berr = (num_scaled_berr / den_scaled_berr) if den_scaled_berr != 0.0 else (0.0 if num_scaled_berr == 0.0 else np.float32(np.inf))

    # berr_to_macheps = berr / macheps
    berr_to_macheps = np.float32(berr / macheps)

    return {
        "pivot_growth": pivot_growth,
        "cond_estimate": cond_estimate,
        "cond_to_tcond": cond_to_tcond,
        "ferr": ferr,
        "ferr_to_terr": ferr_to_terr,
        "terr_to_cond_macheps": terr_to_cond_macheps,
        "scaled_berr": scaled_berr,
        "berr_to_macheps": berr_to_macheps,
    }

def gepp_report(n, A, x, b, c):
    """
    Solve A x' = b with LAPACK SGESVX and report:
      - pivot growth factor
      - estimated condition number 1/RCOND
      - (1/RCOND)/c  where c is the provided exact cond(A, inf)
      - FERR
      - FERR / true_error, where true_error = ||x' - x||_inf / ||x'||_inf
      - true_error / (macheps * cond_estimate)
      - scaled backward error: ||A x' - b||_inf / ((||A||_inf*||x||_inf + ||b||_inf) * macheps)
      - BERR / macheps
    Returns a dict of these values.

    Notes:
      * Single-precision throughout (float32).
      * This function calls gepp_outputs then postprocessing.
    """
    xprime, pivot_growth, cond_estimate, ferr, berr = gepp_outputs(A, b)
    return postprocessing(xprime, pivot_growth, cond_estimate, ferr, berr, A, x, b, c)

def gecp_report(n, A, x, b, c):
    """
    Solve A x' = b with a variant of LAPACK SGESVX but with complete pivoting (GECP) and report:
      - pivot growth factor
      - estimated condition number 1/RCOND (via SGECON on ORIGINAL A, ∞-norm)
      - (1/RCOND)/c  where c = exact cond(A, ∞) provided by caller
      - FERR (normwise forward-error bound using cond(A,∞) and residual; no refinement)
      - FERR / true_error, where true_error = ||x' - x||_inf / ||x'||_inf
      - true_error / (macheps * cond_estimate)
      - scaled backward error: ||A x' - b||_inf / ((||A||_inf*||x||_inf + ||b||_inf) * macheps)
      - BERR / macheps   (componentwise relative backward error)
    Returns a dict of these values.

    Notes
    -----
    * Single-precision throughout (float32).
    * This function calls gecp_outputs then postprocessing.
    """
    xprime, pivot_growth, cond_estimate, ferr, berr = gecp_outputs(A, b)
    return postprocessing(xprime, pivot_growth, cond_estimate, ferr, berr, A, x, b, c)

def timed_report(n, A, x, b, c, m=10000):
    """
    Time and return reports for both GEPP and GECP.

    Returns:
      {
        "gepp": {<gepp_report fields...>, "runtime_sec": t_gepp},
        "gecp": {<gecp_report   fields...>, "runtime_sec": t_gecp},
      }

    Timing scheme (identical for both):
        t1 = time
        repeat m times: set up problem; solve the problem (Fortran only)
        t2 = time
        repeat m times: set up problem
        t3 = time
        t = ((t2 - t1) - (t3 - t2)) / m

    Notes:
      * "solve the problem" calls gepp_outputs(...) or gecp_outputs(...) (Fortran only).
      * Postprocessing is NOT timed, only the Fortran solver calls.
      * "set up problem" mirrors the float32 casting, Fortran order, and reshaping
        used inside gepp_outputs/gecp_outputs, so the subtraction removes setup overhead.
    """
    if m <= 0:
        raise ValueError("m must be a positive integer")

    # One real run of each to warm up and to capture full outputs (with postprocessing)
    out_gepp = gepp_report(n, A, x, b, c)
    out_gecp = gecp_report(n, A, x, b, c)

    def _setup_only():
        # Mirrors input preparation done in gepp_outputs/gecp_outputs
        A32 = np.array(A, dtype=np.float32, order='F')
        b32 = np.array(b, dtype=np.float32).reshape(-1)
        return A32, b32  # prevent optimization-away

    # ----- Time GEPP path (Fortran only) -----
    t1 = time.perf_counter()
    for _ in range(m):
        _ = gepp_outputs(A, b)
    t2 = time.perf_counter()
    for _ in range(m):
        _ = _setup_only()
    t3 = time.perf_counter()
    t_gepp = ((t2 - t1) - (t3 - t2)) / m
    if t_gepp < 0:
        t_gepp = 0.0

    # ----- Time GECP path (Fortran only) -----
    t4 = time.perf_counter()
    for _ in range(m):
        _ = gecp_outputs(A, b)
    t5 = time.perf_counter()
    for _ in range(m):
        _ = _setup_only()
    t6 = time.perf_counter()
    t_gecp = ((t5 - t4) - (t6 - t5)) / m
    if t_gecp < 0:
        t_gecp = 0.0

    # Attach runtimes
    out_gepp_rt = dict(out_gepp)
    out_gepp_rt["runtime_sec"] = t_gepp

    out_gecp_rt = dict(out_gecp)
    out_gecp_rt["runtime_sec"] = t_gecp

    return {"gepp": out_gepp_rt, "gecp": out_gecp_rt}