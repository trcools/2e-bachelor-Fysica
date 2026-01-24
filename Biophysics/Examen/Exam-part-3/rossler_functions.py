"""
Rössler attractor system.

The Rössler system is a three-dimensional autonomous ordinary differential equation:
    x' = -y - z
    y' = x + a*y
    z' = b + z*(x - c)

where (a, b, c) are system parameters.

The system displays chaotic behavior for certain parameter values and has rich dynamics
including bifurcations and transitions to chaos as parameters vary.
"""

"""
Note:

This module was refactored with assistance from an AI tool to improve structure,
naming, and documentation for readability and maintainability.

"""

# =============================================================================
# Imports
# =============================================================================

import numpy as np 

import sympy as sp
from sympy.solvers import solve
from sympy import Symbol, Matrix, symbols, Eq
    
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D 

import scipy
from scipy.signal import find_peaks

from ipywidgets import interact
    

# =============================================================================
# General Parameters
# =============================================================================

def get_parameters():
    
    # Rössler system parameters
    a = 0.2
    b = 0.2
    c = 5.7
    # c will be varied as the bifurcation parameter
        
    # --- Time settings ---
    t_min, t_max = 0.0, 200.0  # Longer time for Rössler to show attractor
    dt = 0.01
        
    # --- Initial conditions ---
    x0, y0, z0 = 1, 1, 1
    initial_condition = (x0, y0, z0)
    
    initial_conditions = [
        (0.1, 0.1, 0.1),
        (5.0, 5.0, 5.0),
        (-5.0, 2.0, 1.0),
        (0.5, 0.5, 0.5)
    ]
    return a, b, c, t_min, t_max, dt, initial_condition, initial_conditions

# ------- parameters ----

a, b, c, t_min, t_max, dt, initial_condition, initial_conditions = get_parameters()

# ============================================================
# 1. Numerical integration 
# ============================================================

    # ========================================
    # 1.1. Time grid and state arrays
    # ========================================

def get_time_array(t_min=t_min, t_max=t_max, dt=dt):
    """
    Construct a time array for numerical integration.

    Parameters
    ----------
    t_min, t_max : float
        Start and end of the time interval.
    dt : float
        Time step.

    Returns
    -------
    t : ndarray
        One-dimensional array of time points in [t_min, t_max).
    """
    t = np.arange(t_min, t_max, dt)
    return t

def get_state_arrays(n, initial_condition):
    """
    Initialize state arrays for numerical integration.
    
    Parameters
    ----------
    n : int
        Number of time points
    initial_condition : tuple of float
        Initial state (x0, y0, z0)
        
    Returns
    -------
    x, y, z : ndarray
        Arrays of size n with initial conditions set
    """
    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)
    x[0], y[0], z[0] = initial_condition
    return x, y, z


# ========================================
# 1.2 Core Rössler Dynamics
# ========================================


def rossler_rhs(state, a=0.2, b=0.2, c=5.7):
    """
    Compute the right-hand side of the Rössler system.
    
    Parameters
    ----------
    state : array-like, shape (3,)
        Current state (x, y, z)
    a, b, c : float
        System parameters
        
    Returns
    -------
    dstate : ndarray, shape (3,)
        Time derivatives (dx/dt, dy/dt, dz/dt)
    """
    x, y, z = state
    dx = -y - z
    dy = x + a * y
    dz = b + z * (x - c)
    return np.array([dx, dy, dz])


# ============================================================
# 1.3. Euler integration (using explicit Euler method)
# # ============================================================


def get_euler_vectors(a, b, c, initial_condition, t_min, t_max, h):
    """Integrate Rössler system using explicit Euler method."""

    t = get_time_array(t_min=t_min, t_max=t_max + 1e-15, dt=h)
    n = len(t)
    x, y, z = get_state_arrays(n, initial_condition)

    for i in range(n - 1):
        dx, dy, dz = rossler_rhs((x[i], y[i], z[i]), a=a, b=b, c=c)
        x[i + 1] = x[i] + h * dx
        y[i + 1] = y[i] + h * dy
        z[i + 1] = z[i] + h * dz

    return t, x, y, z


# ============================================================
# 1.4. Manual Integrator : RK4
# ============================================================

def rk4_integrate_rossler(a, b, c, initial_condition, t_min, t_max, h):
    """
    Classic RK4 integrator for the 3D Rössler system.
    Returns t, x, y, z on a uniform grid with step h.
    """
    # Arrays and time grid
    t = np.arange(t_min, t_max + 1e-15, h)
    n = len(t)
    x, y, z = get_state_arrays(n, initial_condition)

    def f(state):
        "Helper function to reduce code duplication."
        return rossler_rhs(state, a=a, b=b, c=c)

    for i in range(n - 1):
        s = np.array([x[i], y[i], z[i]])
        k1 = f(s)
        k2 = f(s + 0.5 * h * k1)
        k3 = f(s + 0.5 * h * k2)
        k4 = f(s + h * k3)

        s_new = s + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        x[i + 1], y[i + 1], z[i + 1] = s_new[0], s_new[1], s_new[2]

    return t, x, y, z

def get_RK4_vectors(a=a, b=b, c=c, initial_condition=initial_condition,
                    t_min=t_min, t_max=t_max, dt=dt):
    "Wrapper function to get RK4 vectors."
    t, x, y, z = rk4_integrate_rossler(a, b, c, initial_condition, t_min, t_max, dt)
    return t, x, y, z


# ============================================================
# 2. Symbolic equilibria of the Rössler system
# ============================================================


def get_solutions():
    """
    Compute the symbolic equilibria of the Rössler system using SymPy.

    Returns
    -------
    sols : list[dict]
        SymPy solutions as dicts {x: expr, y: expr, z: expr}.
    (x, y, z) : tuple of SymPy symbols
        The state variables used in the symbolic system.
    """
    x, y, z, c, a_sym, b_sym = symbols('x y z c a b', real=True)

    eqs = [
        Eq(-y - z, 0),
        Eq(x + a_sym * y, 0),
        Eq(b_sym + z * (x - c), 0)]
    
    solutions = solve(eqs, (x, y, z), dict=True)

    return solutions, (x, y, z)

def print_solutions():
    """
    Print all symbolic equilibria in a human-readable format.
    """
    solutions, (x, y, z) = get_solutions()
    for i,solution in enumerate(solutions):
        print(f"Solution {i+1}:\n X = {solution[x]} Y = {solution[y]} Z = {solution[z]} \n")


# ============================================================
# 3. Classification and grouping of equilibria
# ============================================================


def equilibrium_classification(eigvals):
    """
    Classify an equilibrium based on the eigenvalues of the Jacobian.

    Uses a small tolerance to decide whether eigenvalues are effectively zero.

    Parameters
    ----------
    eigvals : array-like of complex
        Eigenvalues of the Jacobian at the equilibrium.

    Returns
    -------
    eq_type : str
         One of:
        - "stable node"
        - "stable spiral"
        - "unstable node"
        - "unstable spiral"
        - "saddle point"
        - "saddle-spiral"
        - "pitchfork bifurcation"
        - "Hopf bifurcation"
        - "other equilibrium"
    """
    # Introduce a tolerance to find the eigenvalues of zero:
    tol = 1e-9
    re = np.real(eigvals)
    im = np.imag(eigvals)
    has_complex = np.any(np.abs(im) > tol)
    
    # Count the number 'n' of positive, negative and zero eigenvalues
    n_pos = np.sum(re > tol)
    n_neg = np.sum(re < -tol)
    n_zero = len(eigvals) - n_pos - n_neg
    
    # 1) All real parts < 0  -> stable
    if n_pos == 0 and n_zero == 0:
        return "stable spiral" if has_complex else "stable node"
    
    # 2) At least one > 0, none negative -> purely unstable
    if n_pos > 0 and n_neg == 0 and n_zero == 0:
        return "unstable spiral" if has_complex else "unstable node"
    
    # 3) Both positive and negative real parts -> saddle / saddle-focus
    if n_pos > 0 and n_neg > 0:
        return "saddle spiral" if has_complex else "saddle point"
    
    # 4) At least one eigenvalue ≈ 0, all real -> pitchfork
    if n_zero > 0 and not has_complex and n_pos == 0:
        return "pitchfork bifurcation"
    
    # 5) At least one eigenvalue ≈ 0, complex pair present -> Hopf
    if n_zero > 0 and has_complex and n_pos == 0:
        return "Hopf bifurcation"
    
    # 6) Exotic / other cases
    return "other equilibrium"
 
def get_groups(dim="1D"):
    """
    Create plotting groups for equilibria.

    Parameters
    ----------
    dim : {"1D", "2D"}
        - "1D": (c, val)
        - "2D": (c, v1, v2)

    Returns
    -------
    dict
        keys = equilibrium type
        values = dict with coordinate lists + style dict
    """
    if dim not in {"1D", "2D"}:
        raise ValueError(f"dim must be '1D' or '2D', got {dim!r}")
        
    
    base_styles = styles_for_equilibria()
                  
    groups = {}

    if dim == "2D":
        for key, style in base_styles.items():
            groups[key] = {
                "c": [],
                "v1": [],
                "v2": [],
                "style": style.copy(),
            }
    else: # dim == "1D"
        for key, style in base_styles.items():
            groups[key] = {
                "c": [],
                "val": [],
                "style": style.copy(),
            }

    return groups

def styles_for_equilibria():
    """Return the full style mapping for all equilibrium types."""
    return {
        "stable node": dict(marker="o", color="C2", label="stable node"),
        "stable spiral": dict(marker="o", color="darkgreen", label="stable-spiral"),
        "unstable node": dict(marker="x", color="red", label="unstable node"),
        "unstable spiral": dict(marker="x", color="black", label="unstable-spiral"),
        "saddle point": dict(marker="o", color="orange", label="saddle"),
        "saddle spiral": dict(marker="o", color="green", label="saddle-spiral"),
        "pitchfork bifurcation": dict(marker="o", color="red", label="Pitchfork bifurcation"),
        "Hopf bifurcation": dict(marker="o", color="blue", label="Hopf bifurcation"),
    }


def styles_for_equilibrium(eq_type):
    """Return the style for a single equilibrium type."""
    return styles_for_equilibria().get(eq_type, dict(marker="o", color="gray", label=eq_type))

# ============================================================
# 4. Jacobian, equilibria and eigenvalues as function of c
# ============================================================

def jacobian_rossler(state, a, b, c):
    """
    Compute Jacobian matrix of Rössler system at a given state.
    
    Parameters
    ----------
    state : array-like, shape (3,)
        Current state (x, y, z)
    a, b, c : float
        Rössler parameters
        
    Returns
    -------
    J : ndarray, shape (3, 3)
        Jacobian matrix at state
    """
    x, y, z = state
    J = np.array([
        [0.0,   -1.0,   -1.0    ],
        [1.0,    a,      0.0    ],
        [z,      0.0,   x - c   ]
    ], dtype=float)
    return J 

def get_equilibria_for_c(sol, x, y, z, c_sym, a_sym, b_sym, c_val, a_val=a, b_val=b):
    """
    Compute the equilibrium for a single symbolic solution and a given c.

    Only real equilibria are kept.

    Parameters
    ----------
    sol : dict
        SymPy solution {x: expr, y: expr, z: expr}.
    x, y, z : SymPy symbols
        State variables.
    c_sym, a_sym, b_sym : SymPy symbols
        Parameter symbols (c, a, b).
    c_val : float
        Value at which to evaluate the equilibrium.
    a_val : float
        Rössler parameter a.
    b_val : float
        Rössler parameter b.

    Returns
    -------
    eq : dict or None
        If real, a dict with keys
        {"c", "x_eq", "y_eq", "z_eq", "J", "eigvals", "eq_type"},
        otherwise None.
    """

    x_expr = sol[x]
    y_expr = sol[y]
    z_expr = sol[z]
        
    subs_dict = {c_sym: c_val, a_sym: a_val, b_sym: b_val}

    # Here we substitute c, a and b into the symbolic equations
    x_num = x_expr.subs(subs_dict).evalf()
    y_num = y_expr.subs(subs_dict).evalf()
    z_num = z_expr.subs(subs_dict).evalf()
    # evalf: evaluate to a floating-point number

    #Only keep real solutions
    if not (x_num.is_real and y_num.is_real and z_num.is_real):
        return None
                
    xsol = float(x_num)
    ysol = float(y_num)
    zsol = float(z_num)

    # Jacobian at equilibrium for Rössler system
    J_eq = jacobian_rossler(state=(xsol, ysol, zsol), a=a_val, b=b_val, c=c_val)

    eigvals, _ = np.linalg.eig(J_eq)
    eq_type = equilibrium_classification(eigvals)
            
    return {
        "c": c_val,
        "x_eq": xsol,
        "y_eq": ysol,
        "z_eq": zsol,
        "J": J_eq,
        "eigvals": eigvals,
        "eq_type": eq_type}
    
def get_equilibria_from_jacobian(sol, x, y, z, c_sym, a_sym, b_sym, c_values,
                                 a_val=a, b_val=b):
    """
    Helper function: compute equilibria for many c values for a single solution.

    Parameters
    ----------
    sol : dict
        SymPy solution {x: expr, y: expr, z: expr}.
    x, y, z : SymPy symbols
        State variables.
    c_sym, a_sym, b_sym : SymPy symbols
        Parameter symbols.
    c_values : iterable of float
        Values of c to scan.
    a_val, b_val : float
        Rössler parameters.

    Returns
    -------
    equilibrium : list[dict]
        List of equilibrium dictionaries (see get_equilibria_for_c).
    """
    equilibrium = []
    
    for c_val in c_values:
        eq = get_equilibria_for_c(sol, x, y, z, c_sym, a_sym, b_sym, c_val,
                                    a_val=a_val, b_val=b_val)
        if eq is not None:
                equilibrium.append(eq)
    return equilibrium

def compute_jacobian_and_stability(sols, x, y, z, small_c=False):
    """
    Compute Jacobian matrices, eigenvalues and stability classification
    for equilibria over a grid of c values.

    Parameters
    ----------
    sols : list[dict]
        SymPy solutions {x: expr, y: expr, z: expr}.
    x, y, z : SymPy symbols
        State variables.
    small_c : bool
        If True, use c in [0, 2]; otherwise use c in [0, 18].

    Returns
    -------
    all_equilibria : list[dict]
        Each dict contains keys:
        {"c", "x_eq", "y_eq", "z_eq", "J", "eigvals", "eq_type"}.
    """
    c_sym, a_sym, b_sym = symbols('c a b', real=True)
    all_equilibria = []
    
    if small_c: 
        c_values = np.linspace(0, 2, 51)
    else:
        c_values = np.linspace(0, 18, 181)
        
    for sol in sols:
        all_equilibria.extend(get_equilibria_from_jacobian(sol, x, y, z,
                                                           c_sym, a_sym, b_sym,
                                                           c_values,
                                                           a_val=a, b_val=b))
    return all_equilibria
 
def get_eigvals_for_c(sols, x, y, z, c_sym, a_sym, b_sym, a_val=a, b_val=b, h=dt, small_c=False):
    """
    Collect the Jacobian eigenvalues of all real equilibria as a function of c.

    For small_c = True, c is scanned over [0, 2).
    For small_c = False, c is scanned over [2, 18).

    Parameters
    ----------
    sols : list of SymPy solutions {x: expr, y: expr, z: expr}
    x, y, z : SymPy symbols
        State variables.
    c_sym, a_sym, b_sym : SymPy symbols
        Parameter symbols.
    a_val, b_val : float
        Rössler parameters.
    h : float
        Step in c-space.
    small_c : bool
        Whether to scan [0,2) or [2,18).

    Returns
    -------
    eigvals_by_c : dict
        keys: c values (floats)
        values: list of eigenvalue arrays (one per equilibrium at that c).
    """
    if small_c:
        c_values = np.arange(0, 2, h)
    else:
        c_values = np.arange(2, 18, h)
    eigvals_by_c = {}

    for c_val in c_values:
        eigvals_list = []
        for sol in sols:
            eq = get_equilibria_for_c(sol, x, y, z,
                                        c_sym, a_sym, b_sym, c_val,
                                        a_val=a_val, b_val=b_val)
            if eq is not None:
                eigvals_list.append(eq["eigvals"])
        eigvals_by_c[c_val] = eigvals_list

    return eigvals_by_c

def equilibria_at_c(c_val, a_val=a, b_val=b):
    sols, (x_sym, y_sym, z_sym) = get_solutions()
    c_sym, a_sym, b_sym = symbols("c a b", real=True)

    eqs = []
    for sol in sols:
        eq = get_equilibria_for_c(sol, x_sym, y_sym, z_sym,
                                  c_sym, a_sym, b_sym,
                                  c_val, a_val=a_val, b_val=b_val)
        if eq is not None:
            eqs.append(eq)
    return eqs


# ============================================================
#  Find fitting Euler dt for given RK4 precision 
# ============================================================


def max_abs_error_on_common_grid(t_ref, x_ref, t_test, x_test):
    "Compute max abs error by interpolating test solution onto reference grid."
    x_test_on_ref = np.interp(t_ref, t_test, x_test)
    return np.max(np.abs(x_test_on_ref - x_ref))

def find_euler_dt_for_rk4_precision(a, b, c, ic, t_min, t_max, h_rk4=0.1, h0_euler=0.1, max_halvings=20):
    """
    Find Euler step size that achieves comparable accuracy to RK4 reference.
    
    Uses RK4 with step h_rk4 as reference. Halves Euler step size until 
    error <= RK4 baseline error.
    
    Returns
    -------
    h_best, err_best : tuple[float, float]
        Best Euler step size found and maximum error vs. RK4 reference.
        If no suitable step is found, returns last tested (h, err) pair.
    """
    # Reference: RK4(h=0.1)
    t_ref, x_ref, y_ref, z_ref = rk4_integrate_rossler(a, b, c, ic, t_min, t_max, h_rk4)

    # Baseline "precision target": use RK4 vs RK4 with smaller h as internal reference
    # Avoids needing an exact solution (theory from course python for scientists)
    t_fine, x_fine, y_fine, z_fine = rk4_integrate_rossler(a, b, c, ic, t_min, t_max, h_rk4/2)

    target = max(
        max_abs_error_on_common_grid(t_ref, x_ref, t_fine, x_fine),
        max_abs_error_on_common_grid(t_ref, y_ref, t_fine, y_fine),
        max_abs_error_on_common_grid(t_ref, z_ref, t_fine, z_fine))

    h = h0_euler
    best_h = None
    best_err = None

    for _ in range(max_halvings):
        t_eu, x_eu, y_eu, z_eu = get_euler_vectors(a, b, c, ic, t_min, t_max, h)
        err = max(
            max_abs_error_on_common_grid(t_ref, x_ref, t_eu, x_eu),
            max_abs_error_on_common_grid(t_ref, y_ref, t_eu, y_eu),
            max_abs_error_on_common_grid(t_ref, z_ref, t_eu, z_eu))

        # If error is not finite (overflow/NaN), reduce step and continue
        if not np.isfinite(err):
            h /= 2
            continue

        if err <= target:
            best_h = h
            best_err = err
            break

        h /= 2

    # Return result (found or last attempt)
    if best_h is not None and best_err is not None:
        return best_h, best_err
    else:
        # Recompute last error if needed
        t_eu, x_eu, y_eu, z_eu = get_euler_vectors(a, b, c, ic, t_min, t_max, h)
        err = max(
            max_abs_error_on_common_grid(t_ref, x_ref, t_eu, x_eu),
            max_abs_error_on_common_grid(t_ref, y_ref, t_eu, y_eu),
            max_abs_error_on_common_grid(t_ref, z_ref, t_eu, z_eu))
        return h, err


# ============================================================
# 6. Local maxima of Z(t) + return map
# ============================================================


def get_Z_maxima(a, b, c, initial_condition, skip_first=5, t_min=0, t_max=100, dt=0.01,
                 prominence=None, distance=None):
    """
    Determine the local maxima Z_n of Z(t) for a given c.

    Parameters
    ----------
    a, b, c : float
        Rössler parameters.
    initial_condition : tuple
        Initial condition (x0, y0, z0).
    skip_first : int
        Number of first maxima to discard as transient.
    t_min, t_max : float
        Time interval.
    dt : float
        Time step.

    Additional peak parameters (passed to scipy.signal.find_peaks)
    -------------------------------------------------------------
    prominence : float or None
        Minimum prominence of peaks; helps ignore small noisy bumps.
    distance : float or None
        Minimum horizontal distance between neighboring peaks (in samples).

    Returns
    -------
    t_max : ndarray
        Times at which local maxima occur (after discarding the first ones).
    Z_max : ndarray
        Values of Z at those maxima.
    """
    t, _, _, z = get_RK4_vectors(a=a, b=b, c=c, initial_condition=initial_condition,
                             t_min=t_min, t_max=t_max, dt=dt)

    # Robust peak detection using prominence/distance to suppress noise-driven peaks
    peaks, _ = find_peaks(z, prominence=prominence, distance=distance)

    Z_max = z[peaks]
    t_max = t[peaks]

    # Throw the first maxima away ('skip_first').
    if len(Z_max) > skip_first:
        Z_max = Z_max[skip_first:]
        t_max = t_max[skip_first:]

    return t_max, Z_max, t, z

def rossler_return_map_from_data(Z_n):
    """
    Construct an approximate 1D map Z_{n+1} = F(Z_n) from a sequence of maxima.

    Parameters
    ----------
    Z_n : array-like
        Sequence of maxima Z_0, Z_1, ..., Z_{N-1}.

    Returns
    -------
    F : callable
        Interpolated map F(z) based on the pairs (Z_n, Z_{n+1}).
    """
    Z_n = np.asarray(Z_n)
    x_data = Z_n[:-1]
    y_data = Z_n[1:]

    # Sort by x so that np.interp is well-defined
    idx = np.argsort(x_data)
    x_sorted = x_data[idx]
    y_sorted = y_data[idx]

    def F(z):
        return np.interp(z, x_sorted, y_sorted)

    return F

def local_maxima(values):
    values = np.asarray(values)
    idx = np.where((values[1:-1] > values[:-2]) & (values[1:-1] > values[2:]))[0] + 1
    return idx





def lyapunov_exponent_estimate(a=0.2, b=0.2, c=5.7, initial_condition=initial_condition,
                               t_min=0, t_max=100, dt=dt, eps=1e-8):
    """
    Estimate the largest Lyapunov exponent numerically.
    
    Parameters
    ----------
    x0, y0, z0 : float
        Initial conditions
    a, b, c : float
        System parameters
    t_span : tuple
        Time interval
    eps : float
        Initial perturbation magnitude
        
    Returns
    -------
    t : ndarray
        Time points
    lyap_exp : ndarray
        Lyapunov exponent estimate over time
    """
    # Reference trajectory
    t, x1, y1, z1 = get_RK4_vectors(
        a=a, b=b, c=c,
        initial_condition=initial_condition,
        t_min=t_min, t_max=t_max, dt=dt
        )
    
    # Perturbed trajectory
    t, x2, y2, z2 = get_RK4_vectors(
        a=a, b=b, c=c,
        initial_condition=(initial_condition[0] + eps,
                           initial_condition[1],
                           initial_condition[2]),
        t_min=t_min, t_max=t_max, dt=dt
        )
    
    # Distance between trajectories
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    distance = np.sqrt(dx**2 + dy**2 + dz**2)
    
    # Avoid log(0)
    distance = np.maximum(distance, 1e-12)

    def lyapunov_from_two_trajs(t, distance, eps, t_fit=(20, 60)):
        # Take frame where log-distance approximately linear
        mask = (t >= t_fit[0]) & (t <= t_fit[1])
        y = np.log(distance[mask] / eps)
        x = t[mask]
        lam = np.polyfit(x, y, 1)[0]   # slope
        return lam
    
    lam = lyapunov_from_two_trajs(t, distance, eps=1e-8, t_fit=(20, 60))
    print("Largest Lyapunov exponent ≈", lam)
    return

def bifurcation_maxima_X(param_name, param_values, a, b, c, ic,
                         t_min=0.0, t_max=200.0, h=0.01, skip_frac=0.5):
    """
    Returns arrays (p_plot, x_max_plot) for maxima of X after transients.
    param_name: "b" or "c"
    """
    p_plot = []
    x_plot = []

    for p in param_values:
        aa, bb, cc = a, b, c
        if param_name == "b":
            bb = p
        elif param_name == "c":
            cc = p
        else:
            raise ValueError("param_name must be 'b' or 'c'")

        t, x, _, _ = rk4_integrate_rossler(aa, bb, cc, ic, t_min, t_max, h)
        skip = int(len(t) * skip_frac)
        x_use = x[skip:]

        idx = local_maxima(x_use)
        for j in idx:
            x_plot.append(x_use[j])
            p_plot.append(p)

    return np.array(p_plot), np.array(x_plot)


# ============================================================
# 7. Euler vs RK4 Comparison
# ============================================================


def find_optimal_euler_step(a, b, c, initial_condition, t_min=0.0, t_max=60.0, 
                            h_rk4=0.1, h0_euler=0.1, max_halvings=20):
    """
    Find optimal Euler step size for comparable accuracy to RK4 reference.
    
    Parameters
    ----------
    a, b, c : float
        System parameters
    initial_condition : tuple of float
        Initial state (x0, y0, z0)
    t_min, t_max : float
        Integration time interval
    h_rk4 : float
        Reference RK4 step size
    h0_euler : float
        Initial Euler step size for halving search
    max_halvings : int
        Maximum number of step halvings to attempt
    
    Returns
    -------
    h_best : float
        Optimal Euler step size
    err_best : float
        Maximum error vs RK4 reference
    """
    h_best, err_best = find_euler_dt_for_rk4_precision(
        a, b, c, initial_condition, t_min, t_max, h_rk4, h0_euler, max_halvings
    )
    return h_best, err_best


def compare_euler_vs_rk4(a, b, c, initial_condition, t_min=0.0, t_max=60.0,
                         h_rk4=0.1, h_euler=None, h0_euler=0.1, max_halvings=20):
    """
    Compare Euler and RK4 integration methods for Rössler system.
    
    Integrates with both methods and returns trajectories plus error statistics.
    If h_euler is not specified, finds optimal step size automatically.
    
    Parameters
    ----------
    a, b, c : float
        System parameters
    initial_condition : tuple of float
        Initial state (x0, y0, z0)
    t_min, t_max : float
        Integration time interval
    h_rk4 : float
        RK4 step size (reference)
    h_euler : float, optional
        Euler step size. If None, determined automatically.
    h0_euler : float
        Initial step size for automatic search (if h_euler is None)
    max_halvings : int
        Max halving iterations for automatic search (if h_euler is None)
    
    Returns
    -------
    results : dict
        Dictionary containing:
        - 't_ref', 'x_ref', 'y_ref', 'z_ref': RK4 trajectories
        - 't_eu', 'x_eu', 'y_eu', 'z_eu': Euler trajectories
        - 'h_rk4': RK4 step size used
        - 'h_euler': Euler step size used
        - 'err_x', 'err_y', 'err_z': Component errors
        - 'max_err_x', 'max_err_y', 'max_err_z': Maximum component errors
    """
    # RK4 reference solution
    t_ref, x_ref, y_ref, z_ref = rk4_integrate_rossler(a, b, c, initial_condition, 
                                                         t_min, t_max, h_rk4)
    
    # Determine Euler step size
    if h_euler is None:
        h_euler, _ = find_optimal_euler_step(a, b, c, initial_condition, t_min, t_max, 
                                            h_rk4, h0_euler, max_halvings)
    
    # Euler solution
    t_eu, x_eu, y_eu, z_eu = get_euler_vectors(a, b, c, initial_condition, 
                                                t_min, t_max, h_euler)
    
    # Interpolate Euler results on RK4 grid for error estimation
    x_eu_on_ref = np.interp(t_ref, t_eu, x_eu)
    y_eu_on_ref = np.interp(t_ref, t_eu, y_eu)
    z_eu_on_ref = np.interp(t_ref, t_eu, z_eu)
    
    err_x = np.abs(x_eu_on_ref - x_ref)
    err_y = np.abs(y_eu_on_ref - y_ref)
    err_z = np.abs(z_eu_on_ref - z_ref)
    
    return {
        't_ref': t_ref, 'x_ref': x_ref, 'y_ref': y_ref, 'z_ref': z_ref,
        't_eu': t_eu, 'x_eu': x_eu, 'y_eu': y_eu, 'z_eu': z_eu,
        'h_rk4': h_rk4,
        'h_euler': h_euler,
        'err_x': err_x, 'err_y': err_y, 'err_z': err_z,
        'max_err_x': err_x.max(), 'max_err_y': err_y.max(), 'max_err_z': err_z.max()
    }









