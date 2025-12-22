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
from sympy import Symbol, re, im, Matrix, symbols, Eq
    
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D 

import scipy
from scipy.integrate import solve_ivp

from ipywidgets import interact
    
import tqdm
from tqdm.auto import tqdm

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
    
    initial_conditions = [(1, 1, 1),
                          (2, 1, 0),
                          (-1, 2, 1),]
    return a, b, c, t_min, t_max, dt, initial_condition, initial_conditions

# ------- parameters ----
a, b, c, t_min, t_max, dt, initial_condition, initial_conditions = get_parameters()

# ============================================================
# 1. Time grid and numerical integration 
# ============================================================


    # ========================================
    # 1.1. Time grid
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


    # ========================================
    # 1.2 Core Rössler Dynamics
    # ========================================

def rossler_rhs(t, state, a=0.2, b=0.2, c=5.7):
    """
    Compute the right-hand side of the Rössler system.
    
    Parameters
    ----------
    t : float
        Time (not explicitly used, included for compatibility with solve_ivp)
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

def rossler_chaos(state, t, c, a=a, b=b):
    """
    Right-hand side of the Rössler system, specifically designed for use with `scipy.integrate.odeint`.

    This function is identical to `rossler_rhs`, but is written specifically for use with `odeint`, 
    because `odeint` requires the function signature `f(state, t)` rather than `f(t, state)` as in `solve_ivp`.

    While the functionality is the same,
    the difference lies in the expected parameter order by the numerical solvers.

    Parameters
    ----------
    state : array-like, shape (3,)
        Current state [x, y, z].
    t : float
        Time (not used explicitly, but required by odeint).
    c : float
        Control parameter c.
    a, b : float
        Rössler parameters a and b.

    Returns
    -------
    rhs : list of float
        Time derivatives [dx/dt, dy/dt, dz/dt].
    """
    x, y, z = state
    dx = -y - z
    dy = x + a * y
    dz = b + z * (x - c)
    return np.array([dx, dy, dz])

    # ========================================
    # 1.3. Euler integration
    # ========================================

def get_euler_vectors(c=2.0, initial_condition=initial_condition, t_min=t_min, t_max=t_max, dt=dt):
    """
    Integrate the Rössler system using a simple explicit Euler scheme.

    Parameters
    ----------
    c : float
        Control parameter c in the Rössler system.
    initial_condition : tuple of float
        Initial state (X0, Y0, Z0).
    t_min, t_max : float
        Time interval.
    dt : float
        Time step.

    Returns
    -------
    x, y, z : ndarray
        Arrays containing X(t), Y(t), Z(t) along the trajectory.
    """
    t = get_time_array(t_min=t_min, t_max=t_max , dt=dt)
    # ---Arrays for the results---
    def get_arrays(t=t):
        x,y,z = np.zeros_like(t), np.zeros_like(t), np.zeros_like(t)
        return x,y,z
        
    # ---initial conditions---
    def get_initial_conditions(x,y,z,x0,y0,z0):
        x[0],y[0],z[0] = x0,y0,z0
        return x,y,z
    
    # ---Euler integration---
    def get_euler_integration(x,y,z,c,dt=dt,a=a, b=b):
        for i in range (1,len(t)):
            x[i] = x[i-1] + dt*(-y[i-1] - z[i-1])
            y[i] = y[i-1] + dt*(x[i-1] + a*y[i-1])
            z[i] = z[i-1] + dt*(b + z[i-1]*(x[i-1] - c))
        return x,y,z

    x,y,z = get_arrays()
    x,y,z = get_initial_conditions(x,y,z,*initial_condition)
    x,y,z = get_euler_integration(x,y,z,c,dt=dt)
    
    return x,y,z

    # ============================================================
    # 1.4. Integrator based on solve_ivp (Runge–Kutta)
    # ============================================================

def get_RK_vectors(a=a, b=b, c=c,
                    initial_condition=initial_condition,
                    t_min=t_min, t_max=t_max, dt=dt,
                    method="RK45",
                    rtol=1e-8, atol=1e-10):
    """
    Integrate the Rössler system using scipy.integrate.solve_ivp
    with a higher-order Runge–Kutta method.

    Parameters
    ----------
    c : float
        Control parameter c (bifurcation parameter).
    initial_condition : tuple of float
        Initial state (X0, Y0, Z0).
    t_min, t_max : float
        Time interval.
    dt : float
        Time step used to *sample* the solution (t_eval grid).
        Internal steps are chosen adaptively by the solver.
    method : str
        Integrator used by solve_ivp (e.g. "RK45", "DOP853").
    rtol, atol : float
        Relative and absolute tolerances for the solver.

    Returns
    -------
    x, y, z : ndarray
        Arrays containing X(t), Y(t), Z(t) evaluated on the uniform
        grid defined by get_time_array(t_min, t_max, dt).
    """
    t_eval = get_time_array(t_min=t_min, t_max=t_max, dt=dt)

    sol = solve_ivp(
        fun=rossler_rhs,
        t_span=(t_min, t_max),
        y0=np.array(initial_condition, dtype=float),
        t_eval=t_eval,
        args=(a,b,c),
        method=method,
        rtol=rtol,
        atol=atol)

    if not sol.success:
        raise RuntimeError(f"Rössler integration failed: {sol.message}")

    x, y, z = sol.y  # shape (3, len(t_eval))
    return x, y, z

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
    
    sols = solve(eqs, (x, y, z), dict=True)

    return sols, (x, y, z)



def show_solutions():
    """
    Print all symbolic equilibria in a human-readable format.
    """
    sols, (x, y, z) = get_solutions()
    for i,sol in enumerate(sols):
            print(f"Solution {i+1}:\n X = {sol[x]} Y = {sol[y]} Z = {sol[z]} \n")

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
    real_eigvals = np.real(eigvals)
    imag_eigvals = np.imag(eigvals)
    has_complex = np.any(np.abs(imag_eigvals) > tol)
    
    # Count the number 'n' of positive, negative and zero eigenvalues
    n_pos = np.sum(real_eigvals >  tol)
    n_neg = np.sum(real_eigvals < -tol)
    n_zero = len(eigvals) - n_pos - n_neg
    
    # 1) All real parts < 0  -> stable
    if n_pos == 0 and n_zero == 0:
        if has_complex:
            eq_type = "stable spiral"
        else:
            eq_type = "stable node"
    
    # 2) At least one > 0, none negative -> purely unstable
    elif n_pos > 0 and n_neg == 0 and n_zero == 0:
        if has_complex:
            eq_type = "unstable spiral"
        else:
            eq_type = "unstable node"
    
    # 3) Both positive and negative real parts -> saddle / saddle-focus
    elif n_pos > 0 and n_neg > 0:
        if has_complex:
            eq_type = "saddle spiral"
        else:
            eq_type = "saddle point"
    
    # 4) At least one eigenvalue ≈ 0, all real -> pitchfork
    elif n_pos == 0 and n_zero > 0 and not has_complex:
        eq_type = "pitchfork bifurcation"
    
    # 5) At least one eigenvalue ≈ 0, complex pair present -> Hopf
    elif n_pos == 0 and n_zero > 0 and has_complex:
        eq_type = "Hopf bifurcation"
    
    # 6) Exotic / other cases
    else:
        eq_type = "other equilibrium"
    
    return eq_type

    

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
        
    
    base_styles = {"stable node": dict(marker="o", s=10, label="stable node"),
                   "stable spiral": dict(marker="o", s=10, label="stable spiral"),
                   "unstable node": dict(marker="x", s=10, label="unstable node"),
                   "unstable spiral": dict(marker="x", s=10, label="unstable spiral"),
                   "saddle point": dict(marker="o", s=10, label="saddle point"),
                   "saddle spiral": dict(marker="o", s=10, label="saddle-spiral"),
                   "pitchfork bifurcation": dict(marker="D", s=50, facecolors="blue", edgecolors="blue",label="pitchfork bifurcation"),
                   "Hopf bifurcation": dict(marker="s", s=50, facecolors="purple", edgecolors="purple", label="Hopf bifurcation"),
                   "other equilibrium": dict(marker="^", color="gray", s=20, label="other / unspecified"),
                  }
    
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

# ============================================================
# 4. Jacobian, equilibria and eigenvalues as function of c
# ============================================================
        

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

    #Only keep real solutions
    if not (x_num.is_real and y_num.is_real and z_num.is_real):
        return None
                
    xsol = float(x_num)
    ysol = float(y_num)
    zsol = float(z_num)

    # Jacobian at equilibrium for Rössler system
    J_eq = np.array([
        [0.0,        -1.0,      -1.0],
        [1.0,        a_val,      0.0],
        [zsol,       0.0,    xsol - c_val]
    ], dtype=float)

    eigvals, eigvecs = np.linalg.eig(J_eq)
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

    

def get_eigvals_for_c(sols, x, y, z, c_sym, a_sym, b_sym, a_val=a, b_val=b, dt=dt, small_c=False):
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
    dt : float
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
        c_values = np.arange(0, 2, dt)
    else:
        c_values = np.arange(2, 18, dt)
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


# ============================================================
# 5. Wing classification and derived data
# ============================================================


def ZY_data_and_wings(c, t_min, t_max, dt=dt):
    """
    Helper: calculate Y(t), Z(t) and left/right-wing masks.

    Returned:
        y, z, left_mask, right_mask
    """
    if t_min >= t_max:
        raise ValueError("t_min must be smaller than t_max")

    x, y, z = get_RK_vectors(c, t_min=t_min, t_max=t_max, dt=dt)

    right = x >= 0
    left  = x < 0

    return y, z, left, right


def XYZ_data_and_wings(c, t_min, t_max, dt=dt):
    """
    Helper: calculate X(t), Y(t), Z(t) and left/right-wing masks.

    Returned:
        x, y, z, left_mask, right_mask
    """
    if t_min >= t_max:
        raise ValueError("t_min must be smaller than t_max")

    x, y, z = get_RK_vectors(c, t_min=t_min, t_max=t_max, dt=dt)

    right = x >= 0
    left  = x < 0

    return x, y, z, left, right


# ============================================================
# 6. Local maxima of Z(t)
# ============================================================


def get_Z_maxima(a, b, c, initial_condition, skip_first=5, t_min=0, t_max=100, dt=0.01):
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

    Returns
    -------
    t_max : ndarray
        Times at which local maxima occur (after discarding the first ones).
    Z_max : ndarray
        Values of Z at those maxima.
    """
    t = get_time_array(t_min, t_max, dt)
    _, _, z = get_RK_vectors(a, b, c, initial_condition=initial_condition, 
                             t_min=t_min, t_max=t_max, dt=dt)


    t_max = []
    Z_max = []

    # How we find the maxima: z[i-1] < z[i] > z[i+1]
    for i in range(1, len(z) - 1):
        if z[i] > z[i-1] and z[i] > z[i+1]:
            t_max.append(t[i])
            Z_max.append(z[i])

    Z_max = np.array(Z_max)
    t_max = np.array(t_max)

    # Throw the first maxima away ('skip_first').
    if len(Z_max) > skip_first:
        Z_max = Z_max[skip_first:]
        t_max = t_max[skip_first:]

    return t_max, Z_max




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



# ============================================================
# Figures
# ============================================================



# ====================================
# Visualization: 2D Projections
# ====================================


def plot_2d_projection(ax, t, x, y, z, plane='XY',
                       initial_condition=initial_condition, title=None,
                       alpha=0.7, linewidth=1.0, color='blue'):
    """
    Plot a 2D projection of the Rössler trajectory.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object
    t, x, y, z : ndarray
        Time and solution components
    plane : str
        Projection plane: 'XY', 'XZ', or 'YZ'
    title : str, optional
        Subplot title
    alpha : float
        Line transparency
    linewidth : float
        Line width
    color : str
        Line color
    """
    if plane == 'XY':
        ax.plot(x, y, color=color, alpha=alpha, linewidth=linewidth)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
    elif plane == 'XZ':
        ax.plot(x, z, color=color, alpha=alpha, linewidth=linewidth)
        ax.set_xlabel('X')
        ax.set_ylabel('Z')
    elif plane == 'YZ':
        ax.plot(y, z, color=color, alpha=alpha, linewidth=linewidth)
        ax.set_xlabel('Y')
        ax.set_ylabel('Z')
    
    if title:
        ax.set_title(title)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)


def plot_3d_trajectory(ax, t, x, y, z, initial_condition=initial_condition,
                      title=None, alpha=0.8, linewidth=1.0, color='blue',
                      skip_first_frac=0.2):
    """
    Plot a 3D trajectory of the Rössler system.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes3DSubplot
        3D axes object
    t, x, y, z : ndarray
        Time and solution components
    title : str, optional
        Plot title
    alpha : float
        Line transparency
    linewidth : float
        Line width
    color : str
        Line color
    skip_first_frac : float
        Fraction of transient data to skip (0.0 to 1.0)
    """
    # Skip transient behavior
    skip_idx = int(len(t) * skip_first_frac)
    
    ax.plot(x[skip_idx:], y[skip_idx:], z[skip_idx:],
            color=color, alpha=alpha, linewidth=linewidth,
            label=f"Trajectory IC={initial_condition}")
    
    if title:
        ax.set_title(title)

    # --- Show equilibrium ---
    sols, (x_sym, y_sym, z_sym) = get_solutions()
    all_eq = compute_jacobian_and_stability(sols, x_sym, y_sym, z_sym,
                                            small_c=(c <= 2))
    eq_needed = [eq for eq in all_eq if np.isclose(eq["c"], c)]

    # How we can recognize the different kinds of equilibrium:
    eq_styles = {"stable node": dict(marker="o", color="C2", label="stable node"),
                   "stable spiral": dict(marker="o", color="darkgreen", label="stable-spiral"),
                   "unstable node": dict(marker="x", color="red", label="unstable node"),
                   "unstable spiral": dict(marker="x", color="black", label="unstable-spiral"),
                   "saddle point": dict(marker="o", color="orange", label="saddle"),
                   "saddle spiral": dict(marker="o", color="green", label="saddle-spiral"),
                   "pitchfork bifurcation": dict(marker="o", color="red", label="Pitchfork bifurcation"),
                   "Hopf bifurcation": dict(marker="o", color="blue", label="Hopf bifurcation"),
                  }
    
    used_labels = set() # To not get the same legend twice.
    
    for eq in eq_needed:
        x_eq = eq["x_eq"]
        y_eq = eq["y_eq"]
        z_eq = eq["z_eq"]
        eq_type = eq["eq_type"]
        
        style = eq_styles.get(eq_type, dict(marker="o", color="gray", label=eq_type))
        
        label = style["label"]

        if label in used_labels:
            plot_label = "_nolegend_"
        else:
            plot_label = label
            used_labels.add(label)
            
        ax.scatter(x_eq, y_eq, z_eq, s=40, marker=style["marker"],
            color=style["color"], label=plot_label)

   # --- General plot settings --- 
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend(loc='best',fontsize=6)

# ===================================
# Visualization: Multiple Panels
# ===================================

def plot_time_series(a=0.2, b=0.2, c=5.7, initial_condition=initial_condition,
                    t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                    rtol=1e-8, atol=1e-10, figsize=(14, 8)):
    """
    Plot time series of x(t), y(t), z(t).
    
    Parameters
    ----------
    x0, y0, z0 : float
        Initial conditions
    a, b, c : float
        System parameters
    t_span : tuple
        Time interval
    figsize : tuple
        Figure size
    """
    x, y, z = get_RK_vectors(a, b, c, initial_condition,
                                t_min, t_max, dt,
                                method,rtol, atol)
    t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)
    
    axes[0].plot(t, x, 'b-', linewidth=0.8, alpha=0.8)
    axes[0].set_ylabel('x(t)')
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title('Time Series: Rössler System')
    
    axes[1].plot(t, y, 'g-', linewidth=0.8, alpha=0.8)
    axes[1].set_ylabel('y(t)')
    axes[1].grid(True, alpha=0.3)
    
    axes[2].plot(t, z, 'r-', linewidth=0.8, alpha=0.8)
    axes[2].set_ylabel('z(t)')
    axes[2].set_xlabel('time (t)')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


def plot_projections_2D(a=0.2, b=0.2, c=5.7, initial_condition=initial_condition,
                     t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                     rtol=1e-8, atol=1e-10, figsize=(15, 5)):
    """
    Plot all three 2D projections of a trajectory.
    
    Parameters
    ----------
    x0, y0, z0 : float
        Initial conditions
    a, b, c : float
        System parameters
    t_span : tuple
        Time interval
    figsize : tuple
        Figure size
    """
    x, y, z = get_RK_vectors(a, b, c, initial_condition,
                                t_min, t_max, dt,
                                method,rtol, atol)
    t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # Skip transient
    skip = int(len(t) * 0.2)
    
    plot_2d_projection(axes[0], t[skip:], x[skip:], y[skip:], z[skip:],
                      plane='XY', title='x-y plane (top view)')
    plot_2d_projection(axes[1], t[skip:], x[skip:], y[skip:], z[skip:],
                      plane='XZ', title='x-z plane (side view)')
    plot_2d_projection(axes[2], t[skip:], x[skip:], y[skip:], z[skip:],
                      plane='YZ', title='y-z plane')
    
    plt.tight_layout()
    return fig


def plot_3d_attractor(a=0.2, b=0.2, c=5.7, initial_condition=initial_condition,
                    t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                    rtol=1e-8, atol=1e-10, figsize=(12, 5)):
    """
    Plot 3D trajectory and a 2D projection.
    
    Parameters
    ----------
    x0, y0, z0 : float
        Initial conditions
    a, b, c : float
        System parameters
    t_span : tuple
        Time interval
    figsize : tuple
        Figure size
    """

    't, x, y, z = integrate_trajectory_dense(x0, y0, z0, a, b, c, t_span)'
    x ,y, z = get_RK_vectors(a, b, c, initial_condition,
                                t_min, t_max, dt,
                                method,rtol, atol)
    t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
    skip = int(len(t) * 0.2)
    
    fig = plt.figure(figsize=figsize)
    
    # 3D plot
    ax1 = fig.add_subplot(121, projection='3d')
    plot_3d_trajectory(ax1, t, x, y, z, title='3D Rössler Attractor',
                      skip_first_frac=0.2, color='darkblue', linewidth=0.5)
    
    # 2D projection
    ax2 = fig.add_subplot(122)
    plot_2d_projection(ax2, t[skip:], x[skip:], y[skip:], z[skip:],
                      plane='XY', title='x-y Projection')
    
    plt.tight_layout()
    return fig


# ============================================================
# Parameter Study
# ============================================================

def plot_parameter_sweep(a=0.1, b=0.1, c_values=None,
                        initial_condition=initial_condition,
                        t_min=0, t_max=200, dt=0.01,
                        method='RK45', rtol=1e-8, atol=1e-10,
                        figsize=(15, 10)):
    """
    Show how the attractor changes as parameter 'a' varies.
    
    Parameters
    ----------
    a_values : list, optional
        List of 'a' parameter values to explore
    b, c : float
        Other system parameters
    x0, y0, z0 : float
        Initial conditions
    t_span : tuple
        Time interval
    figsize : tuple
        Figure size
    """
    if c_values is None:
        c_values = [0.1, 0.15, 0.2, 0.3, 0.4]
    
    n_c = len(c_values)
    n_cols = 3
    n_rows = (n_c + n_cols - 1) // n_cols
    
    fig = plt.figure(figsize=figsize)
    
    for idx, c in enumerate(c_values, 1):
        ax = fig.add_subplot(n_rows, n_cols, idx, projection='3d')
        
        x, y, z = get_RK_vectors(a, b, c, initial_condition=initial_condition,
                                t_min=t_min, t_max=t_max, dt=dt,
                                method=method, rtol=rtol, atol=atol)
        t = get_time_array(t_min=t_min, t_max=t_max,dt=dt)
        plot_3d_trajectory(ax, t, x, y, z,
                          title=f'a = {a}, b = {b}, c = {c}',
                          skip_first_frac=0.2, color='darkblue', linewidth=0.3)
    
    plt.tight_layout()
    return fig


def bifurcation_analysis_parameter_c(a=0.1, b=0.1, c_values=None,
                        initial_condition=initial_condition,
                        t_min=0, t_max=200, dt=0.01,
                        method='RK45', rtol=1e-8, atol=1e-10,
                        skip_frac=0.5):
    """
    Create a bifurcation diagram: x(t) values versus parameter 'c'.
    """
    if c_values is None:
        c_values = np.linspace(0.05, 0.5, 300)
    
    x_plot = []
    c_plot = []
    
    for c in c_values:
        # FAST version: fewer points, looser tolerances
        x, y, z = get_RK_vectors(a, b, c, initial_condition=initial_condition,
                                t_min=t_min, t_max=t_max, dt=dt,
                                method=method, rtol=rtol, atol=atol)
        t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
        skip = int(len(t) * skip_frac)
        
        # Collect only local maxima of x (characteristic points)
        x_vals = x[skip:]
        local_max_indices = np.where((x_vals[1:-1] > x_vals[:-2]) & 
                                      (x_vals[1:-1] > x_vals[2:]))[0] + 1
        
        if len(local_max_indices) > 0:
            for idx in local_max_indices:
                x_plot.append(x_vals[idx])
                c_plot.append(c)
    
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.scatter(c_plot, x_plot, s=1, alpha=0.5, color='darkblue')
    ax.set_xlabel('Parameter c', fontsize=12)
    ax.set_ylabel('x (local maxima)', fontsize=12)
    ax.set_title('Bifurcation Diagram: Rössler System (varying c)', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    return fig, np.array(c_plot), np.array(x_plot)


# ============================================================
# Comparison and Analysis
# ============================================================

def compare_initial_conditions(a=0.2, b=0.2, c=5.7, initial_conditions= None,
                    t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                    rtol=1e-8, atol=1e-10, figsize=(15, 5)):
    """
    Compare trajectories from different initial conditions.
    
    Parameters
    ----------
    a, b, c : float
        System parameters
    initial_conditions : list of tuples, optional
        List of (x0, y0, z0) tuples
    t_span : tuple
        Time interval
    figsize : tuple
        Figure size
    """
    if initial_conditions is None:
        initial_conditions = [
            (0.1, 0.1, 0.1),
            (1.0, 1.0, 1.0),
            (-1.0, 0.5, 0.5)
        ]
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    
    for ic_idx, (x0, y0, z0) in enumerate(initial_conditions):
        x, y, z = get_RK_vectors(a, b, c, initial_condition=(x0, y0, z0),
                                t_min=t_min, t_max=t_max, dt=dt,
                                method=method, rtol=rtol, atol=atol)
        t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
        skip = int(len(t) * 0.2)

        color = colors[ic_idx % len(colors)]
        label = f'IC: ({x0}, {y0}, {z0})'
        
        plot_2d_projection(axes[0], t[skip:], x[skip:], y[skip:], z[skip:],
                          plane='XY', color=color, title='x-y Plane')
        plot_2d_projection(axes[1], t[skip:], x[skip:], y[skip:], z[skip:],
                          plane='XZ', color=color, title='x-z Plane')
        plot_2d_projection(axes[2], t[skip:], x[skip:], y[skip:], z[skip:],
                          plane='YZ', color=color, title='y-z Plane')
        
        axes[0].plot([], [], color=color, label=label, linewidth=2)
    
    for ax in axes:
        ax.legend(loc='best', fontsize=9)
    
    plt.tight_layout()
    return fig


def lyapunov_exponent_estimate(a=0.2, b=0.2, c=5.7, initial_condition=initial_condition,
                               t_min=0, t_max=100, dt=dt,
                               method="RK45", rtol=1e-8,
                               atol=1e-10, eps=1e-8):
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
    x1, y1, z1 = get_RK_vectors(
        a=0.2, b=0.2, c=5.7,
        initial_condition=initial_condition,
        t_min=t_min, t_max=t_max, dt=dt,
        method="RK45",rtol=1e-8, atol=1e-10
        )
    
    # Perturbed trajectory
    x2, y2, z2 = get_RK_vectors(
        a=0.2, b=0.2, c=5.7,
        initial_condition=(initial_condition[0] + eps,
                           initial_condition[1],
                           initial_condition[2]),
        t_min=t_min, t_max=t_max, dt=dt,
        method="RK45",rtol=1e-8, atol=1e-10
        )
    
    t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
    
    # Distance between trajectories
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    distance = np.sqrt(dx**2 + dy**2 + dz**2)
    
    # Avoid log(0)
    distance = np.maximum(distance, 1e-12)
    
    # Lyapunov exponent estimate: ln(distance(t)) / t
    lyap_exp = np.log(distance) / (t + 1e-10)
    
    return t, lyap_exp


# ============================================================
#  Analysis of Chaotic Characteristics
# ============================================================

    # ====================================================
    # Sensitivity to initial conditions (Lyapunov flavour)
    # ====================================================

def plot_sensitivity_to_initial_conditions(c=5.7, t_min=t_min, t_max=t_max, dt=dt, delta=1e-6):
    """
    Show sensitivity to initial conditions in the Rössler system.

    For a given c, we take two initial conditions that differ by only `delta`
    in the X-component, integrate both trajectories, and plot how their
    distance in state space evolves in time (linear and log scale).
    """
    t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)

    # Choose reference initial condition from the general parameters
    base_ic = initial_condition
    ic1 = base_ic
    ic2 = (base_ic[0] + delta, base_ic[1], base_ic[2])

    x1, y1, z1 = get_RK_vectors(a=a, b=b, c=c, t_min=t_min, t_max=t_max, dt=dt,
                                initial_condition=ic1)
    x2, y2, z2 = get_RK_vectors(a=a, b=b, c=c, t_min=t_min, t_max=t_max, dt=dt,
                                initial_condition=ic2)
    # Euclidean distance in (X, Y, Z)
    d = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

    plt.close('all')
    fig, axes = plt.subplots(2, 1, figsize=(8, 8))

    # Upper: distance on linear scale
    ax = axes[0]
    ax.plot(t, d)
    ax.set_xlabel('t')
    ax.set_ylabel('distance d(t)')
    ax.set_title(rf'Distance between two nearby trajectories, $c={c:.2f}$')
    ax.grid(True)

    # Lower: log-scale to highlight exponential growth
    ax = axes[1]
    ax.semilogy(t, d)
    ax.set_xlabel('t')
    ax.set_ylabel('d(t) log scale')
    ax.set_title('Approximate exponential divergence (semilogy)')
    ax.grid(True)

    plt.show()

def show_butterfly_effect(a, b, c, initial_condition, delta=1e-6, t_min=0, t_max=100, dt=0.01):
    """
    Illustrate the butterfly effect in the Rössler system.
    
    We integrate two trajectories with almost identical
    initial conditions and compare X(t), Y(t), Z(t) in time.
    """
    plt.close('all')
    
    # Time array
    t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
    
    # Two initial conditions that differ by delta in X
    base_ic = initial_condition
    ic1 = base_ic
    ic2 = (base_ic[0] + delta, base_ic[1], base_ic[2])
    
    # Integrate both trajectories
    x1, y1, z1 = get_RK_vectors(a=a, b=b, c=c, t_min=t_min, t_max=t_max, dt=dt,
                                initial_condition=ic1)
    x2, y2, z2 = get_RK_vectors(a=a, b=b, c=c, t_min=t_min, t_max=t_max, dt=dt,
                                initial_condition=ic2)
    
    series1 = {"X": x1, "Y": y1, "Z": z1}
    series2 = {"X": x2, "Y": y2, "Z": z2}
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), constrained_layout=True)
    
    for ax, comp in zip(axes, ("X", "Y", "Z")):
        ax.plot(t, series1[comp],
                label=rf'{comp}(t), initial condition',
                alpha=0.8)
        ax.plot(t, series2[comp],
                label=rf'{comp}(t), initial condition + $\delta$',
                alpha=0.8, linestyle='--')
    
        ax.set_ylabel(rf'{comp}(t)')
        ax.grid(True)
        ax.legend(loc='best', fontsize=8)
    
    axes[-1].set_xlabel('t')
    fig.suptitle(rf'Butterfly effect in Rössler system for $c={c:.2f}$, $\delta={delta:.1e}$')
    plt.show()  


    # =================================
    #   Return map: maxima of Z(t)
    # =================================

def plot_Z_time_with_maxima(a, b, c, initial_condition, t_min=0, t_max=200, dt=0.01):
    """
    Plot Z(t) over time and mark the detected local maxima Z_n
    for a given c and time interval [t_min, t_max].

    Parameters
    ----------
    a, b, c : float
        Rössler parameters.
    initial_condition : tuple
        Initial condition (x0, y0, z0).
    t_min, t_max : float
        Time interval.
    dt : float
        Time step used in the integration.
    """
    t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
    _, _, z = get_RK_vectors(a=a, b=b, c=c, t_min=t_min, t_max=t_max, dt=dt,
                                initial_condition=initial_condition)
    t_max_vals, Z_max = get_Z_maxima(a, b, c, initial_condition, t_min=t_min, t_max=t_max, dt=dt)

    plt.close('all')
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(t, z, label='Z(t)', linewidth=0.7)
    ax.scatter(t_max_vals, Z_max, color='red', s=20, label='local maxima $Z_n$')
    ax.set_xlabel('t')
    ax.set_ylabel('Z(t)')
    ax.set_title(rf'Local maxima $Z_n$ of $Z(t)$ for $c={c:.2f}$')
    ax.legend(loc='best',fontsize=8)
    ax.grid(True)

    plt.show()


def plot_Zn_vs_Znplus1(a, b, c, initial_condition, skip_first=20, t_min=0, t_max=200, dt=0.01):
    """
    Plot the return map Z_{n+1} versus Z_n, where Z_n are the local maxima of Z(t).

    This visualises how each maximum Z_n
    (almost) uniquely determines the next one Z_{n+1} for a given c.
    """
    # Extract the local maxima Z_n of Z(t)
    _, Z_max = get_Z_maxima(a=a, b=b, c=c, initial_condition=initial_condition,
                            skip_first=skip_first, t_min=t_min, t_max=t_max, dt=dt)

    # Getting couples (Z_n, Z_{n+1})
    Z_n = Z_max[:-1]
    Z_n1 = Z_max[1:]

    plt.close('all')
    fig, ax = plt.subplots(figsize=(6, 5))

    # Check if Z_n and Z_n1 are not empty
    if Z_n.size > 0 and Z_n1.size > 0:
        z_min = min(Z_n.min(), Z_n1.min())
        z_max = max(Z_n.max(), Z_n1.max())
        ax.scatter(Z_n, Z_n1, s=10, alpha=0.7, label=r"$(Z_n, Z_{n+1})$")

        # Identity line Z_{n+1} = Z_n to guide (45°)
        ax.plot([z_min, z_max], [z_min, z_max], 'k--', linewidth=0.8, label='45° slope')

    ax.set_xlabel(r'Max $Z_n$')
    ax.set_ylabel(r'Max $Z_{n+1}$')
    ax.set_title(rf'Return map for max $Z_n$ at $c = {c:.2f}$')
    ax.grid(True)
    ax.legend(loc='best', fontsize=8)
    plt.show()

def cobweb_plot(f, x0, n_iter: int = 30, x_min: float = 0.0, x_max: float = 1.0, ax=None,
                func_label: str = r"$x_{n+1} = f(x_n)$", cob_color: str = "tab:red",
                cob_lw: float = 0.8, mark_start: bool = True):
    """
    Draw a cobweb diagram for a one-dimensional map x_{n+1} = f(x_n).

    Parameters
    ----------
    f : callable
        One-dimensional map f(x).
    x0 : float
        Initial value x_0.
    n_iter : int, optional
        Number of cobweb iterations to draw.
    x_min, x_max : float, optional
        Range for the axes and for plotting the function.
    ax : matplotlib.axes.Axes or None, optional
        Axis to plot on. If None, a new figure and axis are created.
    func_label : str, optional
        Label for the curve y = f(x).
    cob_color : str, optional
        Color of the cobweb lines.
    cob_lw : float, optional
        Line width of the cobweb lines.
    mark_start : bool, optional
        If True, mark the starting point (x0, f(x0)) on the curve.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))

    # Plot f(x) and the diagonal
    xs = np.linspace(x_min, x_max, 400)
    ax.plot(xs, f(xs), color="C0", lw=2.0, label=func_label)
    ax.plot(xs, xs, linestyle="--", color="0.5", lw=1.0, label=r"$x_{n+1} = x_n$")

    # Mark the starting point on the curve: (x0, f(x0))
    if mark_start:
        y0 = f(x0)
        ax.scatter([x0], [y0], s=20, color="k", 
                   zorder=4, label=r"start $(x_0, f(x_0))$")

    # Cobweb iterations
    x = x0
    for _ in range(n_iter):
        x_next = f(x)

        # vertical: (x_n, x_n) -> (x_n, x_{n+1})
        ax.plot([x, x], [x, x_next], color=cob_color, lw=cob_lw)

        # horizontal: (x_n, x_{n+1}) -> (x_{n+1}, x_{n+1})
        ax.plot([x, x_next], [x_next, x_next], color=cob_color, lw=cob_lw)

        x = x_next

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(x_min, x_max)
    ax.set_xlabel(r"$x_n$")
    ax.set_ylabel(r"$x_{n+1}$")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)
    ax.legend(loc="best", fontsize=8)

    return ax


def show_rossler_return_cobweb(Z_n, c: float = 5.7, n_iter: int = 30, skip_first: int = 50):
    """
    Plot the Rössler return map Z_{n+1} vs Z_n with an overlaid cobweb diagram.

    Parameters
    ----------
    Z_n : array-like
        Sequence of local maxima of z(t).
    c : float
        Parameter c for the Rössler system (for title/context only).
    n_iter : int
        Number of cobweb iterations to draw.
    skip_first : int
        Number of initial maxima to skip as transients when constructing the map.

    Returns
    -------
    ax : matplotlib.axes.Axes
        Axes with the return map scatter and cobweb overlay.
    """
    Z_n = np.asarray(Z_n).astype(float)
    if Z_n.size < 2:
        raise ValueError("Z_n must contain at least two maxima to form a return map")

    if skip_first < 0:
        skip_first = 0
    if skip_first >= len(Z_n) - 1:
        skip_first = max(0, len(Z_n) - 2)

    # Build the empirical return map from data after skipping transients
    Z_use = Z_n[skip_first:]
    F = rossler_return_map_from_data(Z_use)

    x_data = Z_use[:-1]
    y_data = Z_use[1:]

    # Plot settings
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x_data, y_data, s=12, alpha=0.7, color="tab:blue", label=r"$(Z_n, Z_{n+1})$")

    # Diagonal
    x_min = float(np.min(x_data))
    x_max = float(np.max(x_data))
    pad = 0.05 * (x_max - x_min if x_max > x_min else 1.0)
    x_min -= pad
    x_max += pad
    xs = np.linspace(x_min, x_max, 400)
    ax.plot(xs, xs, color="k", linestyle="--", linewidth=0.7, label="y = x")

    # Cobweb from first usable point for more representative dynamics
    x0 = 7
    cobweb_plot(F, x0=x0, n_iter=n_iter, x_min=x_min, x_max=x_max, ax=ax,
                func_label=r"$Z_{n+1} = F(Z_n)$", cob_color="tab:red", cob_lw=0.8, mark_start=True)

    ax.set_xlabel(r"$Z_n$")
    ax.set_ylabel(r"$Z_{n+1}$")
    ax.set_title(fr"Rössler return map with cobweb (c={c:.3f})")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    return ax


# ============================================================
# Notebook Convenience Wrappers (produce same results as notebook)
# ============================================================

def _with_defaults(a=None, b=None, c=None, initial_condition=None, t_min=None, t_max=None, dt=None):
    """
    Resolve defaults from get_parameters() while allowing overrides.
    Returns (a, b, c, initial_condition, t_min, t_max, dt).
    """
    a0, b0, c0, tmin0, tmax0, dt0, ic0, _ = get_parameters()
    return (
        a if a is not None else a0,
        b if b is not None else b0,
        c if c is not None else c0,
        initial_condition if initial_condition is not None else ic0,
        t_min if t_min is not None else tmin0,
        t_max if t_max is not None else tmax0,
        dt if dt is not None else dt0,
    )


def nb_print_parameters():
    """Print the default parameters exactly as in the notebook."""
    a, b, c, t_min, t_max, dt, initial_condition, initial_conditions = get_parameters()
    print(f"Parameters: a={a}, b={b}, c={c}")
    print(f"Time: t_min={t_min}, t_max={t_max}, dt={dt}")
    return a, b, c, t_min, t_max, dt, initial_condition, initial_conditions


def nb_show_attractor(a=None, b=None, c=None, initial_condition=None, t_min=None, t_max=None, dt=None):
    """Show 3D attractor (identical to notebook cell)."""
    a, b, c, ic, t_min, t_max, dt = _with_defaults(a, b, c, initial_condition, t_min, t_max, dt)
    plot_3d_attractor(a=a, b=b, c=c, initial_condition=ic,
                      t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                      rtol=1e-8, atol=1e-10, figsize=(12, 5))
    plt.show()


def nb_show_projections(a=None, b=None, c=None, initial_condition=None, t_min=None, t_max=None, dt=None):
    """Show all three 2D projections (identical to notebook cell)."""
    a, b, c, ic, t_min, t_max, dt = _with_defaults(a, b, c, initial_condition, t_min, t_max, dt)
    plot_projections_2D(a=a, b=b, c=c, initial_condition=ic,
                        t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                        rtol=1e-8, atol=1e-10, figsize=(15, 5))
    plt.show()


def nb_time_series(a=None, b=None, c=None, initial_condition=None, t_min=None, t_max=None, dt=None):
    """Show x(t), y(t), z(t) time series."""
    a, b, c, ic, t_min, t_max, dt = _with_defaults(a, b, c, initial_condition, t_min, t_max, dt)
    plot_time_series(a=a, b=b, c=c, initial_condition=ic,
                     t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                     rtol=1e-8, atol=1e-10, figsize=(14, 8))
    plt.show()


def nb_compare_initial_conditions(a=None, b=None, c=None, initial_conditions=None,
                                  t_min=None, t_max=None, dt=None):
    """Compare trajectories from several initial conditions."""
    a, b, c, ic, t_min, t_max, dt = _with_defaults(a, b, c, None, t_min, t_max, dt)
    if initial_conditions is None:
        # Use defaults list from get_parameters
        _, _, _, _, _, _, _, initial_conditions = get_parameters()
    compare_initial_conditions(a=a, b=b, c=c, initial_conditions=initial_conditions,
                               t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                               rtol=1e-8, atol=1e-10, figsize=(15, 5))
    plt.show()


def nb_return_map_and_cobweb(a=None, b=None, c=5.7, initial_condition=None,
                             t_min=None, t_max=None, dt=None,
                             skip_first_maxima=5, cobweb_skip_first=50, n_iter=15):
    """Plot Z(t) with maxima, the return map, and a cobweb plot (same as notebook trio)."""
    a, b, c, ic, t_min, t_max, dt = _with_defaults(a, b, c, initial_condition, t_min, t_max, dt)
    plot_Z_time_with_maxima(a=a, b=b, c=c, initial_condition=ic, t_min=t_min, t_max=200, dt=dt)
    plot_Zn_vs_Znplus1(a=a, b=b, c=c, initial_condition=ic, t_min=0, t_max=2000)
    _, Z_maxima = get_Z_maxima(a=a, b=b, c=c, initial_condition=ic,
                                skip_first=skip_first_maxima, t_min=0, t_max=2000, dt=dt)
    show_rossler_return_cobweb(Z_maxima, c=c, n_iter=n_iter, skip_first=cobweb_skip_first)


def nb_butterfly_effect(a=None, b=None, c=5.7, initial_condition=None, delta=1e-6, t_min=50, t_max=280, dt=None):
    """Show butterfly effect plot (two nearby initial conditions)."""
    a, b, c, ic, _, _, dt = _with_defaults(a, b, c, initial_condition, None, None, dt)
    show_butterfly_effect(a=a, b=b, c=c, initial_condition=ic, delta=delta, t_min=t_min, t_max=t_max, dt=dt)


def nb_plot_sensitivity(c=5.7, t_min=None, t_max=None, dt=None, delta=1e-6):
    """Plot distance growth on lin/log scale as sensitivity visualization."""
    _, _, _, _, t_min, t_max, dt = _with_defaults(None, None, None, None, t_min, t_max, dt)
    plot_sensitivity_to_initial_conditions(c=c, t_min=t_min, t_max=t_max, dt=dt, delta=delta)


def nb_parameter_sweep(c_values=None, a=0.1, b=0.1, initial_condition=None,
                       t_min=None, t_max=200, dt=None):
    """Show attractor panels for a list of c values."""
    a, b, _, ic, t_min, _, dt = _with_defaults(a, b, None, initial_condition, t_min, None, dt)
    if c_values is None:
        c_values = [4, 5, 6, 7, 8, 3]
    plot_parameter_sweep(a=a, b=b, c_values=c_values, initial_condition=ic,
                         t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                         rtol=1e-8, atol=1e-10, figsize=(15, 10))
    plt.show()


def nb_bifurcation(c_range=None, a=None, b=None, initial_condition=None,
                   t_min=0, t_max=60, skip_frac=0.6):
    """Produce the bifurcation diagram for varying c (same defaults as notebook)."""
    a, b, _, ic, _, _, _ = _with_defaults(a, b, None, initial_condition, None, None, None)
    if c_range is None:
        c_range = np.linspace(0.05, 45, 100)
    fig, c_vals, x_vals = bifurcation_analysis_parameter_c(
        a=a, b=b, c_values=c_range, initial_condition=ic,
        t_min=t_min, t_max=t_max, skip_frac=skip_frac
    )
    plt.show()
    return fig, c_vals, x_vals