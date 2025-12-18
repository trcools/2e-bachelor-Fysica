# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
# ---

"""
Core data and analysis functions for the Rössler system:

- time integration (Euler scheme),
- symbolic equilibria,
- Jacobian and eigenvalue analysis,
- grouping equilibria by type,
- helper data for wings and Z-maxima.
"""

"""
Note:

This module was adapted from the Lorenz system to work with the Rössler system.
The Rössler equations are:
    ẋ = -y - z
    ẏ = x + ay
    ż = b + z(x - c)
"""

# ============================================================
# Imports and global parameters
# ============================================================

# ----- Imports -----
import importlib
import imports
from imports import (np, sp, solve, Symbol, re, im, Matrix, symbols, Eq, plt,Axes3D, solve_ivp, interact, tqdm)

import parameters
from parameters import get_parameters

# ------- parameters ----
a, b, t_min, t_max, dt, initial_condition, initial_conditions = get_parameters()


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
    # 1.2. Rössler RHS (new better/faster method)
    # ========================================

def rossler_rhs(t, state, c, a=a, b=b):
    """
    Right-hand side of the Rössler system.
    This function represents the differential equations of the Rössler system.

    Parameters
    ----------
    t : float
        Time (not used explicitly, but required by solve_ivp).
    state : array-like, shape (3,)
        Current state [x, y, z].
    c : float
        Control parameter c (bifurcation parameter).
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
    return [dx, dy, dz]

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
    return [dx, dy, dz]




    # ========================================
    # 1.3. Euler integration (old and slow)
    # ========================================


def get_vectors_old(c=2.0, initial_condition=initial_condition, t_min=t_min, t_max=t_max, dt=dt):
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
    def get_arrays_old(t=t):
        x,y,z = np.zeros_like(t), np.zeros_like(t), np.zeros_like(t)
        return x,y,z
        
    # ---initial conditions---
    def get_initial_conditions_old(x,y,z,x0,y0,z0):
        x[0],y[0],z[0] = x0,y0,z0
        return x,y,z
    
    # ---Euler integration---
    def get_euler_integration_old(x,y,z,c,dt=dt,a=a, b=b):
        for i in range (1,len(t)):
            x[i] = x[i-1] + dt*(-y[i-1] - z[i-1])
            y[i] = y[i-1] + dt*(x[i-1] + a*y[i-1])
            z[i] = z[i-1] + dt*(b + z[i-1]*(x[i-1] - c))
        return x,y,z

    x,y,z = get_arrays_old()
    x,y,z = get_initial_conditions_old(x,y,z,*initial_condition)
    x,y,z = get_euler_integration_old(x,y,z,c,dt=dt)
    
    return x,y,z



    # ============================================================
    # 1.4. New integrator based on solve_ivp (Runge–Kutta)
    # ============================================================


def get_vectors(c=2.0,
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
        args=(c,),
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

    x, y, z = get_vectors(c, t_min=t_min, t_max=t_max, dt=dt)

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

    x, y, z = get_vectors(c, t_min=t_min, t_max=t_max, dt=dt)

    right = x >= 0
    left  = x < 0

    return x, y, z, left, right


# ============================================================
# 6. Local maxima of Z(t)
# ============================================================


def get_Z_maxima(c=5.7, skip_first=5, t_min=t_min, t_max=t_max, dt=dt):
    """
    Determine the local maxima Z_n of Z(t) for a given c.

    Parameters
    ----------
    c : float
        Parameter c.
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
    _, _, z = get_vectors(c, t_min=t_min, t_max=t_max, dt=dt)


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
