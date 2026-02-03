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
Core data and analysis functions for the Lorenz system:

- time integration (Euler scheme),
- symbolic equilibria,
- Jacobian and eigenvalue analysis,
- grouping equilibria by type,
- helper data for wings and Z-maxima.
"""

"""
Note:

This module was refactored with assistance from an AI tool to improve structure,
naming, and documentation for readability and maintainability.

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
sigma, beta, t_min, t_max, dt, initial_condition, initial_conditions = get_parameters()


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
    t = np.arange(t_min, t_max, dt) # or np.linspace(*tspan,times(=150 as example)).
    return t



    # ========================================
    # 1.2. Lorenz RHS (new better/faster method)
    # ========================================

def lorenz_rhs(t, state, rho, sigma=sigma, beta=beta):
    """
    Right-hand side of the Lorenz system.
    This function represents the differential equations of the Lorenz system.

    Parameters
    ----------
    t : float
        Time (not used explicitly, but required by solve_ivp).
    state : array-like, shape (3,)
        Current state [x, y, z].
    rho : float
        Control parameter ρ.
    sigma, beta : float
        Lorenz parameters σ and β.

    Returns
    -------
    rhs : list of float
        Time derivatives [dx/dt, dy/dt, dz/dt].
    """
    x, y, z = state
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]

def lorenz_chaos(state, t, rho, sigma=sigma, beta=beta):
    """
    Right-hand side of the Lorenz system, specifically designed for use with `scipy.integrate.odeint`.

    This function is identical to `lorenz_rhs`, but is written specifically for use with `odeint`, 
    because `odeint` requires the function signature `f(state, t)` rather than `f(t, state)` as in `solve_ivp`.

    While the functionality is the same,
    the difference lies in the expected parameter order by the numerical solvers.

    Parameters
    ----------
    t : float
        Time (not used explicitly, but required by solve_ivp).
    state : array-like, shape (3,)
        Current state [x, y, z].
    rho : float
        Control parameter ρ.
    sigma, beta : float
        Lorenz parameters σ and β.

    Returns
    -------
    rhs : list of float
        Time derivatives [dx/dt, dy/dt, dz/dt].
    """
    x, y, z = state
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]




    # ========================================
    # 1.3. Euler intergation (old and slow)
    # ========================================


def get_vectors_old(rho=0.5, initial_condition=initial_condition, t_min=t_min, t_max=t_max, dt=dt):
    """
    Integrate the Lorenz system using a simple explicit Euler scheme.

    Parameters
    ----------
    rho : float
        Control parameter ρ in the Lorenz system.
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
        #x,y,z,t = arrays(tspan)
        x[0],y[0],z[0] = x0,y0,z0
        return x,y,z
    
    # ---Euler integration---
    def get_euler_integration_old(x,y,z,rho,dt=dt,sigma=sigma, beta =beta):
        #x,y,z,t = arrays(tspan)  
        for i in range (1,len(t)):
            x[i] = x[i-1] + dt*(sigma*(y[i-1] - x[i-1]))
            y[i] = y[i-1] + dt*(x[i-1]*(rho-z[i-1])-y[i-1])
            z[i] = z[i-1] + dt*(x[i-1]*y[i-1]-beta*z[i-1])
        return x,y,z

    x,y,z = get_arrays_old()
    x,y,z = get_initial_conditions_old(x,y,z,*initial_condition)
    x,y,z = get_euler_integration_old(x,y,z,rho,dt=dt)
    
    return x,y,z #The lowercase in the code has the same meaning of the uppercases in theory for x, y and z.



    # ============================================================
    # 1.4. New integrator based on solve_ivp (Runge–Kutta)
    # ============================================================


def get_vectors(rho=0.5,
                    initial_condition=initial_condition,
                    t_min=t_min, t_max=t_max, dt=dt,
                    method="RK45",
                    rtol=1e-8, atol=1e-10):
    """
    Integrate the Lorenz system using scipy.integrate.solve_ivp
    with a higher-order Runge–Kutta method.

    Parameters
    ----------
    rho : float
        Control parameter ρ.
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
        fun=lorenz_rhs,
        t_span=(t_min, t_max),
        y0=np.array(initial_condition, dtype=float),
        t_eval=t_eval,
        args=(rho,),
        method=method,
        rtol=rtol,
        atol=atol)

    if not sol.success:
        raise RuntimeError(f"Lorenz integration failed: {sol.message}")

    x, y, z = sol.y  # shape (3, len(t_eval))
    return x, y, z




# ============================================================
# 2. Symbolic equilibria of the Lorenz system
# ============================================================


def get_solutions():
    """
    Compute the symbolic equilibria of the Lorenz system using SymPy.

    Returns
    -------
    sols : list[dict]
        SymPy solutions as dicts {x: expr, y: expr, z: expr}.
    (x, y, z) : tuple of SymPy symbols
        The state variables used in the symbolic system.
    """
    x, y, z, rho, sigma, beta = symbols('x y z rho sigma beta ', real=True)

    eqs = [
        Eq(sigma * (y - x), 0),
        Eq(x * (rho - z) - y, 0),
        Eq(x * y - beta * z, 0)]
    
    sols = solve(eqs, (x, y, z), dict=True)

    return sols, (x, y, z) #so that get_solutions only gives back print to avoid kinda duplicate prints.



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
    tol = 1e-9 # To not be too strict.
    real_eigvals = np.real(eigvals)
    imag_eigvals = np.imag(eigvals)
    has_complex = np.any(np.abs(imag_eigvals) > tol)
    
    # Count the number 'n' of positive, negative and zero eigenvalues
    n_pos = np.sum(real_eigvals >  tol)
    n_neg = np.sum(real_eigvals < -tol)
    n_zero = len(eigvals) - n_pos - n_neg # Also possibe to use real_eigvals <  tol & real_eigvals > -tol
    
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
        - "1D": (rho, val)
        - "2D": (rho, v1, v2)

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
                "rho": [],
                "v1": [],
                "v2": [],
                "style": style.copy(),
            }
    else: # dim == "1D"
        for key, style in base_styles.items():
            groups[key] = {
                "rho": [],
                "val": [],
                "style": style.copy(),
            }

    return groups


# ============================================================
# 4. Jacobian, equilibria and eigenvalues as function of rho
# ============================================================
        

def get_equilibria_for_rho(sol, x, y, z, rho_sym, beta_sym, rho_val, sigma=sigma, beta_val=beta):
    """
    Compute the equilibrium for a single symbolic solution and a given rho.

    Only real equilibria are kept.

    Parameters
    ----------
    sol : dict
        SymPy solution {x: expr, y: expr, z: expr}.
    x, y, z : SymPy symbols
        State variables.
    rho_sym, beta_sym : SymPy symbols
        Parameter symbols (rho, beta).
    rho_val : float
        Value at which to evaluate the equilibrium.
    sigma : float
        Lorenz parameter σ.
    beta_val : float
        Lorenz parameter β.

    Returns
    -------
    eq : dict or None
        If real, a dict with keys
        {"rho", "x_eq", "y_eq", "z_eq", "J", "eigvals", "eq_type"},
        otherwise None.
    """

    x_expr = sol[x]
    y_expr = sol[y]
    z_expr = sol[z]
        
    subs_dict = {rho_sym: rho_val, beta_sym: beta_val}

    # Here we substitute rho and beta into the symbolic equations
    x_num = x_expr.subs(subs_dict).evalf()
    y_num = y_expr.subs(subs_dict).evalf()
    z_num = z_expr.subs(subs_dict).evalf()

    #Only keep real solutions
    if not (x_num.is_real and y_num.is_real and z_num.is_real):
        return None
                
    xsol = float(x_num)
    ysol = float(y_num)
    zsol = float(z_num)

    # Jacobean at equilibrium
    J_eq = np.array([
        [-sigma,     sigma,  0.0],
        [rho_val - zsol,    -1.0,  -xsol  ],
        [ysol,          xsol,     -beta_val]
    ], dtype=float)

    eigvals, eigvecs = np.linalg.eig(J_eq)
    eq_type = equilibrium_classification(eigvals)
            
    return {
        "rho": rho_val,
        "x_eq": xsol,
        "y_eq": ysol,
        "z_eq": zsol,
        "J": J_eq,
        "eigvals": eigvals,
        "eq_type": eq_type}
    

def get_equilibria_from_jacobian(sol, x, y, z,rho_sym, beta_sym, rho_values,
                                 sigma=sigma,beta_val=beta):
    """
    Helper function: compute equilibria for many rho values for a single solution.

    Parameters
    ----------
    sol : dict
        SymPy solution {x: expr, y: expr, z: expr}.
    x, y, z : SymPy symbols
        State variables.
    rho_sym, beta_sym : SymPy symbols
        Parameter symbols.
    rho_values : iterable of float
        Values of rho to scan.
    sigma, beta_val : float
        Lorenz parameters.

    Returns
    -------
    equilibrium : list[dict]
        List of equilibrium dictionaries (see get_equilibria_for_rho).
    """
    equilibrium = []
    
    for rho_val in rho_values:
        eq = get_equilibria_for_rho(sol, x, y, z, rho_sym, beta_sym, rho_val,
                                    sigma=sigma, beta_val=beta_val)
        if eq is not None:
                equilibrium.append(eq)
    return equilibrium


def compute_jacobian_and_stability(sols, x, y, z, small_rho=False):
    """
    Compute Jacobian matrices, eigenvalues and stability classification
    for equilibria over a grid of rho values.

    Parameters
    ----------
    sols : list[dict]
        SymPy solutions {x: expr, y: expr, z: expr}.
    x, y, z : SymPy symbols
        State variables.
    small_rho : bool
        If True, use rho in [0, 1]; otherwise use rho in [0, 28].

    Returns
    -------
    all_equilibria : list[dict]
        Each dict contains keys:
        {"rho", "x_eq", "y_eq", "z_eq", "J", "eigvals", "eq_type"}.
    """
    rho_sym, beta_sym = symbols('rho beta', real=True)
    all_equilibria = []
    
    if small_rho: 
        """
        If rho < 1 the only real solution is sol[0] (0,0,0).
        """
        sol = sols[0]
        rho_values = np.linspace(0, 1, 51)
        all_equilibria.extend(get_equilibria_from_jacobian(sol, x, y, z,
                                                           rho_sym, beta_sym,
                                                        rho_values,sigma=sigma,
                                                           beta_val=beta))
        
    else:
        """
        The general case where rho > 1 and has 3 solutions.
        """
        rho_values = np.linspace(0,28,113)
        for sol in sols:
            all_equilibria.extend(get_equilibria_from_jacobian(sol, x, y, z,
                                                               rho_sym,
                                                               beta_sym,
                                                               rho_values,
                                                               sigma=sigma,
                                                               beta_val=beta))
    return all_equilibria

    

def get_eigvals_for_rho(sols, x, y, z, rho_sym, beta_sym, sigma=10, beta_val=8/3, dt=dt, small_rho=False):
    """
    Collect the Jacobian eigenvalues of all real equilibria as a function of rho.

    For small_rho = True, rho is scanned over [0, 1).
    For small_rho = False, rho is scanned over [1, 28).

    Parameters
    ----------
    sols : list of SymPy solutions {x: expr, y: expr, z: expr}
    x, y, z : SymPy symbols
        State variables.
    rho_sym, beta_sym : SymPy symbols
        Parameter symbols.
    sigma, beta_val : float
        Lorenz parameters.
    dt : float
        Step in rho-space.
    small_rho : bool
        Whether to scan [0,1) or [1,28).

    Returns
    -------
    eigvals_by_rho : dict
        keys: rho values (floats)
        values: list of eigenvalue arrays (one per equilibrium at that rho).
    """
    if small_rho:
        rho_values = np.arange(0,1,dt=dt)
    else:
        rho_values = np.arange(1,28,dt=dt)
    eigvals_by_rho = {}

    for rho_val in rho_values:
        eigvals_list = []
        for sol in sols:
            eq = get_equilibria_for_rho(sol, x, y, z,
                                        rho_sym, beta_sym, rho_val,
                                        sigma=sigma, beta_val=beta_val)
            if eq is not None:
                eigvals_list.append(eq["eigvals"])
        eigvals_by_rho[rho_val] = eigvals_list

    return eigvals_by_rho


# ============================================================
# 5. Wing classification and derived data
# ============================================================


def ZY_data_and_wings(rho, t_min, t_max, dt=dt):
    """
    Helper: calculate Y(t), Z(t) and left/right-wing masks.

    Returned:
        y, z, left_mask, right_mask
    """
    if t_min >= t_max:
        raise ValueError("t_min must be smaller than t_max")

    x, y, z = get_vectors(rho, t_min=t_min, t_max=t_max, dt=dt)

    right = x >= 0
    left  = x < 0

    return y, z, left, right


# ============================================================
# 6. Local maxima of Z(t)
# ============================================================


def get_Z_maxima(rho=24.06, skip_first=5, t_min=t_min, t_max=t_max, dt=dt):
    """
    Determine the local maxima Z_n of Z(t) for a given rho.

    Parameters
    ----------
    rho : float
        Parameter rho.
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
    _, _, z = get_vectors(rho, t_min=t_min, t_max=t_max, dt=dt)


    t_max = []
    Z_max = []

    # How we find the maxima: z[i-1] < z[i] > z[i+1]
    for i in range(1, len(z) - 1):
        if z[i] > z[i-1] and z[i] > z[i+1]:
            t_max.append(t[i])
            Z_max.append(z[i])

    Z_max = np.array(Z_max)
    t_max = np.array(t_max)

    # Throw the first 20 maxima away ('skip_first').
    if len(Z_max) > skip_first:
        Z_max = Z_max[skip_first:]
        t_max = t_max[skip_first:]

    return t_max, Z_max




def lorenz_return_map_from_data(Z_n):
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













    