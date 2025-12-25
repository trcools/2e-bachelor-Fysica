"""
Note:

This module was refactored with assistance from an AI tool to improve structure,
naming, and documentation for readability and maintainability.

"""

# -----------------------------------------------------------------------------
# imports
# -----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# -----------------------------------------------------------------------------
# parameters
# -----------------------------------------------------------------------------

def default_parameters():
    """
    Return default parameters for the glucose degradation model.
    
    The 2D glucose model is:
        dX/dt = -X + aY + X²Y
        dY/dt = b - aY - X²Y
    
    where X represents enzyme concentration and Y represents protein.
    
    Returns
    -------
    a : float
        Enzyme-protein interaction parameter (default 0.06)
    b : float
        Protein production rate (default 0.4)
    x_lim : tuple
        Domain limits for X coordinate (0, 4)
    y_lim : tuple
        Domain limits for Y coordinate (0, 7)
    """
    a = 0.06
    b = 0.4
    x_lim = (0, 4)
    y_lim = (0, 7)

    return a, b, x_lim, y_lim


# -----------------------------------------------------------------------------
# Glucose system functions part 3
# -----------------------------------------------------------------------------

def glucose_rhs(X, Y, a, b):
    """
    Compute the right-hand side of the glucose degradation model.
    
    The 2D glucose model:
        dX/dt = -X + aY + X²Y
        dY/dt = b - aY - X²Y
    
    Parameters
    ----------
    X : float or ndarray
        Enzyme concentration
    Y : float or ndarray
        Protein concentration
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    
    Returns
    -------
    dX : float or ndarray
        Time derivative of X
    dY : float or ndarray
        Time derivative of Y
    """
    dX = -X + a*Y + (X**2)*Y
    dY =  b - a*Y - (X**2)*Y
    return dX, dY

def get_vector_grid(xlim=(0, 4), ylim=(0, 7), a=0.06, b=0.4, n=250):
    """
    Compute the vector field (dXdt, dYdt) on a regular grid.
    
    Used for streamline plots to visualize the flow of the glucose model.
    
    Parameters
    ----------
    xlim : tuple
        X domain limits (x_min, x_max)
    ylim : tuple
        Y domain limits (y_min, y_max)
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    n : int
        Number of grid points in each direction (default 250)
    
    Returns
    -------
    x : ndarray
        1D array of X coordinates
    y : ndarray
        1D array of Y coordinates
    dXdt : ndarray
        2D grid of dX/dt values
    dYdt : ndarray
        2D grid of dY/dt values
    """
    x = np.linspace(xlim[0], xlim[1], n)
    y = np.linspace(ylim[0], ylim[1], n)
    X, Y = np.meshgrid(x, y)
    dXdt, dYdt = glucose_rhs(X, Y, a, b)
    return x, y, dXdt, dYdt

def nullclines(x, a, b):
    """
    Compute the nullclines (curves where derivatives are zero).
    
    For the glucose model:
        X-nullcline (X'=0): Y = X / (a + X²)
        Y-nullcline (Y'=0): Y = b / (a + X²)
    
    Parameters
    ----------
    x : float or ndarray
        X coordinate value(s)
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    
    Returns
    -------
    Y_xnull : float or ndarray
        Y values on the X-nullcline
    Y_ynull : float or ndarray
        Y values on the Y-nullcline
    """
    return x/(a + x**2), b/(a + x**2)

def equilibrium(a, b):
    """
    Compute the equilibrium point (fixed point) of the glucose model.
    
    At equilibrium: X* = b, Y* = b/(a + b²)
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    
    Returns
    -------
    Xeq : float
        X coordinate of equilibrium
    Yeq : float
        Y coordinate of equilibrium
    """
    return b, b/(a + b**2)

def mask_in_window(x, y, xlim, ylim):
    """
    Create a boolean mask for points inside a rectangular window.
    
    Parameters
    ----------
    x : float or ndarray
        X coordinate(s)
    y : float or ndarray
        Y coordinate(s)
    xlim : tuple
        X domain limits (x_min, x_max)
    ylim : tuple
        Y domain limits (y_min, y_max)
    
    Returns
    -------
    mask : bool or ndarray of bool
        True where points are inside the window, False otherwise
    """
    return ((xlim[0] <= x) & (x <= xlim[1]) &
            (ylim[0] <= y) & (y <= ylim[1]))

def plot_curve_in_window(ax, x, y, xlim, ylim, *, label=None, **plot_kwargs):
    """
    Plot a curve (x, y) only where it lies inside the plot window.
    
    Automatically masks points outside xlim and ylim before plotting.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes for plotting
    x : ndarray
        X coordinates of the curve
    y : ndarray
        Y coordinates of the curve
    xlim : tuple
        X domain limits (x_min, x_max)
    ylim : tuple
        Y domain limits (y_min, y_max)
    label : str, optional
        Legend label for the curve
    **plot_kwargs
        Additional keyword arguments passed to ax.plot()
    """
    mask = mask_in_window(x, y, xlim, ylim)
    if np.any(mask):
        ax.plot(x[mask], y[mask], label=label, **plot_kwargs)


def draw_nullclines_panel(
        ax, a=0.06, b=0.4, xlim=(0, 4), ylim=(0, 7),
        n=250,
        density=0.6, alpha=0.35,
        n_nullcline=800, show_equilibrium=True,
        title=None,
        nullcline_colors=("tab:red", "tab:blue")
        ):
    """
    Draw a complete phase plane analysis panel on an existing matplotlib axis.
    
    Plots the vector field (streamlines), nullclines with direction arrows,
    and the equilibrium point.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes for plotting
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    xlim : tuple
        X domain limits (x_min, x_max)
    ylim : tuple
        Y domain limits (y_min, y_max)
    n : int
        Number of grid points for vector field (default 250)
    density : float
        Density of streamlines (default 0.6)
    alpha : float
        Transparency of vector field (default 0.35)
    n_nullcline : int
        Number of points for nullcline curves (default 800)
    show_equilibrium : bool
        Whether to plot the equilibrium point (default True)
    title : str, optional
        Title for the subplot
    """
    
    # --- Plot setup ---
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    if title is not None:
        ax.set_title(title)

    # --- Vector field ---
    x, y, dXdt, dYdt = get_vector_grid(xlim, ylim, a, b, n) # Get vector field data

    streamplt = ax.streamplot(x, y, dXdt, dYdt, # Vector field (streamlines)
                              density=density, linewidth=0.8, arrowsize=1.0,
                              minlength=0.08, maxlength=4.0, zorder=1)
    # set alpha for better visibility
    streamplt.lines.set_alpha(alpha)
    streamplt.arrows.set_alpha(alpha)

    # --- Nullclines (analytic) ---
    xn = np.linspace(xlim[0], xlim[1], n_nullcline)
    Y_xnull, Y_ynull = nullclines(xn, a, b)

    c_xnull, c_ynull = nullcline_colors
    plot_curve_in_window(ax, xn, Y_xnull, xlim, ylim, label="X' = 0", color=c_xnull, linewidth=2.0)
    plot_curve_in_window(ax, xn, Y_ynull, xlim, ylim, label="Y' = 0", color=c_ynull, linewidth=2.0)

    # --- Equilibrium (intersection of nullclines) ---
    if show_equilibrium:
        Xeq, Yeq = equilibrium(a, b)
        if xlim[0] <= Xeq <= xlim[1] and ylim[0] <= Yeq <= ylim[1]:
            ax.scatter(Xeq, Yeq, label="Equilibrium", color="k", zorder=3)


    ax.grid(True, alpha=0.5)
    

#=============================================================================
# Change vectors and direction arrows
#=============================================================================

def draw_change_vectors(ax, points, a, b, *,
                        len_comp=0.28, len_res=0.35,
                        show_resultant=True,
                        color_x="red", color_y="blue", color_res="k",
                        width=0.006, zorder=4):
    pts = np.asarray(points, dtype=float)
    Xp, Yp = pts[:, 0], pts[:, 1]
    dX, dY = glucose_rhs(Xp, Yp, a, b)

    sx = np.sign(dX)
    sy = np.sign(dY)

    # component arrows (constant length, only sign matters)
    Ux, Vx = sx * len_comp, np.zeros_like(sx)
    Uy, Vy = np.zeros_like(sy), sy * len_comp

    ax.quiver(Xp, Yp, Ux, Vx, angles="xy", scale_units="xy", scale=1,
              color=color_x, width=width, zorder=zorder)
    ax.quiver(Xp, Yp, Uy, Vy, angles="xy", scale_units="xy", scale=1,
              color=color_y, width=width, zorder=zorder)

    if show_resultant:
        # resultant direction (normalize to constant length)
        R = np.sqrt(dX*dX + dY*dY)
        # avoid division by zero
        R = np.where(R == 0, 1.0, R)
        Ur = (dX / R) * len_res
        Vr = (dY / R) * len_res
        ax.quiver(Xp, Yp, Ur, Vr, angles="xy", scale_units="xy", scale=1,
                  color=color_res, width=width, zorder=zorder)


def region_representative_points(a, b, xlim, ylim, nx=50, ny=50, delta=None):
    """
    Select one representative point per region bounded by the nullclines.

    We classify regions by the signs of:
      sX = sign(Y - Y_xnull(X))  -> sign of dX
      sY = sign(Y_ynull(X) - Y)  -> sign of dY

    For clarity, choose points far from both nullclines and away from
    plot edges so arrows are fully visible.

    Returns up to four points (one per sign combination).
    """
    x = np.linspace(xlim[0], xlim[1], nx)
    y = np.linspace(ylim[0], ylim[1], ny)
    X, Y = np.meshgrid(x, y)

    Yx, Yy = nullclines(X, a, b)

    if delta is None:
        delta = 0.06 * (ylim[1] - ylim[0])

    sX = np.sign(Y - Yx)
    sY = np.sign(Yy - Y)

    edge_margin_x = 0.05 * (xlim[1] - xlim[0])
    edge_margin_y = 0.05 * (ylim[1] - ylim[0])
    away_from_edges = ((X > xlim[0] + edge_margin_x) & (X < xlim[1] - edge_margin_x) &
                       (Y > ylim[0] + edge_margin_y) & (Y < ylim[1] - edge_margin_y))

    combos = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    pts = []
    # Distance metric: product of distances to each nullcline for better separation
    dist = np.abs(Y - Yx) * np.abs(Y - Yy)

    for cx, cy in combos:
        # progressively relax delta if region not found
        found = False
        for scale in (1.0, 0.7, 0.5, 0.3):
            far = (np.abs(Y - Yx) > (delta * scale)) & (np.abs(Y - Yy) > (delta * scale))
            valid = far & away_from_edges & (sX != 0) & (sY != 0)
            m = valid & (sX == cx) & (sY == cy)
            if np.any(m):
                iy, ix = np.where(m)
                dsel = dist[m]
                k = int(np.argmax(dsel))
                i, j = iy[k], ix[k]
                pts.append((float(X[i, j]), float(Y[i, j])))
                found = True
                break
        if not found:
            # as last resort, pick mid-grid point in this sign region
            m = (sX == cx) & (sY == cy)
            if np.any(m):
                iy, ix = np.where(m)
                k = len(iy) // 2
                i, j = iy[k], ix[k]
                pts.append((float(X[i, j]), float(Y[i, j])))

    return pts


def _equilibrium_xy(a, b):
    xe = b
    ye = b / (a + b * b)
    return xe, ye


def nudge_points(points, a, b, xlim, ylim, delta=None):
    """Apply small, targeted nudges to improve arrow placement clarity.

    - Move upper-right arrows slightly toward the equilibrium to avoid the legend.
    - Push top-left arrows a bit away from the Y'=0 nullcline.
    - Shift bottom-left arrows slightly to the right to avoid the left edge.
    """
    if not points:
        return points
    if delta is None:
        delta = 0.06 * (ylim[1] - ylim[0])

    xe, ye = _equilibrium_xy(a, b)
    w = (xlim[1] - xlim[0])
    h = (ylim[1] - ylim[0])

    nudged = []
    for (x, y) in points:
        # Rule 1: if near the upper-right corner, move closer to equilibrium
        if (x > xlim[1] - 0.6) and (y > ylim[1] - 0.8):
            x = x - 0.35 * (x - xe)
            y = y - 0.35 * (y - ye)

        # Rule 2: if in the upper-left area and too close to Y'=0, nudge away
        if (x < xlim[0] + 0.8) and (y > ylim[0] + 0.8 * h):
            Yy = b / (a + x * x)
            if abs(y - Yy) < 0.6 * delta:
                sign = 1.0 if (y > Yy) else -1.0
                y = y + sign * 0.35 * delta
            # Slight left shift per request (very small)
            x = x - 0.03 * w

        # Rule 3: if in the lower-left area, nudge a tiny bit right
        if (y < ylim[0] + 0.20 * h) and (x < xlim[0] + 0.35 * w):
            x = x + 0.015 * w

        # Keep within bounds
        x = float(np.clip(x, xlim[0] + 1e-6, xlim[1] - 1e-6))
        y = float(np.clip(y, ylim[0] + 1e-6, ylim[1] - 1e-6))
        nudged.append((x, y))

    return nudged


def plot_nullclines(a=0.06, b=0.4, xlim=(0,4), ylim=(0,7),
                             n=250, density=0.6, alpha=0.35,
                             show_resultant=True):
    """
    Combine the classic nullclines figure (title + legend) with the
    photo-style directional arrows.

    - Keeps the original title and legend from the classic plot.
    - Uses colored nullclines and component/resultant arrows like the photo style.

    Parameters
    ----------
    a, b : float
        Model parameters.
    xlim, ylim : tuple
        Axis limits.
    n : int
        Grid resolution for the streamline plot.
    density : float
        Streamline density.
    alpha : float
        Streamline transparency.
    show_resultant : bool
        Whether to draw the black resultant arrows.
    """
    plt.close("all")
    fig, ax = plt.subplots(figsize=(8, 6))

    # classic panel with colored nullclines
    draw_nullclines_panel(
        ax, a=a, b=b, xlim=xlim, ylim=ylim,
        n=n, density=density, alpha=alpha,
        nullcline_colors=("red", "blue"),
        title=rf"Degradation of Glucose: Nullclines with $a={a}$ and $b={b}$"
    )

    # Overlay photo-style change vectors: one pair per region
    pts = region_representative_points(a, b, xlim, ylim)
    pts = nudge_points(pts, a, b, xlim, ylim)
    draw_change_vectors(ax, pts, a, b, show_resultant=show_resultant,
                        color_x="red", color_y="blue", color_res="k",
                        len_comp=0.32, len_res=0.40, width=0.008, zorder=4)

    # Keep the classic legend (nullclines + equilibrium)
    ax.legend(loc="upper right", framealpha=0.9)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.grid(True, alpha=0.35)
    plt.show()


# -----------------------------------------------------------------------------
# Glucose system functions part 2
# -----------------------------------------------------------------------------

def jacobian(X, Y, a, b):
    """
    Compute the Jacobian matrix of the glucose model at point (X, Y).
    
    For the 2D system:
        dX/dt = -X + aY + X²Y
        dY/dt = b - aY - X²Y
    
    The Jacobian is:
        J(X,Y) = [[-1 + 2XY,    a + X²],
                  [-2XY,        -(a + X²)]]
    
    Parameters
    ----------
    X : float or ndarray
        X coordinate(s)
    Y : float or ndarray
        Y coordinate(s)
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate (not used in Jacobian but included for consistency)
    
    Returns
    -------
    J : ndarray
        Jacobian matrix (2×2) or array of 2×2 matrices
    """
    return np.array([
        [-1 + 2*X*Y,  a + X**2],
        [-2*X*Y,     -(a + X**2)]
    ], dtype=float)


def get_equilibrium_eigvals(a, b):
    """
    Compute the equilibrium point, Jacobian, and its eigenvalues.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    
    Returns
    -------
    eq : tuple
        (Xeq, Yeq) coordinates of equilibrium
    J_eq : ndarray
        2×2 Jacobian matrix at equilibrium
    eig : ndarray
        Eigenvalues of the Jacobian (complex numbers in general)
    """
    Xeq, Yeq = equilibrium(a, b)
    J_eq = jacobian(Xeq, Yeq, a, b)
    eig = np.linalg.eigvals(J_eq)
    return (Xeq, Yeq), J_eq, eig

# -----------------------------------------------------------------------------
# Glucose system functions part 3
# -----------------------------------------------------------------------------

def rhs_state(t, z, a, b):
    """
    Right-hand side function for scipy.integrate.solve_ivp.
    
    Wraps glucose_rhs for compatibility with ODE solver API.
    
    Parameters
    ----------
    t : float
        Time (not used but required by solve_ivp)
    z : list or ndarray
        State vector [X, Y]
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    
    Returns
    -------
    dzdt : list
        Time derivatives [dX/dt, dY/dt]
    """
    X, Y = z
    dX, dY = glucose_rhs(X, Y, a, b)
    return [dX, dY]


def default_initial_conditions(a, b, xlim, ylim, step=0.08):
    """
    Generate sets of initial conditions for trajectory integration.
    
    Creates two groups of initial conditions:
    - Near equilibrium: for studying local stability
    - Near edges: for studying global dynamics
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    xlim : tuple
        X domain limits (x_min, x_max)
    ylim : tuple
        Y domain limits (y_min, y_max)
    step : float
        Relative step size (fraction of domain width/height), default 0.08
    
    Returns
    -------
    near_eq : list of tuples
        Initial conditions near equilibrium
    near_edge : list of tuples
        Initial conditions near domain boundaries
    """
    Xeq, Yeq = equilibrium(a, b)
    dx = step * (xlim[1] - xlim[0])
    dy = step * (ylim[1] - ylim[0])

    # near equilibrium (local)
    near_eq = [
        (Xeq + dx, Yeq),
        (Xeq - dx, Yeq),
        (Xeq, Yeq + dy),
        (Xeq, Yeq - dy),
    ]

    # near edges
    near_edge = [
        (xlim[0] + 0.05*(xlim[1]-xlim[0]), ylim[0] + 0.05*(ylim[1]-ylim[0])),
        (xlim[1] - 0.05*(xlim[1]-xlim[0]), ylim[0] + 0.05*(ylim[1]-ylim[0])),
        (xlim[0] + 0.05*(xlim[1]-xlim[0]), ylim[1] - 0.05*(ylim[1]-ylim[0])),
        (xlim[1] - 0.05*(xlim[1]-xlim[0]), ylim[1] - 0.05*(ylim[1]-ylim[0])),
        (xlim[1] - 0.02*(xlim[1]-xlim[0]), 0.5*(ylim[0]+ylim[1])),
    ]

    return near_eq , near_edge

def draw_trajectories(ax, a, b, xlim, ylim, initials, t_span=(0,60), t_eval_n=2500,
                      line_kw=None, start_kw=None):
    """
    Integrate and plot trajectories from multiple initial conditions.
    
    Solves the glucose model ODE for each initial condition and plots the result
    on the provided axes. Optionally marks starting points.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes for plotting
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    xlim : tuple
        X domain limits (not enforced, used for reference)
    ylim : tuple
        Y domain limits (not enforced, used for reference)
    initials : list of tuples
        Initial conditions [(X0, Y0), ...]
    t_span : tuple
        Integration time interval (t_min, t_max), default (0, 60)
    t_eval_n : int
        Number of evaluation points, default 2500
    line_kw : dict, optional
        Keyword arguments for trajectory lines (passed to ax.plot)
    start_kw : dict, optional
        Keyword arguments for starting point markers (passed to ax.plot)
    """

    t_eval = np.linspace(t_span[0], t_span[1], t_eval_n)
    # Ensure plotting kwargs are mappings
    if line_kw is None:
        line_kw = {}
    if start_kw is None:
        start_kw = {}

    for (X0, Y0) in initials:
        sol = solve_ivp(rhs_state, t_span, [X0, Y0],
                        args=(a, b), t_eval=t_eval,
                        rtol=1e-8, atol=1e-10)
        ax.plot(sol.y[0], sol.y[1], **line_kw)
        ax.plot([X0], [Y0], **start_kw)



def plot_2glucose_trajectories(a=0.06, b=0.4, xlim=(0,4), ylim=(0,7),
                  n=250, step=0.08,
                  t_span=(0,60), t_eval_n=2500):
    """
    Create side-by-side comparison of local and global trajectory behavior.
    
    Left panel: Trajectories starting near the equilibrium (local stability analysis)
    Right panel: Trajectories starting near domain boundaries (global dynamics)
    
    Both panels include nullclines, vector field, and direction information.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter (default 0.06)
    b : float
        Protein production rate (default 0.4)
    xlim : tuple
        X domain limits (default (0, 4))
    ylim : tuple
        Y domain limits (default (0, 7))
    n : int
        Number of grid points for vector field (default 250)
    step : float
        Relative step for generating initial conditions (default 0.08)
    t_span : tuple
        Integration time interval (default (0, 60))
    t_eval_n : int
        Number of evaluation points (default 2500)
    """

    near_eq, near_edge = default_initial_conditions(a, b, xlim, ylim, step=step)

    plt.close("all")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True, sharey=True)

    # --- Local panel ---
    ax = axes[0]
    draw_nullclines_panel(ax, a=a, b=b, xlim=xlim, ylim=ylim, n=n,
                          density=0.8, alpha=0.30, show_equilibrium=True,
                          title="Local trajectories near equilibrium")
    draw_trajectories(ax, a, b, xlim, ylim, near_eq, t_span=t_span, t_eval_n=t_eval_n)

    # --- Global panel ---
    ax = axes[1]
    draw_nullclines_panel(ax, a=a, b=b, xlim=xlim, ylim=ylim, n=n,
                          density=0.8, alpha=0.30, show_equilibrium=True,
                          title="Global trajectories further from equilibrium")
    draw_trajectories(ax, a, b, xlim, ylim, near_edge, t_span=t_span, t_eval_n=t_eval_n)

    # One general legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, framealpha=0.9, loc='upper center')

    fig.tight_layout()
    plt.show()


def plot_zoomed_spiral_convergence(a=0.06, b=0.4, xlim=(0, 1.5), ylim=(0, 3),
                                    t_span=(0, 60), t_eval_n=2500, step=0.08,
                                    single_init=(0.465, 0.6)):
    """
    Plot side-by-side visualization of spiral convergence to stable equilibrium.
    
    Left panel: Multiple trajectories from near_eq points (all marked with green circles)
    Right panel: Single detailed trajectory from single_init point
    
    Parameters
    ----------
    a, b : float
        Model parameters (default: a=0.06, b=0.4)
    xlim, ylim : tuple
        Domain limits for visualization
    t_span, t_eval_n : float, int
        Integration time span and number of evaluation points
    step : float
        Step size for generating initial conditions around equilibrium
    single_init : tuple
        Initial condition (x0, y0) for the right panel
    """
    plt.close('all')
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    ax1, ax2 = axes
    
    # Get equilibrium point
    Xeq, Yeq = equilibrium(a, b)
    
    # --- Left panel: Multiple trajectories ---
    draw_nullclines_panel(
        ax1, a=a, b=b, xlim=xlim, ylim=ylim,
        n=250, density=0.6, alpha=0.35,
        title="Multiple trajectories converging to equilibrium", show_equilibrium=True
    )
    near_eq, _ = default_initial_conditions(a, b, xlim=xlim, ylim=ylim, step=step)
    draw_trajectories(ax1, a, b, xlim=xlim, ylim=ylim, initials=near_eq, 
                      t_span=t_span, t_eval_n=t_eval_n,
                      line_kw=dict(linewidth=1.2, alpha=0.7))
    # Mark all start points
    for x0, y0 in near_eq:
        ax1.plot(x0, y0, 'go', markersize=8, zorder=5)
    # Add one entry to legend for start points
    ax1.plot([], [], 'go', markersize=8, label='Start points', zorder=5)
    ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # --- Right panel: Single detailed trajectory ---
    draw_nullclines_panel(
        ax2, a=a, b=b, xlim=xlim, ylim=ylim,
        n=500, density=2, alpha=0.35,
        title="Single trajectory", show_equilibrium=True
    )
    draw_trajectories(ax2, a, b, xlim=xlim, ylim=ylim, initials=[single_init], 
                      t_span=t_span, t_eval_n=t_eval_n,
                      line_kw=dict(linewidth=2.0, alpha=0.85, color='red'))
    ax2.plot(single_init[0], single_init[1], 'go', markersize=10, label='Start points', zorder=5)
    ax2.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # Main title
    fig.suptitle('Zoomed view: Spiral convergence to stable equilibrium (a={:.2f}, b={:.2f})'.format(a, b), 
                 fontsize=13, fontweight='bold', y=1.00)
    
    # Get all unique labels from both axes and combine
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    # Combine and remove duplicates
    all_handles = handles1 + handles2
    all_labels = labels1 + labels2
    unique_dict = dict(zip(all_labels, all_handles))
    fig.legend(unique_dict.values(), unique_dict.keys(), framealpha=0.95, loc='upper center', 
              bbox_to_anchor=(0.5, -0.02), ncol=4, fontsize=10)
    
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    plt.show()


# -----------------------------------------------------------------------------
# Glucose system functions part 4: Bifurcation Analysis
# -----------------------------------------------------------------------------

def trace_det_at_eq(a, b):
    """
    Compute trace and determinant of Jacobian at equilibrium.
    
    Mathematical derivation: At equilibrium (X*, Y*) = (b, b/(a+b²)),
    the trace is τ = 1 - s - 2a/s where s = a + b².
    
    Args:
        a, b: Parameters of the glucose model
    
    Returns:
        tau: Trace of Jacobian (stability indicator, bifurcation when τ=0)
        Delta: Determinant = a + b² (always > 0 for a > 0)
    """
    s = a + b**2
    tau = 1 - s - 2*a/s  # Formula derived from linearization at equilibrium
    Delta = s
    return tau, Delta

def _max_real_part_from_trace_det_vec(tau, det):
    """
    Vectorized maximum real part of eigenvalues of a 2x2 matrix given trace and determinant.

    For eigenvalues solving λ^2 - τ λ + Δ = 0:
    - If discriminant D = τ^2 - 4Δ >= 0 (real eigenvalues), max real = (τ + sqrt(D))/2
    - If D < 0 (complex pair), real parts are equal τ/2, so max real = τ/2
    """
    tau = np.asarray(tau)
    det = np.asarray(det)
    disc = tau*tau - 4.0*det
    sqrt_disc = np.sqrt(np.maximum(disc, 0.0))
    return np.where(disc >= 0.0, 0.5*(tau + sqrt_disc), 0.5*tau)

def bcrit_values(a):
    """
    Calculate critical b-values where Hopf bifurcations occur (τ(b) = 0).
    
    Mathematical derivation:
    Setting τ = 1 - s - 2a/s = 0 and solving s² - s + 2a = 0 yields:
        s_± = (1 ± √(1 - 8a)) / 2
    
    Then converting back to b using b² = s - a:
        b_± = √(s_± - a)
    
    Real solutions exist only if discriminant 1 - 8a ≥ 0, i.e., a ≤ 1/8.
    
    Args:
        a: Parameter of the glucose model
    
    Returns:
        List of critical b-values [b_-, b_+] where bifurcations occur
        Empty list if no bifurcations exist for this a value
    """
    disc = 1 - 8*a  # Discriminant of quadratic s² - s + 2a = 0
    if disc < 0:
        return []  # No real solutions → no Hopf bifurcation for this a
    
    s1 = (1 + np.sqrt(disc))/2  # s_+ (larger root)
    s2 = (1 - np.sqrt(disc))/2  # s_- (smaller root)
    
    bs = []
    for s in (s2, s1):  # Process in order from smaller to larger b
        b2 = s - a      # Since b² = s - a
        if b2 >= 0:
            bs.append(np.sqrt(b2))  # b_± = √(s_± - a)
    return bs

def scan_bifurcation_in_b(a=0.06, b_min=0.0, b_max=1.2, n=2000):
    """
    Scan the stability of the equilibrium as b varies (with a fixed).
    
    Plots max(Re(eigenvalue)) vs b and marks bifurcation points where
    stability changes (where max Re(λ) crosses zero).
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter (kept fixed)
    b_min : float
        Minimum b value to scan (default 0.0)
    b_max : float
        Maximum b value to scan (default 1.2)
    n : int
        Number of b values to evaluate (default 2000)
    
    Returns
    -------
    b_crit : list
        Approximate critical b values where bifurcations occur
    """
    b_vals = np.linspace(b_min, b_max, n)
    max_real = np.empty_like(b_vals)

    for i, b in enumerate(b_vals):
        (_, _), _, eig = get_equilibrium_eigvals(a, b)
        max_real[i] = np.max(np.real(eig))   # stability indicator

    # detect sign changes in max_real: where it crosses 0
    s = np.sign(max_real)
    idx = np.where(np.diff(s) != 0)[0]  # indices near crossing

    b_crit = []
    for k in idx:
        # linear interpolation for a decent estimate
        b0, b1 = b_vals[k], b_vals[k+1]
        f0, f1 = max_real[k], max_real[k+1]
        bc = b0 - f0*(b1-b0)/(f1-f0)
        b_crit.append(bc)

    # plot max Re(lambda) vs b
    plt.figure(figsize=(7,4))
    plt.plot(b_vals, max_real)
    plt.axhline(0, linewidth=1)
    for bc in b_crit:
        plt.axvline(bc, linestyle="--")
    plt.xlabel("b")
    plt.ylabel("max Re(eigenvalue at equilibrium)")
    plt.title(f"Stability of equilibrium vs b (a={a})")
    plt.grid(True, alpha=0.4)
    plt.show()

    return b_crit

def plot_varying_b(a):
    """
    Plot equilibrium curve with regions of stability vs instability.
    
    Shows how Y* moves with b, colored by stability (trace criterion).
    Vertical dashed lines indicate bifurcation points where stability changes.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter (kept fixed)
    """
    b_vals = np.linspace(0, 1.2, 600)
    Xeq, Yeq = b_vals, b_vals/(a + b_vals**2)
    tau = np.array([trace_det_at_eq(a, bb)[0] for bb in b_vals])

    stable = tau < 0

    plt.figure(figsize=(7,4))
    plt.plot(b_vals[stable],   Yeq[stable],   label="stable eq")
    plt.plot(b_vals[~stable],  Yeq[~stable],  label="unstable eq")
    for bc in bcrit_values(a):
        plt.axvline(bc, linestyle="--")
    plt.xlabel("b")
    plt.ylabel("Y*")
    plt.title("Equilibrium stability vs b (a fixed)")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.show()

def plot_bifurcation_summary_figures(a=0.06, xlim=(0, 1.2), fast=True):
    """
    Create a comprehensive set of bifurcation visuals for varying b with fixed a.

    Figures:
    - Top: max Re(eigenvalue) vs b with shaded stable/unstable regions and critical b lines.
    - Middle: Equilibrium Y*(b) with stability coloring and markers at b_crit.
    - Bottom: Phase portraits for b before, at, and after the first bifurcation.

    Parameters
    ----------
    a : float
        Fixed parameter a.
    xlim : tuple
        Range of b values shown.
    """
    # --- Compute stability curve (vectorized, no eigensolver) ---
    n_b = 30 if fast else 200
    b_vals = np.linspace(xlim[0], xlim[1], n_b)
    s = a + b_vals**2
    tau_vals = 1.0 - s - (2.0*a)/s
    det_vals = s
    max_real = _max_real_part_from_trace_det_vec(tau_vals, det_vals)

    b_crit = bcrit_values(a)

    plt.close('all')
    fig = plt.figure(figsize=(12, 12))

    # Top panel: max Re(lambda) vs b
    ax1 = fig.add_subplot(3,1,1)
    ax1.plot(b_vals, max_real, color='tab:blue', label='max Re(λ) at equilibrium')
    ax1.axhline(0, color='k', linewidth=1)
    for bc in b_crit:
        ax1.axvline(bc, linestyle='--', color='tab:red', label='b_crit' if bc == b_crit[0] else None)
    # Shade stable/unstable regions using tau sign (equivalent here)
    stable_mask = tau_vals < 0
    ax1.fill_between(b_vals, max_real, 0, where=stable_mask, color='tab:green', alpha=0.15, label='stable')
    ax1.fill_between(b_vals, max_real, 0, where=~stable_mask, color='tab:orange', alpha=0.12, label='unstable')
    ax1.set_xlim(*xlim)
    ax1.set_xlabel('b')
    ax1.set_ylabel('max Re(λ)')
    ax1.set_title(f'Stability vs b (a={a})')
    handles, labels = ax1.get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    ax1.legend(uniq.values(), uniq.keys(), framealpha=0.9, loc='upper right')
    ax1.grid(True, alpha=0.3)

    # Middle panel: Y*(b) with stability coloring
    ax2 = fig.add_subplot(3,1,2)
    Xeq = b_vals
    Yeq = b_vals/(a + b_vals**2)
    ax2.plot(b_vals[stable_mask], Yeq[stable_mask], color='tab:green', label='stable equilibrium')
    ax2.plot(b_vals[~stable_mask], Yeq[~stable_mask], color='tab:orange', label='unstable equilibrium')
    for bc in b_crit:
        ax2.axvline(bc, linestyle='--', color='tab:red', label='b_crit' if bc == b_crit[0] else None)
        ax2.scatter(bc, bc/(a+bc**2), color='tab:red', zorder=3)
    ax2.set_xlim(*xlim)
    ax2.set_xlabel('b')
    ax2.set_ylabel('Y*')
    ax2.set_title('Equilibrium position vs b (stability coloring)')
    handles, labels = ax2.get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    ax2.legend(uniq.values(), uniq.keys(), framealpha=0.9, loc='upper right')
    ax2.grid(True, alpha=0.3)

    # Bottom panel: phase portraits around first bifurcation (if exists)
    ax3a = fig.add_subplot(3,3,7)
    ax3b = fig.add_subplot(3,3,8)
    ax3c = fig.add_subplot(3,3,9)
    if len(b_crit) > 0:
        b0 = b_crit[0]
        eps = 0.03
        bs = [max(xlim[0], b0 - eps), b0, min(xlim[1], b0 + eps)]
        titles = [f'before (b={bs[0]:.3f})', f'at bifurcation (b={bs[1]:.3f})', f'after (b={bs[2]:.3f})']
        axes = [ax3a, ax3b, ax3c]
        for axp, bb, tt in zip(axes, bs, titles):
            if fast:
                # Lightweight panel: sparse field, no trajectories
                draw_nullclines_panel(axp, a=a, b=bb, xlim=(0, 1.5), ylim=(0, 3),
                                      n=50, density=0.25, alpha=0.25,
                                      title=tt, show_equilibrium=True)
            else:
                draw_nullclines_panel(axp, a=a, b=bb, xlim=(0, 1.5), ylim=(0, 3),
                                      n=150, density=0.8, alpha=0.3,
                                      title=tt, show_equilibrium=True)
                near_eq, _ = default_initial_conditions(a, bb, xlim=(0,1.5), ylim=(0,3), step=0.1)
                draw_trajectories(axp, a, bb, xlim=(0,1.5), ylim=(0,3), initials=near_eq,
                                  t_span=(0, 40), t_eval_n=1500,
                                  line_kw=dict(linewidth=0.8, alpha=0.85))
    else:
        ax3b.text(0.5, 0.5, 'No bifurcation for this a', ha='center', va='center')
        for axp in (ax3a, ax3b, ax3c):
            axp.axis('off')

    fig.suptitle(f'Bifurcation summary (a={a})', fontsize=14, fontweight='bold', y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

def plot_bifurcation_phase_portraits(a=0.06, b_crit_vals=None, xlim=(0, 1.5), ylim=(0, 3)):
    """
    Plot phase portraits before, at, and after each bifurcation point.
    
    Creates a grid of subplots showing how the qualitative dynamics change
    as a parameter crosses a bifurcation point.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter (kept fixed)
    b_crit_vals : list, optional
        Critical b values where bifurcations occur. If None, computed automatically.
    xlim : tuple
        X domain limits (default (0, 1.5))
    ylim : tuple
        Y domain limits (default (0, 3))
    """
    if b_crit_vals is None:
        b_crit_vals = bcrit_values(a)
    
    n_bifurcations = len(b_crit_vals)
    if n_bifurcations == 0:
        print("No bifurcations found")
        return
    
    # Create subplots for each bifurcation: before, at, and after
    fig, axes = plt.subplots(n_bifurcations, 3, figsize=(14, 4*n_bifurcations))
    if n_bifurcations == 1:
        axes = axes.reshape(1, -1)
    
    for row, b_crit in enumerate(b_crit_vals):
        epsilon = 0.02  # Small perturbation around bifurcation
        b_vals = [b_crit - epsilon, b_crit, b_crit + epsilon]
        titles = [f"Before (b={b_crit-epsilon:.4f})", 
                  f"At bifurcation (b={b_crit:.4f})", 
                  f"After (b={b_crit+epsilon:.4f})"]
        
        for col, (b_val, title) in enumerate(zip(b_vals, titles)):
            ax = axes[row, col]
            
            # Draw nullclines and vector field
            draw_nullclines_panel(ax, a=a, b=b_val, xlim=xlim, ylim=ylim,
                                 n=250, density=0.8, alpha=0.3,
                                 title=title, show_equilibrium=True)
            
            # Add some trajectories
            near_eq, _ = default_initial_conditions(a, b_val, xlim, ylim, step=0.1)
            draw_trajectories(ax, a, b_val, xlim, ylim, initials=near_eq,
                            t_span=(0, 40), t_eval_n=2000,
                            line_kw=dict(linewidth=0.7, zorder=2.5, alpha=0.8))
    
    plt.tight_layout()
    plt.show()

def summarize_bifurcations_in_b(a=0.06):
    """
    Print a detailed summary of bifurcations when varying b.
    
    Displays:
    - Critical b values where bifurcations occur
    - Stability status at various b values
    - Equilibrium positions as b changes
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter (kept fixed)
    """
    print(f"BIFURCATION ANALYSIS: Glucose model with a = {a}")
    print("=" * 60)
    
    b_crit = bcrit_values(a)
    
    if len(b_crit) == 0:
        print(f"No bifurcations exist for a = {a}")
        return
    
    print(f"\nCritical b-values (where τ(b) = 0):")
    for i, b_c in enumerate(b_crit):
        print(f"  b_{i+1} = {b_c:.6f}")
    
    # Evaluate stability at different regions
    print(f"\nStability along the b-axis:")
    test_points = [0.01] + [b_c + (-1)**(i % 2) * 0.01 for i, b_c in enumerate(b_crit)] + [b_crit[-1] + 0.1]
    test_points = sorted(set(test_points))
    
    for b_test in test_points:
        (Xeq, Yeq), _, eig = get_equilibrium_eigvals(a, b_test)
        max_re = np.max(np.real(eig))
        status = "STABLE" if max_re < 0 else "UNSTABLE"
        print(f"  b = {b_test:.4f}: max Re(λ) = {max_re:8.5f}  →  {status}")
    
    print(f"\nEquilibrium displacement:")
    b_test_vals = np.linspace(0.01, b_crit[-1] + 0.1, 5)
    for b_test in b_test_vals:
        Xeq, Yeq = equilibrium(a, b_test)
        print(f"  b = {b_test:.4f}: (X*, Y*) = ({Xeq:.4f}, {Yeq:.4f})")
    
    print("=" * 60)


# -----------------------------------------------------------------------------
# Glucose system functions part 5: 2D Bifurcation Analysis
# -----------------------------------------------------------------------------

def trace_at_equilibrium(a, b):
    """
    Compute the trace of the Jacobian at equilibrium.
    
    The trace is a stability indicator: negative trace (with positive determinant)
    indicates stable equilibrium.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    
    Returns
    -------
    tau : float
        Trace of the Jacobian matrix at (X*, Y*)
    """
    Xeq, Yeq = equilibrium(a, b)
    J = jacobian(Xeq, Yeq, a, b)
    return np.trace(J)

def hopf_curves(a_vals):
    """
    Compute the Hopf bifurcation curves b_-(a) and b_+(a).
    
    These curves define the boundaries of the region where the equilibrium
    transitions from stable to unstable (and vice versa).
    
    Bifurcations exist only for a ≤ 1/8. For larger a, returns NaN.
    
    Parameters
    ----------
    a_vals : float or ndarray
        Enzyme-protein interaction parameter values
    
    Returns
    -------
    b_minus : float or ndarray
        Lower bifurcation curve (smaller b values)
    b_plus : float or ndarray
        Upper bifurcation curve (larger b values)
    """
    a_vals = np.asarray(a_vals)
    disc = 1 - 8*a_vals
    b_minus = np.full_like(a_vals, np.nan, dtype=float)
    b_plus  = np.full_like(a_vals, np.nan, dtype=float)

    ok = disc >= 0
    s_plus  = (1 + np.sqrt(disc[ok]))/2
    s_minus = (1 - np.sqrt(disc[ok]))/2

    b2_plus  = s_plus  - a_vals[ok]
    b2_minus = s_minus - a_vals[ok]

    b_plus[ok]  = np.sqrt(np.maximum(b2_plus,  0.0))
    b_minus[ok] = np.sqrt(np.maximum(b2_minus, 0.0))
    return b_minus, b_plus
        

def stability_map_ab(a_min=0.0, a_max=0.14, b_min=0.0, b_max=1.2,
                     na=250, nb=250, show_legend=True):
    """
    Assignment-style 2D stability map (a,b-plane).

    - Background: stable (τ<0) vs unstable (τ>0) regions using the trace criterion
    - Overlays: Hopf curves b_-(a), b_+(a)
    - Legend: English labels for stable/unstable and the two Hopf branches

    Parameters
    ----------
    show_legend : bool, default True
        Whether to draw the legend (use False for a clean figure without labels).

    Returns
    -------
    a_vals, b_vals : ndarray
        Meshgrid axes used for the plot.
    max_real : ndarray
        Maximum real part of the eigenvalues at each (a,b) (analytic, vectorized).
    """
    a_vals = np.linspace(a_min, a_max, na)
    b_vals = np.linspace(b_min, b_max, nb)
    A, B = np.meshgrid(a_vals, b_vals, indexing="xy")

    # Vectorized stability metrics at equilibrium
    S = A + B*B
    Tau = 1.0 - S - (2.0*A)/S
    Det = S

    # Analytic max real part of eigenvalues from trace/determinant
    discr = Tau*Tau - 4.0*Det
    max_real = np.where(discr >= 0.0,
                       0.5*(Tau + np.sqrt(discr)),  # real eigenvalues
                       0.5*Tau)                      # complex pair: real part = Tau/2

    stable = Tau < 0


    plt.close('all')
    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    cmap = ListedColormap(["#fde2b8", "#c9e9c5"])  # unstable, stable
    ax.pcolormesh(a_vals, b_vals, stable.astype(int), shading="auto", cmap=cmap, vmin=0, vmax=1)

    # Hopf curves
    b_m, b_p = hopf_curves(a_vals)
    line_bm, = ax.plot(a_vals, b_m, color='#1f77b4', linewidth=2.0, linestyle='--', label='Hopf b_-(a)')
    line_bp, = ax.plot(a_vals, b_p, color='#d62728', linewidth=2.0, linestyle=':',  label='Hopf b_+(a)')

    if show_legend:
        stable_patch = Patch(facecolor="#c9e9c5", edgecolor='none', label='Stable (τ < 0)')
        unstable_patch = Patch(facecolor="#fde2b8", edgecolor='none', label='Unstable (τ > 0)')
        handles = [stable_patch, unstable_patch, line_bm, line_bp]
        ax.legend(handles=handles, loc='lower right', framealpha=0.95)

    ax.set_xlim(a_min, a_max)
    ax.set_ylim(b_min, b_max)
    ax.set_xlabel('a')
    ax.set_ylabel('b')
    ax.set_title('Stability of equilibrium in (a,b)-plane (trace criterion)')
    ax.grid(True, alpha=0.25)
    plt.show()

    return a_vals, b_vals, max_real


# -----------------------------------------------------------------------------
# Glucose system functions part 5
# -----------------------------------------------------------------------------