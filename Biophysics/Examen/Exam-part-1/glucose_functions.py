"""
Note:

This module was refactored with assistance from an AI tool to improve structure,
naming, and documentation for readability and maintainability.

"""

# -----------------------------------------------------------------------------
# imports
# -----------------------------------------------------------------------------

import sys
import os
from pathlib import Path

ROOT = Path.cwd()  # or Path(__file__).resolve().parent 
sys.path.insert(0, str(ROOT))

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from ipywidgets import Checkbox

# Import helper functions from centralized module
from personalized_layout import create_interactive_plot, fig_legend_dedup


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
    xlim = (0, 4)
    ylim = (0, 7)

    return a, b, xlim, ylim

a, b, xlim, ylim = default_parameters()
# -----------------------------------------------------------------------------
# Glucose system functions part 1
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
    dX = -X + a*Y + X**2*Y
    dY =  b - a*Y - X**2*Y
    return dX, dY

def get_vector_grid(xlim=xlim, ylim=ylim, a=a, b=b, n=250):
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
        interaction parameter
    b : float
        parameter

    
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
    b : float
    
    Returns
    -------
    Xeq : float
        X coordinate of equilibrium
    Yeq : float
        Y coordinate of equilibrium
    """
    return b, b/(a + b**2)

def mask_in_window(x, y, xlim=xlim, ylim=ylim):
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

def plot_curve_in_window(ax, x, y, xlim=xlim, ylim=ylim, *, label=None, **plot_kwargs):
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
        ax, a, b, xlim=xlim, ylim=ylim,
        n=250,
        density=0.6, alpha=0.35,
        n_nullcline=800, show_equilibrium=True,
        title=None,
        nullcline_colors=("tab:red", "tab:blue"),
        nullcline_linestyle="solid"
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
    nullcline_linestyle : str
        Line style for nullclines: 'solid', 'dashed', 'dotted', or 'dashdot' (default 'solid')
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
                              density, linewidth=0.8, arrowsize=1.0,
                              minlength=0.08, maxlength=4.0, zorder=1)
    # set alpha for better visibility
    streamplt.lines.set_alpha(alpha)
    streamplt.arrows.set_alpha(alpha)

    # --- Nullclines (analytic) ---
    xn = np.linspace(xlim[0], xlim[1], n_nullcline)
    Y_xnull, Y_ynull = nullclines(xn, a, b)

    c_xnull, c_ynull = nullcline_colors
    plot_curve_in_window(ax, xn, Y_xnull, xlim, ylim, label="X' = 0", color=c_xnull, linewidth=2.0, linestyle=nullcline_linestyle)
    plot_curve_in_window(ax, xn, Y_ynull, xlim, ylim, label="Y' = 0", color=c_ynull, linewidth=2.0, linestyle=nullcline_linestyle)

    # --- Equilibrium (intersection of nullclines) ---
    if show_equilibrium:
        Xeq, Yeq = equilibrium(a, b)
        if xlim[0] <= Xeq <= xlim[1] and ylim[0] <= Yeq <= ylim[1]:
            ax.scatter(Xeq, Yeq, label="Equilibrium", color="k", zorder=3)


    ax.grid(True, alpha=0.5)
    

#=============================================================================
# Change vectors and direction arrows
#=============================================================================

def quiver_on_curve(ax, x, y, xlim, ylim,  b, orientation="vertical", length=0.15):
    """
    Plot direction arrows along a curve inside the plot window.
    
    For the glucose model, the direction along both nullclines is determined by sign(b - X).
    
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
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate (determines direction sign)
    orientation : str, default "vertical"
        Direction of arrows ("vertical" for vertical motion on X-nullcline or "horizontal" for Y-nullcline)
    length : float, default 0.15
        Length of arrows in data coordinates
    """
    mask = mask_in_window(x, y, xlim, ylim)
    x, y = x[mask], y[mask]

    if x.size == 0:
        return

    s = np.sign(b - x) * length

    if orientation == "vertical":
        X = np.zeros_like(s)
        Y = s
    elif orientation == "horizontal":
        X = s
        Y = np.zeros_like(s)
    else:
        raise ValueError("orientation must be 'vertical' or 'horizontal'")

    ax.quiver(x, y, X, Y, angles="xy", scale_units="xy", scale=1, width=0.003, color="k", zorder=1.5)

def draw_change_vectors(ax, points, a=a, b=b, *,
                        len_comp=0.32, len_res=0.40,
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

    # resultant direction (always drawn)
    R = np.sqrt(dX*dX + dY*dY)
    R = np.where(R == 0, 1.0, R)
    Ur = (dX / R) * len_res
    Vr = (dY / R) * len_res
    ax.quiver(Xp, Yp, Ur, Vr, angles="xy", scale_units="xy", scale=1,
              color=color_res, width=width, zorder=zorder)

def draw_nullcline_arrows(ax, a=a, b=b, xlim=xlim, ylim=ylim, *, n_points=20, arrow_len=0.15, width=0.003, zorder=1.5):
    xx = np.linspace(xlim[0], xlim[1], n_points)
    Y_xnull, Y_ynull = nullclines(xx, a, b)
    quiver_on_curve(ax, xx, Y_xnull, xlim, ylim, b, orientation="vertical", length=arrow_len)
    quiver_on_curve(ax, xx, Y_ynull, xlim, ylim, b, orientation="horizontal", length=arrow_len)

def region_representative_points():
    """Hardcoded arrow positions for standard plot."""
    return [(0.1, 2.5), (1.2, 2.5), (0.4, 0.6), (1.2, 0.6)]

def nudge_points(points, a=a, b=b, xlim=xlim, ylim=ylim, delta=None):
    """Apply small nudges to improve arrow placement."""
    if not points:
        return points
    if delta is None:
        delta = 0.06 * (ylim[1] - ylim[0])

    xe, ye = equilibrium(a, b)
    w, h = (xlim[1] - xlim[0]), (ylim[1] - ylim[0])
    nudged = []
    
    for (x, y) in points:
        # Upper-right: move toward equilibrium
        if (x > xlim[1] - 0.6) and (y > ylim[1] - 0.8):
            x -= 0.35 * (x - xe)
            y -= 0.35 * (y - ye)
        # Upper-left: nudge away from Y'=0
        elif (x < xlim[0] + 0.8) and (y > ylim[0] + 0.8 * h):
            Yy = b / (a + x * x)
            if abs(y - Yy) < 0.6 * delta:
                y += (1.0 if y > Yy else -1.0) * 0.35 * delta
            x -= 0.03 * w
        # Lower-left: nudge right
        elif (y < ylim[0] + 0.20 * h) and (x < xlim[0] + 0.35 * w):
            x += 0.015 * w

        nudged.append((float(np.clip(x, xlim[0] + 1e-6, xlim[1] - 1e-6)),
                       float(np.clip(y, ylim[0] + 1e-6, ylim[1] - 1e-6))))
    return nudged

def plot_nullclines(a=a, b=b, xlim=xlim, ylim=ylim,
                    n=250, density=0.6, alpha=0.35, nullcline_linestyle="solid",
                    show_arrows=True, show_nullcline_arrows=False):
    """
    Combine the classic nullclines figure (title + legend) with optional directional arrows.

    - Keeps the original title and legend from the classic plot.
    - Uses colored nullclines and component/resultant arrows.
    - Optionally shows direction arrows along the nullclines themselves.

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
    nullcline_linestyle : str
        Line style for nullclines: 'solid', 'dashed', 'dotted', or 'dashdot' (default 'solid')
    show_arrows : bool
        Whether to show the directional change vectors (default True)
    show_nullcline_arrows : bool
        Whether to show arrows along the nullclines indicating flow direction (default False)
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # classic panel with colored nullclines
    draw_nullclines_panel(
        ax, a=a, b=b, xlim=xlim, ylim=ylim,
        n=n, density=density, alpha=alpha,
        nullcline_colors=("red", "blue"),
        title=rf"Degradation of Glucose: Nullclines with $a={a:.2f}$ and $b={b:.2f}$",
        nullcline_linestyle=nullcline_linestyle
    )

    # Optionally show arrows on the nullclines
    if show_nullcline_arrows:
        draw_nullcline_arrows(ax, a=a, b=b, xlim=xlim, ylim=ylim, n_points=25, arrow_len=0.25, zorder=3.5)

    # Optionally show change vectors
    if show_arrows:
        pts = region_representative_points()
        pts = nudge_points(pts, a=a, b=b, xlim=xlim, ylim=ylim)
        draw_change_vectors(ax, pts, a=a, b=b,
                            color_x="red", color_y="blue", color_res="k",
                            len_comp=0.32, len_res=0.40, width=0.006, zorder=4)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.grid(True, alpha=0.35)

    fig_legend_dedup(fig, ax, loc="lower center", ncol=5,
                 framealpha=0.95, bbox_to_anchor=(0.5, 0.0),
                 frameon=True, edgecolor="black")
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    return fig



# -----------------------------------------------------------------------------
# Glucose system functions part 2
# -----------------------------------------------------------------------------

def jacobian(X, Y, a=a):
    """
    Compute the Jacobian matrix of the glucose model at point (X, Y).
    
    For the 2D system:
        dX/dt = -X + aY + X²Y
        dY/dt = b - aY - X²Y
    
    The Jacobian is:
        J = [[∂X'/∂X, ∂X'/∂Y],
             [∂Y'/∂X, ∂Y'/∂Y]]

        J(X,Y) = [[-1 + 2XY,  a + X²],
                  [-2XY,     -(a + X²)]]
    
    Parameters
    ----------
    X : float or ndarray
        X coordinate(s)
    Y : float or ndarray
        Y coordinate(s)
    a : float
        Enzyme-protein interaction parameter
    Returns
    -------
    J : ndarray
        Jacobian matrix (2×2) or array of 2×2 matrices
    """
    return np.array([
        [-1 + 2*X*Y,  a + X**2],
        [-2*X*Y,     -a - X**2]])

def classify_equilibrium(eigenvalues):
    """Classify equilibrium based on eigenvalues"""
    real_parts = np.real(eigenvalues)
    imag_parts = np.imag(eigenvalues)
    
    if np.all(np.abs(imag_parts) < 1e-10):  # Real eigenvalues
        if np.all(real_parts < 0):
            return "Stable node"
        elif np.all(real_parts > 0):
            return "Unstable node"
        else:
            return "Saddle point"
    else:  # Complex eigenvalues
        if real_parts[0] < 0:
            return "Stable spiral (focus)"
        elif real_parts[0] > 0:
            return "Unstable spiral (focus)"
        else:
            return "Center"

def get_equilibrium_eigvals(a=a, b=b):
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
    Xeq, Yeq = equilibrium(a=a, b=b)
    J_eq = jacobian(Xeq, Yeq, a=a)
    eig = np.linalg.eigvals(J_eq)
    return (Xeq, Yeq), J_eq, eig

def show_equilibrium_info(a=a, b=b):
    """
    Print equilibrium coordinates, Jacobian, and eigenvalues.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    """
    (Xeq, Yeq), J_eq, eig = get_equilibrium_eigvals(a=a, b=b)
    print(f"Equilibrium point: Xeq = {Xeq:.4f}, Yeq = {Yeq:.4f}\n")
    print("Jacobian at equilibrium:")
    print(J_eq)
    print(f"\n Eigenvalues: \n λ1 = {eig[0]:.4f} \n λ2 = {eig[1]:.4f}")
    
    # Classify equilibrium
    classification = classify_equilibrium(eig)
    print(f"\n Classification: {classification}")

#-----------------------------------------------------------------------------
# Glucose system functions part 3: Trajectory plotting
#-----------------------------------------------------------------------------


def trajectory_initial_conditions(a=a, b=b, xlim=xlim, ylim=ylim, step=0.08):
    """
    Generate sets of initial conditions for trajectory integration.
    
    Creates two groups of initial conditions:
    - Near equilibrium: for studying local stability
    - Near edges: for studying global dynamics
    
    Parameters
    ----------
    a : float
        interaction parameter
    b : float
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
    Xeq, Yeq = equilibrium(a=a, b=b)
    dx = step * (xlim[1] - xlim[0])
    dy = step * (ylim[1] - ylim[0])

    # near equilibrium
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
    ]

    return near_eq , near_edge

def draw_trajectories(ax, a=a, b=b, initials=None, t_span=(0,60), t_eval_n=2500,
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

    # Base styles for smoother, cleaner trajectories
    base_line = {
        "linewidth": 2.2,
        "alpha": 0.92,
        "antialiased": True,
        "zorder": 2,
        "solid_capstyle": "round",
        "solid_joinstyle": "round",
        "path_effects": [
            pe.Stroke(linewidth=3.2, foreground="white", alpha=0.6),
            pe.Normal(),
        ],
    }
    base_start = {
        "s": 18,
        "linewidth": 0.8,
        "zorder": 3,
    }

    # Use a pleasant colormap across trajectories (avoid extremes)
    cmap = plt.get_cmap("viridis")
    m = max(len(initials), 1)

    for i, (X0, Y0) in enumerate(initials):
        color = cmap(0.15 + 0.7 * (i / max(m - 1, 1)))
        sol = solve_ivp(
            lambda t, z: glucose_rhs(z[0], z[1], a, b),
            t_span, [X0, Y0], t_eval=t_eval,
            rtol=1e-6, atol=1e-9
        )

        # Resample for extra-smooth rendering
        t_dense = np.linspace(t_eval[0], t_eval[-1], max(t_eval_n, 4000))
        fx = interp1d(t_eval, sol.y[0], kind="cubic")
        fy = interp1d(t_eval, sol.y[1], kind="cubic")
        x_smooth = fx(t_dense)
        y_smooth = fy(t_dense)

        style_line = {**base_line, "color": color, **line_kw}
        ax.plot(x_smooth, y_smooth, **style_line)

        style_start = {**base_start, "edgecolor": "white", "facecolor": color}
        style_start.update(start_kw)
        ax.scatter([X0], [Y0], **style_start)

def plot_2glucose_trajectories(a=a, b=b, xlim=xlim, ylim=ylim,
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

    near_eq, near_edge = trajectory_initial_conditions(a=a, b=b, xlim=xlim, ylim=ylim, step=step)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True, sharey=True)

    # --- Local trajectories panel ---
    ax = axes[0]
    draw_nullclines_panel(ax, a=a, b=b, xlim=xlim, ylim=ylim, n=n,
                          density=0.4, alpha=0.30, show_equilibrium=True,
                          title="Local trajectories near equilibrium",
                          nullcline_linestyle="dashed")
    draw_trajectories(ax, a=a, b=b, initials=near_eq, t_span=t_span, t_eval_n=t_eval_n)

    # --- Global trajectories panel ---
    ax = axes[1]
    draw_nullclines_panel(ax, a=a, b=b, xlim=xlim, ylim=ylim, n=n,
                          density=0.4, alpha=0.30, show_equilibrium=True,
                          title="Global trajectories further from equilibrium",
                          nullcline_linestyle="dashed")
    draw_trajectories(ax, a=a, b=b, initials=near_edge, t_span=t_span, t_eval_n=t_eval_n)

   # specify legend at figure level
    fig_legend_dedup(fig, list(axes), loc="lower center", ncol=5,
                 framealpha=0.95, bbox_to_anchor=(0.5, 0.0),
                 frameon=True, edgecolor="black")

    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    plt.show()

def plot_zoomed_spiral_convergence(a=a, b=b, xlim=xlim, ylim=ylim,
                                    t_span=(0, 60), t_eval_n=2500, step=0.08,
                                    single_init=(0.465, 0.6)):
    """
    Plot side-by-side visualization of trajectories spiraling from unstable equilibrium to limit cycle.
    
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
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    ax1, ax2 = axes
    
    # --- Left panel: Multiple trajectories ---
    draw_nullclines_panel(
        ax1, a=a, b=b, xlim=xlim, ylim=ylim,
        n=250, density=0.6, alpha=0.35,
        title="Multiple trajectories spiraling away from unstable equilibrium", show_equilibrium=True,
        nullcline_linestyle="dashed"
    )
    near_eq, _ = trajectory_initial_conditions(a=a, b=b, xlim=xlim, ylim=ylim, step=step)
    draw_trajectories(ax1, a=a, b=b, initials=near_eq, 
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
        title="Single trajectory", show_equilibrium=True,
        nullcline_linestyle="dashed"
    )
    draw_trajectories(ax2, a=a, b=b, initials=[single_init], 
                      t_span=t_span, t_eval_n=t_eval_n,
                      line_kw=dict(linewidth=2.0, alpha=0.85, color='red'))
    ax2.plot(single_init[0], single_init[1], 'go', markersize=10, label='Start points', zorder=5)
    ax2.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # Main title
    fig.suptitle('Zoomed view: Trajectories spiral from unstable equilibrium to limit cycle (a={:.2f}, b={:.2f})'.format(a, b), 
                 fontsize=13, fontweight='bold', y=0.95)
    
  
    for ax in (axes):
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()
    
    fig_legend_dedup(fig, [ax1, ax2], loc="upper center",bbox_to_anchor=(0.5, 0.05), frameon=True,
                     ncol=4,framealpha=0.95,fontsize=10, edgecolor='black')

    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    plt.show()

# -----------------------------------------------------------------------------
# Glucose system functions part 4: Bifurcation Analysis
# -----------------------------------------------------------------------------

def trace_det_at_eq(a=a, b=b):
    """
    Compute trace and determinant of Jacobian at equilibrium.
    
    Calculates the Jacobian matrix at the equilibrium point and computes
    its trace (sum of diagonal elements) and determinant.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter
    b : float
        Protein production rate
    
    Returns
    -------
    trace : float
        Trace of Jacobian (stability indicator; bifurcation when trace=0)
    det : float
        Determinant of Jacobian
    """
    Xeq, Yeq = equilibrium(a=a, b=b)
    J_eq = jacobian(Xeq, Yeq, a=a)
    return np.trace(J_eq), np.linalg.det(J_eq)

def _max_real_part(trace, det):
    """
    Vectorized maximum real part of eigenvalues of a 2x2 matrix given trace and determinant.

    For eigenvalues solving λ^2 - trace*λ + det = 0:
    - If discriminant D = trace^2 - 4*det >= 0 (real eigenvalues), max real = (trace + sqrt(D))/2
    - If D < 0 (complex pair), real parts are equal trace/2, so max real = trace/2
    """
    trace = np.asarray(trace)
    det = np.asarray(det)
    disc = trace**2 - 4.0*det # Discriminant
    sqrt_disc = np.sqrt(np.maximum(disc, 0.0))
    # this is more efficient than np.linalg.eigvals
    # to calculate the eigenvalues for large arrays.
    return np.where(disc >= 0.0, 0.5*(trace + sqrt_disc), 0.5*trace)

def trace_det_from_ab(a, b_vals):
    """
    Vectorized computation of trace and determinant at equilibrium for multiple b values.
    
    Much faster than calling trace_det_at_eq() in a loop.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter
    b_vals : array-like
        Array of b values
    
    Returns
    -------
    trace : ndarray
        Trace of Jacobian at equilibrium for each b
    det : ndarray
        Determinant of Jacobian at equilibrium for each b
    """
    b_vals = np.asarray(b_vals, dtype=float)
    s = a + b_vals*b_vals          # s = a + b²
    trace = 1.0 - s - (2.0*a)/s    # trace at equilibrium
    det = s                         # det at equilibrium
    return trace, det

def stable_mask_from_ab(a, b_vals):
    """
    Vectorized computation of stability mask for multiple b values.
    
    Returns True where equilibrium is stable (max real part of eigenvalues < 0).
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter
    b_vals : array-like
        Array of b values
    
    Returns
    -------
    stable : ndarray of bool
        True where equilibrium is stable, False where unstable
    """
    trace, det = trace_det_from_ab(a, b_vals)
    max_real = _max_real_part(trace, det)
    return max_real < 0

def classify_from_trace_det(trace, det):
    disc = trace**2 - 4.0*det
    if det < 0:
        return -1, "Saddle"
    stable = trace < 0
    if disc < 0:
        return (1 if stable else 0), ("Stable focus" if stable else "Unstable focus")
    else:
        return (1 if stable else 0), ("Stable node" if stable else "Unstable node")

def bcrit_values(a=a):
    """
    Calculate critical b-values where Hopf bifurcations occur (trace(J_eq)=0).
    
    Uses vectorized trace computation to find where trace changes sign.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter
    
    Returns
    -------
    list of float
        Critical b-values where bifurcations occur
    """
    b_m, b_p = hopf_curves(np.array([a])) # the function hopf_curves stands further down this file
    out = [b_m[0], b_p[0]]
    return [float(x) for x in out if np.isfinite(x) and x >= 0]

def plot_equilibrium_vs_b(a=a, b_min=0.0, b_max=1.2, n=2000):
    """
    Plot both X_eq and Y_eq vs b with stability coloring in a single figure.
    
    Left panel: X_eq = b (linear relationship)
    Right panel: Y_eq = b/(a + b²) (nonlinear relationship)
    
    Both are colored by stability (green = stable, orange = unstable).
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter (kept fixed)
    b_min, b_max : float
        Range of b values to plot
    n : int
        Number of points
    """
    b_vals = np.linspace(b_min, b_max, n)
    
    # Compute equilibrium positions
    Xeq, Yeq = equilibrium(a, b_vals) 
    # Compute stability (vectorized)
    stable = stable_mask_from_ab(a, b_vals)

   
    b_crit = bcrit_values(a)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharex=True)

    def panel(ax, yvals, ylabel, y_at_bcrit):
        ax.scatter(b_vals[stable],  yvals[stable],  s=1, color="tab:green",  label="stable equilibrium")
        ax.scatter(b_vals[~stable], yvals[~stable], s=1, color="tab:orange", label="unstable equilibrium")

        for i, bc in enumerate(b_crit):
            ax.axvline(bc, ls="--", color="tab:red", label=r"$b_{crit}$" if i == 0 else None)
            ax.scatter(bc, y_at_bcrit(bc), color="tab:red", zorder=3, s=50, label='Bifurcation points' if i == 0 else None)

        ax.set(xlabel="b", ylabel=ylabel)
        ax.grid(True, alpha=0.4)
    
    # --- Left panel: X_eq vs b ---
    panel(
        axes[0],
        Xeq,
        r"$X_{eq}$",
        y_at_bcrit=lambda bc: bc)

    # --- Right panel: Y_eq vs b ---
    panel(
        axes[1],
        Yeq,
        r"$Y_{eq}$",
        y_at_bcrit=lambda bc: bc / (a + bc**2))

    # --- Shared legend at bottom ---
    fig_legend_dedup(fig, list(axes), loc="lower center", ncol=5,
                 framealpha=0.95, bbox_to_anchor=(0.5, 0.0),
                 frameon=True, edgecolor="black")
    
    fig.suptitle(rf"Equilibrium $X_{{eq}}$ and $Y_{{eq}}$ for variable $b$ and fixed $a={a}$", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0.06, 1, 0.96])
    plt.show()

def plot_equilibrium_3d(a=a, b_min=0.0, b_max=1.2, n=2000):
    """
    Create a 3D visualization of equilibria as a function of b.
    
    Shows X_eq and Y_eq varying with b parameter, with stability coloring.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter (kept fixed)
    b_min, b_max : float
        Range of b values to plot
    n : int
        Number of points
    """
    
    # Compute parameter grid
    b_vals = np.linspace(b_min, b_max, n)
    
    # Compute equilibria
    Xeq, Yeq = equilibrium(a, b_vals)
    
    # Determine stability
    stable = stable_mask_from_ab(a, b_vals)
    b_crit = bcrit_values(a)
    
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot stable points (green)
    if np.any(stable):
        ax.scatter(b_vals[stable], Xeq[stable], Yeq[stable], 
                   c='tab:green', s=3, label='Stable equilibrium', alpha=0.8)
    
    # Plot unstable points (orange)
    if np.any(~stable):
        ax.scatter(b_vals[~stable], Xeq[~stable], Yeq[~stable], 
                   c='tab:orange', s=3, label='Unstable equilibrium', alpha=0.8)
    
    # Mark bifurcation points (red)
    if b_crit:
        bcrit_Xeq = np.array([b for b in b_crit])
        bcrit_Yeq = np.array([b / (a + b**2) for b in b_crit])
        ax.scatter(b_crit, bcrit_Xeq, bcrit_Yeq, 
                   c='tab:red', s=30, label='Bifurcation points', zorder=50)
    
    # Draw lines connecting the equilibrium curve
    ax.plot(b_vals, Xeq, Yeq, 'k-', alpha=0.3, linewidth=0.5, zorder=1)
    
    ax.set_xlabel('b', fontsize=11, fontweight='bold')
    ax.set_ylabel('$X_{eq}$', fontsize=11, fontweight='bold')
    ax.set_zlabel('$Y_{eq}$', fontsize=11, fontweight='bold', labelpad=8)
    ax.set_title(f'3D Equilibrium with fixed a={a}', fontsize=13, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3)
    
    # Legend at bottom
    fig_legend_dedup(fig, ax, loc="lower center", ncol=5,
                 framealpha=0.95, bbox_to_anchor=(0.5, 0.0),
                 frameon=True, edgecolor="black")
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    plt.show()

def plot_trace_and_determinant(a=a, b_min=0.0, b_max=1.2, n=2000):
    """
    Plot trace and determinant of Jacobian at equilibrium vs b.
    
    Combines trace and determinant visualization in a single figure with two panels.
    
    Parameters
    ----------
    a : float
        Enzyme-protein interaction parameter (kept fixed)
    b_min, b_max : float
        Range of b values to plot
    n : int
        Number of points
    """
    b_vals = np.linspace(b_min, b_max, n)
    b_crit = bcrit_values(a)
    
    # Compute trace and determinant (vectorized)
    trace_vals, det_vals = trace_det_from_ab(a, b_vals)
    stable = stable_mask_from_ab(a, b_vals)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharex=True)
    
    # Panel 1: trace vs b
    axes[0].scatter(b_vals[stable],  trace_vals[stable],  s=1, color="tab:green",  label="stable")
    axes[0].scatter(b_vals[~stable], trace_vals[~stable], s=1, color="tab:orange", label="unstable")
    axes[0].axhline(0, color="k", lw=1.5, linestyle='-', alpha=0.5)
    for i, bc in enumerate(b_crit):
        axes[0].axvline(bc, ls="--", color="tab:red", label=r"$b_{crit}$" if i == 0 else None)
        axes[0].scatter(bc, 0, color="tab:red", zorder=3, s=50, label='Bifurcation points' if i == 0 else None)
    axes[0].set(xlabel="b", ylabel=r"$\mathrm{trace}(J)$", title=r"$\mathrm{trace}(J_{eq})$ for parameter $b$")
    axes[0].grid(True, alpha=0.4)
    
    # Panel 2: determinant vs b
    axes[1].scatter(b_vals[stable],  det_vals[stable],  s=1, color="tab:green",  label="stable")
    axes[1].scatter(b_vals[~stable], det_vals[~stable], s=1, color="tab:orange", label="unstable")
    axes[1].axhline(0, color="k", lw=1.5, linestyle='-', alpha=0.5)
    for i, bc in enumerate(b_crit):
        axes[1].axvline(bc, ls="--", color="tab:red", label=r"$b_{crit}$" if i == 0 else None)
        axes[1].scatter(bc, bc**2 + a, color="tab:red", zorder=3, s=50)
    axes[1].set(xlabel="b", ylabel=r"$\det(J)$", title=r"$\det(J_{eq})$ for parameter $b$")
    axes[1].grid(True, alpha=0.4)
    
    fig.suptitle(rf"Jacobian properties at equilibrium with a={a} fixed", fontsize=13, fontweight="bold")
    
    # Shared legend at the bottom
    fig_legend_dedup(fig, list(axes), loc="lower center", ncol=5,
                 framealpha=0.95, bbox_to_anchor=(0.5, 0.0),
                 frameon=True, edgecolor="black")

    fig.tight_layout(rect=[0, 0.06, 1, 0.96])
    plt.show()

def plot_bifurcation_summary_figures(a, b_range,
                                     zoom_xlim=(0, 1.5),zoom_ylim=(0, 3),
                                     n_b=2000, eps=0.03):
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
    b_vals = np.linspace(b_range[0], b_range[1], n_b)
    trace, det = trace_det_from_ab(a, b_vals)
    max_real = _max_real_part(trace, det)

    stable = max_real < 0  # most direct definition
    b_crit = bcrit_values(a)

    fig = plt.figure(figsize=(12, 10))

    # --- Top: stability curve ---
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(b_vals, max_real, label="max Re(λ) at equilibrium")
    ax1.axhline(0, color="k", lw=1)
    for i, bc in enumerate(b_crit):
        ax1.axvline(bc, color="tab:red", label="b_crit" if i == 0 else None)

    ax1.fill_between(b_vals, max_real, 0, where=stable,  alpha=0.15, label="stable")
    ax1.fill_between(b_vals, max_real, 0, where=~stable, alpha=0.12, label="unstable")

    ax1.set(xlim=b_range, xlabel="b", ylabel="max Re(λ)", title=f"Stability for parameter b and a={a} fixed")
    ax1.grid(True, alpha=0.3)
    # Don't add legend here - we'll create a shared one at the bottom

    # --- Bottom: phase portraits around first bifurcation ---
    axes = [fig.add_subplot(2, 3, k) for k in (4, 5, 6)]
    if b_crit:
        b0 = b_crit[0]
        bs = [np.clip(b0 - eps, *b_range), b0, np.clip(b0 + eps, *b_range)]
        titles = [f"before bifurcation b={bs[0]:.3f}", f"at bifurcation b={bs[1]:.3f}", f"after bifurcation b={bs[2]:.3f}"]

        for axp, bb, tt in zip(axes, bs, titles):
            draw_nullclines_panel(axp, a=a, b=bb, xlim=(0.0,0.75), ylim=(1,2.5),
                                  n=150, density=0.6, alpha=0.1,
                                  title=tt, show_equilibrium=True,
                                  nullcline_linestyle="dashed")
            near_eq, _ = trajectory_initial_conditions(a, bb, xlim=zoom_xlim, ylim=zoom_ylim, step=0.1)
            draw_trajectories(axp, a, bb, initials=near_eq[:1],
                              t_span=(0, 150), t_eval_n=5000,
                              line_kw=dict(linewidth=0.8, alpha=0.85))
            axp.plot(near_eq[0][0], near_eq[0][1], 'go', markersize=4, label='Start point', zorder=5)    
    else:
        axes[1].text(0.5, 0.5, "No bifurcation for this a", ha="center", va="center")
        for axp in axes:
            axp.axis("off")

    fig.suptitle(f"Bifurcation summary (a={a})", fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0.04, 1, 0.99])
    
    # --- Shared legend at the bottom ---
    fig_legend_dedup(fig, list(axes), loc="lower center", ncol=5,
                 framealpha=0.95, bbox_to_anchor=(0.5, 0.0),
                 frameon=True, edgecolor="black")

    plt.show()

def plot_bifurcation_phase_portraits(a=a, b_crit_vals=None, xlim=xlim, ylim=ylim):
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
    
    axes_flat = np.ravel(axes)

    for row, b_crit in enumerate(b_crit_vals):
        epsilon = 0.02  # Small perturbation around bifurcation
        b_vals = [b_crit - epsilon, b_crit, b_crit + epsilon]
        titles = [f"Before bifurcation b={b_crit-epsilon:.4f}", 
                  f"At bifurcation b={b_crit:.4f}", 
                  f"After bifurcation b={b_crit+epsilon:.4f}"]
        
        for col, (b_val, title) in enumerate(zip(b_vals, titles)):
            ax = axes[row, col]
            
            # Draw nullclines and vector field
            draw_nullclines_panel(ax, a=a, b=b_val, xlim=xlim, ylim=ylim,
                                 n=250, density=0.8, alpha=0.1,
                                 title=title, show_equilibrium=True,
                                 nullcline_linestyle="dashed")
            
            # Add a single trajectory
            near_eq, _ = trajectory_initial_conditions(a=a, b=b_val, xlim=xlim, ylim=ylim, step=0.1)
            draw_trajectories(ax, a=a, b=b_val, initials=near_eq[:1],
                            t_span=(0, 150), t_eval_n=5000,
                            line_kw=dict(linewidth=0.7, zorder=2.5, alpha=0.9))
            # Mark the start point
            ax.plot(near_eq[0][0], near_eq[0][1], 'go', markersize=8, label='Start point', zorder=5)

    fig_legend_dedup(fig, list(axes_flat), loc="lower center", ncol=5,
                 framealpha=0.95, bbox_to_anchor=(0.5, 0.0),
                 frameon=True, edgecolor="black")

    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    plt.show()


# ==============================================================================
# Glucose system functions part 5: 2D Bifurcation Analysis
# ==============================================================================

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
                     na=250, nb=250):
    """
    Assignment-style 2D stability map (a,b-plane).

    - Background: stable (trace<0) vs unstable (trace>0) regions using the trace criterion
    - Overlays: Hopf curves b_-(a), b_+(a)
    - Legend: English labels for stable/unstable and the two Hopf branches

    Returns
    -------
    Stability map figure
    """
    a_vals = np.linspace(a_min, a_max, na)
    b_vals = np.linspace(b_min, b_max, nb)
    A, B = np.meshgrid(a_vals, b_vals, indexing="xy")

    # Vectorized stability metrics at equilibrium
    S = A + B*B
    trace = 1.0 - S - (2.0*A)/S

    # Analytic max real part of eigenvalues from trace/determinant
    stable = trace < 0

    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    cmap = ListedColormap(["#fde2b8", "#c9e9c5"])  # unstable, stable
    ax.pcolormesh(a_vals, b_vals, stable.astype(int), shading="auto", cmap=cmap, vmin=0, vmax=1)

    # Hopf curves
    b_m, b_p = hopf_curves(a_vals)
    line_bm, = ax.plot(a_vals, b_m, color='#1f77b4', linewidth=2.0, linestyle='--', label=r'Hopf $b_-(a)$')
    line_bp, = ax.plot(a_vals, b_p, color='#d62728', linewidth=2.0, linestyle=':',  label=r'Hopf $b_+(a)$')

    ax.set_xlim(a_min, a_max)
    ax.set_ylim(b_min, b_max)
    ax.set_xlabel('a')
    ax.set_ylabel('b')
    ax.set_title(r'Stability of equilibrium in $ab$-plane')
    ax.grid(True, alpha=0.25)
    
    
    stable_patch = Patch(facecolor="#c9e9c5", edgecolor='none', label='Stable equilibrium')
    unstable_patch = Patch(facecolor="#fde2b8", edgecolor='none', label='Unstable equilibrium')
    handles = [stable_patch, unstable_patch, line_bm, line_bp]
    fig.legend(handles=handles, loc='lower center', ncol=4, framealpha=0.95, 
                bbox_to_anchor=(0.5, 0.0), frameon=True, edgecolor='black').get_frame().set_linewidth(0.5)
    plt.tight_layout(rect=[0, 0.06, 1, 0.96])
    
    plt.show()



# ==============================================================================
# interactive functions for notebook
# ==============================================================================

def interactive_nullclines(xlim=xlim, ylim=ylim):
    """
    Create an interactive nullclines plot with sliders for parameters.
    
    This function provides an interactive widget-based visualization where users
    can dynamically adjust model parameters (a, b), plot resolution (n), 
    streamline density, transparency (alpha), and toggle various arrow displays.
    
    Parameters
    ----------
    xlim : tuple
        X domain limits (default from module settings)
    ylim : tuple
        Y domain limits (default from module settings)
    
    Notes
    -----
    Requires ipywidgets and the personalized_layout module.
    Best used in Jupyter notebooks with interactive widget support.
    
    Example
    -------
    >>> interactive_nullclines()
    # Display interactive plot with sliders for a, b, n, density, alpha, 
    # and checkboxes for show_arrows and show_nullcline_arrows
    """
    
    # Define the interactive plotting function
    def _plot(a, b, n, density, alpha, show_arrows, show_nullcline_arrows):
        plot_nullclines(a=a, b=b, xlim=xlim, ylim=ylim,
                       n=int(n), density=density, alpha=alpha, 
                       show_arrows=show_arrows, show_nullcline_arrows=show_nullcline_arrows)
    
    # Configure sliders for continuous parameters
    slider_configs = {
        'a': {'value': 0.06, 'min': 0.01, 'max': 0.14, 'step': 0.01},
        'b': {'value': 0.4, 'min': 0.0, 'max': 1.2, 'step': 0.05},
        'n': {'value': 250, 'min': 50, 'max': 500, 'step': 50},
        'density': {'value': 0.6, 'min': 0.2, 'max': 1.5, 'step': 0.1},
        'alpha': {'value': 0.35, 'min': 0.1, 'max': 1.0, 'step': 0.05}
    }
    
    # Create checkboxes for arrow display options
    show_arrows_checkbox = Checkbox(
        value=True,
        description='Show change vectors',
        indent=False)
    
    show_nullcline_arrows_checkbox = Checkbox(
        value=False,
        description='Show nullcline arrows',
        indent=False)
    
    # Use centralized helper with extra widgets
    extra_widgets = {
        'show_arrows': show_arrows_checkbox,
        'show_nullcline_arrows': show_nullcline_arrows_checkbox
    }
    create_interactive_plot(_plot, slider_configs, n_cols=2, extra_widgets=extra_widgets, extra_widgets_inline=True)


