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

import sys
import os

# Add parent directory to path to import personalized_layout
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import importlib
import rossler_functions as rossler

import numpy as np 

import sympy as sp
from sympy.solvers import solve
from sympy import Symbol, Matrix, symbols, Eq
    
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D 

import scipy

# Import helper functions from centralized module
from personalized_layout import create_interactive_plot, legend_dedup, fig_legend_dedup


# =============================================================================
# General Parameters
# =============================================================================


a, b, c, t_min, t_max, dt, initial_condition, initial_conditions = rossler.get_parameters()


# ============================================================
# Figures
# ============================================================

# ====================================
# Visualization: 2D Projections
# ====================================


def plot_2d_projection(ax, t, x, y, z, plane='XY', title=None,
                       alpha=0.7, linewidth=1.0, color='blue',
                       skip_first_frac=0.2):
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
    # Skip transient
    skip = int(len(t) * skip_first_frac)

    if plane == 'XY':
        ax.plot(x[skip:], y[skip:], color=color, alpha=alpha, linewidth=linewidth)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
    elif plane == 'XZ':
        ax.plot(x[skip:], z[skip:], color=color, alpha=alpha, linewidth=linewidth)
        ax.set_xlabel('X')
        ax.set_ylabel('Z')
    elif plane == 'YZ':
        ax.plot(y[skip:], z[skip:], color=color, alpha=alpha, linewidth=linewidth)
        ax.set_xlabel('Y')
        ax.set_ylabel('Z')
    
    if title:
        ax.set_title(title)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)


# ====================================
# Visualization: 3D Projections
# ====================================


def plot_3d_trajectory(ax, t, x, y, z, *, c, a=a, b=b, initial_condition=None,
                      title=None, alpha=0.8, linewidth=1.0, color='blue',
                      skip_first_frac=0.2, show_legend=True):
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
    show_legend : bool
        Whether to show legend on this axis
    """
    # Skip transient behavior
    skip_idx = int(len(t) * skip_first_frac)
    
    # Create label for trajectory
    if initial_condition is not None:
        traj_label = f"Trajectory IC=({initial_condition[0]:.2f}, {initial_condition[1]:.2f}, {initial_condition[2]:.2f})"
    else:
        traj_label = "Trajectory"
    
    ax.plot(x[skip_idx:], y[skip_idx:], z[skip_idx:],
            color=color, alpha=alpha, linewidth=linewidth,
            label=traj_label)
    
    if title:
        ax.set_title(title)

    # --- Show equilibrium ---
    eq_needed = rossler.equilibria_at_c(c, a_val=a, b_val=b)
    eq_styles = rossler.styles_for_equilibria()
    
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
    # Increase legend line width for better visibility
    if show_legend:
        legend = ax.legend(loc='best', fontsize=6)
        for line in legend.get_lines():
            line.set_linewidth(2.0)


# ===================================
# Visualization: Multiple Panels
# ===================================


def plot_3d_attractor(a=0.2, b=0.2, c=5.7, initial_condition=initial_condition,
                    t_min=t_min, t_max=t_max, dt=dt, figsize=(12, 5)):
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

    t, x, y, z = rossler.get_RK4_vectors(a=a, b=b, c=c, initial_condition=initial_condition,
                                t_min=t_min, t_max=t_max, dt=dt)
    
    fig = plt.figure(figsize=figsize)
    
    # 3D plot
    ax1 = fig.add_subplot(121, projection='3d')
    plot_3d_trajectory(ax1, t, x, y, z, c=c, a=a, b=b, initial_condition=initial_condition,
                      title='3D Rössler Attractor', skip_first_frac=0.2, color='darkblue', 
                      linewidth=0.5, show_legend=False)
    
    # 2D projection
    ax2 = fig.add_subplot(122)
    plot_2d_projection(ax2, t, x, y, z, plane='XY', title='x-y Projection')
    
    # Collect legend from 3D plot and display below
    fig_legend_dedup(fig, ax1, loc='lower center', ncol=4, bbox_to_anchor=(0.5, 0.1), 
                     frameon=True, fontsize=8, edgecolor='black')
    
    plt.tight_layout(rect=[0, 0.08, 1, 1]) 
    return fig

def plot_projections_2D(a=0.2, b=0.2, c=5.7, initial_condition=initial_condition,
                     t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                     rtol=1e-8, atol=1e-10, figsize=(15, 5), skip_first_frac=0.2):
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
    t, x, y, z = rossler.get_RK4_vectors(a=a, b=b, c=c, initial_condition=initial_condition,
                                t_min=t_min, t_max=t_max, dt=dt)
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    plot_2d_projection(axes[0], t, x, y, z, plane='XY', title='x-y plane (top view)', skip_first_frac=skip_first_frac)
    plot_2d_projection(axes[1], t, x, y, z, plane='XZ', title='x-z plane (side view)', skip_first_frac=skip_first_frac)
    plot_2d_projection(axes[2], t, x, y, z, plane='YZ', title='y-z plane', skip_first_frac=skip_first_frac)
    
    plt.tight_layout()
    return fig

def plot_time_series(a=0.2, b=0.2, c=5.7, initial_condition=initial_condition,
                    t_min=t_min, t_max=t_max, dt=dt, figsize=(14, 8)):
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
    t, x, y, z = rossler.get_RK4_vectors(a=a, b=b, c=c, initial_condition=initial_condition,
                                t_min=t_min, t_max=t_max, dt=dt)
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


# ============================================================
# Comparison and Analysis
# ============================================================


def compare_initial_conditions(a=0.2, b=0.2, c=5.7, initial_conditions= None,
                    t_min=t_min, t_max=t_max, dt=dt, figsize=(15, 5), skip_first_frac=0.2):
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
        t, x, y, z = rossler.get_RK4_vectors(a=a, b=b, c=c, initial_condition=(x0, y0, z0),
                    t_min=t_min, t_max=t_max, dt=dt)

        color = colors[ic_idx % len(colors)]
        label = f'IC: ({x0}, {y0}, {z0})'
        
        plot_2d_projection(axes[0], t, x, y, z, plane='XY', color=color,
                           title='x-y Plane', skip_first_frac=skip_first_frac)
        plot_2d_projection(axes[1], t, x, y, z, plane='XZ', color=color,
                           title='x-z Plane', skip_first_frac=skip_first_frac)
        plot_2d_projection(axes[2], t, x, y, z, plane='YZ', color=color,
                           title='y-z Plane', skip_first_frac=skip_first_frac)
        
        axes[0].plot([], [], color=color, label=label, linewidth=2)
    
    # Place a single shared legend below the axes
    fig_legend_dedup(fig, axes[0], loc='upper center', bbox_to_anchor=(0.5, 0.1),
                     ncol=min(len(initial_conditions), 4), frameon=True, edgecolor='black')
    
    plt.tight_layout(rect=[0.0, 0.08, 1.0, 1.0])
    return fig


# ============================================================
#  Analysis of Chaotic Characteristics
# ============================================================


# --- Sensitivity to initial conditions (Butterfly effect) ---
def show_butterfly_effect(a, b, c, initial_condition, delta=1e-6, t_min=0, t_max=100, dt=0.01):
    """
    Illustrate the butterfly effect in the Rössler system.
    
    We integrate two trajectories with almost identical
    initial conditions and compare X(t), Y(t), Z(t) in time.
    """
    plt.close('all')
    
    # Two initial conditions that differ by delta in X
    base_ic = initial_condition
    ic1 = base_ic
    ic2 = (base_ic[0] + delta, base_ic[1], base_ic[2])
    
    # Integrate both trajectories
    t, x1, y1, z1 = rossler.get_RK4_vectors(a=a, b=b, c=c, t_min=t_min, t_max=t_max, dt=dt,
                                initial_condition=ic1)
    _, x2, y2, z2 = rossler.get_RK4_vectors(a=a, b=b, c=c, t_min=t_min, t_max=t_max, dt=dt,
                                initial_condition=ic2)
    
    series1 = {"X": x1, "Y": y1, "Z": z1}
    series2 = {"X": x2, "Y": y2, "Z": z2}
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 8))
    
    for (ax, comp) in zip(axes, ("X", "Y", "Z")):
        ax.plot(t, series1[comp],
                label=rf'{comp}(t), initial condition',
                alpha=0.8)
        ax.plot(t, series2[comp],
                label=rf'{comp}(t), initial condition + $\delta$',
                alpha=0.8, linestyle='--')
    
        ax.set_ylabel(rf'{comp}(t)')
        ax.grid(True)
    
    axes[-1].set_xlabel('t')
    fig.suptitle(rf'Butterfly effect in Rössler system for $c={c:.2f}$, $\delta={delta:.1e}$')
    
    # Shared legend below all panels
    fig_legend_dedup(fig, axes[0], loc='lower center', bbox_to_anchor=(0.5, 0.06), ncol=2, 
                     frameon=True, fontsize=9, edgecolor='black')
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.show()  

# --- Sensitivity to initial conditions (Lyapunov flavour) ---
def plot_sensitivity_to_initial_conditions(c=5.7, t_min=t_min, t_max=t_max, dt=dt, delta=1e-6):
    """
    Show sensitivity to initial conditions in the Rössler system.

    For a given c, we take two initial conditions that differ by only `delta`
    in the X-component, integrate both trajectories, and plot how their
    distance in state space evolves in time (linear and log scale).
    """

    # Choose reference initial condition from the general parameters
    base_ic = initial_condition
    ic1 = base_ic
    ic2 = (base_ic[0] + delta, base_ic[1], base_ic[2])

    t, x1, y1, z1 = rossler.get_RK4_vectors(a=a, b=b, c=c, t_min=t_min, t_max=t_max, dt=dt,
                                initial_condition=ic1)
    _, x2, y2, z2 = rossler.get_RK4_vectors(a=a, b=b, c=c, t_min=t_min, t_max=t_max, dt=dt,
                                initial_condition=ic2)
    # Euclidean distance in (X, Y, Z)
    d = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

    plt.close('all')
    fig, axes = plt.subplots(2, 1, figsize=(8, 8))

    # Upper: distance on linear scale
    ax = axes[0]
    ax.plot(t, d)
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



    # =================================
    #   Return map: maxima of Z(t)
    # =================================

# --- return map: Z_n versus Z_{n+1} ---
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
    t_max_vals, Z_max, t, z = rossler.get_Z_maxima(a, b, c, initial_condition, t_min=t_min, t_max=t_max, dt=dt)

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
    _, Z_max, _, _ = rossler.get_Z_maxima(a=a, b=b, c=c, initial_condition=initial_condition,
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
                cob_lw: float = 0.8, mark_start: bool = True, show_legend: bool = True):
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
    show_legend : bool, optional
        If True, display legend on this axis. If False, legend entries are still
        added but legend is not shown (useful for shared legends).
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
    if show_legend:
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
    F = rossler.rossler_return_map_from_data(Z_use)

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
# Parameter Study
# ============================================================


def plot_parameter_sweep(a=0.1, b=0.1, c_values=None,
                        initial_condition=None,
                        t_min=0, t_max=200, dt=0.01,
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
    
    if initial_condition is None:
        initial_condition = (0.1, 0.1, 0.1)
    
    n_c = len(c_values)
    n_cols = 3
    n_rows = (n_c + n_cols - 1) // n_cols
    
    fig = plt.figure(figsize=figsize)
    
    # Collect legend info from first plot only
    first_ax = None
    
    for idx, c in enumerate(c_values, 1):
        ax = fig.add_subplot(n_rows, n_cols, idx, projection='3d')
        
        t, x, y, z = rossler.get_RK4_vectors(a=a, b=b, c=c, initial_condition=initial_condition,
                    t_min=t_min, t_max=t_max, dt=dt)
        
        # Only show legend on first subplot to collect handles/labels
        show_legend_on_axis = (idx == 1)
        
        plot_3d_trajectory(ax, t, x, y, z, c=c, a=a, b=b,
                  title=f'a = {a}, b = {b}, c = {c}',
                  skip_first_frac=0.2, color='darkblue', linewidth=0.3,
                  initial_condition=initial_condition,
                  show_legend=show_legend_on_axis)
        
        if idx == 1:
            first_ax = ax
    
    # Create shared legend below the figure
    if first_ax is not None:
        # Remove the temporary legend from first subplot
        if first_ax.get_legend():
            first_ax.get_legend().remove()
        # Add shared deduplicated legend below with increased line width for trajectory
        legend = fig_legend_dedup(fig, first_ax, loc='upper center', bbox_to_anchor=(0.5, 0.02),
                      ncol=4, frameon=True, fontsize=9, edgecolor='black')
        # Make trajectory line thicker in legend for better visibility
        if legend:
            for line in legend.get_lines():
                if line.get_label() and 'Trajectory' in line.get_label():
                    line.set_linewidth(2.0)
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    return fig

def plot_bifurcation(param_name, param_values, a, b, c, ic, **kwargs):
    p_plot, x_plot = rossler.bifurcation_maxima_X(param_name, param_values, a, b, c, ic, **kwargs)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(p_plot, x_plot, s=1, alpha=0.5)
    ax.set_xlabel(param_name)
    ax.set_ylabel("maxima of X")
    ax.set_title(f"Bifurcation diagram: maxima of X vs {param_name}")
    ax.grid(True, alpha=0.3)
    return fig


# ============================================================
# Wrappers used in notebook cells
# ===========================================================

# --- 3.1 Basic Visualization ---
def nb_show_attractor(a, b, c, initial_condition, t_min, t_max, dt):
    """Show 3D attractor with interactive sliders for parameters and initial conditions."""
    x0, y0, z0 = initial_condition
    
    # Define slider configurations
    slider_configs = {
        'a': {'value': a, 'min': 0.05, 'max': 0.4, 'step': 0.05},
        'b': {'value': b, 'min': 0.05, 'max': 0.4, 'step': 0.05},
        'c': {'value': c, 'min': 0.05, 'max': 11.4, 'step': 0.05},
        'x0': {'value': x0, 'min': -30, 'max': 30, 'step': 0.1},
        'y0': {'value': y0, 'min': -30, 'max': 30, 'step': 0.1},
        'z0': {'value': z0, 'min': -30, 'max': 30, 'step': 0.1},
    }
    
    # Define plotting function
    def _plot(a, b, c, x0, y0, z0):
        ic_val = (x0, y0, z0)
        plot_3d_attractor(a=a, b=b, c=c, initial_condition=ic_val,
                          t_min=t_min, t_max=t_max, dt=dt, figsize=(12, 5))
        plt.show()
    
    # Create interactive plot
    create_interactive_plot(_plot, slider_configs, n_cols=3)

# --- All three 2D projections ---
def nb_show_projections(a, b, c, initial_condition, t_min, t_max, dt, skip_first_frac=0.2):
    """Show all three 2D projections with interactive sliders."""
    x0, y0, z0 = initial_condition
    
    # Define slider configurations
    slider_configs = {
        'a': {'value': a, 'min': 0.05, 'max': 0.4, 'step': 0.05},
        'b': {'value': b, 'min': 0.05, 'max': 0.4, 'step': 0.05},
        'c': {'value': c, 'min': 0.05, 'max': 11.4, 'step': 0.05},
        'x0': {'value': x0, 'min': -10, 'max': 10, 'step': 0.1},
        'y0': {'value': y0, 'min': -10, 'max': 10, 'step': 0.1},
        'z0': {'value': z0, 'min': -10, 'max': 10, 'step': 0.1},
    }
    
    # Define plotting function
    def _plot(a, b, c, x0, y0, z0):
        ic_val = (x0, y0, z0)
        plot_projections_2D(a=a, b=b, c=c, initial_condition=ic_val,
                           t_min=t_min, t_max=t_max, dt=dt, method="RK45",
                           figsize=(15, 5), skip_first_frac=skip_first_frac)
        plt.show()
    
    # Create interactive plot
    create_interactive_plot(_plot, slider_configs, n_cols=3)

# --- Time series plots ---
def nb_time_series(a, b, c, initial_condition, t_min, t_max, dt):
    """Show x(t), y(t), z(t) time series with interactive sliders."""
    x0, y0, z0 = initial_condition
    
    # Define slider configurations
    slider_configs = {
        't_min': {'value': t_min, 'min': 0.0, 'max': 100.0, 'step': 10},
        't_max': {'value': t_max, 'min': 10, 'max': 500.0, 'step': 10},
    }
    
    # Define plotting function
    def _plot(t_min, t_max):
        ic_val = (x0, y0, z0)
        plot_time_series(a=a, b=b, c=c, initial_condition=ic_val,
                        t_min=t_min, t_max=t_max, dt=dt, figsize=(14, 8))
        plt.show()
    
    # Create interactive plot
    create_interactive_plot(_plot, slider_configs, n_cols=3)

# --- Part 3.3: Compare Different Initial Conditions ---
def nb_compare_initial_conditions(a, b, c, initial_conditions, t_min, t_max, dt, skip_first_frac=0.2):
    """Compare trajectories from several initial conditions."""
    if initial_conditions is None:
        # Use defaults list from get_parameters
        _, _, _, _, _, _, _, initial_conditions = rossler.get_parameters()
    compare_initial_conditions(a=a, b=b, c=c, initial_conditions=initial_conditions,
                               t_min=t_min, t_max=t_max, dt=dt, figsize=(15, 5), skip_first_frac=skip_first_frac)
    plt.show()

# --- Return map + cobweb in one go ---
def nb_return_map_and_cobweb(a, b, c, initial_condition, t_min, t_max, dt,
                             skip_first_maxima=5, cobweb_skip_first=50, n_iter=15,
                             show_cobweb=True, use_widget=False):
    """Plot Z(t) with maxima plus return map; cobweb overlay is optional."""
    ic = initial_condition

    # Time series (short window) with marked maxima
    t_ts, z_ts_data, t_series, z_series = rossler.get_Z_maxima(a, b, c, ic, t_min=0, t_max=200, dt=dt)

    # Maxima for return map and cobweb (longer window)
    _, Z_maxima_all, _, _ = rossler.get_Z_maxima(a, b, c, ic, t_min=0, t_max=2000, dt=dt)
    if len(Z_maxima_all) < 2:
        raise ValueError("Not enough maxima to build a return map. Increase integration time or reduce skip.")

    Z_use = np.asarray(Z_maxima_all[skip_first_maxima:], dtype=float)
    if Z_use.size < 2:
        Z_use = np.asarray(Z_maxima_all, dtype=float)
    if Z_use.size < 2:
        raise ValueError("Not enough maxima to build a return map after skipping transients.")

    Z_n = Z_use[:-1]
    Z_n1 = Z_use[1:]

    # Build empirical return map function
    F = rossler.rossler_return_map_from_data(Z_use)

    # Helper to draw panels so we can re-use when toggling cobweb
    def _draw(fig, axes, cobweb_flag: bool):
        # Panel 1: Z(t) with maxima
        ax = axes[0]
        ax.cla()
        ax.plot(t_series, z_series, linewidth=0.8, label='Z(t)')
        ax.scatter(t_ts, z_ts_data, color='tab:red', s=14, label='local maxima $Z_n$')
        ax.set_xlabel('t')
        ax.set_ylabel('Z(t)')
        ax.set_title(rf'Local maxima $Z_n$ for $c={c:.2f}$')
        ax.grid(True, alpha=0.3)

        # Panel 2: Return map with optional cobweb
        ax = axes[1]
        ax.cla()
        ax.scatter(Z_n, Z_n1, s=12, alpha=0.7, color='tab:blue', label=r"$(Z_n, Z_{n+1})$")
        z_min = float(min(Z_n.min(), Z_n1.min()))
        z_max = float(max(Z_n.max(), Z_n1.max()))
        pad = 0.05 * (z_max - z_min if z_max > z_min else 1.0)
        z_min -= pad
        z_max += pad
        ax.plot([z_min, z_max], [z_min, z_max], 'k--', linewidth=0.8, label='y = x')
        if cobweb_flag:
            x0 = float(Z_use[max(0, cobweb_skip_first):][0]) if Z_use.size > cobweb_skip_first else float(Z_use[0])
            # Use cobweb_plot function but without showing legend on axis
            cobweb_plot(F, x0=x0, n_iter=n_iter, x_min=z_min, x_max=z_max, ax=ax,
                        func_label=r"$Z_{n+1} = F(Z_n)$", cob_color="tab:red", cob_lw=0.8, 
                        mark_start=True, show_legend=False)
        ax.set_xlim(z_min, z_max)
        ax.set_ylim(z_min, z_max)
        ax.set_xlabel(r"$Z_n$")
        ax.set_ylabel(r"$Z_{n+1}$")
        ax.set_title(rf'Return map (cobweb={"on" if cobweb_flag else "off"}) for $c={c:.2f}$')
        ax.grid(True, alpha=0.3)

    # Figure layout (two panels: time series + return map/cobweb)
    plt.close('all')
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [1.4, 1]})
    _draw(fig, axes, cobweb_flag=show_cobweb)
    
    # Make right panel square
    axes[1].set_aspect('equal', adjustable='box')

    # Shared legend below both panels
    handles, labels = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    if handles:
        fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.05),
                   ncol=4, frameon=True, fontsize=9, edgecolor='black')

    plt.tight_layout(rect=[0, 0.12, 1, 1])

    # Optional interactive toggle (requires ipywidgets)
    if use_widget:
        try:
            import ipywidgets as widgets
            from IPython.display import display

            toggle = widgets.Checkbox(value=show_cobweb, description="Show cobweb overlay", indent=False)

            def _on_change(change):
                if change.get('name') == 'value':
                    _draw(fig, axes, cobweb_flag=change['new'])
                    handles.clear(); labels.clear()
                    for ax in axes:
                        h, l = ax.get_legend_handles_labels()
                        handles.extend(h)
                        labels.extend(l)
                    if handles:
                        fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.1), ncol=3,
                                   frameon=True, fontsize=9, edgecolor='black')
                    fig.canvas.draw_idle()

            toggle.observe(_on_change, names='value')
            display(toggle)
        except Exception:
            pass


# --- Butterfly effect plot --- 
def nb_butterfly_effect(a, b, c, initial_condition, delta=1e-6, t_min=50, t_max=280, dt=0.01):
    """Show butterfly effect plot with interactive sliders."""
    x0, y0, z0 = initial_condition
    
    # Define slider configurations
    slider_configs = {
        't_min': {'value': t_min, 'min': 0.0, 'max': 100.0, 'step': 10},
        't_max': {'value': t_max, 'min': 10, 'max': 500.0, 'step': 10},
    }
    
    # Define plotting function
    def _plot(t_min, t_max):
        ic_val = (x0, y0, z0)
        show_butterfly_effect(a=a, b=b, c=c, initial_condition=ic_val, 
                            delta=delta, t_min=t_min, t_max=t_max, dt=dt)
    
    # Create interactive plot
    create_interactive_plot(_plot, slider_configs, n_cols=3)

# --- Sensitivity to initial conditions plot ---
def nb_plot_sensitivity(c, t_min, t_max, dt, delta=1e-6):
    """Plot distance growth with interactive sliders."""
    # Define slider configurations
    slider_configs = {
        'c': {'value': c, 'min': 0.05, 'max': 11.4, 'step': 0.05},
        'delta': {'value': delta, 'min': 1e-8, 'max': 1e-4, 'step': 1e-7},
    }
    
    # Define plotting function
    def _plot(c, delta):
        plot_sensitivity_to_initial_conditions(c=c, t_min=t_min, t_max=t_max, 
                                              dt=dt, delta=delta)
    
    # Create interactive plot
    create_interactive_plot(_plot, slider_configs, n_cols=2)

# --- Parameter sweep for c values ---
def nb_parameter_sweep(c_values, a, b, initial_condition, t_min, t_max, dt):
    """Show attractor panels for a list of c values."""
    plot_parameter_sweep(a=a, b=b, c_values=c_values,
                        initial_condition=initial_condition, t_min=t_min, t_max=t_max, dt=dt,
                        figsize=(15, 10))
    plt.show()

# --- Part 3.5: Bifurcation Diagram ---
def nb_bifurcation_vs_c(a=0.2, b=0.2, c_min=0.0, c_max=18.0, n=400,
                        ic=initial_condition, t_min=0.0, t_max=200.0, h=0.01):
    c_values = np.linspace(c_min, c_max, n)
    fig = plot_bifurcation("c", c_values, a=a, b=b, c=0.0, ic=ic,
                           t_min=t_min, t_max=t_max, h=h, skip_frac=0.5)
    plt.show()
    return fig

def nb_bifurcation_vs_b(a=0.2, c=5.7, b_min=0.0, b_max=1.0, n=400,
                        ic=initial_condition, t_min=0.0, t_max=200.0, h=0.01):
    b_values = np.linspace(b_min, b_max, n)
    fig = plot_bifurcation("b", b_values, a=a, b=0.0, c=c, ic=ic,
                           t_min=t_min, t_max=t_max, h=h, skip_frac=0.5)
    plt.show()
    return fig
# --- Euler vs RK4 Comparison ---
def nb_compare_euler_vs_rk4(a, b, c, initial_condition, t_min, t_max, 
                            h_rk4=0.1, h_euler=None, h0_euler=0.1, max_halvings=20):
    """Compare Euler vs RK4 trajectories with interactive sliders."""
    x0, y0, z0 = initial_condition
    
    # Define slider configurations
    slider_configs = {
        't_max': {'value': t_max, 'min': 1.0, 'max': 500.0, 'step': 1.0},
    }
    
    # Define plotting function
    def _plot(t_max):
        ic_val = (x0, y0, z0)
        
        results = rossler.compare_euler_vs_rk4(
            a, b, c, ic_val, t_min, t_max, h_rk4, h_euler, h0_euler, max_halvings
        )

        t_ref, x_ref, y_ref, z_ref = results['t_ref'], results['x_ref'], results['y_ref'], results['z_ref']
        t_eu, x_eu, y_eu, z_eu = results['t_eu'], results['x_eu'], results['y_eu'], results['z_eu']
        h_best = results['h_euler']
        
        # Print error statistics
        print(f"Max |Δx|: {results['max_err_x']:.3e} | Max |Δy|: {results['max_err_y']:.3e} | Max |Δz|: {results['max_err_z']:.3e}")
        
        # Create comparison plots
        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

        axes[0].plot(t_ref, x_ref, label=f'RK4 h={h_rk4:.3f}')
        axes[0].plot(t_eu, x_eu, '--', label=f'Euler h={h_best:.4f}')
        axes[0].set_ylabel('x(t)')
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(t_ref, y_ref, label=f'RK4 h={h_rk4:.3f}')
        axes[1].plot(t_eu, y_eu, '--', label=f'Euler h={h_best:.4f}')
        axes[1].set_ylabel('y(t)')
        axes[1].grid(True, alpha=0.3)

        axes[2].plot(t_ref, z_ref, label=f'RK4 h={h_rk4:.3f}')
        axes[2].plot(t_eu, z_eu, '--', label=f'Euler h={h_best:.4f}')
        axes[2].set_ylabel('z(t)')
        axes[2].set_xlabel('t')
        axes[2].grid(True, alpha=0.3)

        fig_legend_dedup(fig, axes[0], loc='lower center', ncol=2, bbox_to_anchor=(0.5, 0.0),
                         frameon=True, fontsize=10, edgecolor='black')

        plt.tight_layout(rect=[0, 0.05, 1, 1])
        plt.show()
    
    # Create interactive plot
    create_interactive_plot(_plot, slider_configs, n_cols=3)