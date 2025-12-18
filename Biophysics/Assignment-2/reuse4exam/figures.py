"""
figures.py implements all figure plotting functions for the Rössler project.
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
from imports import (np, sp, scipy, solve, Symbol, re, im, Matrix, symbols, Eq, Line2D, plt,Axes3D, solve_ivp, interact, tqdm)

import parameters
from parameters import get_parameters

from data import (get_time_array,get_vectors,get_solutions, get_groups ,compute_jacobian_and_stability,get_eigvals_for_c, get_equilibria_for_c,
ZY_data_and_wings, XYZ_data_and_wings, get_Z_maxima, rossler_return_map_from_data, rossler_chaos)

# ------- parameters ----
a, b, t_min, t_max, dt, initial_condition, initial_conditions = get_parameters()



# ============================================================
# 1. 1D bifurcation diagrams: equilibria vs c
# ============================================================


def plot_bifuractions(results, ax=None, component="X", only_small_c=False, show_figure=False):
    """
    Plot the bifurcation with equilibria for the Rössler system.
    (This is a general function for X, Y and Z to reduce duplicate code).

    results dictionary: {c_val, xsol, ysol, zsol, J_eq, eigvals, stable}.

    'show_figure' is to know if we have to make a subplot here or later, for later use.
    """
    
    key_map = {"X": "x_eq", "Y": "y_eq", "Z": "z_eq"}
    key_eq = key_map[component]


    groups = get_groups(dim="1D")

    for entry in results:
        c_val  = entry["c"]
        value_eq = entry[key_eq]
        eq_type  = entry["eq_type"]

        if only_small_c and c_val >= 2:
            continue

        if eq_type in groups:
            groups[eq_type]["c"].append(c_val)
            groups[eq_type]["val"].append(value_eq) 
        else: 
            # Unknowns in "unstable"
            groups["unstable node"]["c"].append(c_val)
            groups["unstable node"]["val"].append(value_eq)


    if show_figure:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)

    # Equilibria
    for cat in groups.values():
        if cat["c"]:  
            ax.scatter(cat["c"], cat["val"], **cat["style"])
    
    ax.set_ylabel(rf"{component}")
    ax.set_title(rf"Equilibria for {component} in function of the parameter $c$")
    
    ax.set_xlabel(r"$c$")
    ax.legend(loc='best',fontsize=6)
    ax.grid(True)

    if show_figure:
        plt.show()


def plot_eigvals_vs_c(small_c=False, part="real", title='Jacobian eigenvalues at equilibria'):
    """
    Plot the Jacobian eigenvalues at the equilibria as a function of c.

    Parameters
    ----------
    small_c : bool
        If True: use the 'small c' grid (0 <= c <= 2).
        If False: use the full grid (0 <= c <= 18).
    part : {"real", "imag"}
        Which part of the eigenvalues to plot on the y-axis.
    """
    # Get symbolic equilibria
    sols, (x_sym, y_sym, z_sym) = get_solutions()

    # Compute all equilibria in terms of different c
    all_eq = compute_jacobian_and_stability(sols, x_sym, y_sym, z_sym, small_c=small_c)


    cs = []
    eigvals_list = []

    for entry in all_eq:
        c_val = float(entry["c"])
        eigvals = np.array(entry["eigvals"], dtype=complex)

        if part == "real":
            vals = np.real(eigvals)
            ylabel = r"$\Re(\lambda)$"
        elif part == "imag":
            vals = np.imag(eigvals)
            ylabel = r"$\Im(\lambda)$"
        else:
            raise ValueError("part must be 'real' or 'imag'")

        cs.append(c_val)
        eigvals_list.append(vals)

    cs = np.array(cs)
    eigvals_arr = np.vstack(eigvals_list)

    order = np.argsort(cs)
    cs = cs[order]
    eigvals_arr = eigvals_arr[order, :]

    # Plot
    plt.close("all")
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)

    n_eig = eigvals_arr.shape[1]
    for j in range(n_eig):
        ax.scatter(cs, eigvals_arr[:, j], label=fr"$\lambda_{j+1}$")

    # Axis styling
    ax.axhline(0.0, color="k", linewidth=0.5)
    ax.set_xlabel(r"$c$")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True)
    ax.legend(loc='best')

    plt.show()

        
def show_2D_bif_figs(results, components=("X","Y","Z"), figsize=(11,9)):
    """
    Create a 2x2 grid of 2D bifurcation diagrams for the Rössler equilibria.

    Parameters
    ----------
    results : list[dict]
        Output of compute_jacobian_and_stability, containing equilibria info.
    components : tuple of str
        State components to plot on the y-axes (default: ("X", "Y", "Z")).
    figsize : tuple
        Figure size passed to plt.subplots.
    """
    plt.close()
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    
    for ax, component in zip(axes.flat, components):
        plot_bifuractions(results, ax=ax, component=component, only_small_c=False)
    axes[1][1].axis('off')
    plt.show()


# ============================================================
# 2. Utility: arrows, axes and 3D helpers
# ============================================================
    

def add_direction_arrows(ax, x, y, z, length=0.5):
    """
    Add a few arrows along a 3D trajectory to indicate the direction of motion.
    The places of the arrows are manually chosen to avoid bad placed arrows.

    Parameters
    ----------
    ax : mpl_toolkits.mplot3d.Axes3D
        The 3D axis to draw on.
    x, y, z : array-like
        Coordinates of the trajectory, all of the same length.
    length : float, optional
        Base length passed to `ax.quiver`.
    """
    # --- Manually entered arrows ---
    arrows = [(0,1), (300,350)]
    for r1, r2 in arrows:
        dx = x[r2] - x[r1]
        dy = y[r2] - y[r1]
        dz = z[r2] - z[r1]
        ax.quiver(
            [x[r1]], [y[r1]], [z[r1]], #start
             [dx], [dy], [dz], #direction
             length=length, normalize=True,
        color='k', linewidth=1,arrow_length_ratio=0.3)
            

def plot_minmax_values_and_axes(ax, xlim=(-10, 10), ylim=(-10, 10), zlim=(0, 20)):
    """
    Draw coordinate axes with arrows and set fixed 3D limits.

    Parameters
    ----------
    ax : mpl_toolkits.mplot3d.Axes3D
        The 3D axes on which to draw.
    xlim, ylim, zlim : tuple of float
        (min, max) limits for the respective axes.
    """
    xmin,xmax = xlim
    ymin,ymax = ylim
    zmin,zmax = zlim

    # --- Plot the coordinate axes as thin black lines ---
    ax.plot([xmin,xmax], [0, 0], [0, 0], color='k', lw=0.5) # x-axes
    ax.plot([0, 0], [ymin,ymax], [0, 0], color='k', lw=0.5) # y-axes
    ax.plot([0, 0], [0, 0], [zmin,zmax], color='k', lw=0.5) # z-axes

    # --- Arrows at the positive ends of the axes ---
    arrow_kw = dict(color='k', linewidths=1, arrow_length_ratio=0.5)

    Lx = xmax - xmin
    Ly = ymax - ymin
    Lz = zmax - zmin

    ax.quiver(xmax - 0.05*Lx, 0, 0,  0.05*Lx, 0, 0, **arrow_kw)  # x-arrow
    ax.quiver(0, ymax - 0.05*Ly, 0,  0, 0.05*Ly, 0, **arrow_kw)  # y-arrow
    ax.quiver(0, 0, zmax - 0.05*Lz,  0, 0, 0.05*Lz, **arrow_kw)  # z-arrow

    # --- Set limit to the axes --- 
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)



# ============================================================
# 3. 3D trajectories for fixed c
# ============================================================


def show_3Dfigure_plot(c=5.0, small_c=False):
    """
    Show Rössler trajectories in 3D for several initial conditions and mark
    the corresponding equilibria for a given c.

    Parameters
    ----------
    c : float
        Value of the control parameter c.
    small_c : bool, optional
        If True, draw axis bounds and arrows adapted to small c regime.
    """
    plt.close('all')
    fig, ax = plt.subplots(1, 1, figsize=(8, 6), subplot_kw={'projection': '3d'})
    ax.set_title(rf'3D trajectories for $c = {c:.2f}$')
    
    if small_c:
        plot_minmax_values_and_axes(ax)
        
    # --- Trajectories ---    
    for initial_c in initial_conditions:
        x,y,z = get_vectors(c,initial_condition=initial_c)
                
        ax.plot3D(x, y, z,label=f"Initial conditions: {initial_c}")
                
        # --- Direction at given point of trajectory ---
        add_direction_arrows(ax, x, y, z, length=0.5)
        

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
    plt.show()

        

def show_3Dfigure_scatter(c=5.0, small_c=False):
    """
    3D scatter plot of Rössler trajectories for a single value of c.

    Intended use:
    - Call this function from an ipywidgets.interact slider for c,
      so you can see how the trajectories change as c varies.

    Parameters
    ----------
    c : float, optional
        Parameter c in the Rössler system.
    small_c : bool, optional
        If True, draw coordinate axes and bounds tailored for small c.
    """
    plt.close('all')
    fig, ax = plt.subplots(1, 1, figsize=(8, 7), subplot_kw={'projection': '3d'}, constrained_layout=True)
    fig.suptitle(r"Rössler trajectories for different values of $c$")
    
    if small_c:
        plot_minmax_values_and_axes(ax)
    ax.set_title(rf"3D scatter of Rössler trajectories, $c = {c:.2f}$")
    
    for initial_c_val in tqdm(initial_conditions,desc=f"Scatter Plot (c={c})"):
        x, y, z = get_vectors(c,initial_condition=initial_c_val)
        ax.scatter3D(x, y, z, s=5,
                     label=f"Initial conditions {initial_c_val}")
    
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend(loc='best',fontsize=6)
    plt.show()



# ============================================================
# 4. 3D bifurcation diagrams in state space
# ============================================================
    

def plot_bifurcation_3D_two_components(results, ax=None, component1='X', component2='Y'):
    """
    Plot a 3D bifurcation diagram of the Rössler equilibria.

    The axes show (c, component1_eq, component2_eq). Points are grouped by
    equilibrium type (stable, unstable, saddle, pitchfork, Hopf, other),
    with marker/colour styles defined in get_groups(dim="2D").

    Parameters
    ----------
    results : list of dict
        Output of compute_jacobian_and_stability. Each entry must contain
        at least the keys "c", "x_eq", "y_eq", "z_eq" and "eq_type".
    ax : mpl_toolkits.mplot3d.Axes3D
        3D axis to draw on. This function does NOT call plt.show().
    component1, component2 : {"X", "Y", "Z"}
        Which equilibrium coordinates to use on the Y- and Z-axes.
    """

    key_map = {"X": "x_eq", "Y": "y_eq", "Z": "z_eq"}
    key1_eq = key_map[component1]
    key2_eq = key_map[component2]

    groups = get_groups(dim="2D")

    for entry in results:
        c_val = entry["c"]
        v1 = entry[key1_eq]
        v2 = entry[key2_eq]
        eq_type = entry["eq_type"] 

        if eq_type not in groups:
            eq_type = "other equilibrium"

        groups[eq_type]["c"].append(c_val)
        groups[eq_type]["v1"].append(v1)
        groups[eq_type]["v2"].append(v2)

        
    # --- Scatter each non-empty type (group) ---
    for group in groups.values():
        c_vals = np.asarray(group["c"])
        if c_vals.size == 0:
            continue  # skip empty groups from the legend
        v1_vals = np.asarray(group["v1"])
        v2_vals = np.asarray(group["v2"])
        ax.scatter(c_vals, v1_vals, v2_vals, **group["style"])

    # ---- General plot stuff ----- 
    ax.set_xlabel(r"$c$")
    ax.set_ylabel(rf"{component1}")
    ax.set_zlabel(rf"{component2}")
    ax.set_title(f"3D bifurcation diagram of equilibria in "
                 f"({component1}_eq, {component2}_eq) space")
    ax.legend()


def three_3D_bifurcation_figures(results, figsize=(11, 9)):
    """
    Plot three 3D bifurcation diagrams for the component pairs (X, Y), (Y, Z) and (Z, X).

    Parameters
    ----------
    results : dict or pandas.DataFrame
        Container with the bifurcation data, passed directly to
        `plot_bifurcation_3D_two_components`.
    figsize : tuple of float, optional
        Figure size passed to `plt.subplots`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure.
    axes : numpy.ndarray of Axes3D
        Array of axes with shape (2, 2). The bottom-right axis is turned off.
    """
    plt.close()
    fig, axes = plt.subplots(2, 2, figsize=(11, 9),subplot_kw={'projection': '3d'})

    # Which component pairs to plot in 3D
    component_pairs = [("X", "Y"), ("Y", "Z"), ("Z", "X")]
    
    flat_axes = axes.ravel()

    for ax, (c1, c2) in zip(flat_axes, component_pairs):
        plot_bifurcation_3D_two_components(results, ax, component1=c1, component2=c2)

    for ax in flat_axes[len(component_pairs):]: # To turn the last plot off.
        ax.axis("off") 

    plt.show()


    
def plot_equilibria_all_parameters_3D(results):
    """
    Plot all equilibria in (X_eq, Y_eq, Z_eq)-space for the Rössler system.

    Colour encodes the parameter c, marker encodes stability
    (circles for stable equilibria, triangles for unstable ones).
    """
    #extract data
    X = np.array([e["x_eq"] for e in results])
    Y = np.array([e["y_eq"] for e in results])
    Z = np.array([e["z_eq"] for e in results])
    C = np.array([e["c"]  for e in results])
    STABLE = np.array([e["eq_type"] == "stable equilibrium" for e in results], dtype=bool,)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")

    # stable equilibria: circles
    sc_stable = ax.scatter(
        X[STABLE], Y[STABLE], Z[STABLE],
        c=C[STABLE],
        cmap="viridis",
        marker="o",
        edgecolors="k",
        linewidths=0.5,
        s=10,
        label="stable",
        depthshade=False)

    # unstable equilibria: triangles
    sc_unstable = ax.scatter(
        X[~STABLE], Y[~STABLE], Z[~STABLE],
        c=C[~STABLE],
        cmap="viridis",
        marker="^",
        s=50,
        linewidths=0.25,
        edgecolors="k",
        label="unstable",
        depthshade=False)

    # Labels and title
    ax.set_xlabel(r"$X_{\mathrm{eq}}$")
    ax.set_ylabel(r"$Y_{\mathrm{eq}}$")
    ax.set_zlabel(r"$Z_{\mathrm{eq}}$")
    ax.set_title(r"Equilibria of Rössler system in state space $(X,Y,Z)$")

    # Colorbar for c values
    cbar = fig.colorbar(sc_stable, ax=ax, pad=0.1)
    cbar.set_label(r"$c$")

    ax.legend()
    plt.show()


# ============================================================
# 5. Time series plots x(t), y(t), z(t)
# ============================================================
    


def plot_xyz_func_of_t(c, ax=None, component='X', y_label=None,
                       t_min=t_min, t_max=t_max, dt=dt):
    """
    Plot one component of the Rössler trajectory as a function of time.

    Parameters
    ----------
    c : float
        Value of the parameter c.
    ax : matplotlib.axes.Axes or None, optional
        Axis to plot on. If None, a new figure and axis are created.
    component : {"X", "Y", "Z"}, optional
        Which component to plot: "X", "Y", or "Z".
    y_label : str or None, optional
        Label for the y-axis. If None, uses "<component>(t)".
    t_min, t_max : float
        Time interval [t_min, t_max] for the simulation.
    dt : float
        Time step used to build the time array.
    """
    
    # Time and trajectories
    t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
    x, y, z = get_vectors(c, t_min=t_min, t_max=t_max, dt=dt)

    
    key_map = {"X": x, "Y": y, "Z": z}
    key_value = key_map[component]
    if y_label is None:
        y_label = f'{component}(t)'

    # --- Plot ---

    ax.plot(t,key_value,label=f'{component}(t)', alpha=0.7)
    ax.set_ylabel(y_label)
    ax.set_xlabel('t')
    ax.grid(True)



def show_4_func_of_t_figs(c=5.0, t_min=t_min, t_max=200, dt=dt):
    """
    Plot X(t), Y(t) and Z(t) in three separate panels and once combined.

    For a given value of c, integrate the Rössler system on [t_min, t_max]
    and show:
        - X(t)
        - Y(t)
        - Z(t)
        - all three components in a single panel for direct comparison.
    """
    plt.close('all')
    fig, axes = plt.subplots(4, 1, figsize=(12, 8), constrained_layout=True)
    fig.suptitle(rf'$X(t)$, $Y(t)$ and $Z(t)$ for $c={c:.2f}$')

    components = ["X", "Y", "Z"]

    # --- Three separate panels ---
    for ax, comp in zip(axes[:3], components):
        plot_xyz_func_of_t(c=c, ax=ax, component=comp, t_min=t_min, t_max=t_max, dt=dt)
        ax.legend(loc='best', fontsize=6)    

    
    # --- Fourth panel: all components together ---
    combined_ax=axes[3]
    y_label = r'$f(t)$'

    for comp in components:
        plot_xyz_func_of_t(c=c, ax=combined_ax,component=comp,
                           t_min=t_min, t_max=t_max, dt=dt)

    combined_ax.set_xlabel("t")
    combined_ax.set_ylabel(y_label)
    combined_ax.legend(loc='best', fontsize=6)  

    plt.show()



def show_diff_c_vals(t_min=t_min, t_max=t_max, dt=dt):
    """
    Plot X(t), Y(t) and Z(t) for several values of c in a vertical stack
    of subplots.

    For each c in `cs`, we plot the three components X(t), Y(t) and Z(t)
    in a single panel to compare their time evolution for that parameter value.
    """
    cs=(4.0, 5.7, 8.5)
    
    plt.close('all')
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), constrained_layout=True)  
    fig.suptitle(rf'$X(t)$, $Y(t)$ and $Z(t)$ for different values of $c$')
    
    y_label = r'$f(t)$'

    for ax, c in zip(axes, cs):
        ax.set_title(fr'$c={c}$')
        for comp in ("X", "Y", "Z"):
            plot_xyz_func_of_t(c=c, ax=ax, component=comp, y_label=y_label,
                               t_min=t_min, t_max=t_max, dt=dt)

        
    ax.legend(loc='best', fontsize=6)  
    plt.show()



# ============================================================
# 6. Stability near equilibria: perturbed trajectories
# ============================================================

    
def show_stability_of_equilibrium(c=3.0, t_min=0.0, t_max=100.0, dt=dt,
                                 show_perturbation1=True, show_perturbation2=True, show_perturbation3=True):
    """
    Illustrate the stability of an equilibrium in the Rössler system.

    For a given c:
    - find the (numerically) stable equilibrium using the Jacobian eigenvalues,
    - start several trajectories with small perturbations around this equilibrium,
    - show that X(t), Y(t) and Z(t) relax back to the equilibrium values.

    Parameters
    ----------
    c : float
        Parameter value in the Rössler system.
    t_min, t_max : float
        Time interval for the simulation.
    dt : float
        Time step for the integration.
    show_perturbation1/2/3 : bool
        Whether to show each of the three perturbed trajectories.
    """
    # Get symbolic solutions and all equilibria
    sols, (x_sym, y_sym, z_sym) = get_solutions()

    c_sym, a_sym, b_sym = sp.symbols("c a b", real=True)

    eqs_here = []
    for sol in sols:
        eq = get_equilibria_for_c(sol, x_sym, y_sym, z_sym, c_sym, a_sym, b_sym, c)
        
        if eq is not None and eq["eq_type"] in ("stable node", "stable spiral"):
            eqs_here.append(eq)

    if not eqs_here:
        raise RuntimeError(f"No stable equilibrium found for c={c:.2f}.")

    eq = eqs_here[0]    
    center = np.array([eq["x_eq"], eq["y_eq"], eq["z_eq"]])

    # Small perturbations around the equilibrium
    perts = np.array([
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]])
    
    show_flags = [show_perturbation1, show_perturbation2, show_perturbation3]

    t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
    
    plt.close('all')
    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
    labels = ["X(t)", "Y(t)", "Z(t)"]
    colors = ["C0", "C1", "C2"]

    for p, show, color in zip(perts, show_flags, colors):
        if not show:
            continue 
        ic = tuple(center + p)
        x, y, z = get_vectors(c=c, initial_condition=ic, t_min=t_min, t_max=t_max, dt=dt)
        for ax, series in zip(axes, [x, y, z]):
            ax.plot(t, series, alpha=0.7, color=color)

    # Draw equilibrium values as dashed lines
    for ax, val, lab in zip(axes, center, labels):
        ax.axhline(val, color="k", linestyle="--", linewidth=1) #Equilibrium
        ax.set_ylabel(lab)
        ax.grid(True)

    axes[-1].set_xlabel("t")
    axes[0].set_title(rf"Convergence towards a stable equilibrium for $c={c:.2f}$ "
        r"with perturbed trajectories")

    # Custom legend
    legend_lines = [ Line2D([0], [0], color="C0", lw=1.5),
                     Line2D([0], [0], color="C1", lw=1.5),
                     Line2D([0], [0], color="C2", lw=1.5),
                     Line2D([0], [0], color="k", lw=1.0, linestyle="--")]
    
    legend_labels = ["perturbed trajectory",
                     "perturbed trajectory",
                     "perturbed trajectory",
                     "equilibrium"]
    
    axes[0].legend(legend_lines, legend_labels, loc="best",fontsize=8)
    
    plt.show()


# ============================================================
# 7. Geometry & wings: Z(Y), 3D orbit, wing colouring
# ============================================================


def plot_ZvY(c=5.7, t_min=0.0, t_max=200.0 , dt=dt):
    """
    Plot Z as a function of Y for a given value of c.
    
    This makes a 2D projection of the Rössler trajectory in the (Y, Z)-plane,
    using the time interval [t_min, t_max].
    """
    y, z, _, _ = ZY_data_and_wings(c, t_min, t_max, dt)

    plt.close('all')
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(y, z, label=r'Z(Y(t))', linewidth=0.7)

    ax.set_xlabel(r'$Y$')
    ax.set_ylabel(r'$Z$')
    ax.set_title(rf'$Z$ as a function of $Y$ for $c = {c:.2f}$,' 
                 rf' for $t \in [{t_min:.1f}, {t_max:.1f}]$')
    
    ax.legend(loc='best')
    ax.grid(True)
    plt.show()

def plot_ZvY_wings(c=5.7, t_min=t_min, t_max=t_max, dt=0.01):
    """
    Plot Z as a function of Y and colour the points according to the wing (sign of X).

    This highlights the 'switching wings' behaviour of the Rössler attractor:
    points with X < 0 and X > 0 are shown in different colours in the (Y, Z)-plane.
    """
    y, z, left, right = ZY_data_and_wings(c, t_min, t_max, dt)
    
    plt.close('all')
    fig, ax = plt.subplots(figsize=(6, 5))

    # Left wing:  X < 0
    ax.scatter(y[left], z[left],s=2,label='left wing (x < 0)',alpha=0.5)
    
    # Right wing:  X > 0
    ax.scatter(y[right], z[right],s=2,label='right wing (x > 0)',alpha=0.5)

    ax.set_xlabel(r'$Y$')
    ax.set_ylabel(r'$Z$')
    ax.set_title(rf'$Z(Y)$ for $c = {c:.2f}$, '
                 rf'$ t \in [{t_min:.1f}, {t_max:.1f}]$, with visual wing switching')
    
    ax.legend(loc='best', fontsize=6)
    ax.grid(True)
    plt.show()



def plot_XYZ_3D(c=5.7, t_min=0.0, t_max=200.0 , dt=dt):
    """
    Plot the Rössler trajectory in 3D state space (X, Y, Z) for a given c
    and time interval [t_min, t_max].
    """
    x ,y, z, _, _ = XYZ_data_and_wings(c, t_min, t_max, dt)

    plt.close('all')
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")
    
    ax.plot(x, y, z, label=r'trajectory in configuration space', linewidth=0.7)
    ax.set_xlabel(r'$X$')
    ax.set_ylabel(r'$Y$')
    ax.set_zlabel(r'$Z$')
    ax.set_title(rf'Rössler attractor for $c = {c:.2f}$,' 
                 rf' for $t \in [{t_min:.1f}, {t_max:.1f}]$')
    
    ax.legend(loc='best', fontsize=6)
    ax.grid(True)
    plt.show()


# ============================================================
# 8. Return map: maxima of Z(t)
# ============================================================

    
def plot_Z_time_with_maxima(c=5.7, t_min=t_min, t_max=200, dt=dt):
    """
    Plot Z(t) over time and mark the detected local maxima Z_n
    for a given c and time interval [t_min, t_max].

    Parameters
    ----------
    c : float
        Rössler parameter c.
    t_min, t_max : float
        Time interval.
    dt : float
        Time step used in the integration.
    """
    t = get_time_array(t_min=t_min, t_max=t_max, dt=dt)
    _, _, z = get_vectors(c, t_min=t_min, t_max=t_max, dt=dt)
    t_max_vals, Z_max = get_Z_maxima(c, t_min=t_min, t_max=t_max, dt=dt)

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


def plot_Zn_vs_Znplus1(c=5.7, skip_first=20, t_min=t_min, t_max=200, dt=dt):
    """
    Plot the return map Z_{n+1} versus Z_n, where Z_n are the local maxima of Z(t).

    This visualises how each maximum Z_n
    (almost) uniquely determines the next one Z_{n+1} for a given c.
    """
    # Extract the local maxima Z_n of Z(t)
    _, Z_max = get_Z_maxima(c=c, skip_first=skip_first,
                            t_min=t_min, t_max=t_max, dt=dt)

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
        Sequence of maxima of Z(t).
    c : float, optional
        Rössler parameter c (only used for the title).
    n_iter : int, optional
        Number of cobweb iterations to draw.
    skip_first : int, optional
        Number of initial maxima to discard as transient.
    """
    # Discard transients
    Z_n = np.asarray(Z_n)[skip_first:]

    F = rossler_return_map_from_data(Z_n)

    # Use percentiles to ignore extreme outliers in the domain
    x_min, x_max = np.percentile(Z_n, [1, 99])

    # Choose a representative starting value, e.g. the median
    x0 = np.median(Z_n)

    plt.close("all")
    fig, ax = plt.subplots(figsize=(5, 5))

    # Scatter plot of the data
    ax.scatter(Z_n[:-1], Z_n[1:], s=8, alpha=0.3, label=r"data $(Z_n, Z_{n+1})$")

    # Cobweb on top of the data
    cobweb_plot(F, x0=x0, n_iter=n_iter, x_min=x_min, x_max=x_max, ax=ax, 
                func_label=r"$Z_{n+1} = F(Z_n)$",
                cob_color="tab:red", cob_lw=0.8)

    ax.set_title(rf"Cobweb for Rössler return map, $c={c:.2f}$")
    ax.set_xlabel(r"$Z_n$")
    ax.set_ylabel(r"$Z_{n+1}$")
    plt.show()
    

# ============================================================
# 9. Sensitivity to initial conditions (Lyapunov flavour)
# ============================================================    


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

    x1, y1, z1 = get_vectors(c,t_min=t_min, t_max=t_max, dt=dt, initial_condition=ic1)
    x2, y2, z2 = get_vectors(c,t_min=t_min, t_max=t_max, dt=dt, initial_condition=ic2)

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


def show_butterfly_effect(c=5.7, delta=1e-6, t_min=t_min, t_max=t_max, dt=dt):
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
    x1, y1, z1 = get_vectors(c, t_min=t_min, t_max=t_max, dt=dt,
                                initial_condition=ic1)
    x2, y2, z2 = get_vectors(c, t_min=t_min, t_max=t_max, dt=dt,
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




# ============================================================
# 10. EXTRA: Bifurcation Diagram
# ============================================================  


def show_bifurcation_diagram():
    """
    Plot a numerical bifurcation diagram for the Rössler system as a function of c.

    The Rössler parameters a and b are fixed to their standard values (a = 0.2, b = 0.2),
    while the control parameter c is varied over a range (0.1–15, geometrically
    spaced). For each value of c, the Rössler equations are integrated in time starting
    from a random initial condition.

    After discarding the initial transient, the long-time values of x(t) are sampled
    and plotted versus c. This produces a bifurcation diagram in the (c, x)-plane that
    visualizes the different dynamical regimes of the system (stable equilibria,
    periodic orbits and chaotic behavior).
    """
    a_val = a
    b_val = b
    T = 500
    dt_val = 0.01
    time_points = np.arange(0., T, dt_val)
    
    c_list = np.geomspace(0.1, 15., 1501)
    sol_list = []

    for c_val in c_list:
        init = np.random.rand(3)*2-1
        sol = scipy.integrate.odeint(rossler_chaos, init, time_points, args=(c_val, a_val, b_val))
        sol_list.append(sol)

    plt.close('all')
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.set_title(rf'Bifurcation Diagram for the Rössler system')
    
    for i in range(len(c_list)):
        c_val = c_list[i]
        sol = sol_list[i]
        y = sol[int(T/dt_val*0.75):,0]
        x = [c_val] * len(y)
        plt.scatter(x, y, s=0.1, c='k', marker='.', edgecolor='none')
    
    ax.set_xlabel(r'$c$')
    ax.set_ylabel(r'$x$')
    ax.legend()
    plt.show()
