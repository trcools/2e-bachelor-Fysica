# --- imports ---
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

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

def quiver_on_curve(ax, x, y, xlim, ylim, a, b, orientation="vertical", length=0.15):
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
        Direction of arrows ("vertical" for vertical nullcline or "horizontal")
    length : float, default 0.15
        Length of arrows in data coordinates
    """
    # keep only points inside the window
    mask = mask_in_window(x, y, xlim, ylim)
    x, y = x[mask], y[mask]

    if x.size == 0:
        return

    s = np.sign(b - x) * length  # direction from remaining component

    if orientation == "vertical":
        X = np.zeros_like(s)
        Y = s
    elif orientation == "horizontal":
        X = s
        Y = np.zeros_like(s)
    else:
        raise ValueError("orientation must be 'vertical' or 'horizontal'")

    ax.quiver(x, y, X, Y, angles="xy", scale_units="xy", scale=1, width=0.003, color="k", zorder=1.5)





def draw_nullclines_panel(
        ax, a=0.06, b=0.4, xlim=(0, 4), ylim=(0, 7),
        n=250, arrows=18,
        density=0.6, alpha=0.35,
        n_nullcline=800, show_equilibrium=True,
        title=None
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
    arrows : int
        Number of direction arrows on each nullcline (default 18)
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
    plot_curve_in_window(ax, xn, Y_xnull, xlim, ylim, label="X-nullcline (X'=0)")
    plot_curve_in_window(ax, xn, Y_ynull, xlim, ylim, label="Y-nullcline (Y'=0)")

    # --- Direction arrows on nullclines (change vectors) ---
    xs = np.linspace(xlim[0], xlim[1], arrows)
    y_on_xnull, y_on_ynull = nullclines(xs, a, b)
    quiver_on_curve(ax, xs, y_on_xnull, xlim, ylim, a, b, orientation="vertical", length=0.35)
    quiver_on_curve(ax, xs, y_on_ynull, xlim, ylim, a, b, orientation="horizontal", length=0.35)

    # --- Equilibrium (intersection of nullclines) ---
    if show_equilibrium:
        Xeq, Yeq = equilibrium(a, b)
        if xlim[0] <= Xeq <= xlim[1] and ylim[0] <= Yeq <= ylim[1]:
            ax.scatter(Xeq, Yeq, label="Equilibrium", color="red", zorder=2)

    ax.grid(True, alpha=0.5)
    
def plot_nullclines(
        a=0.06, b=0.4, xlim=(0, 4), ylim=(0, 7),
        n=250, arrows=18
        ):
    """
    Create and display a complete phase plane portrait.
    
    Shows the nullclines, vector field, direction arrows on nullclines,
    and the equilibrium point for the given parameters.
    
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
    n : int
        Number of grid points for vector field (default 250)
    arrows : int
        Number of direction arrows on each nullcline (default 18)
    """
    plt.close("all")
    fig, ax = plt.subplots(figsize=(8, 6))
    draw_nullclines_panel(
        ax, a=a, b=b, xlim=xlim, ylim=ylim,
        n=n, arrows=arrows, density =0.6, alpha=0.35,
        title=rf"Degradation of Glucose: Nullclines with $a={a}$ and $b={b}$"
        )
    ax.legend(loc="upper right", framealpha=0.9)
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

def classify_equilibrium(a, b):
    """
    Classify the equilibrium point and return its coordinates, Jacobian, and eigenvalues.
    
    Returns:
        (Xeq, Yeq): Equilibrium coordinates
        J_eq: Jacobian matrix at equilibrium
        eig: Eigenvalues of the Jacobian
    """
    return get_equilibrium_eigvals(a, b)

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

    for (X0, Y0) in initials:
        sol = solve_ivp(rhs_state, t_span, [X0, Y0],
                        args=(a, b), t_eval=t_eval,
                        rtol=1e-8, atol=1e-10)
        ax.plot(sol.y[0], sol.y[1], **line_kw)
        ax.plot([X0], [Y0], **start_kw)



def plot_2glucose_trajectories(a=0.06, b=0.4, xlim=(0,4), ylim=(0,7),
                  n=250, arrows=12, step=0.08,
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
    arrows : int
        Number of direction arrows on nullclines (default 12)
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
    draw_nullclines_panel(ax, a=a, b=b, xlim=xlim, ylim=ylim, n=n, arrows=arrows,
                          density=0.8, alpha=0.30, show_equilibrium=True,
                          title="Local trajectories near equilibrium")
    draw_trajectories(ax, a, b, xlim, ylim, near_eq, t_span=t_span, t_eval_n=t_eval_n)

    # --- Global panel ---
    ax = axes[1]
    draw_nullclines_panel(ax, a=a, b=b, xlim=xlim, ylim=ylim, n=n, arrows=arrows,
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
        n=250, arrows=0, density=0.6, alpha=0.35,
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
        n=500, arrows=0, density=2, alpha=0.35,
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
                                 n=250, arrows=0, density=0.8, alpha=0.3,
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

def plot_stability_map(a_max=0.14, b_max=1.2, na=300, nb=400):
    """
    Plot 2D stability map using the trace criterion.
    
    Colors regions where the equilibrium is stable (trace < 0) vs unstable.
    Overlays the Hopf bifurcation curves as black lines.
    
    Parameters
    ----------
    a_max : float
        Maximum a value to display (default 0.14)
    b_max : float
        Maximum b value to display (default 1.2)
    na : int
        Number of a grid points (default 300)
    nb : int
        Number of b grid points (default 400)
    """
    a_vals = np.linspace(0, a_max, na)
    b_vals = np.linspace(0, b_max, nb)
    A, B = np.meshgrid(a_vals, b_vals, indexing="xy")

    # trace on grid
    Tau = np.zeros_like(A)
    for i in range(nb):
        for j in range(na):
            Tau[i, j] = trace_at_equilibrium(A[i, j], B[i, j])

    stable = Tau < 0  # stable if trace < 0 (det>0 for a>0)

    plt.figure(figsize=(8, 5))
    plt.pcolormesh(a_vals, b_vals, stable, shading="auto")

    # overlay Hopf curves
    b_m, b_p = hopf_curves(a_vals)
    plt.plot(a_vals, b_m, linewidth=2.0)
    plt.plot(a_vals, b_p, linewidth=2.0)

    plt.xlim(0, a_max)
    plt.ylim(0, b_max)
    plt.xlabel("a")
    plt.ylabel("b")
    plt.title("Stability of equilibrium in (a,b)-plane (trace criterion)")
    plt.grid(True, alpha=0.3)
    plt.show()

def stability_map_ab(a_min=0.0, a_max=0.14, b_min=0.0, b_max=1.2,
                     na=250, nb=250):
    """
    Plot 2D stability map showing max real eigenvalue at equilibrium.
    
    Color intensity represents max Re(λ):
    - Negative (blue): stable equilibrium
    - Positive (red): unstable equilibrium
    
    A contour line at max Re(λ) = 0 shows the bifurcation set.
    
    Parameters
    ----------
    a_min, a_max : float
        Range of a parameter to scan (default 0.0 to 0.14)
    b_min, b_max : float
        Range of b parameter to scan (default 0.0 to 1.2)
    na : int
        Number of a grid points (default 250)
    nb : int
        Number of b grid points (default 250)
    
    Returns
    -------
    a_vals : ndarray
        Grid of a values
    b_vals : ndarray
        Grid of b values
    max_real : ndarray
        Maximum real eigenvalue at each (a, b) point
    """
    a_vals = np.linspace(a_min, a_max, na)
    b_vals = np.linspace(b_min, b_max, nb)
    A, B = np.meshgrid(a_vals, b_vals, indexing="xy")

    max_real = np.zeros_like(A, dtype=float)

    for i in range(nb):
        for j in range(na):
            a = A[i, j]
            b = B[i, j]
            if a == 0 and b == 0:
                max_real[i, j] = np.nan
                continue
            (_, _), _, eig = get_equilibrium_eigvals(a, b)
            max_real[i, j] = np.max(np.real(eig))

    plt.figure(figsize=(8,5))
    plt.pcolormesh(a_vals, b_vals, max_real, shading="auto")
    plt.colorbar(label="max Re(eigenvalue at equilibrium)")
    plt.contour(a_vals, b_vals, max_real, levels=[0.0], linewidths=2)
    plt.xlabel("a")
    plt.ylabel("b")
    plt.title("2D stability map of equilibrium in (a,b)")
    plt.grid(True, alpha=0.25)
    plt.show()

    return a_vals, b_vals, max_real


# -----------------------------------------------------------------------------
# Glucose system functions part 5
# -----------------------------------------------------------------------------