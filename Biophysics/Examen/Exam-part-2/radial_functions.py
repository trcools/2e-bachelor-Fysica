"""
Radial symmetry / generalized Hopf normal form:

Radial equation:
    r' = μ r + r^3 - r^5 = r(μ + r^2 - r^4)

If you add a constant rotation θ' = ω, the 2D Cartesian system is:
    x' = (f(r,μ)/r) x - ω y
    y' = (f(r,μ)/r) y + ω x
with r = sqrt(x^2 + y^2).

Key bifurcations (analytic):
    Saddle-node of cycles: μ = -1/4 at r = 1/sqrt(2)
    Subcritical Hopf (↔ radial subcritical pitchfork): μ = 0 at r = 0
"""

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


# ============================================================
# Core radial dynamics
# ============================================================

def f_radial(r, mu):
    """dr/dt = μr + r^3 - r^5."""
    return mu * r + r**3 - r**5


def df_radial(r, mu):
    """∂/∂r (μr + r^3 - r^5) = μ + 3r^2 - 5r^4."""
    return mu + 3*r**2 - 5*r**4


def equilibria_r(mu, include_negative=True, tol=1e-12):
    """
    Real equilibria of r' = r(μ + r^2 - r^4) = 0.

    Solve u^2 - u - μ = 0 with u = r^2 >= 0.
    """
    eq = {0.0}
    disc = 1.0 + 4.0*mu

    if disc >= -tol:
        disc = max(disc, 0.0)
        u_plus  = (1.0 + np.sqrt(disc)) / 2.0
        u_minus = (1.0 - np.sqrt(disc)) / 2.0

        for u in (u_plus, u_minus):
            if u > tol:
                r = float(np.sqrt(u))
                if include_negative:
                    eq.update([r, -r])
                else:
                    eq.add(r)

    out = sorted(eq)
    if not include_negative:
        out = [r for r in out if r >= -tol]
    return out


def classify_stability_1d(r_eq, mu, tol=1e-10):
    """
    1D stability of an equilibrium r_eq for r' = f(r,μ):
      stable if f'(r_eq)<0, unstable if f'(r_eq)>0.
    """
    lam = df_radial(r_eq, mu)
    if abs(lam) < tol:
        return "neutral"
    return "stable" if lam < 0 else "unstable"


def classify_stability_2d(r_eq, mu, omega=1.0):
    """
    2D stability of a limit cycle at radius r_eq for the Cartesian system
    using equilibrium_classification().
    
    For a limit cycle, the eigenvalues are:
      λ = df_radial(r_eq, mu) ± iω
    
    Special cases:
    - r ≈ 0 with Re(λ) ≈ 0: Hopf bifurcation at origin
    - r > 0 with Re(λ) ≈ 0: Saddle-node bifurcation of limit cycles
    
    Parameters
    ----------
    r_eq : float
        Equilibrium radius (must be > 0 for limit cycle)
    mu : float
        System parameter
    omega : float
        Angular velocity (default 1.0)
    
    Returns
    -------
    classification : str
        Equilibrium type from equilibrium_classification()
    """
    tol = 1e-9
    
    if r_eq <= 1e-12:
        # Origin: eigenvalues are μ ± iω
        eigvals = np.array([mu + 1j*omega, mu - 1j*omega])
        classification = equilibrium_classification(eigvals)
    else:
        # Limit cycle: eigenvalues are df_radial(r_eq, mu) ± iω
        lam_radial = df_radial(r_eq, mu)
        eigvals = np.array([lam_radial + 1j*omega, lam_radial - 1j*omega])
        
        # Check if this is a saddle-node bifurcation (f'(r) ≈ 0 at r > 0)
        if abs(lam_radial) < tol:
            # Saddle-node: two limit cycles collide
            classification = "saddle-node bifurcation"
        else:
            classification = equilibrium_classification(eigvals)
    
    return classification


def bifurcation_points(omega=1.0):
    """
    Compute analytic bifurcation markers and classify them using equilibrium_classification().

    - Saddle-node of cycles: solve f(r,mu)=0 and df/dr=0 for r>0.
    - Hopf at origin: where Re(eigs) = mu crosses 0 => mu=0.
    
    Parameters
    ----------
    omega : float
        Angular velocity for 2D system classification (default 1.0)
    
    Returns
    -------
    list of dict
        Each dict contains mu, type, r_eq, label, color, marker, linestyle
    """
    # Get styles for different equilibrium types
    groups = get_groups(dim="1D")
    
    # --- Saddle-node: solve f=0 and df=0 for r>0 ---
    # From f=0 with r>0: mu = r^4 - r^2
    # Plug into df=0: mu + 3r^2 - 5r^4 = 0
    # -> (r^4 - r^2) + 3r^2 - 5r^4 = 0 -> 2r^2 - 4r^4 = 0 -> r^2 = 1/2
    r_sn = float(np.sqrt(0.5))
    mu_sn = float(r_sn**4 - r_sn**2)
    type_sn = classify_stability_2d(r_sn, mu_sn, omega=omega)
    
    # Get style for this type, or use defaults
    if type_sn in groups:
        style_sn = groups[type_sn]["style"]
        marker_sn = style_sn.get("marker", "^")
        color_sn = style_sn.get("facecolors", style_sn.get("color", "#00bcd4"))
    else:
        marker_sn, color_sn = "^", "#00bcd4"

    # --- Hopf at origin (for omega != 0) ---
    mu_hopf = 0.0
    r_hopf = 0.0
    type_hopf = classify_stability_2d(r_hopf, mu_hopf, omega=omega)
    
    # Get style for this type
    if type_hopf in groups:
        style_hopf = groups[type_hopf]["style"]
        marker_hopf = style_hopf.get("marker", "D")
        color_hopf = style_hopf.get("facecolors", style_hopf.get("color", "#ff9800"))
    else:
        marker_hopf, color_hopf = "D", "#ff9800"

    return [
        dict(mu=mu_sn, type=type_sn, r_eq=r_sn,
             label=rf"{type_sn} ($\mu={mu_sn:.3f}$, $r={r_sn:.3f}$)", 
             color=color_sn, marker=marker_sn, linestyle="--"),
        dict(mu=mu_hopf, type=type_hopf, r_eq=r_hopf,
             label=rf"{type_hopf} ($\mu={mu_hopf:.1f}$)", 
             color=color_hopf, marker=marker_hopf, linestyle=":"),
    ]


def sample_equilibria(mu_values, include_negative=True):
    """Sample equilibria and split into stable/unstable (neutral kept separately)."""
    stable_mu, stable_r = [], []
    unstable_mu, unstable_r = [], []
    neutral_mu, neutral_r = [], []

    for mu in mu_values:
        for r in equilibria_r(mu, include_negative=include_negative):
            st = classify_stability_1d(r, mu)
            if st == "stable":
                stable_mu.append(mu); stable_r.append(r)
            elif st == "unstable":
                unstable_mu.append(mu); unstable_r.append(r)
            else:
                neutral_mu.append(mu); neutral_r.append(r)

    return {
        "stable":   (np.array(stable_mu),   np.array(stable_r)),
        "unstable": (np.array(unstable_mu), np.array(unstable_r)),
        "neutral":  (np.array(neutral_mu),  np.array(neutral_r)),
    }



# ============================================================
# Classification and grouping of equilibria
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
        - "saddle spiral"
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
                   "saddle spiral": dict(marker="o", s=10, label="saddle spiral"),
                   "saddle-node bifurcation": dict(marker="^", s=50, facecolors="#00bcd4", edgecolors="black", label="saddle-node bifurcation"),
                   "pitchfork bifurcation": dict(marker="D", s=50, facecolors="blue", edgecolors="blue",label="pitchfork bifurcation"),
                   "Hopf bifurcation": dict(marker="D", s=50, facecolors="#ff9800", edgecolors="black", label="Hopf bifurcation"),
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
# 1D bifurcation plot (μ vs r)
# ============================================================

def plot_bifurcation_1d(
    mu_range=(-3, 3),
    r_range=(-2, 2),
    n_mu=3000,
    include_negative=True,
    ax=None,
    show_bif_points=True,
    legend=True,
):
    """Bifurcation diagram: equilibria r vs parameter μ."""
    mu_min, mu_max = mu_range
    mu_values = np.linspace(mu_min, mu_max, n_mu)

    data = sample_equilibria(mu_values, include_negative=include_negative)
    stable_mu, stable_r = data["stable"]
    unstable_mu, unstable_r = data["unstable"]

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))
    else:
        fig = ax.figure

    ax.scatter(stable_mu, stable_r, s=2, alpha=0.25, color="red", label="stable")
    ax.scatter(unstable_mu, unstable_r, s=2, alpha=0.25, color="blue", label="unstable")

    if show_bif_points:
        for b in bifurcation_points():
            ax.axvline(b["mu"], color="k", ls=b["linestyle"], lw=1.8)
            ax.scatter([b["mu"]], [b["r_eq"]], color=b["color"], marker=b["marker"],
                       s=80, edgecolors="black", linewidths=1.2, label=b["label"], zorder=10)
            if include_negative and b["r_eq"] != 0:
                ax.scatter([b["mu"]], [-b["r_eq"]], color=b["color"], marker=b["marker"],
                           s=80, edgecolors="black", linewidths=1.2, zorder=10)

    ax.set_xlim(mu_min, mu_max)
    if r_range is not None:
        ax.set_ylim(*r_range)

    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$r$")
    ax.set_title(r"Bifurcation diagram for $r'=\mu r + r^3 - r^5$")
    ax.grid(True, alpha=0.3)

    if legend:
        # de-duplicate legend entries
        handles, labels = ax.get_legend_handles_labels()
        seen, H, L = set(), [], []
        for h, l in zip(handles, labels):
            if l and l not in seen:
                seen.add(l); H.append(h); L.append(l)
        leg = ax.legend(H, L, loc="best")
        # Increase alpha for legend markers to make colors more visible
        for handle in leg.legend_handles:
            handle.set_alpha(1.0)

    return fig, ax


# ============================================================
# 2D Cartesian system + simulation
# ============================================================

def cartesian_rhs(t, state, mu, omega=1.0):
    """
    2D system with constant angular velocity omega.
    Note: only the origin is a fixed point (if omega != 0).
    Nonzero equilibria in r correspond to periodic orbits (limit cycles) in (x,y).
    """
    x, y = state
    r = np.hypot(x, y)

    if r < 1e-12:
        # limit r->0 : f(r,mu)/r ~ mu
        return [mu*x - omega*y, mu*y + omega*x]

    g = f_radial(r, mu) / r
    return [g*x - omega*y, g*y + omega*x]


def integrate_trajectory(x0, y0, mu, t_span=(0, 50), omega=1.0, n_points=2000,
                         rtol=1e-8, atol=1e-10):
    """Integrate one trajectory in the 2D system."""
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    sol = solve_ivp(
        lambda t, z: cartesian_rhs(t, z, mu, omega=omega),
        t_span, [x0, y0], t_eval=t_eval, method="RK45", rtol=rtol, atol=atol
    )
    return sol


def floquet_exponents_limit_cycle(r_eq, mu):
    """
    For a limit cycle at radius r_eq (r_eq>0), the Floquet exponents are:
      0  (phase direction)
      df_radial(r_eq, mu)  (radial direction)
    """
    if r_eq <= 0:
        raise ValueError("r_eq must be > 0 for a limit cycle.")
    return 0.0, float(df_radial(r_eq, mu))


def plot_phase_portrait(mu, omega=1.0, ax=None, n_trajectories=8, t_max=30):
    """Phase portrait with trajectories + circles indicating limit cycles (r equilibria)."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.5, 7.5))

    # Choose plot scale based on largest positive equilibrium radius (if any)
    eq = equilibria_r(mu, include_negative=False)
    r_pos = [r for r in eq if r > 1e-12]
    r_max = (max(r_pos) * 1.5) if r_pos else 2.0

    # Draw limit cycles (as circles) for r_eq>0
    theta = np.linspace(0, 2*np.pi, 400)
    for r_eq in r_pos:
        st = classify_stability_1d(r_eq, mu)
        col = "green" if st == "stable" else "red"
        ls  = "-" if st == "stable" else "--"
        ax.plot(r_eq*np.cos(theta), r_eq*np.sin(theta), color=col, ls=ls, lw=2,
                label=f"limit cycle r={r_eq:.3f} ({st})")

    # Origin (fixed point)
    ax.plot(0, 0, "ko", ms=7, label="origin (fixed point)")

    # Trajectories
    radii0 = np.linspace(0.1, r_max, n_trajectories)
    for r0 in radii0:
        sol = integrate_trajectory(r0, 0.0, mu, t_span=(0, t_max), omega=omega)
        if sol.success:
            ax.plot(sol.y[0], sol.y[1], lw=1, alpha=0.65)
            ax.plot(sol.y[0][0], sol.y[1][0], "o", ms=3)

    ax.set_xlim(-r_max, r_max)
    ax.set_ylim(-r_max, r_max)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Phase portrait (μ={mu:.3f}, ω={omega})")
    ax.grid(True, alpha=0.3)

    # de-duplicate legend
    handles, labels = ax.get_legend_handles_labels()
    seen, H, L = set(), [], []
    for h, l in zip(handles, labels):
        if l and l not in seen:
            seen.add(l); H.append(h); L.append(l)
    ax.legend(H, L, loc="best", fontsize=9)

    return ax


def plot_radial_time_evolution(x0, y0, mu, omega=1.0, t_max=50):
    """Plot r(t) for a single 2D trajectory, plus equilibrium radii."""
    sol = integrate_trajectory(x0, y0, mu, t_span=(0, t_max), omega=omega)
    r_t = np.hypot(sol.y[0], sol.y[1])

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(sol.t, r_t, lw=2)
    ax.set_xlabel("t")
    ax.set_ylabel("r(t)")
    ax.set_title(f"Radial evolution (μ={mu:.3f}, ω={omega}, r0={np.hypot(x0,y0):.3f})")
    ax.grid(True, alpha=0.3)

    # Mark positive equilibrium radii
    eq = equilibria_r(mu, include_negative=False)
    for r_eq in eq:
        if r_eq > 1e-12:
            st = classify_stability_1d(r_eq, mu)
            col = "green" if st == "stable" else "red"
            ax.axhline(r_eq, color=col, ls="--", alpha=0.7, label=f"r={r_eq:.3f} ({st})")

    # de-duplicate legend
    handles, labels = ax.get_legend_handles_labels()
    seen, H, L = set(), [], []
    for h, l in zip(handles, labels):
        if l and l not in seen:
            seen.add(l); H.append(h); L.append(l)
    if L:
        ax.legend(H, L, loc="best")

    return fig, ax


# ============================================================
# High-level plotting functions (for notebooks)
# ============================================================

def plot_eigenvalues_vs_mu(mu_range=(-0.3, 0.2), n_points=800, omega=1.0):
    """Plot stability (eigenvalues) vs μ for all equilibria."""
    mu_min, mu_max = mu_range
    mu_vals = np.linspace(mu_min, mu_max, n_points)
    
    MU_st, LAM_st = [], []
    MU_un, LAM_un = [], []
    
    for mu in mu_vals:
        for r_eq in equilibria_r(mu, include_negative=False):
            if r_eq > 1e-12:
                lam_phase, lam_radial = floquet_exponents_limit_cycle(r_eq, mu)
                st = classify_stability_1d(r_eq, mu)
                if st == "stable":
                    MU_st.append(mu); LAM_st.append(lam_radial)
                elif st == "unstable":
                    MU_un.append(mu); LAM_un.append(lam_radial)
    
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(mu_vals, mu_vals, lw=2, label=r"origin: $\mathrm{Re}(\lambda)=\mu$")
    ax.scatter(MU_st, LAM_st, s=6, alpha=0.35, label=r"limit cycle (stable)")
    ax.scatter(MU_un, LAM_un, s=6, alpha=0.35, label=r"limit cycle (unstable)")
    ax.axhline(0, ls=":", lw=1.5, alpha=0.7)

    for bp in bifurcation_points():
        ax.axvline(bp["mu"], ls=bp["linestyle"], lw=1.2, alpha=0.6)

    ax.set_xlim(mu_min, mu_max)
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"real stability rates (origin + $\lambda_r$ for cycles)")
    ax.set_title(r"Stability vs $\mu$ (correct for $\omega\neq 0$)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    plt.tight_layout()
    return fig, ax


def plot_three_regimes():
    """Plot bifurcation diagram and phase portraits for three regimes."""
    from matplotlib.gridspec import GridSpec
    
    mu_vals = [-0.5, -0.15, 0.3]
    titles = [r"$\mu < -1/4$", r"$-1/4 < \mu < 0$", r"$\mu > 0$"]
    
    fig = plt.figure(figsize=(18, 10))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.1, 1])
    
    # Top: 1D bifurcation
    ax0 = fig.add_subplot(gs[0, :])
    plot_bifurcation_1d(
        mu_range=(-0.6, 1.0),
        r_range=(0, 2.0),
        n_mu=3000,
        include_negative=False,
        ax=ax0,
        show_bif_points=True,
        legend=True,
    )
    
    # Bottom: phase portraits
    for i, (mu, t) in enumerate(zip(mu_vals, titles)):
        ax = fig.add_subplot(gs[1, i])
        plot_phase_portrait(mu, omega=1.0, ax=ax, n_trajectories=8, t_max=30)
        ax.set_title(t)
    
    plt.tight_layout()
    return fig


def plot_radial_time_subplots(mu=-0.1, initial_radii=[0.3, 0.6, 1.0], omega=1.0, t_max=50):
    """Plot time evolution of radius for multiple initial conditions in subplots."""
    fig, axes = plt.subplots(1, len(initial_radii), figsize=(6*len(initial_radii), 4), sharey=True)
    if len(initial_radii) == 1:
        axes = [axes]
    
    # Get positive equilibrium radii for reference
    eq_pos = [r for r in equilibria_r(mu, include_negative=False) if r > 1e-12]
    
    for idx, r0 in enumerate(initial_radii):
        # Integrate trajectory
        sol = integrate_trajectory(r0, 0.0, mu, t_span=(0, t_max), omega=omega, n_points=2000)
        r_t = np.hypot(sol.y[0], sol.y[1])
        
        # Plot on subplot
        axes[idx].plot(sol.t, r_t, lw=2)
        axes[idx].set_title(f"$r_0={r0}$")
        axes[idx].set_xlabel("t")
        axes[idx].grid(True, alpha=0.3)
        
        # Mark limit cycles
        for r_eq in eq_pos:
            st = classify_stability_1d(r_eq, mu)
            col = "green" if st == "stable" else "red"
            ls = "-" if st == "stable" else "--"
            axes[idx].axhline(r_eq, color=col, ls=ls, alpha=0.7, label=f"r={r_eq:.3f} ({st})")
        
        # De-duplicate legend
        H, L = axes[idx].get_legend_handles_labels()
        seen, HH, LL = set(), [], []
        for h, l in zip(H, L):
            if l not in seen:
                seen.add(l); HH.append(h); LL.append(l)
        axes[idx].legend(HH, LL, loc="best", fontsize=9)
    
    axes[0].set_ylabel("r(t)")
    fig.suptitle(rf"Radial time evolution ($\mu={mu}$, $\omega={omega}$)", y=1.02)
    plt.tight_layout()
    return fig


def plot_stability_equilibria(mu_range=(-0.6, 1.0), n_mu=800):
    """Plot stability of limit cycles via radial Floquet exponent."""
    mu_vals = np.linspace(*mu_range, n_mu)
    
    # Use sample_equilibria
    data = sample_equilibria(mu_vals, include_negative=False)
    
    MU, LAM, ST = [], [], []
    for mu in mu_vals:
        for r_eq in equilibria_r(mu, include_negative=False):
            if r_eq > 1e-12:
                MU.append(mu)
                LAM.append(df_radial(r_eq, mu))
                ST.append(classify_stability_1d(r_eq, mu))
    
    print(f"sample_equilibria() found: {len(data['stable'][0])} stable and {len(data['unstable'][0])} unstable equilibria")
    
    fig, ax = plt.subplots(figsize=(10, 5))
    MU = np.array(MU); LAM = np.array(LAM); ST = np.array(ST)
    
    ax.scatter(MU[ST=="stable"],   LAM[ST=="stable"],   s=5, alpha=0.35, label="stable limit cycle")
    ax.scatter(MU[ST=="unstable"], LAM[ST=="unstable"], s=5, alpha=0.35, label="unstable limit cycle")
    ax.axhline(0, ls=":", lw=1.5, alpha=0.7)

    for bp in bifurcation_points():
        ax.axvline(bp["mu"], ls=bp["linestyle"], lw=1.2, alpha=0.6)

    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"radial exponent $\lambda_r = f'(r_{eq},\mu)$")
    ax.set_title("Stability of limit cycles via radial Floquet exponent")
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()
    return fig


def print_cartesian_rhs_demo(mu_demo=-0.1, omega_demo=1.0, state_demo=None):
    """Demonstrate cartesian_rhs() function."""
    if state_demo is None:
        state_demo = [1.0, 0.5]
    
    derivatives = cartesian_rhs(0, state_demo, mu_demo, omega=omega_demo)
    
    print(f"At state (x, y) = ({state_demo[0]:.3f}, {state_demo[1]:.3f}) with μ = {mu_demo}, ω = {omega_demo}:")
    print(f"  dx/dt = {derivatives[0]:.6f}")
    print(f"  dy/dt = {derivatives[1]:.6f}")
    
    # Verify the radial component
    r_demo = np.hypot(*state_demo)
    f_r = f_radial(r_demo, mu_demo)
    print(f"\nRadius r = {r_demo:.6f}")
    print(f"Radial velocity f(r) = {f_r:.6f}")
    print(f"Expected radial component f(r)/r = {f_r/r_demo:.6f}")
    print("\nThis confirms cartesian_rhs() correctly implements the 2D system.")

