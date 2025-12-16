"""
Helper functions for analyzing the radial symmetry system:
dr/dt = μr + r³ - r⁵ = r(μ + r² - r⁴)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ====== 1D Radial Dynamics ======

def f_radial(r, mu):
    """
    Radial equation: dr/dt = μr + r³ - r⁵
    
    Parameters
    ----------
    r : float or array
        Radial coordinate
    mu : float
        System parameter
        
    Returns
    -------
    float or array
        dr/dt
    """
    return mu * r + r**3 - r**5


def df_radial(r, mu):
    """
    Derivative of radial equation: d/dr[f(r,μ)]
    Used for stability analysis.
    
    Parameters
    ----------
    r : float or array
        Radial coordinate
    mu : float
        System parameter
        
    Returns
    -------
    float or array
        df/dr
    """
    return mu + 3*r**2 - 5*r**4


def find_equilibria(mu):
    """
    Find equilibria of the radial equation: r(μ + r² - r⁴) = 0
    
    Equilibria: r = 0 or r² = (1 + sqrt(1 + 4μ))/2 or r² = (1 - sqrt(1 + 4μ))/2
    
    Parameters
    ----------
    mu : float
        System parameter
        
    Returns
    -------
    list
        List of real non-negative equilibrium values
    """
    equilibria = [0.0]  # r = 0 is always an equilibrium
    
    # Solve r² - r⁴ + μ = 0, or r⁴ - r² - μ = 0
    # Let u = r², then u² - u - μ = 0
    discriminant = 1 + 4*mu
    
    if discriminant >= 0:
        u1 = (1 + np.sqrt(discriminant)) / 2
        u2 = (1 - np.sqrt(discriminant)) / 2
        
        # r² = u1 (always positive when discriminant >= 0)
        if u1 > 0:
            equilibria.extend([np.sqrt(u1), -np.sqrt(u1)])
        
        # r² = u2 (only positive when u2 > 0)
        if u2 > 0:
            equilibria.extend([np.sqrt(u2), -np.sqrt(u2)])
    
    # Return only non-negative r values, sorted
    return sorted([r for r in equilibria if r >= 0])


def classify_stability(r_eq, mu):
    """
    Classify stability of an equilibrium point.
    
    Parameters
    ----------
    r_eq : float
        Equilibrium point
    mu : float
        System parameter
        
    Returns
    -------
    str
        'stable' if df/dr < 0, 'unstable' if df/dr > 0, 'neutral' if df/dr = 0
    """
    derivative = df_radial(r_eq, mu)
    
    if abs(derivative) < 1e-10:
        return 'neutral'
    elif derivative < 0:
        return 'stable'
    else:
        return 'unstable'


# ====== 2D Cartesian Dynamics ======

def polar_to_cartesian_system(omega=1.0):
    """
    Convert polar system (r, θ) to Cartesian system (x, y).
    
    With dr/dt = f(r,μ) and dθ/dt = ω, we have:
    dx/dt = dr/dt * cos(θ) - r * sin(θ) * dθ/dt
    dy/dt = dr/dt * sin(θ) + r * cos(θ) * dθ/dt
    
    Substituting r = sqrt(x² + y²), cos(θ) = x/r, sin(θ) = y/r:
    dx/dt = f(r,μ) * x/r - ω * y
    dy/dt = f(r,μ) * y/r + ω * x
    
    Parameters
    ----------
    omega : float
        Angular velocity (default: 1.0)
        
    Returns
    -------
    function
        Function that takes (t, state, mu) and returns [dx/dt, dy/dt]
    """
    def cartesian_system(t, state, mu):
        x, y = state
        r = np.sqrt(x**2 + y**2)
        
        if r < 1e-10:  # Near origin, use limit
            dxdt = mu * x - omega * y
            dydt = mu * y + omega * x
        else:
            f_r = f_radial(r, mu)
            dxdt = f_r * x / r - omega * y
            dydt = f_r * y / r + omega * x
        
        return [dxdt, dydt]
    
    return cartesian_system


def integrate_trajectory(x0, y0, mu, t_span=(0, 50), omega=1.0, n_points=1000):
    """
    Integrate a trajectory in the 2D system.
    
    Parameters
    ----------
    x0, y0 : float
        Initial conditions
    mu : float
        System parameter
    t_span : tuple
        Time span for integration (t_start, t_end)
    omega : float
        Angular velocity
    n_points : int
        Number of time points
        
    Returns
    -------
    solution : OdeResult
        Solution object from solve_ivp
    """
    system = polar_to_cartesian_system(omega)
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    
    sol = solve_ivp(
        system, 
        t_span, 
        [x0, y0], 
        args=(mu,),
        t_eval=t_eval,
        method='RK45',
        rtol=1e-8,
        atol=1e-10
    )
    
    return sol


# ====== Visualization Functions ======

def plot_phase_portrait(mu, omega=1.0, ax=None, n_trajectories=8, t_max=30):
    """
    Plot phase portrait for given μ value.
    
    Parameters
    ----------
    mu : float
        System parameter
    omega : float
        Angular velocity
    ax : matplotlib axis
        Axis to plot on (creates new if None)
    n_trajectories : int
        Number of trajectories to plot
    t_max : float
        Maximum integration time
        
    Returns
    -------
    ax : matplotlib axis
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    
    # Find equilibria to determine plot range
    equilibria = find_equilibria(mu)
    if len(equilibria) > 1:
        r_max = max(equilibria) * 1.5
    else:
        r_max = 2.0
    
    # Plot equilibrium circles
    theta = np.linspace(0, 2*np.pi, 100)
    for r_eq in equilibria:
        if r_eq > 0:
            stability = classify_stability(r_eq, mu)
            color = 'green' if stability == 'stable' else 'red'
            linestyle = '-' if stability == 'stable' else '--'
            ax.plot(r_eq * np.cos(theta), r_eq * np.sin(theta), 
                   color=color, linestyle=linestyle, linewidth=2,
                   label=f'r={r_eq:.3f} ({stability})')
    
    # Plot origin
    ax.plot(0, 0, 'ko', markersize=8, label='r=0 (origin)')
    
    # Plot sample trajectories
    initial_radii = np.linspace(0.1, r_max, n_trajectories)
    
    for r0 in initial_radii:
        x0, y0 = r0, 0
        sol = integrate_trajectory(x0, y0, mu, t_span=(0, t_max), omega=omega)
        
        if sol.success:
            ax.plot(sol.y[0], sol.y[1], 'b-', alpha=0.6, linewidth=1)
            ax.plot(sol.y[0][0], sol.y[1][0], 'bo', markersize=4)
    
    ax.set_xlim(-r_max, r_max)
    ax.set_ylim(-r_max, r_max)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(f'Phase Portrait: μ = {mu:.3f}, ω = {omega}', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    ax.legend(loc='best', fontsize=10)
    
    return ax


def plot_multiple_portraits(mu_values, omega=1.0, figsize=(15, 5)):
    """
    Plot phase portraits for multiple μ values.
    
    Parameters
    ----------
    mu_values : list
        List of μ values to plot
    omega : float
        Angular velocity
    figsize : tuple
        Figure size
        
    Returns
    -------
    fig, axes
    """
    n = len(mu_values)
    fig, axes = plt.subplots(1, n, figsize=figsize)
    
    if n == 1:
        axes = [axes]
    
    for ax, mu in zip(axes, mu_values):
        plot_phase_portrait(mu, omega=omega, ax=ax)
    
    plt.tight_layout()
    return fig, axes


def plot_radial_time_evolution(x0, y0, mu, omega=1.0, t_max=50):
    """
    Plot r(t) time evolution for a trajectory.
    
    Parameters
    ----------
    x0, y0 : float
        Initial conditions
    mu : float
        System parameter
    omega : float
        Angular velocity
    t_max : float
        Maximum time
        
    Returns
    -------
    fig, ax
    """
    sol = integrate_trajectory(x0, y0, mu, t_span=(0, t_max), omega=omega)
    
    # Calculate r(t)
    r = np.sqrt(sol.y[0]**2 + sol.y[1]**2)
    
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(sol.t, r, 'b-', linewidth=2)
    ax.set_xlabel('Time t', fontsize=12)
    ax.set_ylabel('Radius r(t)', fontsize=12)
    ax.set_title(f'Radial Evolution: μ={mu:.3f}, r₀={np.sqrt(x0**2+y0**2):.3f}', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Mark equilibria
    equilibria = find_equilibria(mu)
    for r_eq in equilibria:
        if r_eq > 0:
            stability = classify_stability(r_eq, mu)
            color = 'green' if stability == 'stable' else 'red'
            ax.axhline(r_eq, color=color, linestyle='--', alpha=0.7,
                      label=f'r_eq={r_eq:.3f} ({stability})')
    
    ax.legend()
    
    return fig, ax


# ====== Analysis Functions ======

def analyze_equilibria(mu, verbose=True):
    """
    Analyze equilibria for a given μ value.
    
    Parameters
    ----------
    mu : float
        System parameter
    verbose : bool
        If True, print detailed analysis
        
    Returns
    -------
    dict
        Dictionary with equilibria and their properties
    """
    equilibria = find_equilibria(mu)
    results = []
    
    for r_eq in equilibria:
        stability = classify_stability(r_eq, mu)
        derivative = df_radial(r_eq, mu)
        
        results.append({
            'r': r_eq,
            'stability': stability,
            'derivative': derivative
        })
        
        if verbose:
            if r_eq == 0:
                print(f"r_eq = {r_eq:.4f} (origin): {stability} (f'(r_eq) = {derivative:.4f})")
            else:
                print(f"r_eq = {r_eq:.4f}: {stability} (f'(r_eq) = {derivative:.4f})")
    
    return results


def print_bifurcation_analysis():
    """
    Print analysis of the two critical bifurcation points.
    """
    mu_bifurcation_1 = -1/4
    mu_bifurcation_2 = 0
    
    print("=" * 60)
    print("BIFURCATION ANALYSIS")
    print("=" * 60)
    
    # At μ = -1/4
    print(f"\n1. At μ = {mu_bifurcation_1} (Saddle-node bifurcation):")
    eq_1 = find_equilibria(mu_bifurcation_1)
    print(f"   Equilibria: {eq_1}")
    for r_eq in eq_1:
        if r_eq > 0:
            deriv = df_radial(r_eq, mu_bifurcation_1)
            print(f"   r = {r_eq:.4f}: f'(r) = {deriv:.6f}")
    
    # At μ = 0
    print(f"\n2. At μ = {mu_bifurcation_2} (Transcritical/Hopf bifurcation):")
    eq_2 = find_equilibria(mu_bifurcation_2)
    print(f"   Equilibria: {eq_2}")
    for r_eq in eq_2:
        deriv = df_radial(r_eq, mu_bifurcation_2)
        print(f"   r = {r_eq:.4f}: f'(r) = {deriv:.6f}")
    
    print("\n" + "=" * 60)


def print_regime_analysis():
    """
    Print detailed analysis of all three regimes.
    """
    regimes = [
        {"mu": -0.5, "name": "Regime 1: μ < -1/4"},
        {"mu": -0.15, "name": "Regime 2: -1/4 < μ < 0"},
        {"mu": 0.3, "name": "Regime 3: μ > 0"}
    ]
    
    print("=" * 70)
    print("DETAILED REGIME ANALYSIS")
    print("=" * 70)
    
    for regime in regimes:
        mu = regime["mu"]
        name = regime["name"]
        
        print(f"\n{name}")
        print("-" * 70)
        
        equilibria = find_equilibria(mu)
        print(f"Equilibria: {[f'{r:.4f}' for r in equilibria]}")
        print()
        
        for i, r_eq in enumerate(equilibria):
            stability = classify_stability(r_eq, mu)
            deriv = df_radial(r_eq, mu)
            
            if r_eq == 0:
                label = "Origin"
            else:
                label = f"Circle {i}"
            
            print(f"  {label}: r = {r_eq:.4f}")
            print(f"    Stability: {stability}")
            print(f"    f'(r) = {deriv:.4f}")
            print()
    
    print("=" * 70)


def plot_time_evolution_comparison(mu, initial_radii, omega=1.0, t_max=50, figsize=(18, 4)):
    """
    Plot time evolution of radius for multiple initial conditions.
    
    Parameters
    ----------
    mu : float
        System parameter
    initial_radii : list
        List of initial radii to compare
    omega : float
        Angular velocity
    t_max : float
        Maximum integration time
    figsize : tuple
        Figure size
        
    Returns
    -------
    fig, axes
    """
    fig, axes = plt.subplots(1, len(initial_radii), figsize=figsize)
    
    if len(initial_radii) == 1:
        axes = [axes]
    
    for i, r0 in enumerate(initial_radii):
        x0, y0 = r0, 0
        sol = integrate_trajectory(x0, y0, mu, t_span=(0, t_max), omega=omega)
        
        # Calculate r(t)
        r_t = np.sqrt(sol.y[0]**2 + sol.y[1]**2)
        
        axes[i].plot(sol.t, r_t, 'b-', linewidth=2)
        axes[i].set_xlabel('Time t', fontsize=11)
        axes[i].set_ylabel('r(t)', fontsize=11)
        axes[i].set_title(f'r₀ = {r0:.2f}', fontsize=12)
        axes[i].grid(True, alpha=0.3)
        
        # Mark equilibria
        equilibria = find_equilibria(mu)
        for r_eq in equilibria:
            if r_eq > 0:
                stability = classify_stability(r_eq, mu)
                color = 'green' if stability == 'stable' else 'red'
                linestyle = '-' if stability == 'stable' else '--'
                axes[i].axhline(r_eq, color=color, linestyle=linestyle, alpha=0.7, linewidth=2)
    
    plt.suptitle(f'Radial Time Evolution for μ = {mu}', fontsize=14, y=1.02)
    plt.tight_layout()
    
    return fig, axes


def plot_bifurcation_diagram(mu_min=-0.6, mu_max=1.0, n_points=500, figsize=(12, 8)):
    """
    Create a bifurcation diagram showing equilibria vs μ.
    
    Parameters
    ----------
    mu_min, mu_max : float
        Range of μ values
    n_points : int
        Number of points to evaluate
    figsize : tuple
        Figure size
        
    Returns
    -------
    fig, ax
    """
    mu_range = np.linspace(mu_min, mu_max, n_points)
    
    # Store equilibria for each mu
    inner_unstable = []
    outer_stable = []
    
    for mu in mu_range:
        equilibria = find_equilibria(mu)
        
        # Classify nonzero equilibria
        nonzero_eq = [r for r in equilibria if r > 1e-6]
        
        if len(nonzero_eq) == 1:
            # Only one nonzero equilibrium (Regime 3: μ > 0)
            outer_stable.append(nonzero_eq[0])
            inner_unstable.append(np.nan)
        elif len(nonzero_eq) == 2:
            # Two nonzero equilibria (Regime 2: -1/4 < μ < 0)
            r1, r2 = sorted(nonzero_eq)
            
            # Check which is stable (inner is unstable, outer is stable)
            if classify_stability(r1, mu) == 'unstable':
                inner_unstable.append(r1)
                outer_stable.append(r2)
            else:
                inner_unstable.append(np.nan)
                outer_stable.append(np.nan)
        else:
            # No nonzero equilibria (Regime 1: μ < -1/4)
            inner_unstable.append(np.nan)
            outer_stable.append(np.nan)
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot origin - color depends on stability
    mu_neg = mu_range[mu_range < 0]
    mu_pos = mu_range[mu_range >= 0]
    ax.plot(mu_neg, np.zeros_like(mu_neg), 'g-', linewidth=3, label='Stable origin')
    ax.plot(mu_pos, np.zeros_like(mu_pos), 'r--', linewidth=3, label='Unstable origin')
    
    # Plot nonzero equilibria
    inner_unstable_arr = np.array(inner_unstable)
    outer_stable_arr = np.array(outer_stable)
    
    ax.plot(mu_range, outer_stable_arr, 'g-', linewidth=2.5, label='Stable limit cycle')
    ax.plot(mu_range, inner_unstable_arr, 'r--', linewidth=2.5, label='Unstable limit cycle')
    
    # Mark bifurcation points
    ax.axvline(-0.25, color='purple', linestyle=':', linewidth=2, alpha=0.7, 
               label='Saddle-node (μ=-1/4)')
    ax.axvline(0, color='orange', linestyle=':', linewidth=2, alpha=0.7, 
               label='Transcritical (μ=0)')
    
    # Add regime labels
    ax.text(-0.4, 1.4, 'Regime 1\nStable origin only', fontsize=11, ha='center', 
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.text(-0.12, 1.4, 'Regime 2\nBistability', fontsize=11, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    ax.text(0.5, 1.4, 'Regime 3\nGlobal limit cycle', fontsize=11, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    
    ax.set_xlabel('Parameter μ', fontsize=13)
    ax.set_ylabel('Equilibrium radius r', fontsize=13)
    ax.set_title('Bifurcation Diagram: Equilibria vs μ', fontsize=15, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=10)
    ax.set_ylim(-0.1, 1.6)
    ax.set_xlim(mu_min, mu_max)
    
    plt.tight_layout()
    
    return fig, ax
