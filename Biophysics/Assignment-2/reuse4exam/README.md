# Rössler System Analysis

This folder contains adapted code from the Lorenz system analysis to work with the Rössler system.

## The Rössler Equations

The Rössler system is described by the following differential equations:

$$
\begin{align}
\dot{x} &= -y - z \\
\dot{y} &= x + ay \\
\dot{z} &= b + z(x - c)
\end{align}
$$

## Parameters

- **a**: Fixed parameter (default: 0.2)
- **b**: Fixed parameter (default: 0.2)
- **c**: Bifurcation parameter (varied to study different dynamical regimes)

## Main Differences from Lorenz System

| Aspect | Lorenz | Rössler |
|--------|---------|---------|
| Equations | $\dot{X} = \sigma(Y-X)$<br>$\dot{Y} = X(\rho-Z)-Y$<br>$\dot{Z} = XY-\beta Z$ | $\dot{x} = -y-z$<br>$\dot{y} = x+ay$<br>$\dot{z} = b+z(x-c)$ |
| Bifurcation parameter | $\rho$ | $c$ |
| Fixed parameters | $\sigma$, $\beta$ | $a$, $b$ |
| Typical range | $\rho \in [0, 28]$ | $c \in [0, 18]$ |
| Attractor shape | Two-lobed butterfly | Single-lobed scroll |
| Jacobian | $\begin{pmatrix}-\sigma & \sigma & 0 \\ \rho-Z & -1 & -X \\ Y & X & -\beta\end{pmatrix}$ | $\begin{pmatrix}0 & -1 & -1 \\ 1 & a & 0 \\ z & 0 & x-c\end{pmatrix}$ |

## Files

- **parameters.py**: Defines system parameters (a, b, initial conditions, time settings)
- **imports.py**: All necessary Python library imports
- **data.py**: Core functions for numerical integration, equilibrium finding, and stability analysis
- **figures.py**: All plotting and visualization functions
- **Rossler_Analysis.ipynb**: Main notebook demonstrating the Rössler system analysis

## Key Features

### Dynamical Regimes

The Rössler system exhibits different behaviors depending on the value of c:

- **c < 2**: Simple fixed point behavior
- **c ≈ 2-4**: Periodic orbits
- **c ≈ 4-6**: Period-doubling route to chaos
- **c > 5.7**: Chaotic behavior (Rössler attractor)

### Analyses Included

1. **Equilibrium Analysis**: Finding and classifying fixed points
2. **Bifurcation Diagrams**: Tracking equilibria as c varies
3. **3D Phase Space Visualization**: The characteristic scroll attractor
4. **Time Series**: Evolution of x(t), y(t), z(t)
5. **Return Maps**: Z_n vs Z_{n+1} for local maxima
6. **Sensitivity to Initial Conditions**: Butterfly effect demonstration
7. **Stability Analysis**: Perturbed trajectories near equilibria

## Usage

```python
# Import the modules
from parameters import get_parameters
from data import get_vectors, get_solutions
from figures import show_3Dfigure_plot, show_bifurcation_diagram

# Get parameters
a, b, t_min, t_max, dt, initial_condition, initial_conditions = get_parameters()

# Generate trajectories for c = 5.7
x, y, z = get_vectors(c=5.7, t_min=0, t_max=200)

# Visualize the attractor
show_3Dfigure_plot(c=5.7)

# Show bifurcation diagram
show_bifurcation_diagram()
```

## Comparison with Original Lorenz Code

All function names have been adapted:
- `lorenz_rhs` → `rossler_rhs`
- `lorenz_chaos` → `rossler_chaos`
- `rho` → `c` (bifurcation parameter)
- `sigma, beta` → `a, b` (fixed parameters)
- References to "Lorenz" in titles/labels → "Rössler"

The code structure remains the same to facilitate comparison and understanding.

## References

- Rössler, O.E. (1976). "An equation for continuous chaos". Physics Letters A. 57 (5): 397–398.
- The Lorenz system code structure from `../` (parent directory)
