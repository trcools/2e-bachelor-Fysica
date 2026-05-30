# AI Coding Agent Instructions

## Project Overview

This is an **academic physics and computational science workspace** containing:
- **Py4sci/**: Python for Scientific Computing course materials (primary computational content)
- **Biophysics/**: Biophysics assignments with modular code structure
- **Kwantum/**, **Chemie/**: Quantum mechanics and chemistry course materials
- **1e-semester-1e-bachelor/Lam/**: Linear algebra materials from first-year bachelor
- Study notes in Dutch (samenvattingen) and Jupyter notebooks for demonstrations

This is NOT a software project—it's educational material combining theory (markdown), computation (notebooks), and assignments.

## Architecture & Organization

### Py4sci Structure (Primary Computational Content)
The Py4sci folder follows a pedagogical structure organized by topics in numerical methods:

```
Py4sci/
├── notebooks/              # Core teaching notebooks (self-contained)
│   ├── numerical_limitations.ipynb
│   ├── linear_systems.ipynb
│   ├── nonlinear_equations.ipynb
│   ├── eigenvalue_problems.ipynb
│   ├── integration_differentiation.ipynb
│   ├── monte_carlo.ipynb
│   └── ...
├── Samenvatting-cursus/    # Study summaries (Dutch) organized by chapter
│   ├── ALGEMENE-SAMENVATTING.md  # Master summary document
│   ├── README.md                  # Chapter-by-chapter summaries
│   └── Hoofdstuk-XX/              # Individual chapter summaries
├── practical-sessions-main/
│   ├── basic-exercises/           # Basic exercises with solutions
│   └── advanced-exercises/        # Advanced exercises
└── Cursus/                 # Course materials, slides, datasets
```

**Key Pattern**: Notebooks in `Py4sci/notebooks/` are **self-contained teaching documents** that combine:
- Theory with LaTeX math (KaTeX rendering)
- Illustrative code examples with inline function definitions
- Matplotlib visualizations (often with `plt.close("figurename")` to manage figure windows)
- Direct algorithm implementations showing pedagogical steps

### Biophysics Modular Code Pattern
Biophysics assignments (e.g., `Assignment-2/`) demonstrate a **modular separation** pattern:

```python
# Structure:
imports.py      # All library imports consolidated
parameters.py   # Global parameters via get_parameters()
data.py         # Core analysis functions (integration, symbolic math, eigenvalues)
figures.py      # All plotting functions
*.ipynb         # Clean notebooks that import from the above modules
```

**Pattern Example** (from `Assignment-2/`):
```python
# figures.py imports from parameters.py and data.py
from parameters import get_parameters
from data import get_equilibria_for_rho, compute_jacobian_and_stability

# Notebooks then cleanly import specific functions
from figures import plot_bifuractions, show_3Dfigure_plot
```

**Why this matters**: When working on similar assignments, follow this modular pattern rather than inline everything in notebooks.

## Python & Scientific Computing Conventions

### Standard Import Pattern
Notebooks consistently use this import block:
```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize, integrate  # or other scipy modules
```

Additional common imports:
- `matplotlib.animation.FuncAnimation` for animations
- `scipy.optimize.root`, `.root_scalar`, `.brentq` for root finding
- `ipywidgets.interact` for interactive demonstrations

### Matplotlib Visualization Conventions

**Figure Management Pattern** (used throughout `Py4sci/notebooks/`):
```python
def plot_function_name():
    plt.close("unique_name")  # Close any existing figure with this name
    fig, ax = plt.subplots(num="unique_name", figsize=(8, 5))
    # ... plotting code ...
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()
```

**Why**: The `plt.close("name")` + `num="name"` pattern prevents duplicate figures when re-running cells.

**Animation Pattern**:
```python
from matplotlib.animation import FuncAnimation

def plot_and_animate():
    fig, ax = plt.subplots(num="animate_name", clear=True)
    # Setup plot
    (line,) = plt.plot([], [], 'o', markersize=10)
    
    def animate(i):
        # Update data
        line.set_data([x], [y])
        return (line,)
    
    return FuncAnimation(fig, animate, frames=20, interval=1000, repeat=False)
```

**Common styling**:
- Black lines for axes: `ax.axhline(0, color="black")`
- Default color cycle: `"blue"`, `"red"`, or `"C0"`, `"C1"` for automatic cycling
- Marker sizes: typically `markersize=10` for scatter points

### Numerical Methods Code Style

**Teaching-First Approach**: Code prioritizes **clarity over optimization**:
```python
def newton_method(f, fp, x, niter):
    """Illustrative implementation of Newton's method.
    
    Parameters
    ----------
    f : callable
        Function to be rooted.
    fp : callable
        Derivative of f.
    x : float
        Initial guess.
    niter : int
        Number of iterations.
    
    Returns
    -------
    float
        Approximation of the root.
    """
    for _ in range(niter):
        x = x - f(x) / fp(x)
    return x
```

**Pattern Notes**:
- Docstrings use NumPy-style format
- Simple iteration loops (not vectorized) when showing algorithmic steps
- Return intermediate values for visualization (e.g., `history` lists)

### SciPy Integration
After showing educational implementations, notebooks demonstrate production tools:
```python
# After custom implementation, show the SciPy equivalent:
optimize.root_scalar(func_cube, method="newton", x0=25, 
                    fprime=func_cube_prime, xtol=1e-16, maxiter=500)
```

## Language & Documentation

- **Primary Language**: English for code/docstrings, Dutch for study notes (samenvattingen)
- **Math Notation**: Use LaTeX inline `$...$` or block `$$...$$` (renders as KaTeX in notebooks)
- **Comments**: Prefer clear variable names over excessive comments
- **Dutch Terms to Know**:
  - "samenvatting" = summary
  - "hoofdstuk" = chapter
  - "oefeningen" = exercises

## Working with Notebooks

### Cell Execution Pattern
- **Setup cells**: Import statements at top
- **Theory cells**: Markdown with LaTeX math
- **Demo cells**: Function definition + immediate execution
- **Visualization cells**: Often return animation objects or display figures

### Key Notebook Features
- Use `interact` from `ipywidgets` for interactive parameter exploration
- Figures often include multiple subplots: `fig, axs = plt.subplots(2, 3, figsize=(...))`
- Cell execution creates inline outputs (matplotlib figures auto-display)

## File Naming & Structure

- **Notebooks**: `snake_case.ipynb` (e.g., `nonlinear_equations.ipynb`)
- **Python modules**: `snake_case.py` 
- **Markdown summaries**: Title case with hyphens (e.g., `Algemene-samenvatting-chemie.md`)
- **Folders**: Either `Title-Case/` or `snake_case/` depending on context

## Common Tasks

### Adding New Visualizations
Follow the established pattern:
```python
def plot_new_concept():
    plt.close("new_concept")
    fig, ax = plt.subplots(num="new_concept", figsize=(8, 5))
    # Implementation
    ax.legend()

plot_new_concept()
```

### Creating Modular Assignment Code
Use Biophysics pattern:
1. `imports.py` - consolidate all imports
2. `parameters.py` - define `get_parameters()` function
3. `data.py` - core computational functions
4. `figures.py` - all plotting functions (imports from above)
5. Notebook imports final functions for clean presentation

### Working with Scientific Equations
- Always show the mathematical formulation first in markdown
- Then implement clearly with matching variable names
- Example from notebooks:
  ```markdown
  $$\mathbf{f}(\mathbf{x})=\mathbf{0}$$
  ```
  ```python
  def func_2d(x):
      return [x[0] + 2*x[1] - 2, x[0]**2 + 4*x[1]**2 - 4]
  ```

## Development Environment

- **IDE**: VS Code (current environment)
- **Python Version**: Likely 3.8+ (uses modern f-strings, type hints occasionally)
- **Key Libraries**: NumPy, SciPy, Matplotlib, SymPy (for symbolic math), ipywidgets
- **No build system**: Direct notebook execution, no compilation steps

## What NOT to Do

- ❌ Don't optimize teaching code—keep it pedagogically clear
- ❌ Don't remove `plt.close()` calls—they prevent figure duplication
- ❌ Don't consolidate notebooks—they're meant to be standalone teaching units
- ❌ Don't expect production code patterns—this is educational material
- ❌ Don't assume English everywhere—study notes are in Dutch

## Quick Reference

**Running notebooks**: Execute cells sequentially in Jupyter/VS Code
**Reusing code**: Import from modular `.py` files in assignment folders
**Math rendering**: Use `$...$` for inline, `$$...$$` for block equations
**Figures**: Always use named figures with `plt.close("name")` first
**SciPy methods**: After custom implementation, show production equivalent





Some of the above content may be a bit outdated or not perfectly aligned with the current state of the repository, but it should give you a solid understanding of the overall structure, patterns, and conventions used in this workspace. If you have any specific questions about certain files or need clarification on any part of the codebase, feel free to ask!
