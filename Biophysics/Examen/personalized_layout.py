"""
Personalized Layout Utilities

This module provides reusable helper functions for creating interactive plots
and managing legends across different exam parts of the Biophysics course.

Author: Tristan Cools
Date: January 2026

Note:
This module was refactored with assistance from an AI tool to improve structure,
naming, and documentation for readability and maintainability.

"""

# =============================================================================
# Imports
# =============================================================================

from ipywidgets import interact, FloatSlider, HBox, interactive_output, Checkbox, Output, VBox
from IPython.display import display
from collections import OrderedDict
import matplotlib.pyplot as plt
import time

# =============================================================================
# Interactive Plot Helper Function
# =============================================================================

def create_interactive_plot(plot_func, slider_configs, n_cols=3, extra_widgets=None, extra_widgets_inline=False):
    """
    Create an interactive plot with dynamically generated sliders and optional extra widgets.
    
    This helper function eliminates code duplication for interactive plotting
    by handling widget creation, layout, and output display.
    
    Parameters
    ----------
    plot_func : callable
        Function that takes widget values as keyword arguments and creates the plot.
        Example: def my_plot(a_val, b_val, c_val, x0_val, ...): ...
    
    slider_configs : dict
        Dictionary mapping slider names to their configurations.
        Each value should be a dict with keys: 'value', 'min', 'max', 'step'
        Example: {'a': {'value': 0.2, 'min': 0.05, 'max': 0.4, 'step': 0.05}, ...}
    
    n_cols : int, optional
        Number of sliders to display per row (default: 3)
    
    extra_widgets : dict, optional
        Dictionary mapping widget names to widget objects (e.g., Checkbox, ToggleButton).
        These widgets will be displayed after the sliders and included in the interactive output.
        Example: {'show_arrows': Checkbox(value=True, description='Show arrows')}
    
    extra_widgets_inline : bool, optional
        When True, display all extra widgets in a single horizontal row (HBox).
        When False (default), display each widget on its own line.
    
    Returns
    -------
    None
        Displays the interactive plot in the notebook cell.
    
    Example
    -------
    >>> slider_configs = {
    ...     'a': {'value': 0.2, 'min': 0.05, 'max': 0.4, 'step': 0.05},
    ...     'b': {'value': 0.2, 'min': 0.05, 'max': 0.4, 'step': 0.05},
    ... }
    >>> checkbox = Checkbox(value=True, description='Show feature')
    >>> extra_widgets = {'show_feature': checkbox}
    >>> def my_plot(a, b, show_feature):
    ...     # plotting code here
    ...     pass
    >>> create_interactive_plot(my_plot, slider_configs, extra_widgets=extra_widgets)
    """

    # 1) create sliders
    sliders = {}
    for name, cfg in slider_configs.items():
        sliders[name] = FloatSlider(
            value=cfg["value"],
            min=cfg["min"],
            max=cfg["max"],
            step=cfg["step"],
            description=f"{name}:",
            continuous_update=False,
        )

    # 2) layout controls
    slider_list = list(sliders.values())
    rows = []
    for i in range(0, len(slider_list), n_cols):
        rows.append(HBox(slider_list[i:i+n_cols]))

    if extra_widgets:
        if extra_widgets_inline:
            rows.append(HBox(list(extra_widgets.values())))
        else:
            rows.extend(list(extra_widgets.values()))

    controls = VBox(rows)
    display(controls)

    # 3) output area
    output = Output()
    display(output)

    last_draw = 0.0
    THROTTLE_SEC = 0.25  # 250 ms: init-burst

    all_widgets = dict(sliders)
    if extra_widgets:
        all_widgets.update(extra_widgets)


    def draw():
        nonlocal last_draw
        now = time.perf_counter()
        if now - last_draw < THROTTLE_SEC:
            return
        last_draw = now

        with output:
            output.clear_output(wait=True)
            plt.close("all")  # extra defensief in VS Code
            fig = plot_func(**{k: w.value for k, w in all_widgets.items()})
            if fig is not None:
                display(fig)
                plt.close(fig)

    def on_change(change):
        draw()

    for w in all_widgets.values():
        w.observe(on_change, names="value")

    draw() 





# =============================================================================
# Legend Helper Functions
# =============================================================================

def legend_dedup(ax, *, loc="best", **kwargs):
    """
    Create a legend on ax, deduplicating entries by (label, handle_type).
    
    Removes "_nolegend_" entries and empty labels.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object on which to place the legend.
    loc : str, optional
        Legend location (passed to ax.legend).
    **kwargs : dict
        Additional keyword arguments passed to ax.legend.
    """
    handles, labels = ax.get_legend_handles_labels()
    uniq = OrderedDict()
    for h, lab in zip(handles, labels):
        if not lab or lab == "_nolegend_":
            continue
        key = (lab, type(h).__name__)
        uniq.setdefault(key, (h, lab))
    if uniq:
        legend = ax.legend([h for h, _ in uniq.values()],
                           [lab for _, lab in uniq.values()],
                           loc=loc, **kwargs)
        # Ensure a thin border when frameon=True
        if legend and legend.get_frame():
            legend.get_frame().set_linewidth(0.5)


def fig_legend_dedup(fig, ax_list, *, loc="best", **kwargs):
    """
    Create a shared legend on figure, collecting and deduplicating from multiple axes.
    
    Removes "_nolegend_" entries and empty labels.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object on which to place the shared legend.
    ax_list : list or Axes
        Axes object or list of Axes objects to collect legend entries from.
    loc : str, optional
        Legend location (passed to fig.legend).
    **kwargs : dict
        Additional keyword arguments passed to fig.legend.
    
    Returns
    -------
    legend : matplotlib.legend.Legend
        The created legend object.
    """
    if not isinstance(ax_list, list):
        ax_list = [ax_list]
    
    uniq = OrderedDict()
    for ax in ax_list:
        handles, labels = ax.get_legend_handles_labels()
        for h, lab in zip(handles, labels):
            if not lab or lab == "_nolegend_":
                continue
            key = (lab, type(h).__name__)
            uniq.setdefault(key, (h, lab))
    
    if uniq:
        legend = fig.legend([h for h, _ in uniq.values()],
                            [lab for _, lab in uniq.values()],
                            loc=loc, **kwargs)
        # Ensure a thin border when frameon=True
        if legend and legend.get_frame():
            legend.get_frame().set_linewidth(0.5)
        return legend
    return None
