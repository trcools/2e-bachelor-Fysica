"""
Personalized Layout Utilities

This module provides reusable helper functions for creating interactive plots
and managing legends across different exam parts of the Biophysics course.


Note:
This module was refactored with assistance from an AI tool to improve structure,
naming, and documentation for readability and maintainability.

"""

# =============================================================================
# Imports
# =============================================================================

from ipywidgets import interact, FloatSlider, HBox, interactive_output, Checkbox, Output, VBox
from IPython.display import display, update_display
from collections import OrderedDict
import matplotlib.pyplot as plt
import time

# =============================================================================
# Interactive Plot Helper Function
# =============================================================================

"""
Note:

I have a lot of problems with my interactive plots giving the figure back to sometimes 6 times, I was working on this problem for at least 2 days,
and haven't been able to solve it. Even when I asked AI to help me it couldn't give me an answer why I have this bug.
Maybe it is because I swithed to VSCode and something in my environment or settings is not right, but i haven't been able to find out what the problem was.

I tried some code what from AI recommanded to fix my problem, but this only made it messy and did not fix the code.
Now I am in some time issues, so I will just leave it like this for now.
"""

def create_interactive_plot(plot_func, slider_configs, n_cols=3, extra_widgets=None, extra_widgets_inline=False):
    # --- build widgets ---
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

    slider_list = list(sliders.values())
    rows = [HBox(slider_list[i:i+n_cols]) for i in range(0, len(slider_list), n_cols)]

    all_widgets = dict(sliders)
    if extra_widgets:
        if extra_widgets_inline:
            rows.append(HBox(list(extra_widgets.values())))
        else:
            rows.extend(list(extra_widgets.values()))
        all_widgets.update(extra_widgets)

    display(VBox(rows))

    out = Output()
    display(out)

    
    last_call = 0.0
    MIN_DT = 0.15  

    display_id = None

    def render():
        nonlocal last_call, display_id
        now = time.perf_counter()
        if now - last_call < MIN_DT:
            return
        last_call = now

        plt.close("all")

        kwargs = {k: w.value for k, w in all_widgets.items()}
        fig = plot_func(**kwargs)

        if fig is None:
            return

        with out:
        
            if display_id is None:
                h = display(fig, display_id=True)
                display_id = h.display_id
            else:
                update_display(fig, display_id=display_id)

        plt.close(fig)

    def on_change(change):
        # ignore pure sync no-ops
        if change.get("old") == change.get("new"):
            return
        render()

  
    for w in all_widgets.values():
        w.observe(on_change, names="value")

    render()





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
