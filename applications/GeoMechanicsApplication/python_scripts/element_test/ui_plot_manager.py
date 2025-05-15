import numpy as np
import plotly.graph_objects as go


def render_plots(figs, axes, canvas):
    for i, ax in enumerate(axes):
        ax.clear()

        if i >= len(figs):
            continue

        fig = figs[i]
        x_label = _get_axis_label(fig.layout.xaxis)
        y_label = _get_axis_label(fig.layout.yaxis)

        x, y = _plot_traces(fig, ax)

        _configure_axis_labels(ax, x_label, y_label)
        _configure_plot_layout(ax, i, x, y)

    canvas.draw()
    canvas.get_tk_widget().config(width=2, height=7)


def _get_axis_label(axis_layout):
    return axis_layout.title.text if axis_layout.title else None


def _plot_traces(fig, ax):
    x, y = [], []
    for trace in fig.data:
        if isinstance(trace, go.Scatter):
            x = trace.x
            y = trace.y
            label = getattr(trace, 'name', None)
            ax.plot(x, y, label=label)
    ax.legend()
    return x, y


def _configure_axis_labels(ax, x_label, y_label):
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)


def _configure_plot_layout(ax, index, x, y):
    ax.set_title(_get_title(index))
    ax.invert_xaxis()

    if index in [1, 2, 4]:
        ax.invert_yaxis()
        ax.set_ylim(0, np.min(y) * 1.2)
    else:
        ax.set_ylim(0, np.max(y) * 1.2)

    ax.set_xlim(0, np.min(x) * 1.2)

def _get_title(index):
    titles = [
        "Delta Sigma",
        "Volumetric Strain",
        "Sigma Plot",
        "p-q Plot",
        "Mohr-Coulomb Circle"
    ]
    return titles[index] if index < len(titles) else f"Plot {index + 1}"
