import numpy as np
import plotly.graph_objects as go


def render_plots(figs, axes, canvas):
    for i, ax in enumerate(axes):
        ax.clear()

        if i < len(figs):
            fig = figs[i]
            x_label = fig.layout.xaxis.title.text if fig.layout.xaxis.title else None
            y_label = fig.layout.yaxis.title.text if fig.layout.yaxis.title else None

            for trace in fig.data:
                if isinstance(trace, go.Scatter):
                    x = trace.x
                    y = trace.y
                    label = getattr(trace, 'name', None)
                    ax.plot(x, y, label=label)

            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title(_get_title(i))
            ax.legend()

            ax.invert_xaxis()
            if i in [1, 2, 4]:
                ax.invert_yaxis()
                ax.set_ylim(0, np.min(y) * 1.2)
            else:
                ax.set_ylim(0, np.max(y) * 1.2)

            ax.set_xlim(0, np.min(x) * 1.2)

    canvas.draw()
    canvas.get_tk_widget().config(width=2, height=7)


def _get_title(index):
    titles = [
        "Delta Sigma",
        "Volumetric Strain",
        "Sigma Plot",
        "p-q Plot",
        "Mohr-Coulomb Circle"
    ]
    return titles[index] if index < len(titles) else f"Plot {index + 1}"
