import matplotlib.pyplot as plt
import pathlib


class DataSeries:
    def __init__(self, data_points, label="", line_style="-", marker=None):
        self.data_points = data_points
        self.label = label
        self.line_style = line_style
        self.marker = marker


def _make_plot(
    data_series_collection,
    plot_file_path,
    xlabel=None,
    ylabel=None,
    yaxis_inverted=False,
    xscale=None,
):
    figure, axes = plt.subplots(layout="constrained")
    if xscale is not None:
        axes.set_xscale(xscale)

    for series in data_series_collection:
        axes.plot(
            *zip(*series.data_points),
            label=series.label,
            linestyle=series.line_style,
            marker=series.marker,
        )
    axes.grid()
    axes.grid(which="minor", color="0.9")
    axes.yaxis.set_inverted(yaxis_inverted)
    if xlabel is not None:
        axes.set_xlabel(xlabel)
    if ylabel is not None:
        axes.set_ylabel(ylabel)
    figure.legend(loc="outside center right")

    if isinstance(plot_file_path, pathlib.Path):
        plot_file_path = str(plot_file_path.resolve())
    print(f"Saving plot to {plot_file_path}")
    plt.savefig(plot_file_path)


def make_settlement_history_plot(data_series_collection, plot_file_path):
    _make_plot(
        data_series_collection,
        plot_file_path,
        xlabel="Time [day]",
        ylabel="Settlement [m]",
        yaxis_inverted=True,
        xscale="log",
    )


def make_stress_over_y_plot(data_series_collection, plot_file_path):
    _make_plot(
        data_series_collection, plot_file_path, xlabel="Stress [kPa]", ylabel=r"$y$ [m]"
    )
