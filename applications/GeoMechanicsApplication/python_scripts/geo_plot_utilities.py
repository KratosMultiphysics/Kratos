import matplotlib.pyplot as plt
import pathlib


class DataSeries():
    def __init__(self, x_values, y_values, label='', line_style='-', marker=None):
        self.x_values = x_values
        self.y_values = y_values
        self.label = label
        self.line_style = line_style
        self.marker = marker


def make_data_series(data_points, label='', line_style='-', marker=None):
    x_values = [point[0] for point in data_points]
    y_values = [point[1] for point in data_points]
    return DataSeries(x_values, y_values, label=label, line_style=line_style, marker=marker)


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
            series.x_values,
            series.y_values,
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


def plot_settlement_results(data_series_collection, plot_file_path):
    _make_plot(
        data_series_collection,
        plot_file_path,
        xlabel="Time [day]",
        ylabel="Settlement [m]",
        yaxis_inverted=True,
        xscale="log",
    )


def make_stress_plot(data_series_collection, plot_file_path):
    _make_plot(
        data_series_collection, plot_file_path, xlabel="Stress [kPa]", ylabel="y [m]"
    )
