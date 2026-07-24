import matplotlib.pyplot as plt
import pathlib
import numpy as np


class DataSeries:
    def __init__(self, data_points, label="", line_style="-", marker=None, color=None):
        self.data_points = data_points
        self.label = label
        self.line_style = line_style
        self.marker = marker
        self.color = color


def _make_plot(
    data_series_collection,
    plot_file_path,
    xlabel=None,
    ylabel=None,
    yaxis_inverted=False,
    xscale=None,
    title=None,
):
    figure, axes = plt.subplots(layout="constrained")
    if xscale is not None:
        axes.set_xscale(xscale)

    _plot_data_series_on_axis(axes, data_series_collection)
    axes.yaxis.set_inverted(yaxis_inverted)
    if xlabel is not None:
        axes.set_xlabel(xlabel)
    if ylabel is not None:
        axes.set_ylabel(ylabel)
    figure.legend(loc="outside center right")

    if title is not None:
        plt.title(title)

    if isinstance(plot_file_path, pathlib.Path):
        plot_file_path = str(plot_file_path.resolve())
    print(f"Saving plot to {plot_file_path}")
    plt.savefig(plot_file_path)
    plt.close(figure)


def _plot_data_series_on_axis(axes, data_series_collection):
    result = []
    for series in data_series_collection:
        # Unpack the data from pairs into two lists. See
        # https://stackoverflow.com/questions/21519203/plotting-a-list-of-x-y-coordinates for details.
        result.extend(axes.plot(
            *zip(*series.data_points),
            label=series.label,
            linestyle=series.line_style,
            marker=series.marker,
            color=series.color
        ))
    axes.grid()
    axes.grid(which="minor", color="0.9")

    return result


def make_sub_plots(
    data_series_collections,
    plot_file_path,
    titles,
    xlabel=None,
    ylabel=None,
    yaxis_inverted=False,
    xscale=None,
):
    figure, axes = plt.subplots(1, len(data_series_collections), figsize=(20, 6))
    axes = np.atleast_1d(axes)

    if ylabel is not None:
        axes[0].set_ylabel(ylabel)

    lines = []
    for ax, collection, title in zip(axes, data_series_collections, titles):
        lines.extend(_plot_data_series_on_axis(ax, collection))
        ax.yaxis.set_inverted(yaxis_inverted)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if xscale is not None:
            ax.set_xscale(xscale)
        ax.set_title(title)

    lines_with_unique_labels = []
    unique_labels = set()
    for line in lines:
        if line.get_label() not in unique_labels:
            lines_with_unique_labels.append(line)
            unique_labels.add(line.get_label())

    figure.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), handles=lines_with_unique_labels)

    if isinstance(plot_file_path, pathlib.Path):
        plot_file_path = str(plot_file_path.resolve())
    print(f"Saving plot to {plot_file_path}")
    plt.savefig(plot_file_path, bbox_inches="tight")
    plt.close(figure)


def make_settlement_history_plot(data_series_collection, plot_file_path):
    _make_plot(
        data_series_collection,
        plot_file_path,
        xlabel="Time [day]",
        ylabel="Settlement [m]",
        yaxis_inverted=True,
        xscale="log",
    )


def make_stress_over_y_plot(data_series_collection, plot_file_path, title=None):
    _make_plot(
        data_series_collection,
        plot_file_path,
        xlabel="Stress [kPa]",
        ylabel=r"$y$ [m]",
        title=title,
    )
