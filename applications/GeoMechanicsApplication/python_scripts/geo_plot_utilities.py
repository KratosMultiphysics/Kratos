import matplotlib.pyplot as plt
import pathlib


def _make_plot(series_collection, figure_filename, xlabel=None, ylabel=None, yaxis_inverted=False, xscale=None):
    figure, axes = plt.subplots(layout="constrained")
    if xscale is not None:
        axes.set_xscale(xscale)

    for series in series_collection:
        axes.plot(series.x_values, series.y_values, label=series.label, linestyle=series.linestyle, marker=series.marker)
    axes.grid()
    axes.grid(which="minor", color="0.9")
    axes.yaxis.set_inverted(yaxis_inverted)
    if xlabel is not None:
        axes.set_xlabel(xlabel)
    if ylabel is not None:
        axes.set_ylabel(ylabel)
    figure.legend(loc='outside center right')

    if isinstance(figure_filename, pathlib.Path):
        figure_filename = str(figure_filename.resolve())
    print(f"Saving plot to {figure_filename}")
    plt.savefig(figure_filename)


def plot_settlement_results(series_collection, figure_filename):
    _make_plot(series_collection, figure_filename, xlabel='Time [day]', ylabel='Settlement [m]', yaxis_inverted=True, xscale='log')


def make_stress_plot(series_collection, figure_filename):
    _make_plot(series_collection, figure_filename, xlabel='Stress [kPa]', ylabel='y [m]')
