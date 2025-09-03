import matplotlib.pyplot as plt
import pathlib


def plot_settlement_results(series_collection, figure_filename):
    figure, axes = plt.subplots(layout="constrained")
    axes.set_xscale('log')
    for series in series_collection:
        axes.plot(series.x_values, series.y_values, label=series.label, linestyle=series.linestyle, marker=series.marker)
    axes.grid()
    axes.grid(which="minor", color="0.9")
    axes.yaxis.set_inverted(True)
    axes.set_xlabel('Time [day]')
    axes.set_ylabel('Settlement [m]')
    figure.legend(loc='outside center right')

    if isinstance(figure_filename, pathlib.Path):
        figure_filename = str(figure_filename.resolve())
    print(f"Save plot to {figure_filename}")
    plt.savefig(figure_filename)


def make_stress_plot(series_collection, figure_filename):
    figure, axes = plt.subplots(layout="constrained")
    for series in series_collection:
        axes.plot(series.x_values, series.y_values, label=series.label, linestyle=series.linestyle, marker=series.marker)
    axes.grid()
    axes.grid(which="minor", color="0.9")
    axes.set_xlabel('Stress [kPa]')
    axes.set_ylabel('y [m]')
    figure.legend(loc='outside center right')

    if isinstance(figure_filename, pathlib.Path):
        figure_filename = str(figure_filename.resolve())
    print(f"Save plot to {figure_filename}")
    plt.savefig(figure_filename)