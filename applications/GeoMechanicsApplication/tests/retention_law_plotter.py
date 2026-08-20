import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import Pa_to_kPa
import os
from pathlib import Path


def plot_van_genuchten_retention_law_characteristics(
    properties: Kratos.Properties, plot_file_path: Path
) -> None:

    law = KratosGeo.VanGenuchtenLaw()

    num_points = 250
    pressures = []
    for i in range(0, num_points):
        pressures.append(-1000 + i * 21000 / num_points)

    parameters = KratosGeo.RetentionLawParameters(properties)
    relative_permeabilities = []
    degrees_of_saturation = []
    diffusivities = []
    capacities = []
    for pressure in pressures:
        parameters.SetFluidPressure(pressure)

        relative_permeability = law.CalculateRelativePermeability(parameters)
        relative_permeabilities.append(relative_permeability)
        degrees_of_saturation.append(law.CalculateSaturation(parameters))
        capacity = -1.0 * law.CalculateDerivativeOfSaturation(parameters)
        capacities.append(capacity)
        capacity = max(capacity, 1e-6)
        diffusivities.append(relative_permeability / capacity)

    data_points = zip(
        [Pa_to_kPa(pressure) for pressure in pressures], relative_permeabilities
    )
    data_series_collections = []
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_points,
            line_style="-",
            marker="",
        )
    )
    data_series_collections.append(data_series_collection)

    # plot_utils._make_plot(
    #     data_series_collection=data_series_collection,
    #     plot_file_path=plot_file_path / "relative_permeability_vs_pressure_plot.svg",
    #     xlabel="water pressure [kPa]",
    #     ylabel="Relative Permeability [-]",
    #     logy=True,
    # )

    data_points = zip(
        [Pa_to_kPa(pressure) for pressure in pressures], degrees_of_saturation
    )
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_points,
            line_style="-",
            marker="",
        )
    )

    data_series_collections.append(data_series_collection)
    # plot_utils._make_plot(
    #     data_series_collection=data_series_collection,
    #     plot_file_path=plot_file_path / "degree_of_saturation_vs_pressure_plot.svg",
    #     xlabel="water pressure [kPa]",
    #     ylabel="Saturation [-]",
    # )

    data_points = zip([Pa_to_kPa(pressure) for pressure in pressures], diffusivities)
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_points,
            line_style="-",
            marker="",
        )
    )

    data_series_collections.append(data_series_collection)
    # plot_utils._make_plot(
    #     data_series_collection=data_series_collection,
    #     plot_file_path=plot_file_path / "diffusivity_vs_pressure_plot.svg",
    #     xlabel="water pressure [kPa]",
    #     ylabel="Diffusivity [-]",
    #     logy=True,
    # )

    data_points = zip([Pa_to_kPa(pressure) for pressure in pressures], capacities)
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_points,
            line_style="-",
            marker="",
        )
    )

    data_series_collections.append(data_series_collection)
    # plot_utils._make_plot(
    #     data_series_collection=data_series_collection,
    #     plot_file_path=plot_file_path / "capacity_vs_pressure_plot.svg",
    #     xlabel="water pressure [kPa]",
    #     ylabel="Capacity [-]",
    # )
    plot_utils.make_separate_sub_plots(
        data_series_collections=data_series_collections,
        plot_file_path=plot_file_path / "van_genuchten_characteristics.svg",
        subplot_options=[
            plot_utils.SubPlotOptions(
                title="Relative Permeability vs Pressure",
                xlabel="Water pressure [kPa]",
                ylabel="Relative Permeability [-]",
                log_y_plot=True,
            ),
            plot_utils.SubPlotOptions(
                title="Saturation vs Pressure",
                xlabel="Water pressure [kPa]",
                ylabel="Saturation [-]",
            ),
            plot_utils.SubPlotOptions(
                title="Relative Diffusivity vs Pressure",
                xlabel="Water pressure [kPa]",
                ylabel="Relative Diffusivity [-]",
                log_y_plot=True,
            ),
            plot_utils.SubPlotOptions(
                title="Capacity vs Pressure",
                xlabel="Water pressure [kPa]",
                ylabel="Capacity [-]",
            ),
        ],
        max_plots_per_row=2,
    )
