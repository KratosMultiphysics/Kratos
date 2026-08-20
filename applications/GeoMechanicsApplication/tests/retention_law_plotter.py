from dataclasses import dataclass, field

import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import Pa_to_kPa
import os
from pathlib import Path
import csv


@dataclass
class Characteristics:
    pressures: list[float] = field(default_factory=list)
    relative_permeabilities: list[float] = field(default_factory=list)
    degrees_of_saturation: list[float] = field(default_factory=list)
    diffusivities: list[float] = field(default_factory=list)
    capacities: list[float] = field(default_factory=list)


def read_characteristics(comparison_data_path: Path) -> Characteristics:
    with open(
        comparison_data_path,
        newline="",
    ) as csv_file:
        reader = csv.DictReader(csv_file, skipinitialspace=True)
        result = Characteristics()
        for row in reader:
            result.pressures.append(-1.0 * float(row["pressure"]))
            result.relative_permeabilities.append(float(row["relative_permeability"]))
            result.degrees_of_saturation.append(float(row["saturation"]))
            result.diffusivities.append(float(row["diffusivity"]))
            result.capacities.append(float(row["capacity"]))
        return result


def plot_van_genuchten_retention_law_characteristics(
    properties: Kratos.Properties,
    plot_file_path: Path,
    comparison_data_path: Path | None = None,
) -> None:
    comparison_characteristics: Characteristics | None = None
    if comparison_data_path:
        comparison_characteristics = read_characteristics(comparison_data_path)

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
        degrees_of_saturation.append(law.CalculateEffectiveSaturation(parameters))
        capacity = -1.0 * law.CalculateDerivativeOfSaturation(parameters)
        capacities.append(capacity)
        diffusivities.append(relative_permeability / capacity)

    data_points = zip(
        [Pa_to_kPa(pressure) for pressure in pressures], relative_permeabilities
    )
    data_series_collections = []
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_points,
            "Kratos",
        )
    )
    if comparison_characteristics:
        data_series_collection.append(
            plot_utils.DataSeries(
                zip(
                    comparison_characteristics.pressures,
                    comparison_characteristics.relative_permeabilities,
                ),
                "DG-Flow",
            )
        )
    data_series_collections.append(data_series_collection)
    data_points = zip(
        [Pa_to_kPa(pressure) for pressure in pressures], degrees_of_saturation
    )
    data_series_collection = []
    data_series_collection.append(plot_utils.DataSeries(data_points, "Kratos"))

    if comparison_characteristics:
        data_series_collection.append(
            plot_utils.DataSeries(
                zip(
                    comparison_characteristics.pressures,
                    comparison_characteristics.degrees_of_saturation,
                ),
                "DG-Flow",
            )
        )

    data_series_collections.append(data_series_collection)

    data_points = zip([Pa_to_kPa(pressure) for pressure in pressures], diffusivities)
    data_series_collection = []
    data_series_collection.append(plot_utils.DataSeries(data_points, "Kratos"))

    if comparison_characteristics:
        data_series_collection.append(
            plot_utils.DataSeries(
                zip(
                    comparison_characteristics.pressures,
                    comparison_characteristics.diffusivities,
                ),
                "DG-Flow",
            )
        )
    data_series_collections.append(data_series_collection)

    data_points = zip([Pa_to_kPa(pressure) for pressure in pressures], capacities)
    data_series_collection = []
    data_series_collection.append(plot_utils.DataSeries(data_points, "Kratos"))
    if comparison_characteristics:
        data_series_collection.append(
            plot_utils.DataSeries(
                zip(
                    comparison_characteristics.pressures,
                    comparison_characteristics.capacities,
                ),
                "DG-Flow",
            )
        )

    data_series_collections.append(data_series_collection)
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
