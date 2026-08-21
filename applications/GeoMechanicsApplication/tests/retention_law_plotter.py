import csv
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import \
    Pa_to_kPa


@dataclass
class RetentionCharacteristics:
    pressures: list[float] = field(default_factory=list)
    relative_permeabilities: list[float] = field(default_factory=list)
    degrees_of_saturation: list[float] = field(default_factory=list)
    diffusivities: list[float] = field(default_factory=list)
    capacities: list[float] = field(default_factory=list)

    @property
    def relative_permeability_vs_pressure_data(self) -> list[tuple[float, float]]:
        return list(zip(self.pressures, self.relative_permeabilities))

    @property
    def degrees_of_saturation_vs_pressure_data(self) -> list[tuple[float, float]]:
        return list(zip(self.pressures, self.degrees_of_saturation))

    @property
    def diffusivities_vs_pressure_data(self) -> list[tuple[float, float]]:
        return list(zip(self.pressures, self.diffusivities))

    @property
    def capacities_vs_pressure_data(self) -> list[tuple[float, float]]:
        return list(zip(self.pressures, self.capacities))


def _read_reference_characteristics(
    comparison_data_path: Path,
) -> RetentionCharacteristics:
    with open(
        comparison_data_path,
        newline="",
    ) as csv_file:
        reader = csv.DictReader(csv_file, skipinitialspace=True)
        result = RetentionCharacteristics()
        for row in reader:
            result.pressures.append(-1.0 * float(row["pressure"]))
            result.relative_permeabilities.append(float(row["relative_permeability"]))
            result.degrees_of_saturation.append(float(row["saturation"]))
            result.diffusivities.append(float(row["diffusivity"]))
            result.capacities.append(float(row["capacity"]))
        return result


def _obtain_kratos_characteristics(
    retention_law, parameters
) -> RetentionCharacteristics:
    result = RetentionCharacteristics()
    num_points = 250
    pressures = []
    for i in range(0, num_points):
        pressures.append(-1000 + i * 21000 / num_points)

    relative_permeabilities = []
    degrees_of_saturation = []
    diffusivities = []
    capacities = []
    for pressure in pressures:
        parameters.SetFluidPressure(pressure)

        relative_permeability = retention_law.CalculateRelativePermeability(parameters)
        result.relative_permeabilities.append(relative_permeability)
        result.degrees_of_saturation.append(
            retention_law.CalculateEffectiveSaturation(parameters)
        )
        capacity = -1.0 * retention_law.CalculateDerivativeOfSaturation(parameters)
        result.capacities.append(capacity)
        if capacity > 0:
            result.diffusivities.append(relative_permeability / capacity)

    result.pressures = [Pa_to_kPa(pressure) for pressure in pressures]
    return result


def _create_characteristic_data_series(
    kratos_characteristics: RetentionCharacteristics,
    comparison_characteristics: RetentionCharacteristics | None,
    data_getter: Callable[[RetentionCharacteristics], list[tuple[float, float]]],
):
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_getter(kratos_characteristics),
            "Kratos",
        )
    )
    if comparison_characteristics:
        data_series_collection.append(
            plot_utils.DataSeries(
                data_getter(comparison_characteristics),
                "DG-Flow",
            )
        )

    return data_series_collection


def plot_van_genuchten_retention_law_characteristics(
    properties: Kratos.Properties,
    plot_file_path: Path,
    comparison_data_path: Path | None = None,
) -> None:
    comparison_characteristics: RetentionCharacteristics | None = None
    if comparison_data_path:
        comparison_characteristics = _read_reference_characteristics(
            comparison_data_path
        )

    kratos_characteristics = _obtain_kratos_characteristics(
        KratosGeo.VanGenuchtenLaw(), KratosGeo.RetentionLawParameters(properties)
    )
    data_series_collections = []
    data_series_collections.append(
        _create_characteristic_data_series(
            kratos_characteristics=kratos_characteristics,
            comparison_characteristics=comparison_characteristics,
            data_getter=lambda p: p.relative_permeability_vs_pressure_data,
        )
    )
    data_series_collections.append(
        _create_characteristic_data_series(
            kratos_characteristics=kratos_characteristics,
            comparison_characteristics=comparison_characteristics,
            data_getter=lambda p: p.degrees_of_saturation_vs_pressure_data,
        )
    )
    data_series_collections.append(
        _create_characteristic_data_series(
            kratos_characteristics=kratos_characteristics,
            comparison_characteristics=comparison_characteristics,
            data_getter=lambda p: p.diffusivities_vs_pressure_data,
        )
    )
    data_series_collections.append(
        _create_characteristic_data_series(
            kratos_characteristics=kratos_characteristics,
            comparison_characteristics=comparison_characteristics,
            data_getter=lambda p: p.capacities_vs_pressure_data,
        )
    )
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
