import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import Pa_to_kPa
import os
def plot_retention_law(properties : Kratos.Properties, plot_file_path: Path):
    parameters = KratosGeo.RetentionLawParameters(properties)
    parameters.SetFluidPressure(7.0)
    print(f"Fluid pressure = {parameters.GetFluidPressure()}")

    law = KratosGeo.VanGenuchtenLaw()

    num_points = 250
    pressures =[]
    for i in range(0, num_points):
        pressures.append(-1000 + i * 21000 / num_points)

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

    data_points = zip([Pa_to_kPa(pressure) for pressure in pressures], relative_permeabilities)
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_points,
            line_style="-",
            marker="",
            color='blue',
        )
    )

    plot_utils._make_plot(
        data_series_collection= data_series_collection,
        plot_file_path = os.path.join(plot_file_path, "relative_permeability_vs_pressure_plot.svg"),
        xlabel="water pressure [kPa]",
        ylabel="Relative Permeability [-]",
    )


    data_points = zip([Pa_to_kPa(pressure) for pressure in pressures], degrees_of_saturation)
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_points,
            line_style="-",
            marker="",
            color='blue',
        )
    )

    plot_utils._make_plot(
        data_series_collection= data_series_collection,
        plot_file_path = os.path.join(plot_file_path, "saturation_vs_pressure_plot.svg"),
        xlabel="water pressure [kPa]",
        ylabel="Saturation [-]",
    )

    data_points = zip([Pa_to_kPa(pressure) for pressure in pressures], diffusivities)
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_points,
            line_style="-",
            marker="",
            color='blue',
        )
    )

    plot_utils._make_plot(
        data_series_collection= data_series_collection,
        plot_file_path = os.path.join(plot_file_path, "diffusivity_vs_pressure_plot.svg"),
        xlabel="water pressure [kPa]",
        ylabel="Diffusivity [?]",
    )

    data_points = zip([Pa_to_kPa(pressure) for pressure in pressures], capacities)
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_points,
            line_style="-",
            marker="",
            color='blue',
        )
    )

    plot_utils._make_plot(
        data_series_collection= data_series_collection,
        plot_file_path = os.path.join(plot_file_path, "capacity_vs_pressure_plot.svg"),
        xlabel="water pressure [kPa]",
        ylabel="Capacity [-]",
    )

