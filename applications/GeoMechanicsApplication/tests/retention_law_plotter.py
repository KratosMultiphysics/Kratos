import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import Pa_to_kPa

def plot_retention_law(properties : Kratos.Properties, plot_file_path: Path):
    parameters = KratosGeo.RetentionLawParameters(properties)
    parameters.SetFluidPressure(7.0)
    print(f"Fluid pressure = {parameters.GetFluidPressure()}")

    law = KratosGeo.VanGenuchtenLaw()

    degree_of_saturation = law.CalculateSaturation(parameters)
    num_points = 250
    pressures =[]
    for i in range(0, num_points):
        pressures.append(-1000 + i * 21000 / num_points)

    degrees_of_saturation = []
    for pressure in pressures:
        parameters.SetFluidPressure(pressure)
        degrees_of_saturation.append(law.CalculateSaturation(parameters))

    print("Pressures:", pressures)
    print("Degrees of Saturation:", degrees_of_saturation)

    data_points = zip([Pa_to_kPa(pressure) for pressure in pressures], degrees_of_saturation)
    data_series_collection = []
    data_series_collection.append(
        plot_utils.DataSeries(
            data_points,
            label=f"My very unique label",
            line_style="-",
            marker="",
            color='blue',
        )
    )

    plot_utils._make_plot(
        data_series_collection= data_series_collection,
        plot_file_path = plot_file_path,
        xlabel="water pressure [kPa]",
        ylabel="Saturation [-]",
    )

    print(f"Degree of Saturation: {degree_of_saturation}")

