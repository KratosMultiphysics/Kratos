import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics as Kratos
def plot_retention_law():
    properties = Kratos.Properties(0)
    properties[KratosGeo.RETENTION_LAW] = "VanGenuchtenLaw"

    parameters = KratosGeo.RetentionLawParameters(properties)
    parameters.SetFluidPressure(7.0)
    print(f"Fluid pressure = {parameters.GetFluidPressure()}")

    law = KratosGeo.VanGenuchtenLaw()

    degree_of_saturation = law.CalculateSaturation(parameters)
    print(f"Degree of Saturation: {degree_of_saturation}")




if __name__ == "__main__":
    plot_retention_law()

