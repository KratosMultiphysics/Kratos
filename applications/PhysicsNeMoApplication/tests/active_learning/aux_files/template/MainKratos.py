import KratosMultiphysics as Kratos

# The dummy_analysis fixture module is made importable through PYTHONPATH by
# the test driving this template.
import dummy_analysis

if __name__ == "__main__":
    with open("ProjectParameters.json", "r") as f:
        parameters = Kratos.Parameters(f.read())
    model = Kratos.Model()
    dummy_analysis.Create(model, parameters).Run()
