import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics.testing.utilities import ReadModelPart


def create_test_model() -> KM.Model:
    model = KM.Model()
    model_part = model.CreateModelPart("test")
    model_part.AddNodalSolutionStepVariable(KM.DENSITY)
    model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
    model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
    model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
    model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
    ReadModelPart("quads", model_part)

    for node in model_part.Nodes:
        id = node.Id
        node.SetSolutionStepValue(KM.VELOCITY, KM.Array3([id+3, id+4, id+5]))
        node.SetSolutionStepValue(KM.PRESSURE, id+3)
        node.SetValue(KM.PRESSURE, id+3)
        node.SetValue(KM.VELOCITY, KM.Array3([id+3, id+4, id+5]))

    KOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(
        model_part, model_part.Conditions)
    for condition in model_part.Conditions:
        id = condition.Id
        condition.Properties[KM.PRESSURE] = id+400
        condition.Properties[KM.VELOCITY] = KM.Array3([id+500, id+600, id+700])
        condition.SetValue(KM.PRESSURE, id+4)
        condition.SetValue(KM.VELOCITY, KM.Array3([id+5, id+6, id+7]))

    KOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(
        model_part, model_part.Elements)
    for element in model_part.Elements:
        id = element.Id
        element.Properties[KM.PRESSURE] = id+500
        element.Properties[KM.VELOCITY] = KM.Array3([id+600, id+700, id+800])
        element.SetValue(KM.PRESSURE, id+5)
        element.SetValue(KM.VELOCITY, KM.Array3([id+6, id+7, id+8]))

    return model


model = create_test_model()

print(model)
print(model.GetModelPart("test"))

test_model_part = model.GetModelPart("test")
container_non_hist = KOA.NodalNonHistoricalVariableData(test_model_part)
container_hist = KOA.HistoricalVariableData(test_model_part)

print(container_non_hist)
print(container_hist)

# proper_container = KOA.ContainerVariableDataUtils.MapNodalVariableToContainerVariable(
#     container_hist, test_model_part.Nodes)

vtk_output = KOA.ContainerVariableDataVtkOutput(
    test_model_part, KM.Parameters())
vtk_output.TestFunction()
# vtk_output.WriteContainerDataToFile()

# KOA.ContainerVariableDataVtkOutput.WriteContainerDataToFile()

vtk_output.printNumberType(int(20))
vtk_output.printNumberType(float(10.5))

generic_vtk_output = KOA.GenericVtkOutput()
pos = np.linspace(0, 1, 10)
cell = pos**2
point = pos**0.5

generic_vtk_output.outputStructuredGrid("./test_generic_output.vtk", pos, pos, pos, cell, point)
