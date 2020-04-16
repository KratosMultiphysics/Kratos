import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

def ReadModelPart(model_part, mdpa_file_name):
    # adding variables used for mapping
    historical_vars = [KM.PRESSURE, KM.FORCE, KM.TEMPERATURE, KM.VELOCITY]
    for var in historical_vars:
        model_part.AddNodalSolutionStepVariable(var)

    import_flags = KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER
    KM.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)


current_model = KM.Model()
model_part_origin = current_model.CreateModelPart("origin")
model_part_destination = current_model.CreateModelPart("destination")

model_part_coupling = current_model.CreateModelPart("coupling")
model_part_coupling_quadrature_points = current_model.CreateModelPart("coupling_quadrature_points")

ReadModelPart(model_part_origin, "coupled_cantilever/domainA")
ReadModelPart(model_part_destination, "coupled_cantilever/domainB")

# Create submodelparts
model_part_origin.CreateSubModelPart("interface")
originInterface = model_part_origin.GetSubModelPart("interface")
for nodeIndex in range(328,334):
    originInterface.AddNodes([nodeIndex])
model_part_destination.CreateSubModelPart("interface")
destinationInterface = model_part_destination.GetSubModelPart("interface")
destinationInterface.AddNodes([10])
destinationInterface.AddNodes([6])
destinationInterface.AddNodes([3])
destinationInterface.AddNodes([1])


# Create dummy line conditions on submodel parts to couple together
originProps = originInterface.GetProperties()[1]
nodeOffset = 327
for conditionIndex in range(1,6):
    cNode = nodeOffset + conditionIndex
    originInterface.CreateNewCondition("LineLoadCondition2D2N", conditionIndex+100, [cNode,cNode+1], originProps)

destProps = destinationInterface.GetProperties()[1]
destinationInterface.CreateNewCondition("LineLoadCondition2D2N", 101, [1,3], destProps)
destinationInterface.CreateNewCondition("LineLoadCondition2D2N", 102, [3,6], destProps)
destinationInterface.CreateNewCondition("LineLoadCondition2D2N", 103, [6,10], destProps)


# TODO, add some logic to call the correct intersection functions. that would depend on the dimension of the coupling interfaces
KratosMapping.FindIntersection1DGeometries2D(originInterface, destinationInterface, model_part_coupling, 1e-6)
KratosMapping.CreateQuadraturePointsCoupling1DGeometries2D(model_part_coupling, model_part_coupling_quadrature_points, 1e-6)

# TODO remove (for information only)
print(model_part_coupling)
print(model_part_coupling_quadrature_points)

mapper_params = KM.Parameters("""{
    "mapper_type": "coupling_geometry",
    "echo_level" : 0
}""")

mapper = KratosMapping.MapperFactory.CreateMapper(model_part_coupling_quadrature_points, model_part_destination, mapper_params)

for node in model_part_origin.Nodes:
    node.SetSolutionStepValue(KM.DISPLACEMENT_X, 1.0, model_part_origin.ProcessInfo)

mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)