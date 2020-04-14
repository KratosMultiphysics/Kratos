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

# TODO in the future we should just submit the interface submodel parts to FindIntersection1DGeometries2D
KratosMapping.FindIntersection1DGeometries2D(model_part_origin, model_part_destination, model_part_coupling, 1e-6)
KratosMapping.CreateQuadraturePointsCoupling1DGeometries2D(model_part_coupling, model_part_coupling_quadrature_points, 1e-6)

mapper_params = KM.Parameters("""{
    "mapper_type": "coupling_geometry",
    "echo_level" : 0
}""")

mapper = KratosMapping.MapperFactory.CreateMapper(model_part_coupling_quadrature_points, model_part_destination, mapper_params)

for node in model_part_origin.Nodes:
    node.SetSolutionStepValue(KM.DISPLACEMENT_X, 1.0, model_part_origin.ProcessInfo)

mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)