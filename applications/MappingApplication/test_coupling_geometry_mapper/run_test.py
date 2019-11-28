import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping

def ReadModelPart(model_part, mdpa_file_name):
    # adding varibables used for mapping
    historical_vars = [KM.PRESSURE, KM.FORCE, KM.TEMPERATURE, KM.VELOCITY]
    for var in historical_vars:
        model_part.AddNodalSolutionStepVariable(var)

    import_flags = KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER
    KM.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)


current_model = KM.Model()
model_part_origin = current_model.CreateModelPart("origin")
model_part_destination = current_model.CreateModelPart("destination")

ReadModelPart(model_part_origin, "../tests/mdpa_files/blade_quad")
ReadModelPart(model_part_destination, "../tests/mdpa_files/blade_tri")

mapper_params = KM.Parameters("""{
    "mapper_type": "coupling_geometry",
    "echo_level" : 0
}""")

mapper = KratosMapping.MapperFactory.CreateMapper(model_part_origin, model_part_destination, mapper_params)

# mapper.Map(KM.PRESSURE, KM.TEMPERATURE)