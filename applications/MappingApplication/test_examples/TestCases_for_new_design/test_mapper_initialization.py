import KratosMultiphysics

try:
    import KratosMultiphysics.mpi as KratosMPI
except ImportError:
    pass
import KratosMultiphysics.MappingApplication as KratosMapping

model_part_origin = KratosMultiphysics.ModelPart()
model_part_destination = KratosMultiphysics.ModelPart()

model_part_origin.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
model_part_origin.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
model_part_destination.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
model_part_destination.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

model_part_origin.CreateNewNode(1,1.0,2.0,3.0)
model_part_destination.CreateNewNode(1,2.0,2.0,3.0)

mapper_settings = KratosMultiphysics.Parameters(
    """
    {
        "mapper_type" : "nearest_neighbor",
        "echo_level"  : 2
    }
    """
);

# try:
mapper = KratosMapping.MapperFactory.CreateMapper(model_part_origin, model_part_destination, mapper_settings)
# except AttributeError:
#     mapper = KratosMapping.MapperFactory.CreateMapper(model_part_origin, model_part_destination, mapper_settings)


print(mapper)

mapper.Map(KratosMultiphysics.PRESSURE, KratosMultiphysics.TEMPERATURE)

mapper.UpdateInterface(KratosMapping.Mapper.REMESHED)
mapper.UpdateInterface()