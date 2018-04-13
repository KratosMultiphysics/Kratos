import KratosMultiphysics

# import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MappingApplication as KratosMapping

model_part_origin = KratosMultiphysics.ModelPart()
model_part_destination = KratosMultiphysics.ModelPart()

mapper_settings = KratosMultiphysics.Parameters(
    """
    {
        "mapper_type" : "nearest_neighbor",
        "echo_level"  : 1
    }
    """
);

mapper = KratosMapping.MapperFactory.CreateMapper(model_part_origin, model_part_destination, mapper_settings)


print(mapper)