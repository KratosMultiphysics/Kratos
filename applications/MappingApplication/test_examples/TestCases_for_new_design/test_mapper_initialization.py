import KratosMultiphysics

try:
    import KratosMultiphysics.mpi as KratosMPI
except ImportError:
    pass
import KratosMultiphysics.MappingApplication as KratosMapping

model_part_origin = KratosMultiphysics.ModelPart()
model_part_destination = KratosMultiphysics.ModelPart()

mapper_settings = KratosMultiphysics.Parameters(
    """
    {
        "mapper_type" : "nearest_element",
        "echo_level"  : 1
    }
    """
);

# try:
mapper = KratosMapping.MapperFactory.CreateMapper(model_part_origin, model_part_destination, mapper_settings)
# except AttributeError:
#     mapper = KratosMapping.MapperFactory.CreateMapper(model_part_origin, model_part_destination, mapper_settings)


print(mapper)