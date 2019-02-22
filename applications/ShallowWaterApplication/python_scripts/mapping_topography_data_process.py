import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow
import KratosMultiphysics.MappingApplication as Mapping

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MappingTopographyDataProcess(Model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class MappingTopographyDataProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "topographic_model_part_name" : "",
                "computing_model_part_name"   : "",
                "mapper_settings"             : {
                    "mapper_type"         : "nearest_neighbor"
                }
            }
            """
            )
        settings.RecursivelyValidateAndAssignDefaults(default_settings)

        # Data process
        self.topography_process = Mapping.MapperFactory.CreateMapper(
            Model[settings["topographic_model_part_name"].GetString()],
            Model[settings["computing_model_part_name"].GetString()],
            settings["mapper_settings"]
        )
        # Z-coordinate process
        self.coordinate_process = Shallow.ShallowWaterVariablesUtility(Model[settings["computing_model_part_name"].GetString()])

    def ExecuteInitialize(self):
        # Data process
        self.topography_process.Map(Shallow.BATHYMETRY, Shallow.BATHYMETRY)
        # Z-coordinate process: set the z coordinate to zero
        self.coordinate_process.ResetMeshPosition(True)
