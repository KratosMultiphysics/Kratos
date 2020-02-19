from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ProjectVectorOnSurfaceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ProjectVectorOnSurfaceProcess(KM.Process):

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KM.Process.__init__(self)

        default_settings = KM.Parameters("""{
            "model_part_name"  : "Structure",
            "echo_level"       : 0,
            "projection_type"  : "planar",
            "global_direction" : [1,0,0],
            "variable_name"    : "PLEASE_SPECIFY",
            "visualize_in_vtk" : false,
            "method_specific_settings" : { },
            "check_local_space_dimension" : true
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        model_part = Model[settings["model_part_name"].GetString()]
        StructuralMechanicsApplication.ProjectVectorOnSurfaceUtility.Execute(model_part, settings)

        if settings["visualize_in_vtk"].GetBool():
            vtk_parameters = KM.Parameters("""{
                "file_format"                  : "binary",
                "output_precision"             : 7,
                "output_control_type"          : "step",
                "output_sub_model_parts"       : false,
                "save_output_files_in_folder"  : false,
                "custom_name_prefix"           : "ProjectVectorOnSurface_",
                "custom_name_postfix"          : "",
                "element_data_value_variables" : []
            }""")

            variable_name = settings["variable_name"].GetString()
            vtk_parameters["custom_name_postfix"].SetString("_"+variable_name)
            vtk_parameters["element_data_value_variables"].Append(variable_name)

            vtk_io = KM.VtkOutput(model_part, vtk_parameters)
            vtk_io.PrintOutput()
