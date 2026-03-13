# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication
import KratosMultiphysics.kratos_utilities

VTK_QUAD = 9

try:
    import h5py as h5
except ImportError as e:
    raise ImportError(
        "h5py is required for this functionality but it is not installed. "
        "Please install it with: pip install h5py"
    ) from e

def Factory(settings, model):
    if not (isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return IgaVTKOutputProcess(model, settings["Parameters"])

class IgaVTKOutputProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""{
            "output_file_name"          : "",
            "brep_surface_ids"          : [],                              
            "model_part_name"           : "",
            "file_label"                : "step",
            "output_control_type"       : "step",
            "output_frequency"          : 1.0
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        ## Get the model part
        self.model_part = model[self.params["model_part_name"].GetString()]
        print('---------------MODEL PART------------------')
        print(self.model_part)


        ## Get the geometry
        self.brep_surface_ids = [
            self.params["brep_surface_ids"][i].GetInt()
            for i in range(self.params["brep_surface_ids"].size())
        ]

        self.output_file_name = self.params["output_file_name"].GetString()
        if not self.output_file_name.endswith(".vtkhdf"):
            self.output_file_name += ".vtkhdf"
        
        # The geometries are not stored in the model part, but they are in the model
        self.brep_surface = self.model_part.GetGeometry(self.brep_surface_ids[0])
        print(self.brep_surface)

        # self.nodal_results_scalar, self.nodal_results_vector = \
        #     CreateVariablesListFromInput(self.params["nodal_results"])

        # self.integration_point_results_scalar, self.integration_point_results_vector = \
        #     CreateVariablesListFromInput(self.params["integration_point_results"])

        # Set up output frequency and format
        output_file_label = self.params["file_label"].GetString()
        if output_file_label == "time":
            self.output_label_is_time = True
        elif output_file_label == "step":
            self.output_label_is_time = False
        else:
            msg = '{} Error: Unknown value "{}" read for parameter "file_label"'.format(self.__class__.__name__,output_file_label)
            raise Exception(msg)

        output_control_type = self.params["output_control_type"].GetString()
        if output_control_type == "time":
            self.output_control_is_time = True
        elif output_control_type == "step":
            self.output_control_is_time = False
        else:
            err_msg  = 'The requested "output_control_type" "' + output_control_type
            err_msg += '" is not available!\nAvailable options are: "time", "step"'
            raise Exception(err_msg)

        self.output_frequency = self.params["output_frequency"].GetDouble()

        self.step_count = 0
        self.printed_step_count = 0
        self.next_output = 0.0

    def ExecuteBeforeSolutionLoop(self):
        with open(self.output_file_name, 'w') as output_file:
            output_file.write("Post Results File 1.0\n")

    def PrintOutput(self):
        time = GetPrettyTime(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        self.printed_step_count += 1
        self.model_part.ProcessInfo[KratosMultiphysics.PRINTED_STEP] = self.printed_step_count
        if self.output_label_is_time:
            label = time
        else:
            label = self.printed_step_count

        #TODO:
        file = h5.File(self.output_file_name, "w")
        root = file.create_group("VTKHDF")
        root.attrs["Version"] = (2, 4)
        root.attrs["Type"] = "PartotopmedDataSetCollection"

        # for brep_id in self.brep_surface_ids:
        #     geometry = self.model_part.GetGeometry(brep_id)
        #     knots_u = geometry.KnotsU()
        #     knots_v = geometry.KnotsV()




        # Schedule next output
        if self.output_frequency > 0.0: # Note: if == 0, we'll just always print
            if self.output_control_is_time:
                while GetPrettyTime(self.next_output) <= time:
                    self.next_output += self.output_frequency
            else:
                while self.next_output <= self.step_count:
                    self.next_output += self.output_frequency

    def IsOutputStep(self):
        if self.output_control_is_time:
            time = GetPrettyTime(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return (time >= GetPrettyTime(self.next_output))
        else:
            return ( self.step_count >= self.next_output )

    def __write_out_surface(self, grid_points) -> h5.Group:
        None



def GetPrettyTime(time):
    pretty_time = "{0:.12g}".format(time)
    pretty_time = float(pretty_time)
    return pretty_time

def CreateVariablesListFromInput(param):
    '''Parse a list of variables from input.'''
    scalar_variables = []
    vector_variables = []
    admissible_scalar_types = ["Bool", "Integer", "Unsigned Integer", "Double"]
    admissible_vector_types = ["Array", "Vector"]

    variable_list = KratosMultiphysics.kratos_utilities.GenerateVariableListFromInput(param)

    # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
    for variable in variable_list:
        if KratosMultiphysics.KratosGlobals.GetVariableType(variable.Name()) in admissible_scalar_types:
            scalar_variables.append(variable)
        elif KratosMultiphysics.KratosGlobals.GetVariableType(variable.Name()) in admissible_vector_types:
            vector_variables.append(variable)
        else:
            raise Exception("unsupported variables type: " + str(type(variable)))

    return scalar_variables, vector_variables
