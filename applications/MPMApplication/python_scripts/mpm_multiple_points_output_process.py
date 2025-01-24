# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.multiple_points_output_process import MultiplePointsOutputProcess
from KratosMultiphysics.MPMApplication.mpm_point_output_process import MPMPointOutputProcess

def Factory(settings, Model):
    if not isinstance(Model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MPMMultiplePointsOutputProcess(Model, settings["Parameters"])

class MPMMultiplePointsOutputProcess(MultiplePointsOutputProcess):
    """This process writes several points to a file
    Internally it holds objects of type "MPMPointOutputProcess"
    """
    def __init__(self, model, params):
        KratosMultiphysics.OutputProcess.__init__(self)

        params.ValidateAndAssignDefaults(self.GetDefaultParameters())

        positions = params["positions"].GetMatrix()

        num_points = positions.Size1()
        if num_points == 0:
            raise Exception('No positions were specified"')
        if positions.Size2() != 3:
            raise Exception('The positions have to be provided with 3 coordinates!')

        self.point_output_processes = []
        params.RemoveValue("positions")
        params.AddEmptyValue("position")
        position_vec = KratosMultiphysics.Vector(3)

        if params["output_file_settings"]["file_name"].GetString().endswith(".dat"):
            params["output_file_settings"]["file_name"].SetString(params["output_file_settings"]["file_name"].GetString()[:-4])

        # Create the individual point_output_processes
        for i in range(num_points):
            point_proc_params = params.Clone()

            for j in range(3): position_vec[j] = positions[i,j]

            point_proc_params["position"].SetVector(position_vec)
            point_proc_params["output_file_settings"]["file_name"].SetString(params["output_file_settings"]["file_name"].GetString() + "_" + str(i+1))

            self.point_output_processes.append(MPMPointOutputProcess(model, point_proc_params))

    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters('''{
            "model_part_name"   : "",
            "entity_type"       : "element",
            "interval"          : [0.0, 1e30],
            "positions"         : [[]],
            "output_variables"  : [],
            "search_tolerance"  : 1e-6,
            "print_format"      : "",
            "output_file_settings": {}
        }''')
