import KratosMultiphysics
from KratosMultiphysics.check_scalar_base_process import CheckScalarBaseProcess

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckScalarToNodesProcess(Model, settings["Parameters"])


class CheckScalarToNodesProcess(CheckScalarBaseProcess, KratosUnittest.TestCase):
    """This process checks analytically from a function the solution (scalar) in a set of nodes belonging a certain submodelpart

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"            : "This process checks analytically from a function the solution (scalar) in a set of nodes belonging a certain submodelpart",
            "mesh_id"         : 0,
            "model_part_name" : "please_specify_model_part_name",
            "variable_name"   : "SPECIFY_VARIABLE_NAME",
            "interval"        : [0.0, 1e30],
            "value"           : 0.0,
            "tolerance_rank"  : 3,
            "reference_conf"  : false,
            "local_axes"      : {}
        }
        """
        )

        # We transfer the parameters to the base class
        base_settings = KratosMultiphysics.Parameters("""{}""")
        if settings.Has("model_part_name"):
            base_settings.AddValue("model_part_name", settings["model_part_name"])
        if settings.Has("variable_name"):
            base_settings.AddValue("variable_name", settings["variable_name"])
        if settings.Has("interval"):
            base_settings.AddValue("interval", settings["interval"])
        if settings.Has("value"):
            base_settings.AddValue("value", settings["value"])
        if settings.Has("tolerance_rank"):
            base_settings.AddValue("tolerance_rank", settings["tolerance_rank"])

        # Construct the base process
        super(CheckScalarToNodesProcess, self).__init__(Model, base_settings)

        # Copy from input
        self.settings = settings

        # Here i do a trick, since i want to allow "value" to be a string or a double value
        if(self.settings.Has("value")):
            if(self.settings["value"].IsString()):
                default_settings["value"].SetString("0.0")

        # Validate and assign
        self.settings.ValidateAndAssignDefaults(default_settings)

        # Admissible values for local axes, are "empty" or
        #"local_axes"               :{
        #    "origin" : [0.0, 0.0, 0.0]
        #    "axes"  : [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
        #    }
        self.non_trivial_local_system = False
        if(self.settings["local_axes"].Has("origin")): # not empty
            self.non_trivial_local_system = True
            self.R = KratosMultiphysics.Matrix(3,3)
            self.R[0,0] = self.settings["local_axes"]["axes"][0][0].GetDouble()
            self.R[0,1] = self.settings["local_axes"]["axes"][0][1].GetDouble()
            self.R[0,2] = self.settings["local_axes"]["axes"][0][2].GetDouble()
            self.R[1,0] = self.settings["local_axes"]["axes"][1][0].GetDouble()
            self.R[1,1] = self.settings["local_axes"]["axes"][1][1].GetDouble()
            self.R[1,2] = self.settings["local_axes"]["axes"][1][2].GetDouble()
            self.R[2,0] = self.settings["local_axes"]["axes"][2][0].GetDouble()
            self.R[2,1] = self.settings["local_axes"]["axes"][2][1].GetDouble()
            self.R[2,2] = self.settings["local_axes"]["axes"][2][2].GetDouble()

            self.x0 = KratosMultiphysics.Vector(3)
            self.x0[0] = self.settings["local_axes"]["origin"][0].GetDouble()
            self.x0[1] = self.settings["local_axes"]["origin"][1].GetDouble()
            self.x0[2] = self.settings["local_axes"]["origin"][2].GetDouble()

        # Getting the mesh
        self.mesh = self.model_part.GetMesh(self.settings["mesh_id"].GetInt())

        # Getting reference configuration
        self.reference_configuration = self.settings["reference_conf"].GetBool()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        process_info = self.model_part.ProcessInfo
        current_time = process_info[KratosMultiphysics.TIME]

        if(current_time >= self.interval[0] and  current_time<self.interval[1]):

            if self.value_is_numeric:
                for node in self.mesh.Nodes:
                    value = node.GetSolutionStepValue(self.variable, 0)
                    self.assertAlmostEqual(self.value, value, self.tolerance_rank)
            else:
                if self.is_time_function:
                    self.value = self.aux_function.f(0.0,0.0,0.0,current_time)
                    for node in self.mesh.Nodes:
                        value = node.GetSolutionStepValue(self.variable, 0)
                        self.assertAlmostEqual(self.value, value, self.tolerance_rank)
                else: #most general case - space varying function (possibly also time varying)
                    if self.non_trivial_local_system is False:
                        for node in self.model_part.Nodes:
                            value = node.GetSolutionStepValue(self.variable, 0)
                            if (self.reference_configuration is False):
                                self.assertAlmostEqual(self.aux_function.f(node.X,node.Y,node.Z,current_time), value, self.tolerance_rank)
                            else:
                                self.assertAlmostEqual(self.aux_function.f(node.X0,node.Y0,node.Z0,current_time), value, self.tolerance_rank)
                    else: #TODO: OPTIMIZE!!
                        for node in self.model_part.Nodes:
                            dx = node - self.x0
                            x = self.R[0,0]*dx[0] + self.R[0,1]*dx[1] + self.R[0,2]*dx[2]
                            y = self.R[1,0]*dx[0] + self.R[1,1]*dx[1] + self.R[1,2]*dx[2]
                            z = self.R[2,0]*dx[0] + self.R[2,1]*dx[1] + self.R[2,2]*dx[2]
                            value = node.GetSolutionStepValue(self.variable, 0)
                            self.assertAlmostEqual(self.aux_function.f(x,y,z,current_time), value, self.tolerance_rank)
