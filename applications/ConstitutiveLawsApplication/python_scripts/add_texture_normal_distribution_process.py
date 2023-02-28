import KratosMultiphysics as KM
import numpy as np

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return AddTextureNormalDistributionProcess(Model, settings["Parameters"])

class AddTextureNormalDistributionProcess(KM.Process):

    """
    This process perturbates the coordinates of the nodes of the designated submodelpart according to a normal distribution with [mu,sigma] by using numpy capabilities (https://numpy.org/doc/stable/reference/random/generated/numpy.random.normal.html)

    Note: the process assumes that the submodelpart includes a set of conditions to be able to compute the normal field.


    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KM.Process.__init__(self)

        # The value can be a double or a string (function)
        default_settings = KM.Parameters("""{
            "help"                       : "This process perturbates the coordinates of the nodes according to a normal distribution",
            "model_part_name"            : "please_specify_model_part_name",
            "normal_mu"                  : 0.0,
            "normal_sigma"               : 0.1,
            "print_modified_coordinates" : true,
            "echo_level"                 : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.mu    = settings["normal_mu"].GetDouble()    # average
        self.sigma = settings["normal_sigma"].GetDouble() # standard deviation
        self.print_coords = settings["print_modified_coordinates"].GetBool()
        self.echo = settings["echo_level"].GetInt()


    def ExecuteInitialize(self):
        """This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # Here we compute the normal field
        normal_calculation_utils = KM.NormalCalculationUtils()
        normal_calculation_utils.CalculateUnitNormalsNonHistorical(self.model_part, 3)
        number_nodes = len(self.model_part.Nodes)
        normal = KM.Array3([0,0,0])

        # We generate a random normal distribution
        norm_distr = np.random.normal(self.mu, self.sigma, number_nodes)
        counter = 0

        for node in self.model_part.Nodes:
            normal = node.GetValue(KM.NORMAL)
            node.X0 += norm_distr[counter] * normal[0]
            node.Y0 += norm_distr[counter] * normal[1]
            node.Z0 += norm_distr[counter] * normal[2]
            node.X  += norm_distr[counter] * normal[0]
            node.Y  += norm_distr[counter] * normal[1]
            node.Z  += norm_distr[counter] * normal[2]
            counter += 1
        
        # Now we print a backup of the nodes coordinates used in the simulation
        


