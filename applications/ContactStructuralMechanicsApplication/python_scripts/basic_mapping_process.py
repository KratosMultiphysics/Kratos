# Importing the Kratos Library
import KratosMultiphysics

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return BasicMappingProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"

class BasicMappingProcess(KratosMultiphysics.Process):
    """This process allows to do a simple mapping using the SimpleMortarMapperProcess

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

        KratosMultiphysics.Process.__init__(self)

        #The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                             : "This process allows to do a simple mapping using the SimpleMortarMapperProcess",
            "mesh_id"                          : 0,
            "origin_model_part_name"           : "please_specify_model_part_name",
            "destination_model_part_name"      : "please_specify_model_part_name",
            "interval"                         : [0.0, 1e30],
            "echo_level"                       : 0,
            "consider_tessellation"            : false,
            "using_average_nodal_normal"       : true,
            "discontinuous_interface"          : false,
            "discontinous_interface_factor"    : 1.0e-4,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2,
            "distance_threshold"               : 1.0e24,
            "remove_isolated_conditions"       : false,
            "mapping_coefficient"              : 1.0e0,
            "origin_variable"                  : "TEMPERATURE",
            "destination_variable"             : "",
            "origin_variable_historical"       : true,
            "origin_are_conditions"            : true,
            "destination_variable_historical"  : true,
            "destination_are_conditions"       : true,
            "update_interface"                 : true,
            "search_parameters"                : {
                "allocation_size"                  : 1000,
                "bucket_size"                      : 4,
                "search_factor"                    : 3.5
            },
            "linear_solver_settings": { }
        }
        """
        )

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # Validating the settings
        settings.ValidateAndAssignDefaults(default_settings)

        self.origin_model_part = Model[settings["origin_model_part_name"].GetString()]
        self.destination_model_part = Model[settings["destination_model_part_name"].GetString()]

        linear_solver_configuration = settings["linear_solver_settings"]
        if linear_solver_configuration.Has("solver_type"): # user specified a linear solver
            linear_solver = linear_solver_factory.ConstructSolver(linear_solver_configuration)
        else:
            linear_solver = None

        # Removing unused parameters
        settings.RemoveValue("help")
        settings.RemoveValue("mesh_id")
        settings.RemoveValue("origin_model_part_name")
        settings.RemoveValue("destination_model_part_name")
        settings.RemoveValue("interval")
        settings.RemoveValue("linear_solver_settings")

        # Creating the mapper
        self.mapper = KratosMultiphysics.SimpleMortarMapperProcess(self.origin_model_part, self.destination_model_part, settings, linear_solver)

        self.step_is_active = False

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        current_time = self.origin_model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):
            self.mapper.ExecuteInitializeSolutionStep()
