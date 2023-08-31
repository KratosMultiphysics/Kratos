import KratosMultiphysics as KM
import numpy as np

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignMPCsToNeighboursProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignMPCsToNeighboursProcess(KM.Process):
    """The process facilitates the discovery of neighboring nodes in a 
    master model part within a designated radius for each node in the 
    slave model part. Following this, it establishes a multipoint constraint (MPC) 
    and calculates its weights using a spatial function that employs radial basis functions.

    Public member variables:
    model -- the container of the different model parts.
    settings -- Kratos parameters containing process settings.
    """

    def __init__(self, model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        KM.Process.__init__(self) # calling the baseclass constructor

        default_settings = KM.Parameters("""{
            "computing_model_part_name": "",
            "slave_model_part_name": "",
            "master_model_part_name": "",
            "variable_names": [],
            "search_radius": 1.0,
            "minimum_number_of_neighbouring_nodes": 3,
            "assign_mpcs_every_time_step": false
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        settings.ValidateAndAssignDefaults(default_settings)
        # Alternative:
        # settings.RecursivelyValidateAndAssignDefaults(default_settings)

        #Assign settings
        self.model = model

        # Get the model part on which the mpcs are going to be applied
        if not settings["computing_model_part_name"].GetString():
            raise Exception("\'computing_model_part_name\' not provided. Please specify the model part to apply to MPCs to.")
        self.computing_model_part_name = settings["computing_model_part_name"].GetString() #mpcs are applied to computing model part 
        
        # Get the slave model part 
        if not settings["slave_model_part_name"].GetString():
            raise Exception("\'slave_model_part_name\' not provided. Please specify the slave model part.")
        self.slave_model_part_name = settings["slave_model_part_name"].GetString()
        
        # Get the master model part 
        if not settings["master_model_part_name"].GetString():
            raise Exception("\'master_model_part_name\' not provided. Please specify the master model part.")
        self.master_model_part_name = settings["master_model_part_name"].GetString()

        # Search radius for the mpcs
        self.search_radius = settings["search_radius"].GetDouble()

        # Minimum number of neighboring nodes retrieved from the search radius
        self.minimum_number_of_neighbouring_nodes = settings["minimum_number_of_neighbouring_nodes"].GetInt()

        # Apply mpcs at each time step (True) or only once (False)
        self.assign_mpcs_every_time_step = settings["assign_mpcs_every_time_step"].GetBool()

        # Retrieve and check if variables exist
        variable_names = settings["variable_names"].GetStringArray()
        if len(variable_names) == 0:
            err_msg = "The variable names need to be specified by the user in the \'variable_names\' string array."
            raise Exception(err_msg)
        if any(variable_names.count(var_name) > 1 for var_name in variable_names):
            err_msg = "There are repeated variables in the \'variable_names\' string array."
            raise Exception(err_msg)
        variable_names.sort()

        self.variables_list = []
        for var_name in variable_names:
            if not KM.KratosGlobals.HasVariable(var_name):
                err_msg = "\'{}\' variable in \'variable_names\' is not in KratosGlobals. Please check provided value.".format(var_name)
            if not KM.KratosGlobals.GetVariableType(var_name):
                err_msg = "\'{}\' variable in \'variable_names\' is not double type. Please check provide double type variables (e.g. [\"DISPLACEMENT_X\",\"DISPLACEMENT_Y\"]).".format(var_name)
            self.variables_list.append(KM.KratosGlobals.GetVariable(var_name))
        debug = True

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        ## Assign master and slave model parts
        self.computing_model_part = self.model.GetModelPart(self.computing_model_part_name)
        self.slave_model_part = self.model.GetModelPart(self.slave_model_part_name)
        self.master_model_part = self.model.GetModelPart(self.master_model_part_name)

        # Initialize MPCs to neighbours utility
        self.assign_mpcs_utility = KM.AssignMPCsToNeighboursUtility(self.master_model_part.Nodes)

        # The user may only need to set up the mpcs only once
        if not self.assign_mpcs_every_time_step:
            for variable in self.variables_list:
                self.assign_mpcs_utility.AssignMPCsToNodes(self.slave_model_part.Nodes,self.search_radius,self.computing_model_part, variable, self.minimum_number_of_neighbouring_nodes)

    def Check(self):
        """ This method verifies that the input is correct

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed just before the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # If the user want the mpcs to be updated at each time step, this is usefull for moving meshes.
        if self.assign_mpcs_every_time_step:
            for variable in self.variables_list:
                self.assign_mpcs_utility.AssignMPCsToNodes(self.slave_model_part.Nodes,self.search_radius,self.computing_model_part, variable, self.minimum_number_of_neighbouring_nodes)

    def ExecuteBeforeOutputStep(self):
        """ This method is executed before writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterOutputStep(self):
        """ This method is executed after writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # If mpcs are updated every time step, these are to me removed before being re-assigned.
        if self.assign_mpcs_every_time_step:
            self.RemoveConstraints()
        
    def RemoveConstraints(self):
        #Remove master-slave constraints
        KM.VariableUtils().SetFlag(KM.TO_ERASE, True, self.computing_model_part.MasterSlaveConstraints)
        self.computing_model_part.RemoveMasterSlaveConstraintsFromAllLevels(KM.TO_ERASE)

    def ExecuteFinalize(self):
        """ This method is executed after the computations, at the end of the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass