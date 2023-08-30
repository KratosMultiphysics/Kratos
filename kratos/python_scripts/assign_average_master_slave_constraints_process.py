import KratosMultiphysics as KM

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignAverageMasterSlaveConstraintsProcess(model, settings["Parameters"])

class AssignAverageMasterSlaveConstraintsProcess(KM.Process):
    """This class defines a process to assign averaged multi-point constraints (MPC) 
    between a slave and a master model part in a computing model part.
    Depending on the settings, these constraints can be reformed at each step or kept constant."""

    def __init__(self, model, settings):
        KM.Process.__init__(self)

        # Default settings for the process.
        default_settings = KM.Parameters("""{
            "computing_model_part_name" : "ComputingModelPartName",
            "slave_model_part_name" : "SlaveModelPartName",
            "master_model_part_name" : "MasterModelPartName",
            "variable_name" : "",
            "reform_constraints_at_each_step": false
        }""")

        # Overwrite the default settings with user-provided ones.
        settings.ValidateAndAssignDefaults(default_settings)
        
        self.model = model
        self.settings = settings

    def ExecuteInitialize(self):
        # Store model part names from the settings
        computing_model_part_name = self.settings["computing_model_part_name"].GetString()
        slave_model_part_name = self.settings["slave_model_part_name"].GetString()
        master_model_part_name = self.settings["master_model_part_name"].GetString()
        
        # Get the actual model parts from the model
        self.computing_model_part = self.model.GetModelPart(computing_model_part_name)
        self.slave_model_part = self.model.GetModelPart(slave_model_part_name)
        self.master_model_part = self.model.GetModelPart(master_model_part_name)
        
        # Decide whether to reform constraints at each step based on settings
        self.reform_constraints_at_each_step = self.settings["reform_constraints_at_each_step"].GetBool()
        
        # Get the variable from the settings
        self.variable = KM.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        
        # If constraints are not reformed at each step, we only assign them once here
        if not self.reform_constraints_at_each_step:
            self.__AverageMultiPointConstraint(self.variable)

    def ExecuteInitializeSolutionStep(self):
        # If constraints are reformed at each step, assign them at the start of each step
        if self.reform_constraints_at_each_step:
            self.__AverageMultiPointConstraint(self.variable)

    def ExecuteFinalizeSolutionStep(self):
        # If constraints are reformed at each step, remove them at the end of each step
        if self.reform_constraints_at_each_step:
            self.__RemoveConstraints()
    
    def __AverageMultiPointConstraint(self, variable):
        """Assigns constraints that make the value of a variable in the slave part
        to be the average of the values in the master part."""
        
        # Prepare containers for the constraint data
        number_of_master_nodes = self.master_model_part.NumberOfNodes()
        weights_vector = KM.Matrix(1, number_of_master_nodes)
        weight = 1/number_of_master_nodes
        weights_vector.fill(weight)
        constant_vector = KM.Vector(number_of_master_nodes)
        constant_vector.fill(0.0)
        
        # Loop over all master nodes to fill the containers
        master_dofs_container = []
        for node in self.master_model_part.Nodes:
            master_dofs_container.append(node.GetDof(variable))
        
        # Loop over all slave nodes to assign the constraints
        counter = 1
        for node in self.slave_model_part.Nodes:
            slave_dof_container = [node.GetDof(variable)]
            self.computing_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", counter, master_dofs_container, slave_dof_container, weights_vector, constant_vector)
            counter += 1

    def __RemoveConstraints(self):
        """This function removes the master-slave constraints from the computing model part"""
        KM.VariableUtils().SetFlag(KM.TO_ERASE, True, self.computing_model_part.MasterSlaveConstraints)
        self.computing_model_part.RemoveMasterSlaveConstraintsFromAllLevels(KM.TO_ERASE)
