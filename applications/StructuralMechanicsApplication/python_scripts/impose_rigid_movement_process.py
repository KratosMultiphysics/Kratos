import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeRigidMovementProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeRigidMovementProcess(KratosMultiphysics.Process):
    """This class is used in order to impose a rigid body movement in a certain region of the problem

    This class constructs the model parts containing the constrains that enforce the rigid body movement
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

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                        : "This process uses LinearMasterSlaveConstraint in order to impose an unified movement in the given submodelpart. The process takes the first node from the submodelpart if no node's ID is provided. The default variable is DISPLACEMENT, and in case no variable is considered for the slave the same variable will be considered",
            "computing_model_part_name"   : "computing_domain",
            "model_part_name"             : "please_specify_model_part_name",
            "new_model_part_name"         : "Rigid_Movement_ModelPart",
            "interval"                    : [0.0, 1e30],
            "master_variable_name"        : "DISPLACEMENT",
            "slave_variable_name"         : "",
            "master_node_id"              : 0
        }
        """)

        # Detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString()):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        # The computing model part
        computing_model_part_name = settings["computing_model_part_name"].GetString()
        self.computing_model_part = Model[computing_model_part_name]

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # The variable to be enforced (master)
        self.var_list_master = []
        master_variable_name = settings["master_variable_name"].GetString()
        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(master_variable_name)
        if (variable_type == "Double"):
            self.var_list_master.append(KratosMultiphysics.KratosGlobals.GetVariable(master_variable_name))
        elif (variable_type == "Component"):
            self.var_list_master.append(KratosMultiphysics.KratosGlobals.GetVariable(master_variable_name))
        elif (variable_type == "Array"):
            self.var_list_master.append(KratosMultiphysics.KratosGlobals.GetVariable(master_variable_name + "_X"))
            self.var_list_master.append(KratosMultiphysics.KratosGlobals.GetVariable(master_variable_name + "_Y"))
            self.var_list_master.append(KratosMultiphysics.KratosGlobals.GetVariable(master_variable_name + "_Z"))
        else:
            raise NameError("Only components variables and double variables can be considered")

        # The variable to be enforced (slave)
        self.var_list_slave = []
        slave_variable_name = settings["slave_variable_name"].GetString()
        if (slave_variable_name != ""):
            variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(slave_variable_name)
            if (variable_type == "Double"):
                self.var_list_slave.append(KratosMultiphysics.KratosGlobals.GetVariable(slave_variable_name))
            elif (variable_type == "Component"):
                self.var_list_slave.append(KratosMultiphysics.KratosGlobals.GetVariable(slave_variable_name))
            elif (variable_type == "Array"):
                self.var_list_slave.append(KratosMultiphysics.KratosGlobals.GetVariable(slave_variable_name + "_X"))
                self.var_list_slave.append(KratosMultiphysics.KratosGlobals.GetVariable(slave_variable_name + "_Y"))
                self.var_list_slave.append(KratosMultiphysics.KratosGlobals.GetVariable(slave_variable_name + "_Z"))
            else:
                raise NameError("Only components variables and double variables can be considered")
        else: # If no variable is assigned to the slave the same considered for the master will be taken
            self.var_list_slave = self.var_list_master

        # The reference node id
        self.master_node_id = settings["master_node_id"].GetInt()

        # We get the corresponding model parts
        self.model_part = Model[settings["model_part_name"].GetString()]
        if (self.model_part.HasSubModelPart(settings["new_model_part_name"].GetString())):
            self.rigid_model_part = self.model_part.GetSubModelPart(settings["new_model_part_name"].GetString())
        else:
            self.rigid_model_part = self.model_part.CreateSubModelPart(settings["new_model_part_name"].GetString())

        # Trasfering the
        transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.rigid_model_part, self.model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODES)
        transfer_process.Execute()

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # TODO: Move a cpp
        # If reference id is 0 we get the first node of the model part
        count = self.model_part.GetRootModelPart().NumberOfMasterSlaveConstraints()
        if (self.master_node_id == 0):
            for master_node in self.rigid_model_part.Nodes:
                break
            for slave_node in self.rigid_model_part.Nodes:
                if (slave_node.Id is not master_node.Id):
                    for var_master, var_slave in zip(self.var_list_master, self.var_list_slave):
                        count += 1
                        constraint = self.rigid_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", count, master_node, var_master, slave_node, var_slave, 1.0, 0.0)
                        self.computing_model_part.AddMasterSlaveConstraint(constraint)
        else:
            master_node = self.model_part.GetRootModelPart().GetNode(self.master_node_id)
            for slave_node in self.rigid_model_part.Nodes:
                if (slave_node.Id is not master_node.Id):
                    for var_master, var_slave in zip(self.var_list_master, self.var_list_slave):
                        count += 1
                        constraint = self.rigid_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", count, master_node, var_master, slave_node, var_slave, 1.0, 0.0)
                        self.computing_model_part.AddMasterSlaveConstraint(constraint)

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        # We activate/deactivate conditions dependeding of interval
        if (self.interval.IsInInterval(current_time)):
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, self.rigid_model_part.MasterSlaveConstraints)
        else:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, False, self.rigid_model_part.MasterSlaveConstraints)
