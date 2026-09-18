import KratosMultiphysics as Kratos

def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> Kratos.Process:
    return ImposeRBE3Process(model, settings["Parameters"])

class ImposeRBE3Process(Kratos.Process):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__()

        default_paramters = Kratos.Parameters("""{
            "main_model_part_name"        : "Structure",
            "model_part_name"             : "please_specify_model_part_name",
            "interval"                    : [0.0, 1e30],
            "master_variable_name"        : "DISPLACEMENT",
            "slave_variable_name"         : "",
            "relation"                    : 1.0,
            "constant"                    : 0.0,
            "slave_node_id"               : 0
        }""")

        parameters.ValidateAndAssignDefaults(default_paramters)

        self.main_model_part = model[parameters["main_model_part_name"].GetString()]
        self.rigid_body_model_part = self.main_model_part.GetSubModelPart(parameters["model_part_name"].GetString())
        slave_node = self.main_model_part.GetNode(parameters["slave_node_id"].GetInt())

        # get master variables
        master_dofs_list = self.__GetDofs(parameters["master_variable_name"].GetString())

        # get slave variables
        slave_variable_name = parameters["slave_variable_name"].GetString()
        if slave_variable_name == "":
            slave_dofs_list = master_dofs_list
        else:
            slave_dofs_list = self.__GetDofs()

        # get the global number of master nodes. Works in MPI too.
        number_of_constraints = self.rigid_body_model_part.GetCommunicator().GlobalNumberOfNodes()
        if self.rigid_body_model_part.HasNode(slave_node.Id):
            number_of_constraints -= 1

        # reading user input for relation and constant matrix.
        # the MPCs are formulated as following:
        # u_slave = relation * u_master + constant
        relation = parameters["relation"].GetDouble() / number_of_constraints
        constant = parameters["constant"].GetDouble()

        master_slave_constraint_id = self.main_model_part.GetCommunicator().GlobalNumberOfMasterSlaveConstraints() + 1
        for master_node in self.rigid_body_model_part.Nodes:
            # creating master slave constraint with same slave, multiple masters
            if master_node.Id != slave_node.Id:
                # adding for each dof a master slave constraint
                for master_dof, slave_dof in zip(master_dofs_list, slave_dofs_list):
                    # constraint needs to be added to the main model part. Otherwise,
                    # it will be ignored when constructing the system matrix.
                    constraint = self.main_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", master_slave_constraint_id, master_node, master_dof, slave_node, slave_dof, relation, constant)

                    # adding same to our rigid body model part. So we can ACTIVATE or
                    # deactivate as we wish later based on the interval.
                    self.rigid_body_model_part.AddMasterSlaveConstraint(constraint)
                    master_slave_constraint_id += 1

        # Assign this here since it will change the "interval" prior to validation
        # this reads the "interval" settings
        self.interval = Kratos.IntervalUtility(parameters)

    def ExecuteInitializeSolutionStep(self):
        current_time = self.main_model_part.ProcessInfo[Kratos.TIME]

        # We activate/deactivate conditions dependeding of interval
        if self.interval.IsInInterval(current_time):
            Kratos.VariableUtils().SetFlag(Kratos.ACTIVE, True, self.rigid_body_model_part.MasterSlaveConstraints)
        else:
            Kratos.VariableUtils().SetFlag(Kratos.ACTIVE, False, self.rigid_body_model_part.MasterSlaveConstraints)

    @staticmethod
    def __GetDofs(variable_name: str) -> 'list[Kratos.DoubleVariable]':
        # get the Kratos variable type from a string variable name
        variable_type = Kratos.KratosGlobals.GetVariableType(variable_name)

        # use to get the kratos variable given the string name
        var_name_getter = Kratos.KratosGlobals.GetVariable
        if variable_type == "Array":
            return [var_name_getter(f"{variable_name}_X"), var_name_getter(f"{variable_name}_Y"), var_name_getter(f"{variable_name}_Z")]
        elif variable_type == "Double":
            return [var_name_getter(variable_name)]
        else:
            raise RuntimeError(f"Unsupported variable = \"{variable_name}\".")



