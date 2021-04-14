import KratosMultiphysics as KM

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SubModelPartEntitiesBooleanOperationProcess(model, settings["Parameters"])

class SubModelPartEntitiesBooleanOperationProcess(KM.Process):
    """Apply a boolean operation to the entities of the given sub model parts.
    The sub model parts must have the same root model part and the entities
    are compared by their Id.
    The possible boolean operations are:
        - 'Union'
        - 'Intersection'
        - 'Difference'
    The possible entities are:
        - 'Nodes'
        - 'Elements'
        - 'Conditions'
        - 'MasterSlaveConstraints'
    """

    def __init__(self, model, settings):
        """Constructor of SubModelPartEntitiesBooleanOperationProcess.

        The boolean operation is executed in the constructor. In this way,
        the sub_model_parts where the BC's are applied, will be affected
        by the boolean operation.
        """
        KM.Process.__init__(self)
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.first_model_part = model.GetModelPart(settings["first_model_part_name"].GetString())
        self.second_model_part = model.GetModelPart(settings["second_model_part_name"].GetString())
        self.result_model_part = self._GetOrCreateSubModelPart(model, settings["result_model_part_name"].GetString())
        self.boolean_operation = settings["boolean_operation"].GetString()
        self.entity_type = settings["entity_type"].GetString()

        self._CheckInput()
        self._ExecuteBooleanOperation()

    @staticmethod
    def GetDefaultParameters():
        """Return the default parameters."""
        return KM.Parameters("""{
            "first_model_part_name"  : "MODEL_PART_NAME",
            "second_model_part_name" : "MODEL_PART_NAME",
            "result_model_part_name" : "MODEL_PART_NAME",
            "boolean_operation"      : "Difference",
            "entity_type"            : "Nodes"
        }""")

    def _CheckInput(self):
        possible_boolean_operations = ["Union", "Intersection", "Difference"]
        if self.boolean_operation not in possible_boolean_operations:
            raise Exception("The possible boolean operations are:\n\t- '{}'\n\t- '{}'\n\t- '{}'".format(*possible_boolean_operations))
        possible_entities_type = ["Nodes", "Elements", "Conditions", "MasterSlaveConstraints"]
        if self.entity_type not in possible_entities_type:
            raise Exception("The possible entities type are:\n\t- '{}'\n\t- '{}'\n\t- '{}'".format(*possible_entities_type))

    def _ExecuteBooleanOperation(self):
        if self.entity_type == "Nodes":
            utility = KM.SubModelPartNodesBooleanOperationUtility
        elif self.entity_type == "Elements":
            utility = KM.SubModelPartElementsBooleanOperationUtility
        elif self.entity_type == "Conditions":
            utility = KM.SubModelPartConditionsBooleanOperationUtility
        elif self.entity_type == "MasterSlaveConstraints":
            utility = KM.SubModelPartConstraintsBooleanOperationUtility
        if self.boolean_operation == "Union":
            utility.Union(self.first_model_part, self.second_model_part, self.result_model_part)
        if self.boolean_operation == "Intersection":
            utility.Intersection(self.first_model_part, self.second_model_part, self.result_model_part)
        if self.boolean_operation == "Difference":
            utility.Difference(self.first_model_part, self.second_model_part, self.result_model_part)

    @staticmethod
    def _GetOrCreateSubModelPart(model, full_name):
        if model.HasModelPart(full_name):
            return model.GetModelPart(full_name)
        else:
            if full_name.find('.') == -1:
                raise Exception("The result model part is a root model part. Only sub model parts are allowed and the input name is '{}'".format(full_name))
            parent_name, name = full_name.rsplit('.', 1)
            parent_model_part = model.GetModelPart(parent_name)
            return parent_model_part.CreateSubModelPart(name)
