import KratosMultiphysics as KM

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SubModelPartSetOperationProcess(model, settings["Parameters"])

class SubModelPartSetOperationProcess(KM.Process):
    """SubModelPartSetOperationProcess

    Apply a set operation to the entities of the given sub model parts.
    The possible set operations are:
        - 'Union'
        - 'Intersection'
        - 'Difference'
    The possible entities are:
        - 'Nodes'
        - 'Elements'
        - 'Conditions'
    """

    def __init__(self, model, settings):
        """Constructor of SubModelPartSetOperationProcess.

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The model to be used
        settings -- The ProjectParameters used

        The set operation is executed in the constructor. In this way,
        the sub_model_parts where the BC's are applied, will be affected
        by the set operation.
        """
        KM.Process.__init__(self)
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.first_model_part = model.GetModelPart(settings["first_model_part_name"].GetString())
        self.second_model_part = model.GetModelPart(settings["second_model_part_name"].GetString())
        self.result_model_part = self._GetOrCreateSubModelPart(model, settings["result_model_part_name"].GetString())
        self.set_operation = settings["set_operation"].GetString()
        self.entity_type = settings["entity_type"].GetString()

        self._CheckInput()
        self._ExecuteSetOperation()

    @staticmethod
    def GetDefaultParameters():
        """Return the default parameters"""
        return KM.Parameters("""{
            "first_model_part_name"  : "MODEL_PART_NAME",
            "second_model_part_name" : "MODEL_PART_NAME",
            "result_model_part_name" : "MODEL_PART_NAME",
            "set_operation"          : "Difference",
            "entity_type"            : "Nodes"
        }""")

    def _CheckInput(self):
        possible_set_operations = ["Union", "Intersection", "Difference"]
        if self.set_operation not in possible_set_operations:
            raise Exception("The possible set operations are:\n\t- '{}'\n\t- '{}'\n\t- '{}'".format(*possible_set_operations))
        possible_entities_type = ["Nodes", "Elements", "Conditions"]
        if self.entity_type not in possible_entities_type:
            raise Exception("The possible entities type are:\n\t- '{}'\n\t- '{}'\n\t- '{}'".format(*possible_entities_type))

    def _ExecuteSetOperation(self):
        if self.entity_type == "Nodes":
            if self.set_operation == "Union":
                KM.SubModelPartOperationsUtility.NodesUnion(self.first_model_part, self.second_model_part, self.result_model_part)
            elif self.set_operation == "Intersection":
                KM.SubModelPartOperationsUtility.NodesIntersection(self.first_model_part, self.second_model_part, self.result_model_part)
            elif self.set_operation == "Difference":
                KM.SubModelPartOperationsUtility.NodesDifference(self.first_model_part, self.second_model_part, self.result_model_part)
        elif self.entity_type == "Elements":
            if self.set_operation == "Union":
                KM.SubModelPartOperationsUtility.ElementsUnion(self.first_model_part, self.second_model_part, self.result_model_part)
            elif self.set_operation == "Intersection":
                KM.SubModelPartOperationsUtility.ElementsIntersection(self.first_model_part, self.second_model_part, self.result_model_part)
            elif self.set_operation == "Difference":
                KM.SubModelPartOperationsUtility.ElementsDifference(self.first_model_part, self.second_model_part, self.result_model_part)
        elif self.entity_type == "Conditions":
            if self.set_operation == "Union":
                KM.SubModelPartOperationsUtility.ConditionsUnion(self.first_model_part, self.second_model_part, self.result_model_part)
            elif self.set_operation == "Intersection":
                KM.SubModelPartOperationsUtility.ConditionsIntersection(self.first_model_part, self.second_model_part, self.result_model_part)
            elif self.set_operation == "Difference":
                KM.SubModelPartOperationsUtility.ConditionsDifference(self.first_model_part, self.second_model_part, self.result_model_part)

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
