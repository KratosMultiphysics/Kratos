# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# other imports
from importlib import import_module


def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CreatePointBasedEntitiesProcess(Model, settings["Parameters"])


class CreatePointBasedEntitiesProcess(KM.Process):
    """This process sets a given value for a certain flag in all the nodes of a submodelpart
    """

    def __init__(self, Model, settings ):
        KM.Process.__init__(self)

        default_settings = KM.Parameters("""{
            "root_model_part_name"       : "PLEASE_SPECIFY",
            "new_sub_model_part_name"    : "PLEASE_SPECIFY",
            "sub_model_part_names"       : [],
            "entity_name"                : "PointLoadCondition3D1N",
            "entity_type"                : "condition",
            "properties_id"              : 0,
            "kratos_application"         : ""
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        root_model_part = Model[settings["root_model_part_name"].GetString()]
        if root_model_part.ProcessInfo[KM.IS_RESTARTED]:
            # Do nothing in case of restart
            return

        root_model_part_name = settings["root_model_part_name"].GetString()
        entity_name = settings["entity_name"].GetString()
        entity_type = settings["entity_type"].GetString()
        properties_id = settings["properties_id"].GetInt()

        # importing the Application where the entites are registered (optional)
        kratos_application = settings["kratos_application"].GetString()
        if kratos_application != "":
            if not CheckIfApplicationsAvailable(kratos_application):
                raise Exception('Application "{}" is not available!'.format(kratos_application))
            import_module("KratosMultiphysics." + kratos_application) # this registers the entities

        if settings["sub_model_part_names"].size() == 0:
            raise Exception('"sub_model_part_names" cannot be empty!')

        model_parts = [Model[root_model_part_name+"."+model_part_name] for model_part_name in settings["sub_model_part_names"].GetStringArray()]

        new_model_part = RecursiveCreateModelParts(root_model_part, settings["new_sub_model_part_name"].GetString())

        node_ids = [node.Id
            for mp in model_parts
                for node in mp.GetCommunicator().LocalMesh().Nodes
        ]
        node_ids = list(set(node_ids)) # make sure the node-Ids are unique. This is needed to not create multiple entities for the same node

        new_model_part.AddNodes(node_ids)

        props = root_model_part.GetProperties(properties_id, 0) # 0 is mesh-id # maybe check and then create the props
        mp_comm = root_model_part.GetCommunicator()

        if entity_type == "element":
            max_id_entities = max([elem.Id for elem in root_model_part.Elements]+[0]) # "+[0]" in case there are no local entities
            creation_fct_ptr = new_model_part.CreateNewElement
        elif entity_type == "condition":
            max_id_entities = max([cond.Id for cond in root_model_part.Conditions]+[0]) # "+[0]" in case there are no local entities
            creation_fct_ptr = new_model_part.CreateNewCondition
        else:
            raise Exception('"entity_type" "{}" is not valid, only "element" or "condition" are possible!'.format(entity_type))

        # using ScanSum to compute the local Id-start. Otherwise the Ids would start with the same value on every rank
        data_comm = mp_comm.GetDataCommunicator()
        scan_sum_nodes = data_comm.ScanSum(len(node_ids))
        max_id_entities_global = data_comm.MaxAll(max_id_entities)
        local_id_start = scan_sum_nodes + 1 + max_id_entities_global - len(node_ids)

        for i, node_id in enumerate(node_ids):
            creation_fct_ptr(entity_name, i+local_id_start, [node_id], props)


def RecursiveCreateModelParts(model_part, model_part_name):
    '''This function creates a hierarchy of SubModelParts on a given ModelPart
    '''
    model_part_name, *sub_model_part_names = model_part_name.split(".")
    model_part = model_part.CreateSubModelPart(model_part_name)
    if len(sub_model_part_names) > 0:
        model_part = RecursiveCreateModelParts(model_part, ".".join(sub_model_part_names))
    return model_part
