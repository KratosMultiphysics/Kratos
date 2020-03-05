import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera
import numpy as np
import os
import json
from copy import deepcopy

def ImportChimeraModelparts(main_modelpart, chimera_mp_import_settings_list, material_file="", parallel_type="OpenMP"):
    '''
        This function extends the functionalities of the
        mpda_manipulator from: https://github.com/philbucher/mdpa-manipulator

        main_modelpart      : The modelpart to which the new modelparts are appended to.
        chimera_mp_import_settings_list   : The list of import setting for all chimera modelparts. each entry has the following format:
        {
            "model_import_settings":{
                    "input_type": "mdpa",
                    "input_filename": "SOME"
            },
            "echo_level":1
        }
    '''
    if parallel_type == "OpenMP":
        for mp_import_setting in chimera_mp_import_settings_list:
            mdpa_file_name = mp_import_setting["input_filename"].GetString()
            if mdpa_file_name.endswith('.mdpa'):
                mdpa_file_name = mdpa_file_name[:-5]

            model = KratosMultiphysics.Model()
            model_part = model.CreateModelPart("new_modelpart")
            KratosChimera.TransferSolutionStepData(main_modelpart, model_part)

            ReadModelPart(mdpa_file_name, model_part, material_file)
            AddModelPart(main_modelpart, model_part)
    elif(parallel_type == "MPI"):
        input_settings = KratosMultiphysics.Parameters("""{
        "model_import_settings":{
            "input_type": "mdpa",
            "input_filename": "SOME"
        },
        "echo_level":1
        }""")

        for mp_import_setting in chimera_mp_import_settings_list:
            model = KratosMultiphysics.Model()
            model_part = model.CreateModelPart("new_modelpart")
            KratosChimera.TransferSolutionStepData(main_modelpart, model_part)
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
            mdpa_file_name = mp_import_setting["input_filename"].GetString()
            if mdpa_file_name.endswith('.mdpa'):
                mdpa_file_name = mdpa_file_name[:-5]
            mp_import_setting["input_filename"].SetString(mdpa_file_name)
            input_settings["model_import_settings"] = mp_import_setting

            from KratosMultiphysics.mpi import distributed_import_model_part_utility

            mpi_import_utility = distributed_import_model_part_utility.DistributedImportModelPartUtility(model_part, input_settings)
            mpi_import_utility.ImportModelPart()
            #mpi_import_utility.CreateCommunicators()
            AddModelPart(main_modelpart, model_part, is_mpi=True)

            ## Construct and execute the Parallel fill communicator (also sets the MPICommunicator)
            import KratosMultiphysics.mpi as KratosMPI
            ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(main_modelpart.GetRootModelPart())
            ParallelFillCommunicator.Execute()


def ReadModelPart(mdpa_file_name, model_part, materials_file_name=""):
    '''
    Read and return a ModelPart from a mdpa file
    '''
    # We reorder because otherwise the numbering might be screwed up when we combine the ModelParts later
    KratosMultiphysics.ReorderConsecutiveModelPartIO(mdpa_file_name, KratosMultiphysics.IO.SKIP_TIMER).ReadModelPart(model_part)

    if materials_file_name != "":
        # in case a materials-file is to be combined, it is read and saved as a string
        # for this the ProcessInfo is used => bcs it is shared among (Sub-)ModelParts
        with open(materials_file_name,'r') as materials_file:
            materials_string = json.dumps(json.load(materials_file))
        model_part.ProcessInfo[KratosMultiphysics.IDENTIFIER] = materials_string
        model_part[KratosMultiphysics.IDENTIFIER] = materials_string

    __RemoveAuxFiles()

def AddModelPart(model_part_1,
                 model_part_2,
                 add_as_submodelpart=False, is_mpi=False):
    '''
    Adding the model_part_2 to model_part_1 (appending)
    '''
    if ( not isinstance(model_part_1, KratosMultiphysics.ModelPart) ):
            raise Exception("input is expected to be provided as a Kratos ModelPart object")
    if ( not isinstance(model_part_2, KratosMultiphysics.ModelPart) ):
            raise Exception("input is expected to be provided as a Kratos ModelPart object")

    comm = model_part_1.GetCommunicator().GetDataCommunicator()
    num_nodes_self = model_part_1.NumberOfNodes()
    total_num_nodes_self = comm.SumAll(num_nodes_self)
    num_elements_self = model_part_1.NumberOfElements()
    total_num_elements_self = comm.SumAll(num_elements_self)
    num_conditions_self = model_part_1.NumberOfConditions()
    total_num_conditions_self = comm.SumAll(num_conditions_self)

    node_id_pi_map ={}

    for node in model_part_2.Nodes:
        node.Id += total_num_nodes_self
        if(is_mpi):
            node_id_pi_map[node.Id] = node.GetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX)
    for element in model_part_2.Elements:
        element.Id += total_num_elements_self
    for condition in model_part_2.Conditions:
        condition.Id += total_num_conditions_self

    KratosMultiphysics.FastTransferBetweenModelPartsProcess(model_part_1, model_part_2,
        KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.ALL).Execute()

    if(is_mpi):
        for node_id, pi in node_id_pi_map.items():
            node = model_part_1.Nodes[node_id]
            node.SetSolutionStepValue(KratosMultiphysics.PARTITION_INDEX, pi)

    if add_as_submodelpart: # add one one lovel lower
        # adding model_part_2 as submodel part to model_part_1 (called recursively)
        __AddAsSubModelPart(model_part_1, model_part_2)
        if model_part_2.ProcessInfo.Has(KratosMultiphysics.IDENTIFIER):
            model_part_name = model_part_1.Name + "." + model_part_2.Name
            CombineMaterialProperties(model_part_1, model_part_2, model_part_name)

    else: # add on same level
        # adding submodel parts of model_part_2 to model_part_1 (called recursively)
        __AddSubModelPart(model_part_1, model_part_2)
        if model_part_2.ProcessInfo.Has(KratosMultiphysics.IDENTIFIER):
            model_part_name = model_part_1.Name
            CombineMaterialProperties(model_part_1,model_part_2,model_part_name)


def __AddEntitiesToSubModelPart(original_sub_model_part,
                                other_sub_model_part):
    '''
    Adds the entities of (nodes, elements and conditions) from
    one SubModelPart to another
    '''
    # making list containing node IDs of particular submodel part
    num_nodes_other = other_sub_model_part.NumberOfNodes()
    smp_node_id_array = np.zeros(num_nodes_other, dtype=np.int)
    for node_i, node in enumerate(other_sub_model_part.Nodes):
        smp_node_id_array[node_i] = node.Id

    # making list containing element IDs of particular submodel part
    num_elements_other = other_sub_model_part.NumberOfElements()
    smp_element_id_array = np.zeros(num_elements_other, dtype=np.int)
    for element_i, element in enumerate(other_sub_model_part.Elements):
        smp_element_id_array[element_i] = element.Id

    # making list containing condition IDs of particular submodel part
    num_conditions_other = other_sub_model_part.NumberOfConditions()
    smp_condition_id_array = np.zeros(num_conditions_other, dtype=np.int)
    for condition_i, condition in enumerate(other_sub_model_part.Conditions):
        smp_condition_id_array[condition_i] = condition.Id

    original_sub_model_part.AddNodes(smp_node_id_array.tolist())
    original_sub_model_part.AddElements(smp_element_id_array.tolist())
    original_sub_model_part.AddConditions(smp_condition_id_array.tolist())

def __AddSubModelPart(original_model_part,
                      other_model_part):
    '''
    Adds the SubModelParts of one ModelPart to the other one
    If the original ModelPart already contains a SMP with the same name,
    the entities are added to it
    '''
    for smp_other in other_model_part.SubModelParts:
        if original_model_part.HasSubModelPart(smp_other.Name):
            smp_original = original_model_part.GetSubModelPart(smp_other.Name)

            # in case we add sth to an existing SubModelPart, we have to make sure that the materials are the same!
            smp_orig_has_materials = smp_original.ProcessInfo.Has(KratosMultiphysics.IDENTIFIER)
            other_mp_has_materials = smp_other.ProcessInfo.Has(KratosMultiphysics.IDENTIFIER)

            if smp_orig_has_materials and other_mp_has_materials: # both MPs have materials, checking if they are the same
                orig_material = json.loads(original_model_part.ProcessInfo[KratosMultiphysics.IDENTIFIER])
                other_material = json.loads(smp_other.ProcessInfo[KratosMultiphysics.IDENTIFIER])

                if not __MaterialsListsAreEqual(orig_material["properties"], other_material["properties"]):
                    err_msg  = 'Trying to add "' + smp_other.GetRootModelPart().Name + '" to "'
                    err_msg += original_model_part.GetRootModelPart().Name + '" but their materials are different!'
                    raise Exception(err_msg)

            elif smp_orig_has_materials and not other_mp_has_materials:
                err_msg  = 'Trying to add "' + smp_other.GetRootModelPart().Name + '" (has NO materials) to "'
                err_msg += original_model_part.GetRootModelPart().Name + '" (has materials)'
                raise Exception(err_msg)
            elif not smp_orig_has_materials and other_mp_has_materials:
                err_msg  = 'Trying to add "' + smp_other.GetRootModelPart().Name + '" (has materials) to "'
                err_msg += original_model_part.GetRootModelPart().Name + '" (has NO materials)'
                raise Exception(err_msg)
            else:
                pass # => none has materials, no checking required

        else:
            smp_original = original_model_part.CreateSubModelPart(smp_other.Name)

        __AddEntitiesToSubModelPart(smp_original, smp_other)

        __AddSubModelPart(smp_original, smp_other) # call recursively to transfer nested SubModelParts

def __AddAsSubModelPart(original_model_part,
                        other_model_part):
    '''
    Adds the SubModelParts of one ModelPart to the other one
    If the original ModelPart already contains a SMP with the same name,
    the entities are added to it
    '''
    smp_original = original_model_part.CreateSubModelPart(other_model_part.Name)

    __AddEntitiesToSubModelPart(smp_original, other_model_part)

    for smp_other in other_model_part.SubModelParts:
        __AddAsSubModelPart(smp_original, smp_other)   #call recursively to transfer nested SubModelParts


def __RemoveAuxFiles():
    '''
    Removes auxiliary files from the directory
    '''
    current_path = os.getcwd()
    files = os.listdir(current_path)
    for file in files:
        if file.endswith(".time") or file.endswith(".lst"):
            os.remove(os.path.join(current_path, file))

def CombineMaterialProperties(to_model_part, from_model_part, model_part_name):

    if not to_model_part.ProcessInfo.Has(KratosMultiphysics.IDENTIFIER):
        # in case the ModelPart does not have materials specified, add empty ones
        empty_props = {"properties" : []}
        to_model_part.ProcessInfo[KratosMultiphysics.IDENTIFIER] = json.dumps(empty_props)

    existing_props = json.loads(to_model_part.ProcessInfo[KratosMultiphysics.IDENTIFIER])
    num_existing_props = len(existing_props["properties"])
    props_by_name_dict = {mat["model_part_name"] : mat for mat in existing_props["properties"]}

    new_props = json.loads(from_model_part.ProcessInfo[KratosMultiphysics.IDENTIFIER])
    new_props_to_add = []

    for i, props in enumerate(new_props["properties"]):
        current_model_part_name = props["model_part_name"].split(".")

        if current_model_part_name[0] == to_model_part.Name: # this means that the name of the MainModelPart is added
            current_model_part_name.pop(0) # remove the MainModelPart-Name

        new_model_part_name = model_part_name + "." + ".".join([name for name in current_model_part_name])
        props["model_part_name"] = new_model_part_name # done here bcs also needed for check
        if new_model_part_name in props_by_name_dict: # properties for this ModelPart exist already
            # check (again, this should have been checked before and should not fail here!) if the props are the same
            if not __MaterialsListsAreEqual([props_by_name_dict[new_model_part_name]], [props]):
                err_msg  = 'Different properties for ModelPart "' + new_model_part_name + '" exist!\n'
                err_msg += 'This should not happen here, the error should have been thrown earlier when adding the ModelParts'
                raise Exception(err_msg)
        else:
            props["new_properties_id"] = num_existing_props + i + 1
            new_props_to_add.append(props)

    existing_props["properties"].extend(new_props_to_add)
    to_model_part.ProcessInfo[KratosMultiphysics.IDENTIFIER] = json.dumps(existing_props)


def __MaterialsListsAreEqual(original_materials,
                             other_materials):
    '''
    In order to compare the materials the "new_properties_id" is removed
    '''
    copy_original_materials = deepcopy(original_materials) # make a copy to not modify the original list
    for mat in copy_original_materials:
        mat.pop("new_properties_id")

    return copy_original_materials == other_materials
