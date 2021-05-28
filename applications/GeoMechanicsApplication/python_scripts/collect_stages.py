
class Sub_model_part():
    def __init__(self):
        self.name =''
        self.tables =[]
        self.nodes =[]
        self.elements =[]
        self.conditions =[]

    def read_submodel_part(self ,mdpa_data ,line_nr=0):
        # while "End SubModelPart\n" not in mdpa_data[line_nr] and line_nr < len(mdpa_data):
        self.name = mdpa_data[line_nr].replace("Begin SubModelPart ", "")
        line_nr += 2

        while "End SubModelPartTables\n" not in mdpa_data[line_nr] and line_nr < len(mdpa_data):
            self.tables.append(mdpa_data[line_nr])
            line_nr += 1

        line_nr += 2
        if line_nr < len(mdpa_data):
            while "End SubModelPartNodes\n" not in mdpa_data[line_nr] and line_nr < len(mdpa_data):
                self.nodes.append(int(mdpa_data[line_nr]))
                line_nr += 1

        line_nr += 2
        if line_nr < len(mdpa_data):
            while "End SubModelPartElements\n" not in mdpa_data[line_nr] and line_nr < len(mdpa_data):
                self.elements.append(int(mdpa_data[line_nr]))
                line_nr += 1

        line_nr += 2
        if line_nr < len(mdpa_data):
            while "End SubModelPartConditions\n" not in mdpa_data[line_nr] and line_nr < len(mdpa_data):
                self.conditions.append(int(mdpa_data[line_nr]))
                line_nr += 1


def clear_node_data_from_mdpa(mdpa_data, line=0):
    """

    :param mdpa_data:
    :param line:
    :return:
    """
    while ('Begin Nodes' not in mdpa_data[line]) and line < len(mdpa_data):
        line += 1

    line += 1

    while 'End Nodes' not in mdpa_data[line] and line < len(mdpa_data):
        del(mdpa_data[line])


def clear_element_data_from_mdpa(mdpa_data,nr_element_types, line=0):
    """

    :param mdpa_data:
    :param line:
    :return:
    """

    for i in range(nr_element_types):
        while ('Begin Elements' not in mdpa_data[line]) and line < len(mdpa_data):
            line += 1

        line += 1

        while 'End Elements' not in mdpa_data[line] and line < len(mdpa_data):
            del(mdpa_data[line])




def count_element_types(mdpa_data):
    nr_element_types = 0
    for line in mdpa_data:
        if "Begin Elements" in line:
            nr_element_types+=1
    return nr_element_types

def collect_submodel_parts_from_mdpa(mdpa_data):
    """

    :param mdpa_data:
    :return:
    """
    all_submodel_parts = []
    for line_nr ,line in enumerate(mdpa_data):
        if "Begin SubModelPart " in line:
            sub_model_part = Sub_model_part()
            sub_model_part.read_submodel_part(mdpa_data, line_nr=line_nr)

            all_submodel_parts.append(sub_model_part)
    return  all_submodel_parts

def remove_sub_model_part(mdpa_data ,sub_model_part ,line=0):
    """

    :param mdpa_data:
    :param sub_model_part:
    :param line:
    :return:
    """
    while (sub_model_part not in mdpa_data[line]) and line < len(mdpa_data):
        line += 1
    while "End SubModelPart\n" not in mdpa_data[line] and line < len(mdpa_data):
        del(mdpa_data[line])
    if line < len(mdpa_data):
        del(mdpa_data[line])

def compare_and_edit_sub_model_parts(duplicated_submodel_part):
    """

    :param duplicated_submodel_parts:
    :return:
    """

    duplicated_nodes = [node_new for node in duplicated_submodel_part[1].nodes
                                 for node_new in duplicated_submodel_part[0].nodes
                                 if node == node_new]

    [duplicated_submodel_part[0].nodes.remove(duplicated_node) for duplicated_node in duplicated_nodes]
    duplicated_submodel_part[1].nodes.extend(duplicated_submodel_part[0].nodes)

    duplicated_elements = [element_new for element in duplicated_submodel_part[1].elements
                                 for element_new in duplicated_submodel_part[0].elements
                                 if element == element_new]

    [duplicated_submodel_part[0].elements.remove(duplicated_element) for duplicated_element in duplicated_elements]
    duplicated_submodel_part[1].elements.extend(duplicated_submodel_part[0].elements)

    duplicated_conditions = [condition_new for condition in duplicated_submodel_part[1].conditions
                                 for condition_new in duplicated_submodel_part[0].conditions
                                 if condition == condition_new]

    [duplicated_submodel_part[0].conditions.remove(duplicated_condition) for duplicated_condition in duplicated_conditions]
    duplicated_submodel_part[1].conditions.extend(duplicated_submodel_part[0].conditions)

    if (not duplicated_submodel_part[0].nodes and not duplicated_submodel_part[0].elements and not
        duplicated_submodel_part[0].conditions):
        submodel_part_new_is_empty = True
    else:
        submodel_part_new_is_empty = False

    return submodel_part_new_is_empty


def add_new_sub_model_part(mdpa_data,submodel_part):
    """

    :param mdpa_data:
    :param submodel_part:
    :return:
    """
    mdpa_data.append('Begin SubModelPart ' + submodel_part.name)

    mdpa_data.append('  Begin SubModelPartTables\n')
    [mdpa_data.append('    ' + str(table) + '\n') for table in submodel_part.tables]
    mdpa_data.append('  End SubModelPartTables\n')

    mdpa_data.append('  Begin SubModelPartNodes\n')
    [mdpa_data.append('    ' + str(node)+'\n') for node in submodel_part.nodes]
    mdpa_data.append('  End SubModelPartNodes\n')

    mdpa_data.append('  Begin SubModelPartElements\n')
    [mdpa_data.append('    ' + str(element) + '\n') for element in submodel_part.elements]
    mdpa_data.append('  End SubModelPartElements\n')

    mdpa_data.append('  Begin SubModelPartConditions\n')
    [mdpa_data.append('    ' + str(condition) + '\n') for condition in submodel_part.conditions]
    mdpa_data.append('  End SubModelPartConditions\n')

    mdpa_data.append('End SubModelPart\n')
    mdpa_data.append('\n')


def read_all_projects(all_project_names):
    """

    :param all_project_names:
    :return:
    """
    import json
    import os

    all_project_parameters = []
    all_material_parameters = []
    all_mdpa_data = []

    for stage ,project_name in enumerate(all_project_names):
        with open(os.path.join(main_folder, project_name, 'ProjectParameters.json'), 'r') as f:
            project_parameters = json.load(f)
        with open(os.path.join(main_folder, project_name, 'MaterialParameters.json'), 'r') as f:
            material_parameters = json.load(f)
        with open(os.path.join(main_folder, project_name, project_name.replace(".gid", ".mdpa")), 'r') as f:
            mdpa_data = f.readlines()

        all_project_parameters.append(project_parameters)
        all_material_parameters.append(material_parameters)
        all_mdpa_data.append(mdpa_data)

    return all_project_parameters, all_material_parameters, all_mdpa_data

def clear_duplications_in_mdpa_files(all_project_names, all_mdpa_data):
    submodel_parts_all = collect_submodel_parts_from_mdpa(all_mdpa_data[0])

    # all_project_parameters[0]["solver_settings"]["material_import_settings"][
    #     "materials_filename"] = 'MaterialParameters_stage1.json'

    for project_nr in range(1, len(all_project_names)):

        clear_node_data_from_mdpa(all_mdpa_data[project_nr])
        nr_element_types = count_element_types(all_mdpa_data[project_nr])
        clear_element_data_from_mdpa(all_mdpa_data[project_nr],nr_element_types)

        submodel_parts = collect_submodel_parts_from_mdpa(all_mdpa_data[project_nr])

        duplicated_submodel_parts = [[submodel_part_new, submodel_part] for submodel_part in submodel_parts_all
                                     for submodel_part_new in submodel_parts
                                     if submodel_part.name == submodel_part_new.name]

        for duplicated_submodel_part in duplicated_submodel_parts:

            submodel_part_new_is_empty = compare_and_edit_sub_model_parts(duplicated_submodel_part)

            remove_sub_model_part(all_mdpa_data[project_nr], duplicated_submodel_part[0].name)
            if not submodel_part_new_is_empty:
                add_new_sub_model_part(all_mdpa_data[project_nr], duplicated_submodel_part[0])

        all_names = [submodel_part.name for submodel_part in submodel_parts_all]
        unique_submodel_parts = [submodel_part_new for submodel_part_new in submodel_parts if
                                 submodel_part_new.name not in all_names]
        submodel_parts_all.extend(unique_submodel_parts)


def update_project_parameters_files(all_project_parameters):
    """

    :param all_project_parameter:
    :return:
    """
    previous_end_time = 0
    for project_nr, project_parameters in enumerate(all_project_parameters):
        project_parameters["solver_settings"]["material_import_settings"][
            "materials_filename"] = 'MaterialParameters_stage' + str(project_nr + 1) + '.json'
        if project_parameters["problem_data"]["start_time"] < previous_end_time:
            dt = project_parameters["problem_data"]["end_time"] - project_parameters["problem_data"]["start_time"]
            project_parameters["problem_data"]["start_time"] = previous_end_time
            project_parameters["solver_settings"]["start_time"] = previous_end_time
            project_parameters["problem_data"]["end_time"] = previous_end_time + dt

        previous_end_time = project_parameters["problem_data"]["end_time"]


def write_all_files_in_project_dir(main_folder,new_project_name,all_project_names,all_project_parameters,
                                   all_material_parameters, all_mdpa_data ):
    import os
    import json

    new_dir_name = os.path.join(main_folder,new_project_name)
    if not os.path.exists(new_dir_name):
        os.makedirs(new_dir_name)

    for stage, project_name in enumerate(all_project_names):

        new_project_parameters_file = os.path.join(new_dir_name,'ProjectParameters_stage' + str(stage+1) + '.json')
        new_material_parameters_file = os.path.join(new_dir_name,'MaterialParameters_stage' + str(stage+1) + '.json')
        new_mdpa_file = os.path.join(new_dir_name, project_name.replace(".gid", ".mdpa"))

        with open(new_project_parameters_file,'w') as f:
            json.dump(all_project_parameters[stage],f,indent=4)

        with open(new_material_parameters_file, 'w') as f:
            json.dump(all_material_parameters[stage], f,indent=4)

        with open(new_mdpa_file, 'w') as f:
            f.writelines(all_mdpa_data[stage])


def run_main(main_folder,all_project_names,new_project_name):
    """

    :param main_folder:         overlaying folder of all stages
    :param all_project_names:   all gid project names
    :param new_project_name:    new project name
    :return:
    """

    all_project_parameters, all_material_parameters, all_mdpa_data = read_all_projects(all_project_names)

    clear_duplications_in_mdpa_files(all_project_names, all_mdpa_data)
    update_project_parameters_files(all_project_parameters)

    write_all_files_in_project_dir(main_folder, new_project_name, all_project_names, all_project_parameters,
                                   all_material_parameters, all_mdpa_data)


if __name__ == '__main__':
    import json
    import os

    main_folder = r"C:\Users\noordam\Documenten\Kratos\applications\GeoMechanicsApplication\test_examples"

    all_project_names = [r"simple_dike_test_stage1.gid", r"simple_dike_test_stage2.gid"] # all project names in correct order of stages
    new_project_name = r'simple_dike_test_all_stages'

    run_main(main_folder, all_project_names, new_project_name)

