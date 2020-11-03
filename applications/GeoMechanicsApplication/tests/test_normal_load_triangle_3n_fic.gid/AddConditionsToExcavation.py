

def write_conditions_to_excavation(filename,excavation_list):
    """

    :param filename:
    :param excavation_list:
    :return:
    """
    if len(excavation_list)>0:
        with open(filename,'r') as f:
            lines = f.readlines()

        conditions=[]
        i = 0
        while i < (len(lines)):
            # find all the conditions in the input file
            if "Begin Conditions" in lines[i]:
                i += 1
                while "End Conditions" not in lines[i] and i < len(lines):
                    conditions.append(list(map(int, lines[i].split())))
                    i += 1

            # find all the excavations in the input file
            for j in range(len(excavation_list)):
                if excavation_list[j] in lines[i]:
                    i += 2

                    excavation_nodes = []
                    excavation_conditions =[]

                    # find the nodes corresponding with the excavation
                    while "End SubModelPartNodes" not in lines[i] and i < len(lines):
                        excavation_nodes.append(lines[i])
                        i += 1

                    # find the condition which is assigned to the same nodes as the the excavation nodes
                    excavation_nodes = list(map(int,excavation_nodes))
                    for condition in conditions:
                        if condition[2] in excavation_nodes and condition[3] in excavation_nodes:
                            excavation_conditions.append(str(condition[0])+'\n')

                    # add the conditions corresponding to the excavation to the file
                    while "End SubModelPartConditions" not in lines[i] and i < len(lines):
                        i += 1
                    lines[i:i] = excavation_conditions
            i += 1

        with open(filename, 'w') as f:
            f.writelines(lines)

if __name__ == '__main__':
    import json
    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = json.load(parameter_file)
	
    input_filename = parameters["solver_settings"]["model_import_settings"]["input_filename"]
    input_type = parameters["solver_settings"]["model_import_settings"]["input_type"]
    constraint_processes = parameters["processes"]["constraints_process_list"]

    excavation_list = []
    for process in constraint_processes:
        settings = process["Parameters"]
        if "variable_name" in settings:
            if settings["variable_name"] == "EXCAVATION":
                excavation_list.append(settings["model_part_name"].split('.')[-1])

    write_conditions_to_excavation(input_filename+"."+input_type,excavation_list)