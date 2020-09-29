def alter_beam_element_strings(lines, new_beam_elements, i):
    """
    Alter the strings in the .mdpa file concerning the beam elements
    :param lines:
    :param new_beam_elements:
    :param i:
    :return:
    """
    lines[i] = "Begin Elements GeoCrBeamElement2D2N\n"
    i += 1
    property_id = int(lines[i].split()[1])
    while "End Elements" not in lines[i]:
        del (lines[i])

    new_lines = ["  " + str(new_beam_element["ID"]) + "  " + str(new_beam_element["property_ID"]) + "  " +
                 " ".join([str(node) for node in new_beam_element["node_connectivities"]]) + "\n"
                 for new_beam_element in new_beam_elements if new_beam_element["property_ID"] == property_id]

    lines[i:i] = new_lines
    return i

def alter_sub_model_strings(lines, sub_model_parts, altered_sub_model_parts, i):
    """
    Alter the lines in the .mdpa file concerning the sub model parts
    :param lines:
    :param sub_model_parts:
    :param altered_sub_model_parts:
    :param i:
    :return:
    """
    sub_model_part_name = lines[i].split()[-1]
    i += 1

    while "End SubModelPart\n" not in lines[i]:
        if sub_model_part_name in altered_sub_model_parts:
            if "Begin SubModelPartElements" in lines[i]:
                i += 1
                while "End SubModelPartElements" not in lines[i]:
                    del (lines[i])

                new_lines = ["    " + str(element_ID) + '\n' for sub_model_part in sub_model_parts
                             for element_ID in sub_model_part['element_IDs']
                             if sub_model_part["name"] == sub_model_part_name]
                lines[i:i] = new_lines
        i += 1

    return i

def get_altered_sub_model_parts(sub_model_parts, new_beam_elements):
    """
    Finds the sub model parts which are influenced by the new beam elements
    :param sub_model_parts:
    :param new_beam_elements:
    :return:
    """
    altered_sub_model_parts = []
    for sub_model_part in sub_model_parts:
        new_element_IDs = [new_beam_elements[idx + 1]["ID"] for
                           idx, new_beam_element in enumerate(new_beam_elements) for element_ID in sub_model_part['element_IDs'] if
                           element_ID == new_beam_element["ID"]]

        if new_element_IDs:
            altered_sub_model_parts.append(sub_model_part["name"])
            sub_model_part['element_IDs'].extend(new_element_IDs)

    return altered_sub_model_parts

def write_mdpa_file(filename, lines, new_beam_elements, sub_model_parts, altered_sub_model_parts):
    """
    Writes the mdpa file with splitted beam elements.
    :param filename:
    :param lines:
    :param new_beam_elements:
    :param sub_model_parts:
    :param altered_sub_model_parts:
    :return:
    """
    i = 0
    while i < len(lines):
        if "GeoCrBeamElementLinear2D3N" in lines[i]:
            i = alter_beam_element_strings(lines, new_beam_elements, i)

        if "Begin SubModelPart " in lines[i]:
            i = alter_sub_model_strings(lines, sub_model_parts, altered_sub_model_parts, i)
        i += 1

    with open(filename, 'w') as f:
        f.writelines(lines)

def read_mdpa_file(filename):
    """
    Reads the .mdpa file
    :param filename:
    :return:
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    elements = []
    sub_model_parts = []
    while i < len(lines):
        if "Begin Elements" in lines[i]:
            i = read_elements(lines, elements, i)

        if "Begin SubModelPart " in lines[i]:
            i = read_sub_model_parts(lines, sub_model_parts, i)
        i += 1
    return elements, sub_model_parts, lines

def read_sub_model_parts(lines, sub_model_parts, i):
    """
    Reads all the sub model parts from the .mdpa file
    :param lines:
    :param sub_model_parts:
    :param i:
    :return:
    """
    sub_model_part_name = lines[i].split()[-1]
    i += 1
    sub_model_part_elements = []
    while "End SubModelPart\n" not in lines[i]:
        if "Begin SubModelPartElements" in lines[i]:
            i += 1
            while "End SubModelPartElements" not in lines[i]:
                sub_model_part_elements.append(int(lines[i]))
                i += 1
        i += 1

    sub_model_part = {"name": sub_model_part_name,
                      "element_IDs": sub_model_part_elements}

    sub_model_parts.append(sub_model_part)
    return i

def read_elements(lines, elements, i):
    """
    Reads all the elements from the .mdpa file
    :param lines:
    :param elements:
    :param i:
    :return:
    """
    element_type = lines[i].split()[-1]
    i += 1
    while "End Elements" not in lines[i]:
        id_s = lines[i].split()
        element = {"type": element_type,
                   "ID": int(id_s[0]),
                   "property_ID": int(id_s[1]),
                   "node_connectivities": [int(id_s[idx]) for idx in range(2, len(id_s))]}

        elements.append(element)
        i += 1
    return i

def define_new_elements(element,n_elements):
    """
    defines the new two-noded line elements, where the same nodes are used as the 3 noded line element
    :param element:
    :param n_elements:
    :return:
    """
    new_element_1 = {"type": "GeoCrBeamElementLinear2D2N",
                     "ID": element["ID"],
                     "property_ID": element["property_ID"],
                     "node_connectivities": [element["node_connectivities"][0], element["node_connectivities"][2]]}

    new_element_2 = {"type": "GeoCrBeamElementLinear2D2N",
                     "ID": n_elements + 1,
                     "property_ID": element["property_ID"],
                     "node_connectivities": [element["node_connectivities"][2], element["node_connectivities"][1]]}

    return [new_element_1, new_element_2]

def main(filename):
    """
    Reads mdpa file.
    Splits 3 noded beam elements in 2 noded beam elements.
    Writes mdpa file.
    :param filename:
    :return:
    """

    elements, sub_model_parts, lines = read_mdpa_file(filename)

    max_elem_n = max([element["ID"] for element in elements])

    beam_elements = [element for element in elements if element["type"] == "GeoCrBeamElementLinear2D3N"]

    new_beam_elements = [new_element for n_added, element in enumerate(beam_elements)
                         for new_element in define_new_elements(element, max_elem_n+n_added)]

    altered_sub_model_parts = get_altered_sub_model_parts(sub_model_parts, new_beam_elements)

    write_mdpa_file(filename, lines, new_beam_elements, sub_model_parts, altered_sub_model_parts)


if __name__ == '__main__':
    import json

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = json.load(parameter_file)

    input_filename = parameters["solver_settings"]["model_import_settings"]["input_filename"]
    input_type = parameters["solver_settings"]["model_import_settings"]["input_type"]

    main(input_filename + "." + input_type)