from typing import Dict, Any
import sys,os
import math

sys.path.append(os.path.join('..','..','..'))

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

sys.path.append(os.path.join('..', 'python_scripts'))
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

from KratosMultiphysics.GeoMechanicsApplication import unit_conversions


def get_file_path(filename):
    import os
    return os.path.join(os.path.dirname(__file__), filename)


def make_geomechanics_analysis(model, project_parameters_file_path):
    with open(project_parameters_file_path, 'r') as f:
        project_parameters = Kratos.Parameters(f.read())

    return analysis.GeoMechanicsAnalysis(model, project_parameters)


def run_kratos(file_path, model=None):
    """
    Runs 1 stage in kratos
    :param file_path:
    :return:
    """
    cwd = os.getcwd()

    parameter_file_name = os.path.join(file_path, 'ProjectParameters.json')
    os.chdir(file_path)

    with open(parameter_file_name, 'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    if model is None:
        model = Kratos.Model()
    simulation = analysis.GeoMechanicsAnalysis(model, parameters)
    simulation.Run()

    os.chdir(cwd)
    return simulation


def get_displacement(simulation):
    """
    Gets displacements from kratos simulation
    :param simulation:
    :return displacements:
    """

    return get_nodal_variable(simulation, Kratos.DISPLACEMENT)


def get_velocity(simulation):
    """
    Gets velocities from kratos simulation
    :param simulation:
    :return velocities:
    """

    return get_nodal_variable(simulation, Kratos.VELOCITY)

def get_temperature(simulation):
    """
    Gets the temperature from kratos simulation
 
    :param simulation:
    :return:
    """
    return get_nodal_variable(simulation, Kratos.TEMPERATURE)

def get_water_pressure(simulation):
    """
    Gets the water pressure from kratos simulation

    :param simulation:
    :return:
    """
    return get_nodal_variable(simulation, Kratos.WATER_PRESSURE)


def get_hydraulic_discharge(simulation):
    """
    Gets displacements from kratos simulation
    :param simulation:
    :return hydraulic discharge:
    """

    return get_nodal_variable(simulation, KratosGeo.HYDRAULIC_DISCHARGE)


def get_nodal_variable(simulation, variable, node_ids=None):
    """
    Gets values of a give nodal variable from kratos simulation
    :param simulation:
    :return values of a variable:
    """

    nodes = simulation._list_of_output_processes[0].model_part.Nodes
    if node_ids:
        nodes = [node for node in nodes if node.Id in node_ids]
    return [node.GetSolutionStepValue(variable) for node in nodes]

def get_nodal_variable_from_ascii(filename: str, variable: str):
    """
    Reads data of Kratos variable from ascii output GID file

   :param filename: ascii output file
   :param variable: variable name in GID output

    :return values of a variable:
    """

    # read data
    with open(filename, "r") as f:
        all_lines = f.readlines()

    time_step = None
    add_var = False
    res = {}

    # read all data at each time step of variable
    for line in all_lines:

        if "End Values" in line and add_var:
            add_var = False

        if add_var:
            if line.startswith("Values"):
                continue
            if time_step is not None:
                add_line_data_to_dictionary(line, res, time_step)

        if r'"' + variable + r'"' in line:
            time_step = float(line.split()[3])
            res[time_step] = {}
            add_var=True

    return res

def add_line_data_to_dictionary(line, dictionary, main_index):
    """
    Adds the data from a GiD Ascii line to a dictionary structure
    :param line: line with data from a GiD Ascii File
    :param dictionary: dictionary to input data (main index already initialized)
    :param main_index: this is the main index under which to store the data
                       (usually time series), i.e. dictionary[main_index] = {}
    """
    if main_index not in dictionary.keys():
        raise KeyError(f"The key '{main_index}' is not in the dictionary.")
    if not isinstance(dictionary[main_index], dict):
        raise TypeError(f"The value for key '{main_index}' is not a dictionary.")
    line_split = line.split()
    line_split[0] = int(line_split[0])
    for ind, str_value in enumerate(line_split[1:]):
        line_split[ind + 1] = float(str_value)
    if (len(line_split[1:]) == 1):
        dictionary[main_index][line_split[0]] = line_split[1]
    else:
        dictionary[main_index][line_split[0]] = line_split[1:]


def get_gauss_coordinates(simulation):
    """
    Gets the coordinates of the gauss points
    :param simulation:
    :return:
    """
    elements = simulation._list_of_output_processes[0].model_part.Elements
    gauss_coordinates = [element.GetIntegrationPoints() for element in elements]
    return gauss_coordinates


def get_nodal_coordinates(simulation):
    """
    Gets original coordinates of all nodes
    :param simulation:
    :return:
    """
    nodes = simulation._list_of_output_processes[0].model_part.Nodes
    coordinates = [[node.X0,node.Y0,node.Z0] for node in nodes]

    return coordinates

def get_green_lagrange_strain_tensor(simulation):
    """
    Gets green lagrange strain tensor from kratos simulation
    :param simulation:
    :return: green lagrange strain tensor
    """
    model_part= simulation._list_of_output_processes[0].model_part
    elements = model_part.Elements
    green_lagrange_strain_tensors = [element.CalculateOnIntegrationPoints(
        Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, model_part.ProcessInfo) for element in elements]

    return green_lagrange_strain_tensors


def get_cauchy_stress_tensor(simulation):
    """
    Gets cauchy stress tensor from kratos simulation
    :param simulation:
    :return: cauchy stress tensor
    """
    model_part = simulation._list_of_output_processes[0].model_part
    elements = model_part.Elements
    cauchy_stress_tensors = [element.CalculateOnIntegrationPoints(
        Kratos.CAUCHY_STRESS_TENSOR, model_part.ProcessInfo) for element in elements]

    return cauchy_stress_tensors

def get_total_stress_tensor(simulation):
    """
    Gets total stress tensor from kratos simulation
    :param simulation:
    :return: total stress tensor
    """
    model_part = simulation._list_of_output_processes[0].model_part
    elements = model_part.Elements
    total_stress_tensors = [element.CalculateOnIntegrationPoints(
        KratosGeo.TOTAL_STRESS_TENSOR, model_part.ProcessInfo) for element in elements]

    return total_stress_tensors


def get_on_integration_points(simulation, kratos_variable):
    """
    Gets the values of a Kratos variables on all integration points within a model part

    :param simulation:
    :return: local stress vector
    """
    model_part = simulation._list_of_output_processes[0].model_part
    elements = model_part.Elements
    results = [element.CalculateOnIntegrationPoints(
        kratos_variable, model_part.ProcessInfo) for element in elements]
    return results


def get_local_stress_vector(simulation):
    """
    Gets local stress vector on all integration points from Kratos simulation
    :param simulation:
    :return: local stress vector
    """

    model_part = simulation._list_of_output_processes[0].model_part
    elements = model_part.Elements

    local_stress_vector = [element.CalculateOnIntegrationPoints(
        KratosGeo.LOCAL_STRESS_VECTOR, model_part.ProcessInfo) for element in elements]
    return local_stress_vector

def get_hydraylic_head_with_intergration_points(simulation):
    """
    Gets hydraylic head on all nodal points from Kratos simulation
    :param simulation:
    :return: force
    """
    model_part = simulation._list_of_output_processes[0].model_part
    elements = model_part.Elements

    x = []
    y = []
    head = []
    for element in elements:
        values = element.CalculateOnIntegrationPoints(KratosGeo.HYDRAULIC_HEAD, model_part.ProcessInfo)
        points = element.GetIntegrationPoints()
        for counter, head_value in enumerate(values):
            x.append(points[counter][0])
            y.append(points[counter][1])
            head.append(head_value)
    return x, y , head

def get_pipe_active_in_elements(simulation):
    """
    Gets the pipe active value on all elements from Kratos simulation
    :param simulation:
    :return: pipe_active : list of booleans determine whether pipe element is active or not
    """
    model_part = simulation._list_of_output_processes[0].model_part
    pipe_elements = [element for element in model_part.Elements if element.Has(KratosGeo.PIPE_ELEMENT_LENGTH)]
    return [element.GetValue(KratosGeo.PIPE_ACTIVE) for element in pipe_elements]

def get_pipe_length(simulation):
    """
    Gets the length of all active pipe elemnets
    :param simulation:
    :return: pipe_length :
    """
    model_part = simulation._list_of_output_processes[0].model_part
    elements = model_part.Elements
    return sum([element.GetValue(KratosGeo.PIPE_ELEMENT_LENGTH) for element in elements if element.GetValue(KratosGeo.PIPE_ACTIVE)])

def get_force(simulation):
    """
    Gets force on all integration points from Kratos simulation
    :param simulation:
    :return: force
    """

    model_part = simulation._list_of_output_processes[0].model_part
    elements = model_part.Elements

    return [element.CalculateOnIntegrationPoints(
        Kratos.FORCE, model_part.ProcessInfo)[0] for element in elements]


def get_moment(simulation):
    """
    Gets bending moment on all integration points from Kratos simulation
    :param simulation:
    :return: bending moment
    """
    model_part = simulation._list_of_output_processes[0].model_part
    elements = model_part.Elements

    moment = [element.CalculateOnIntegrationPoints(
        Kratos.MOMENT, model_part.ProcessInfo)[0] for element in elements]
    return moment


def compute_distance(point1, point2):
    """
    Computes distance between 2 points in a 2D space
    :param point1:
    :param point2:
    :return: distance
    """
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)


def compute_mean_list(list):
    return sum(list)/len(list)


def find_closest_index_greater_than_value(input_list, value):
    """
    Finds closest value in list which is greater than the input value. This method assumes a sorted list from
    low to high

    :param input_list: sorted list
    :param value: value to be checked
    :return:

    """

    for index, list_value in enumerate(input_list):
        if value < list_value:
            return index
    return None


def are_values_almost_equal(expected: Any, actual: Any, abs_tolerance: float = 1e-7) -> bool:
    """
    Checks whether two values are almost equal.

    Args:
        - expected (Any): Expected value.
        - actual (Any): Actual value.

    Returns:
        - True if the values are almost equal, False otherwise.

    """
    # check if the value is a dictionary and check the dictionary
    if isinstance(expected, dict):
        return are_dictionaries_almost_equal(expected, actual)
    elif isinstance(expected, str):
        return expected == actual
    elif isinstance(expected, (list, tuple, set)):
        return are_iterables_almost_equal(expected, actual)
    elif expected is None:
        return actual is None
    elif isinstance(expected, (float, int, complex)):
        return math.isclose(expected, actual, abs_tol=abs_tolerance)
    else:
        raise TypeError(f"Unsupported type {type(expected)}")


def are_iterables_almost_equal(expected: (list, tuple, set), actual: (list, tuple, set),
                               abs_tolerance: float = 1e-7) -> bool:
    """
    Checks whether two iterables are almost equal.

    Args:
        - expected (list, tuple, set): Expected iterable.
        - actual (list, tuple, set): Actual iterable.

    Returns:
        - True if the iterables are almost equal, False otherwise.

    """
    # check if the value is a list, tuple or set and compare the values
    if len(expected) != len(actual):
        return False

    for v_i, actual_i in zip(expected, actual):
        if not are_values_almost_equal(v_i, actual_i, abs_tolerance):
            return False

    return True


def are_dictionaries_almost_equal(expected: Dict[Any, Any],
                                  actual: Dict[Any, Any],
                                  abs_tolerance: float = 1e-7) -> bool:
    """
    Checks whether two dictionaries are equal.

    Args:
        - expected: Expected dictionary.
        - actual: Actual dictionary.

    Returns:
        - True if the dictionaries are equal, False otherwise.

    """
    if len(expected) != len(actual):
        return False

    for k, v in expected.items():

        # check if key is present in both dictionaries
        if k not in actual:
            return False

        # check if values are almost equal
        if not are_values_almost_equal(v, actual[k], abs_tolerance):
            return False

    # all checks passed
    return True


def want_test_plots() -> bool:
    return os.environ.get("KRATOS_GEO_MAKE_TEST_PLOTS", "off").lower() == "on"


def get_data_points_from_file(file_path, data_point_extractor):
    result = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            result.append(data_point_extractor(line))

    return result


def read_coordinates_from_post_msh_file(file_path, node_ids=None):
    node_map = {}

    with open(file_path, "r") as post_msh_file:
        reading_coordinates = False
        for line in post_msh_file:
            line = line.strip()
            if line == "Coordinates":
                reading_coordinates = True
                continue

            if line == "End Coordinates":
                reading_coordinates = False

            if reading_coordinates:
                numbers = line.split()  # [node ID, x, y, z]
                node_map[int(numbers[0])] = tuple([float(number) for number in numbers[1:]])

    if node_ids is None:
        return list(node_map.values())

    result = []
    for id in node_ids:
        result.append(node_map[id])
    return result