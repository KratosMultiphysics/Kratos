from typing import Dict, Any
import sys,os
import math

sys.path.append(os.path.join('..','..','..'))

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

sys.path.append(os.path.join('..', 'python_scripts'))
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

def get_file_path(fileName):
    import os
    return os.path.join(os.path.dirname(__file__), fileName)

def assert_relative_close(test_case, a, b, rel_tol, abs_tol, msg):
    # Calculate the percentage difference
    difference = abs(a - b)
    average = (a + b) / 2
    # Check if the difference is less than the tolerance percentage of the average
    if average != 0:
        test_case.assertTrue((difference / average) < rel_tol if average != 0 else difference < rel_tol, f"{msg} - REL Error {rel_tol}, Expected: {a}, Actual: {b}, error: {difference / average}")
    test_case.assertTrue(difference < abs_tol, f"{msg} - ABS Error {abs_tol}, Expected: {a}, Actual: {b}, error: {difference}")

def compare_stress_results_with_plaxis_table_file(test_case, simulation, plaxis_table_file, int_points, rel_tol=1e-6, abs_tol=1e-6):
    """
    Compares the stress results with the results from a plaxis table file
    :param plaxis_table_file (3 integration points for element, i.e. 6 noded element)
    :return:
    """
    cauchy_stresses = get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_TENSOR)

    # read data
    with open(plaxis_table_file, "r") as f:
        all_data = f.readlines()[1:]

    element = 0
    if int_points == 3:
        map_plx_int_to_kratos_int = [2, 0, 1]
    else:
        raise NotImplementedError("Only 3 integration points per element are supported currently")

    for i in range(0, len(all_data), int_points):
        element_data = all_data[i:i+int_points]
        sig_ele = cauchy_stresses[element]
        for plx_int_no in range(0, int_points):
            kratos_int_no = map_plx_int_to_kratos_int[plx_int_no]
            integration_pt_data = element_data[plx_int_no].split("\t")
            stress_xx_plx = float(integration_pt_data[5]) * 1000 # in kN/m^2
            stress_xx_kratos = sig_ele[kratos_int_no][0,0]
            submsg = f"element: {element+1}, plx int point: {plx_int_no+1}, kratos int point: {kratos_int_no+1}"
            assert_relative_close(test_case, stress_xx_plx, stress_xx_kratos, rel_tol, abs_tol, f"Stress_xx - {submsg}")
            stress_yy_plx = float(integration_pt_data[6]) * 1000 # in kN/m^2
            stress_yy_kratos = sig_ele[kratos_int_no][1,1]
            assert_relative_close(test_case, stress_yy_plx, stress_yy_kratos, rel_tol, abs_tol, f"Stress_yy - {submsg}")
            stress_xy_plx = float(integration_pt_data[8]) * 1000 # in kN/m^2
            stress_xy_kratos = sig_ele[kratos_int_no][0,1]
            assert_relative_close(test_case, stress_xy_plx, stress_xy_kratos, rel_tol, abs_tol, f"Stress_xy - {submsg}")
        element += 1

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

def run_stages(project_path,n_stages):
    """
    Run all construction stages

    :param project_path:
    :param n_stages:
    :return:
    """
    cwd = os.getcwd()
    stages = get_stages(project_path,n_stages)
    [stage.Run() for stage in stages]
    os.chdir(cwd)
    return stages

def get_stages(project_path,n_stages):
    """
    Gets all construction stages

    :param project_path:
    :param n_stages:
    :return:
    """

    parameter_file_names = [os.path.join(project_path, 'ProjectParameters_stage' + str(i + 1) + '.json') for i in
                            range(n_stages)]

    # set stage parameters
    parameters_stages = [None] * n_stages
    os.chdir(project_path)
    for idx, parameter_file_name in enumerate(parameter_file_names):
        with open(parameter_file_name, 'r') as parameter_file:
            parameters_stages[idx] = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    stages = [analysis.GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]
    return stages

def get_separated_directory_names(project_path, n_stages):
    """
    Gets directory names for all construction stages in seperated directories as Stage_0, Stage_1, ...

    :param project_path:
    :param n_stages:
    :return:
    """
    directory_names = [os.path.join(project_path, 'Stage_' +  str(i + 1)) for i in range(n_stages)]

    return directory_names

def get_separated_stages(directory_names):
    """
    Gets all construction stages in seperated directories as Stage_0, Stage_1, ...

    :param project_path:
    :param n_stages:
    :return:
    """
    n_stages = len(directory_names)
    # set stage parameters
    parameters_stages = [None] * n_stages
    for idx, directory_name in enumerate(directory_names):
        parameter_file_name = directory_name + '/ProjectParameters.json'
        with open(parameter_file_name, 'r') as parameter_file:
            parameters_stages[idx] = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    stages = [analysis.GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

    return stages


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

    Force = [element.CalculateOnIntegrationPoints(
        Kratos.FORCE, model_part.ProcessInfo)[0] for element in elements]
    return Force


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


class GiDOutputFileReader:
    def __init__(self):
        self._reset_internal_state()

    def read_output_from(self, gid_output_file_path):
        self._reset_internal_state()

        with open(gid_output_file_path, "r") as result_file:
            for line in result_file:
                line = line.strip()
                if line.startswith("GaussPoints"):
                    self._process_begin_of_gauss_points(line)
                elif line == "End GaussPoints":
                    self._process_end_of_gauss_points(line)
                elif line.startswith("Result"):
                    self._process_result_header(line)
                elif line == "Values":
                    self._process_begin_of_block(line)
                elif line == "End Values":
                    self._process_end_of_block(line)
                elif self.current_block_name == "GaussPoints":
                    self._process_gauss_point_data(line)
                elif self.current_block_name == "Values":
                    self._process_value_data(line)

        return self.output_data

    def _reset_internal_state(self):
        self.output_data = {}
        self.current_block_name = None
        self.result_name = None
        self.result_type = None
        self.result_location = None
        self.gauss_points_name = None
        self.current_integration_point = 0

    def _process_begin_of_gauss_points(self, line):
        self._process_begin_of_block(line)
        self.gauss_points_name = self._strip_off_quotes(line.split()[1])
        if self.current_block_name not in self.output_data:
            self.output_data[self.current_block_name] = {}
        self.output_data[self.current_block_name][self.gauss_points_name] = {}

    def _process_end_of_gauss_points(self, line):
        self._process_end_of_block(line)
        self.gauss_points_name = None

    def _process_result_header(self, line):
        if "results" not in self.output_data:
            self.output_data["results"] = {}
        words = line.split()
        self.result_name = self._strip_off_quotes(words[1])
        if self.result_name not in self.output_data["results"]:
            self.output_data["results"][self.result_name] = []
        self.result_type = words[4]
        self.result_location = words[5]
        this_result = {"time": float(words[3]),
                       "location": self.result_location,
                       "values": []}
        self.output_data["results"][self.result_name].append(this_result)
        if self.result_location == "OnGaussPoints":
            self.current_integration_point = 0
            self.gauss_points_name = self._strip_off_quotes(words[6])

    def _process_gauss_point_data(self, line):
        if line.startswith("Number Of Gauss Points:"):
            pos = line.index(":")
            num_gauss_points = int(line[pos+1:].strip())
            self.output_data[self.current_block_name][self.gauss_points_name]["size"] = num_gauss_points

    def _process_value_data(self, line):
        words = line.split()
        if self.result_location == "OnNodes":
            self._process_nodal_result(words)
        elif self.result_location == "OnGaussPoints":
            self._process_gauss_point_result(words)

    def _process_nodal_result(self, words):
        value = {"node": int(words[0])}
        if self.result_type == "Scalar":
            value["value"] = float(words[1])
        elif self.result_type == "Vector" or self.result_type == "Matrix":
            value["value"] = [float(x) for x in words[1:]]
        self.output_data["results"][self.result_name][-1]["values"].append(value)

    def _process_gauss_point_result(self, words):
        self.current_integration_point %= self.output_data["GaussPoints"][self.gauss_points_name]["size"]
        self.current_integration_point += 1
        if self.current_integration_point == 1:
            value = {"element": int(words[0]),
                     "value": []}
            self.output_data["results"][self.result_name][-1]["values"].append(value)
            words.pop(0)

        value = self.output_data["results"][self.result_name][-1]["values"][-1]["value"]
        if self.result_type == "Scalar":
            value.append(float(words[0]))
        elif self.result_type == "Matrix" or self.result_type == "Vector":
            value.append([float(x) for x in words])
        else:
            raise RuntimeError(f'Unsupported result type "{self.result_type}"')

    def _process_begin_of_block(self, line):
        assert(self.current_block_name is None)  # nested blocks are not supported
        self.current_block_name = line.split()[0]

    def _process_end_of_block(self, line):
        words = line.split()
        assert(words[0] == "End")
        assert(self.current_block_name == words[1])
        self.current_block_name = None

    def _strip_off_quotes(self, quoted_string):
        assert(quoted_string[0] == '"')
        assert(quoted_string[-1] == '"')
        return quoted_string[1:-1]

    @staticmethod
    def get_values_at_time(time, property_results):
        for time_results in property_results:
            if math.isclose(time_results["time"], time):
                return time_results["values"]

    @staticmethod
    def get_value_at_node(node, time_results):
        for node_results in time_results:
            if node_results["node"] == node:
                return node_results["value"]

    @staticmethod
    def nodal_values_at_time(result_item_name, time, output_data, node_ids=None):
        if node_ids and node_ids != sorted(node_ids):
            raise RuntimeError("Node IDs must be sorted")

        matching_item = None
        for item in output_data["results"][result_item_name]:
            if math.isclose(item["time"], time):
                matching_item = item
                break
        if matching_item is None:
            raise RuntimeError(f"'{result_item_name}' does not have results at time {time}")

        if matching_item["location"] != "OnNodes":
            raise RuntimeError(f"'{result_item_name}' is not a nodal result")

        if not node_ids: # return all values
            return [item["value"] for item in matching_item["values"]]

        return [item["value"] for item in matching_item["values"] if item["node"] in node_ids]

    @staticmethod
    def element_integration_point_values_at_time(result_item_name, time, output_data, element_ids=None, integration_point_indices=None):
        if element_ids and element_ids != sorted(element_ids):
            raise RuntimeError("Element IDs must be sorted")

        matching_item = None
        for item in output_data["results"][result_item_name]:
            if math.isclose(item["time"], time):
                matching_item = item
                break
        if matching_item is None:
            raise RuntimeError(f"'{result_item_name}' does not have results at time {time}")

        if matching_item["location"] != "OnGaussPoints":
            raise RuntimeError(f"'{result_item_name}' is not an integration point result")

        if element_ids:
            element_results = [item["value"] for item in matching_item["values"] if item["element"] in element_ids]
        else:
            element_results = [item["value"] for item in matching_item["values"]]

        if integration_point_indices:
            result = []
            for element_result in element_results:
                result.append([item for index, item in enumerate(element_result) if index in integration_point_indices])
            return result
        else:
            return element_results
