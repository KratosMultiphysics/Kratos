import sys,os

sys.path.append(os.path.join('..','..','..'))

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

sys.path.append(os.path.join('..', 'python_scripts'))
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis


def get_file_path(fileName):
    import os
    return os.path.dirname(__file__) + "/" + fileName


def run_kratos(file_path, model=None):
    """
    Runs 1 stage in kratos
    :param file_path:
    :return:
    """
    currentWorking = os.getcwd()

    parameter_file_name = os.path.join(file_path, 'ProjectParameters.json')
    os.chdir(file_path)

    with open(parameter_file_name, 'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    if model is None:
        model = Kratos.Model()
    simulation = analysis.GeoMechanicsAnalysis(model, parameters)
    simulation.Run()

    os.chdir(currentWorking)
    return simulation

def run_stages(project_path,n_stages):
    """
    Run all construction stages

    :param project_path:
    :param n_stages:
    :return:
    """
    currentWorking = os.getcwd()
    stages = get_stages(project_path,n_stages)
    [stage.Run() for stage in stages]
    os.chdir(currentWorking)
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


def get_nodal_variable(simulation, variable):
    """
    Gets values of a give nodal variable from kratos simulation
    :param simulation:
    :return values of a variable:
    """

    nodes = simulation._list_of_output_processes[0].model_part.Nodes
    values = [node.GetSolutionStepValue(variable) for node in nodes]

    return values


def get_nodal_variable_from_ascii(filename: str, variable: str):
    """
    Reads data of Kratos variable from ascii output GID file

   :param filename: ascii output file
   :param variable: variable name in GID output

    :return values of a variable:
    """

    # read data
    with open(filename, "r") as f:
        all_data = f.readlines()

    add_var = False

    data = []
    time_steps = []
    all_var_data = []

    # read all data at each time step of variable
    for line in all_data:

        if "End Values" in line and add_var:
            add_var = False
            all_var_data.append(data)
            data = []
        if add_var:
            data.append(line)

        if r'"' + variable + r'"' in line:
            time_step = float(line.split()[3])

            time_steps.append(time_step)
            add_var=True

    # initialise results dictionary
    res = {"time": time_steps}

    for var_data in all_var_data:
        var_data.pop(0)

        # convert var data to floats
        for i, _ in enumerate(var_data):
            line = var_data[i].split()
            line[1] = float(line[1])
            line[2] = float(line[2])
            line[3] = float(line[3])
            var_data[i] = line

    # add node numbers as dict keys
    for line in var_data:
        res[line[0]] = []

    for var_data in all_var_data:
        for line in var_data:
            res[line[0]].append(line[1:])

    return res


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
    import math
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
