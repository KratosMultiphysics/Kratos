import sys,os

sys.path.append(os.path.join('..','..','..'))

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

sys.path.append(os.path.join('..', 'python_scripts'))
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

def get_file_path(fileName):
    import os
    return os.path.dirname(__file__) + "/" + fileName


def run_kratos(file_path):
    """
    Runs 1 stage in kratos
    :param file_path:
    :return:
    """

    parameter_file_name = os.path.join(file_path, 'ProjectParameters.json')
    os.chdir(file_path)

    with open(parameter_file_name, 'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = analysis.GeoMechanicsAnalysis(model, parameters)
    simulation.Run()
    return simulation

def run_stages(project_path,n_stages):
    """
    Run all construction stages

    :param project_path:
    :param n_stages:
    :return:
    """

    stages = get_stages(project_path,n_stages)
    [stage.Run() for stage in stages]

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


def get_displacement(simulation):
    """
    Gets displacements from kratos simulation
    :param simulation:
    :return displacements:
    """

    return get_nodal_variable(simulation, Kratos.DISPLACEMENT)

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

if __name__ == "__main__":
    file_path = r"C:\Users\noordam\Documenten\Kratos\applications\GeoMechanicsApplication\test_examples\simple_dike_test.gid"
    run_kratos(file_path)