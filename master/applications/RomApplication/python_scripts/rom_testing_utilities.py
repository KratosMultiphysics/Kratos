import importlib

import KratosMultiphysics
import KratosMultiphysics.RomApplication.rom_analysis

def SetUpSimulationInstance(model, parameters, nn_rom_interface=None):
    """ Creates and returns a ROM simulation instance """

    # Get the parent simulation class
    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

    # Set up simulation
    model = KratosMultiphysics.Model()
    instance_factory = KratosMultiphysics.RomApplication.rom_analysis.CreateRomAnalysisInstance

    simulation = instance_factory(
        analysis_stage_class,
        model,
        parameters,
        nn_rom_interface=nn_rom_interface)

    return simulation

def GetNodalResults(model_part, variables_list):
    # Set and return an array containing the values of the variables in the variable list
    # Note that the array type variables need to be specified componentwise
    results_array = []
    for node in model_part.Nodes:
        for variable in variables_list:
            results_array.append(node.GetSolutionStepValue(variable))

    return results_array

def GetScalarNodalResults(model_part, variable):
    # Set and return an array containing the scalar variable values
    results_array = []
    for node in model_part.Nodes:
        results_array.append(node.GetSolutionStepValue(variable))

    return results_array

def GetVectorNodalResults(model_part, variable):
    # Set and return an array containing the vector variable values
    # Note that only the meaningful ones (2 in 2D and 3 in 3D) are returned
    results_array = []
    dim = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
    for node in model_part.Nodes:
        vector_value = node.GetSolutionStepValue(variable)
        for d in range(dim):
            results_array.append(vector_value[d])

    return results_array

def GetNodalAreaVector(model_part):
    # Calculate the NODAL_AREA and save it in the non-historical database
    nodal_area_calculator = KratosMultiphysics.CalculateNonHistoricalNodalAreaProcess(model_part)
    nodal_area_calculator.Execute()

    # Set and return an array containing the NODAL_AREA values
    nodal_area_array = []
    for node in model_part.Nodes:
        nodal_area_array.append(node.GetValue(KratosMultiphysics.NODAL_AREA))

    return nodal_area_array