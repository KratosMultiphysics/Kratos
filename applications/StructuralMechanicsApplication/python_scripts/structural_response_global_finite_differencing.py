"""This module contains the global finite differencing response functions that can be used for any existing response function"""
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
import KratosMultiphysics
import structural_response_function_factory

def _CreateZeroOrderSettings(settings):
    zero_order_settings = KratosMultiphysics.Parameters(settings)
    zero_order_settings.RemoveValue("gradient_settings")
    return zero_order_settings

def _CalculateShapeSensitivities(response_function):
    main_response_settings = response_function.response_settings
    main_gradient_settings = main_response_settings["gradient_settings"]
    main_model_part = response_function.model_part

    if main_gradient_settings["sensitivity_model_part_name"].GetString() == main_model_part.Name:
        sensitivity_model_part = main_model_part
    elif main_gradient_settings["sensitivity_model_part_name"].GetString() == "":
        sensitivity_model_part = main_model_part
    else:
        sensitivity_model_part = main_model_part.GetSubModelPart(main_gradient_settings["sensitivity_model_part_name"].GetString())
    perturbation = main_gradient_settings["step_size"].GetDouble()

    unperturbed_value = response_function.GetValue()

    zero_order_settings = _CreateZeroOrderSettings(main_response_settings)
    for node in sensitivity_model_part.Nodes:
        for direction in range(3):

            perturbed_model = KratosMultiphysics.Model()
            identifier = "Node_{}_direction_{}".format(node.Id, direction)
            perturbed_response = structural_response_function_factory.CreateResponseFunction(identifier, zero_order_settings, perturbed_model)

            perturbed_response.Initialize()

            # perturb
            perturbed_model_part = perturbed_response.model_part
            tmp_node = perturbed_model_part.Nodes[node.Id]
            if direction == 0:
                tmp_node.X += perturbation
                tmp_node.X0 += perturbation
            if direction == 1:
                tmp_node.Y += perturbation
                tmp_node.Y0 += perturbation
            if direction == 2:
                tmp_node.Z += perturbation
                tmp_node.Z0 += perturbation

            perturbed_response.InitializeSolutionStep()
            perturbed_response.CalculateValue()
            perturbed_value = perturbed_response.GetValue()

            sensitivity = node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY)
            sensitivity[direction] = (perturbed_value - unperturbed_value) / perturbation
            node.SetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY, sensitivity)

            perturbed_response.FinalizeSolutionStep()
            perturbed_response.Finalize()

def CalculateResponseGradientWithFiniteDifferencing(response_function):
    """
    For each sensitivity variable:
    - Create a new temporary response with a new model part and perturb it.
    - Calculate the value of the response and use it to calculate the finite difference.
    - Store the sensitivity value in the main model part of this class
    """
    response_settings = response_function.response_settings
    gradient_settings = response_settings["gradient_settings"]
    # nodal_sensitivity_variables
    for i in range(gradient_settings["nodal_sensitivity_variables"].size()):
        nodal_variable = gradient_settings["nodal_sensitivity_variables"][i].GetString()
        if nodal_variable == "SHAPE":
            _CalculateShapeSensitivities(response_function)
        else:
            raise NotImplementedError("FD for {} nodal_sensitivity_variables not implemented!".format(nodal_variable.GetString()) )

    # element_sensitivity_variables
    for i in range(gradient_settings["element_sensitivity_variables"].size()):
        element_variable = gradient_settings["element_sensitivity_variables"][i].GetString()
        raise NotImplementedError("FD for element_sensitivity_variables '{}' not implemented!".format(element_variable))

    # condition_sensitivity_variables
    for i in range(gradient_settings["condition_sensitivity_variables"].size()):
        condition_variable = gradient_settings["condition_sensitivity_variables"][i].GetString()
        raise NotImplementedError("FD for condition_sensitivity_variables '{}' not implemented!".format(condition_variable))
