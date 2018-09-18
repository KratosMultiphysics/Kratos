"""This module contains the global finite differencing response functions that can be used for any existing response function"""
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
import KratosMultiphysics
import structural_mechanics_analysis
import structural_response
import structural_response_function_factory

# ==============================================================================
class GlobalFiniteDifferencingResponseFunction(structural_response.ResponseFunctionBase):

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings
        self.sub_response_settings = response_settings["sub_response_settings"]
        self.model = model

        input_type = response_settings["model_import_settings"]["input_type"].GetString()
        model_part_name = response_settings["model_import_settings"]["input_filename"].GetString()
        model_part_name = model_part_name.split("/")[-1] # to extract the name from the path
        if input_type == "mdpa":
            self.model_part = KratosMultiphysics.ModelPart(model_part_name)
            self.model.AddModelPart(self.model_part)
        elif input_type == "use_input_model_part":
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            raise Exception("Other model part input options are not yet implemented.")

        # this one works with the model/model part of this class
        self.unperturbed_response_function = structural_response_function_factory.CreateResponseFunction(identifier+"_unperturbed", self.sub_response_settings, self.model)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)

    def Initialize(self):
        self.unperturbed_response_function.Initialize()

    def InitializeSolutionStep(self):
        self.unperturbed_response_function.InitializeSolutionStep()

    def CalculateValue(self):
        self.unperturbed_response_function.CalculateValue()

    def CalculateGradient(self):
        """
        For each sensitivity variable:
        - Create a new temporary response with a new model part and perturb it.
        - Calculate the value of the response and use it to calculate the finite difference.
        - Store the sensitivity value in the main model part of this class
        """

        # nodal_sensitivity_variables
        for i in range(self.response_settings["nodal_sensitivity_variables"].size()):
            nodal_variable = self.response_settings["nodal_sensitivity_variables"][i].GetString()
            if nodal_variable == "SHAPE":
                self._CalculateShapeSensitivities()
            else:
                raise NotImplementedError("FD for {} nodal_sensitivity_variables not implemented!".format(nodal_variable.GetString()) )

        # element_sensitivity_variables
        for i in range(self.response_settings["element_sensitivity_variables"].size()):
            element_variable = self.response_settings["element_sensitivity_variables"][i].GetString()
            raise NotImplementedError("FD for element_sensitivity_variables '{}' not implemented!".format(element_variable))

        # condition_sensitivity_variables
        for i in range(self.response_settings["condition_sensitivity_variables"].size()):
            condition_variable = self.response_settings["condition_sensitivity_variables"][i].GetString()
            raise NotImplementedError("FD for condition_sensitivity_variables '{}' not implemented!".format(condition_variable))

    def _CalculateShapeSensitivities(self):
        if self.response_settings["sensitivity_model_part_name"].GetString() == self.model_part.Name:
            sensitivity_model_part = self.model_part
        elif self.response_settings["sensitivity_model_part_name"].GetString() == "":
            sensitivity_model_part = self.model_part
        else:
            sensitivity_model_part = self.model_part.GetSubModelPart(self.response_settings["sensitivity_model_part_name"].GetString())
        perturbation = self.response_settings["step_size"].GetDouble()
        unperturbed_value = self.GetValue()

        for node in sensitivity_model_part.Nodes:
            for direction in range(3):
                perturbed_model = KratosMultiphysics.Model()
                identifier = "Node_{}_direction_{}".format(node.Id, direction)
                perturbed_response = self._CreateNewTmpResponse(identifier, perturbed_model)
                perturbed_response.Initialize()
                # perturb
                model_part = perturbed_response.primal_model_part #TODO find a better way to get the model part
                tmp_node = model_part.Nodes[node.Id]
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


    def _CreateNewTmpResponse(self, identifier, model):
        return structural_response_function_factory.CreateResponseFunction(identifier, self.sub_response_settings, model)

    def FinalizeSolutionStep(self):
        self.unperturbed_response_function.FinalizeSolutionStep()

    def Finalize(self):
        self.unperturbed_response_function.Finalize()

    def GetValue(self):
        return self.unperturbed_response_function.GetValue()

    def GetShapeGradient(self):
        gradient = {}
        for node in self.model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY)
        return gradient
