from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.MeshingApplication as KratosMeshing
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
from KratosMultiphysics.MultilevelMonteCarloApplication.tools import ParametersWrapper

# Import packages
import numpy as np

"""
References:
F. Alauzet, Metric-based anisotropic mesh adaptation, CEA-EDF-INRIA schools: Numerical Analysis Summer School. CEA, Cadarache, France
Kratos wiki: https://github.com/KratosMultiphysics/Kratos/wiki/MMG-Process
"""


class AdaptiveRefinement(object):
    """
    input:  model_coarse       : Kratos model class before refinement
            parameters_coarse  : Kratos parameters class before refinement
            minimal_size_value : minimal size after remeshing
            maximal_size_value : maximal size after remeshing
            metric_param       : Kratos parameters class containing metric custom settings
            remesh_param       : Kratos parameters class containing remeshing custom settings
    """
    def __init__(self,current_level,model_coarse,parameters_coarse,metric_param,remesh_param,metric_name="hessian"):
        self.model_coarse = model_coarse
        self.parameters_coarse = parameters_coarse
        self.metric_param = metric_param
        self.remesh_param = remesh_param
        self.problem_type = self.parameters_coarse["solver_settings"]["solver_type"].GetString()
        self.metric = metric_name
        self.current_level = current_level
        self.wrapper = ParametersWrapper(self.parameters_coarse)

    """
    function computing the refinement of the model based on the solution on the coarse mesh,
    exploiting the hessian metric of the solution
    input:  self: an instance of the class
    output: current_model_refined      : Kratos model class after refinement
            current_parameters_refined : Kratos parameters class after refinement
    """
    def ComputeAdaptiveRefinement(self):
        parameters_coarse = self.parameters_coarse
        model_coarse = self.model_coarse
        metric_param = self.metric_param
        remesh_param = self.remesh_param
        problem_type = self.problem_type
        current_level = self.current_level

        if (self.metric is "hessian"):
            original_interp_error = metric_param["hessian_strategy_parameters"]["interpolation_error"].GetDouble()

            # problem dependent section
            if (problem_type == "potential_flow"):
                model_part_name = parameters_coarse["solver_settings"]["model_part_name"].GetString()
                # set NODAL_AREA and NODAL_H as non historical variables
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, model_coarse.GetModelPart(model_part_name).Nodes)
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, model_coarse.GetModelPart(model_part_name).Nodes)
                # Setting Metric Tensor to 0
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,model_coarse.GetModelPart(model_part_name).Nodes)
                # calculate NODAL_H
                find_nodal_h = KratosMultiphysics.FindNodalHProcess(model_coarse.GetModelPart(model_part_name))
                find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_coarse.GetModelPart(model_part_name))
                find_nodal_h.Execute()
                custom_gradient = KratosMultiphysics.CompressiblePotentialFlowApplication.ComputeCustomNodalGradientProcess(model_coarse.GetModelPart(model_part_name),KratosMultiphysics.VELOCITY, KratosMultiphysics.NODAL_AREA)
                custom_gradient.Execute()

                if metric_param.Has("local_gradient_variable"):
                    metric_param.RemoveValue("local_gradient_variable")
                if current_level > 0:
                    coefficient_interp_error =  metric_param["hessian_strategy_parameters"]["coefficient_interpolation_error"].GetDouble()
                    metric_param["hessian_strategy_parameters"].RemoveValue("coefficient_interpolation_error")
                    # interp_error = original_interp_error*(coefficient_interp_error)**(-current_level)
                    interp_error = original_interp_error/(coefficient_interp_error*current_level)
                    metric_param["hessian_strategy_parameters"]["interpolation_error"].SetDouble(interp_error)

                local_gradient = KratosMeshing.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),KratosMultiphysics.VELOCITY_X,metric_param)
                local_gradient.Execute()
                local_gradient = KratosMeshing.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),KratosMultiphysics.VELOCITY_Y,metric_param)
                local_gradient.Execute()

                # #### OLD APPROACH
                # for node in model_coarse.GetModelPart(model_part_name).Nodes:
                #     vector=node.GetValue(KratosMultiphysics.VELOCITY)
                #     norm=np.linalg.norm(vector)
                #     node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,norm)

                # prepare parameters to calculate the gradient of the designed variable
                # local_gradient_variable_string = metric_param["local_gradient_variable"].GetString()
                # local_gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable(metric_param["local_gradient_variable"].GetString())
                # set interpolation error value (level dependent)
                # calculate the gradient of the variable
                # local_gradient = KratosMeshing.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),local_gradient_variable,metric_param)
                # local_gradient.Execute()

                # # add again the removed variable parameter
                # metric_param.AddEmptyValue("local_gradient_variable")
                # metric_param["local_gradient_variable"].SetString(local_gradient_variable_string)

                #### OLD APPROACH


            elif (problem_type == "monolithic"):
                if metric_param.Has("local_gradient_variable"):
                    metric_param.RemoveValue("local_gradient_variable")
                if current_level > 0:
                    coefficient_interp_error =  metric_param["hessian_strategy_parameters"]["coefficient_interpolation_error"].GetDouble()
                    metric_param["hessian_strategy_parameters"].RemoveValue("coefficient_interpolation_error")
                    # interp_error = original_interp_error*(coefficient_interp_error)**(-current_level)
                    interp_error = original_interp_error/(coefficient_interp_error*current_level)
                    metric_param["hessian_strategy_parameters"]["interpolation_error"].SetDouble(interp_error)
                model_part_name = parameters_coarse["solver_settings"]["model_part_name"].GetString()

                # Setting Metric Tensor to 0
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,model_coarse.GetModelPart(model_part_name).Nodes)

                # calculate NODAL_H
                find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_coarse.GetModelPart(model_part_name))
                find_nodal_h.Execute()
                local_gradient = KratosMeshing.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),KratosFluid.AVERAGE_VELOCITY_X,metric_param)
                local_gradient.Execute()
                local_gradient = KratosMeshing.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),KratosFluid.AVERAGE_VELOCITY_Y,metric_param)
                local_gradient.Execute()

            elif (problem_type == "stationary"):
                if current_level > 0:
                    coefficient_interp_error =  metric_param["hessian_strategy_parameters"]["coefficient_interpolation_error"].GetDouble()
                    metric_param["hessian_strategy_parameters"].RemoveValue("coefficient_interpolation_error")
                    interp_error = original_interp_error*(coefficient_interp_error)**(-current_level)
                    metric_param["hessian_strategy_parameters"]["interpolation_error"].SetDouble(interp_error)
                model_part_name = parameters_coarse["solver_settings"]["model_part_name"].GetString()
                # Setting Metric Tensor to 0
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,model_coarse.GetModelPart(model_part_name).Nodes)
                # calculate NODAL_H
                find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_coarse.GetModelPart(model_part_name))
                find_nodal_h.Execute()
                local_gradient = KratosMeshing.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),KratosMultiphysics.TEMPERATURE,metric_param)
                local_gradient.Execute()

            # create the remeshing process
            MmgProcess = KratosMeshing.MmgProcess2D(model_coarse.GetModelPart(model_part_name),remesh_param)
            MmgProcess.Execute()

            # reset variables if needed
            model_coarse.GetModelPart(model_part_name).ProcessInfo.SetValue(KratosMultiphysics.TIME , 0.0)
            model_coarse.GetModelPart(model_part_name).ProcessInfo.SetValue(KratosMultiphysics.STEP , 0)

            """
            the refinement process empties the coarse model part object and fill it with the refined model part
            the solution on the refined grid is obtained from the interpolation of the coarse solution
            there are not other operations, therefore to build the new model we just need to take the updated coarse model
            """
            current_model_refined = model_coarse
            current_parameters_refined = parameters_coarse
            return current_model_refined,current_parameters_refined

    """
    method computing the mesh size of coarsest level, estimated as minimum nodal_h
    input:  self : an instance of the class
    """
    def ComputeMeshSizeCoarsestLevel(self):
        model_coarse = self.model_coarse
        parameters_coarse = self.parameters_coarse
        model_part_name = self.wrapper.GetModelPartName()
        # set NODAL_AREA and NODAL_H as non historical variables
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, model_coarse.GetModelPart(model_part_name).Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, model_coarse.GetModelPart(model_part_name).Nodes)
        # calculate NODAL_H
        find_nodal_h = KratosMultiphysics.FindNodalHProcess(model_coarse.GetModelPart(model_part_name))
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_coarse.GetModelPart(model_part_name))
        find_nodal_h.Execute()
        # compute average mesh size
        mesh_size = 10.0
        for node in model_coarse.GetModelPart(model_part_name).Nodes:
            if (node.GetValue(KratosMultiphysics.NODAL_H) < mesh_size):
                mesh_size = node.GetValue(KratosMultiphysics.NODAL_H)
        self.mesh_size_coarsest_level = mesh_size

    """
    method estimating the mesh size of current level
    input:  self : an instance of the class
    """
    def EstimateMeshSizeCurrentLevel(self):
        self.ComputeMeshSizeCoarsestLevel()
        current_level = self.current_level
        if (self.metric is "hessian"):
            original_interp_error = self.metric_param["hessian_strategy_parameters"]["interpolation_error"].GetDouble()
            domain_size = self.wrapper.GetDomainSize()
            if (domain_size == 2):
                coefficient = 2/9 # 2d
            elif (domain_size == 3):
                coefficient = 9/32 # 3d
            # TODO: compute below interp error level more automatically
            coefficient_interp_error =  self.metric_param["hessian_strategy_parameters"]["coefficient_interpolation_error"].GetDouble()
            # interp_error_level = original_interp_error*(coefficient_interp_error)**(-current_level)
            if (current_level > 0):
                interp_error_level = original_interp_error/(coefficient_interp_error*current_level)
            else: # current_level == 0
                interp_error_level = original_interp_error
            mesh_size_level = self.mesh_size_coarsest_level*np.sqrt(interp_error_level/original_interp_error) # relation from [Alauzet] eqs. pag 34 and 35
            self.mesh_size = mesh_size_level
