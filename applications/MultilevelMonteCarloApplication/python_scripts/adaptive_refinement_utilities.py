from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.MeshingApplication as KratosMeshing

"""
References:
F. Alauzet, Metric-based anisotropic mesh adaptation, CEA-EDF-INRIA schools: Numerical Analysis Summer School. CEA, Cadarache, France
Kratos wiki: https://github.com/KratosMultiphysics/Kratos/wiki/MMG-Process
"""


class AdaptiveRefinement(object):
    """
    input:  model_coarse:       Kratos model class before refinement
            parameters_coarse:  Kratos parameters class before refinement
            minimal_size_value: minimal size after remeshing
            maximal_size_value: maximal size after remeshing
            metric_param:       Kratos parameters class containing metric custom settings
            remesh_param:       Kratos parameters class containing remeshing custom settings
    """
    def __init__(self,model_coarse,parameters_coarse,metric_param,remesh_param,minimal_size_value=None,maximal_size_value=None,metric_name="hessian"):
        self.model_coarse = model_coarse
        self.parameters_coarse = parameters_coarse
        self.minimal_size = minimal_size_value
        self.maximal_size = maximal_size_value
        self.metric_param = metric_param
        self.remesh_param = remesh_param
        self.problem_type = self.parameters_coarse["problem_data"]["problem_name"].GetString()
        self.metric = metric_name

    """
    function computing the refinement of the model based on the solution on the coarse mesh,
    exploiting the hessian metric of the solution
    input:  self: an instance of the class
    output: current_model_refined:      Kratos model class after refinement
            current_parameters_refined: Kratos parameters class after refinement
    """
    def ComputeAdaptiveRefinement(self):
        parameters_coarse = self.parameters_coarse
        model_coarse = self.model_coarse
        minimal_size_value = self.minimal_size
        maximal_size_value = self.maximal_size
        metric_param = self.metric_param
        remesh_param = self.remesh_param
        problem_type = self.problem_type

        if (self.metric is "hessian"):
            model_part_name = parameters_coarse["problem_data"]["model_part_name"].GetString()
            # set NODAL_AREA and NODAL_H as non historical variables
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, model_coarse.GetModelPart(model_part_name).Nodes)
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, model_coarse.GetModelPart(model_part_name).Nodes)
            # calculate NODAL_H
            find_nodal_h = KratosMultiphysics.FindNodalHProcess(model_coarse.GetModelPart(model_part_name))
            find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_coarse.GetModelPart(model_part_name))
            find_nodal_h.Execute()
            # prepare parameters to calculate the gradient of the designed variable
            local_gradient_variable_string = metric_param["local_gradient_variable"].GetString()
            local_gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable(metric_param["local_gradient_variable"].GetString())
            metric_param.RemoveValue("local_gradient_variable")
            metric_param.AddEmptyValue("minimal_size")
            metric_param["minimal_size"].SetDouble(minimal_size_value)
            metric_param.AddEmptyValue("maximal_size")
            metric_param["maximal_size"].SetDouble(maximal_size_value)
            # calculate the gradient of the variable
            local_gradient = KratosMeshing.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),local_gradient_variable,metric_param)
            local_gradient.Execute()
            # add again the removed variable parameter
            metric_param.AddEmptyValue("local_gradient_variable")
            metric_param["local_gradient_variable"].SetString(local_gradient_variable_string)
            # create the remeshing process
            MmgProcess = KratosMeshing.MmgProcess2D(model_coarse.GetModelPart(model_part_name),remesh_param)
            MmgProcess.Execute()
            """
            the refinement process empties the coarse model part object and fill it with the refined model part
            the solution on the refined grid is obtained from the interpolation of the coarse solution
            there are not other operations, therefore to build the new model we just need to take the updated coarse model
            """
            current_model_refined = model_coarse
            current_parameters_refined = parameters_coarse
            return current_model_refined,current_parameters_refined

