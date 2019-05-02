from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.MeshingApplication as KratosMeshing


'''
References:
F. Alauzet, Metric-based anisotropic mesh adaptation, CEA-EDF-INRIA schools: Numerical Analysis Summer School. CEA, Cadarache, France
Kratos wiki: https://github.com/KratosMultiphysics/Kratos/wiki/MMG-Process
'''

'''
function computing the refinement of the model based on the solution on the coarse mesh,
exploiting the hessian metric of the solution
'''
def compute_refinement_hessian_metric(simulation_coarse,minimal_size_value,maximal_size_value,metric_param,remesh_param):

    simulation_coarse._GetSolver().print_on_rank_zero("::[compute_refinement]:: ", "refinement started")
    '''set NODAL_AREA and NODAL_H as non historical variables'''
    KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, simulation_coarse._GetSolver().main_model_part.Nodes)
    KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, simulation_coarse._GetSolver().main_model_part.Nodes)

    '''calculate NODAL_H'''
    find_nodal_h = KratosMultiphysics.FindNodalHProcess(simulation_coarse._GetSolver().main_model_part)
    find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(simulation_coarse._GetSolver().main_model_part)
    find_nodal_h.Execute()

    '''prepare parameters to calculate the gradient of the designed variable'''
    local_gradient_variable_string = metric_param["local_gradient_variable"].GetString()
    local_gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable(metric_param["local_gradient_variable"].GetString())
    metric_param.RemoveValue("local_gradient_variable")
    metric_param.AddEmptyValue("minimal_size")
    metric_param["minimal_size"].SetDouble(minimal_size_value)
    metric_param.AddEmptyValue("maximal_size")
    metric_param["maximal_size"].SetDouble(maximal_size_value)
    '''calculate the gradient of the variable'''
    local_gradient = KratosMeshing.ComputeHessianSolMetricProcess(simulation_coarse._GetSolver().main_model_part,local_gradient_variable,metric_param)
    local_gradient.Execute()
    '''add again the removed variable parameter'''
    metric_param.AddEmptyValue("local_gradient_variable")
    metric_param["local_gradient_variable"].SetString(local_gradient_variable_string)

    '''create the remeshing process'''
    MmgProcess = KratosMeshing.MmgProcess2D(simulation_coarse._GetSolver().main_model_part,remesh_param)
    MmgProcess.Execute()

    '''the refinement process empties the coarse model part object and fill it with the refined model part
    the solution on the refined grid is obtained from the interpolation of the coarse solution
    there are not other operations, so to build the new model we just need to take the updated coarse model'''

    simulation_coarse._GetSolver().print_on_rank_zero("::[compute_refinement]:: ", "start saving refined model and parameters")

    current_model_refined = simulation_coarse.model

    return current_model_refined