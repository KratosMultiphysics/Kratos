from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from copy import copy

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.MeshingApplication as MeshingApplication

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle


'''see /Examples/mmg_remeshing_examples/validation/hessian2D/source/test_hessian.py for details'''
def compute_refinement_hessian_metric(simulation_coarse,minimal_size_value,maximal_size_value):

    simulation_coarse._GetSolver().print_on_rank_zero("::[compute_refinement]:: ", "refinement started")
    '''set NODAL_AREA and NODAL_H as non historical variables'''
    KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, simulation_coarse._GetSolver().main_model_part.Nodes)
    KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, simulation_coarse._GetSolver().main_model_part.Nodes)

    '''calculate NODAL_H'''
    find_nodal_h = KratosMultiphysics.FindNodalHProcess(simulation_coarse._GetSolver().main_model_part)
    find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(simulation_coarse._GetSolver().main_model_part)
    find_nodal_h.Execute()

    '''calculate the gradient of the distance variable'''
    metric_param = KratosMultiphysics.Parameters(
        """{
            "hessian_strategy_parameters"              :{
                    "estimate_interpolation_error"     : false,
                    "interpolation_error"              : 0.004,
                    "mesh_dependent_constant"          : 0.28125
            },
            "enforce_current"                   : false,
            "anisotropy_remeshing"              : true,
            "anisotropy_parameters":{
                "reference_variable_name"          : "TEMPERATURE",
                "hmin_over_hmax_anisotropic_ratio" : 0.15,
                "boundary_layer_max_distance"      : 1.0,
                "interpolation"                    : "Linear"
            }
        }"""
        )
    metric_param.AddEmptyValue("minimal_size")
    metric_param["minimal_size"].SetDouble(minimal_size_value)
    metric_param.AddEmptyValue("maximal_size")
    metric_param["maximal_size"].SetDouble(maximal_size_value)

    local_gradient = MeshingApplication.ComputeHessianSolMetricProcess2D(simulation_coarse._GetSolver().main_model_part, KratosMultiphysics.TEMPERATURE, metric_param)
    local_gradient.Execute()
    
    '''create the remeshing process: echo_level: 0 for no output at all, 3 for standard output'''
    remesh_param = KratosMultiphysics.Parameters(
        """{
            "echo_level"                       : 0}"""
            )
    # remesh_param.AddEmptyValue("echo_level")
    # remesh_param["echo_level"].SetInt(mmg_remeshing_info)
    MmgProcess = MeshingApplication.MmgProcess2D(simulation_coarse._GetSolver().main_model_part, remesh_param)
    MmgProcess.Execute()
    
    '''the refinement process empties the coarse model part object and fill it with the refined model part
    the solution on the refined grid is obtained from the interpolation of the coarse solution
    there are not other operations, so to build the new model we just need to take the updated coarse model'''

    simulation_coarse._GetSolver().print_on_rank_zero("::[compute_refinement]:: ", "start saving refined model and parameters")

    current_model_refined = simulation_coarse.model

    return current_model_refined