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
def compute_refinement(serialized_model_coarse,minimal_size_value,maximal_size_value,serialized_parameters_coarse,name_model_part):

    #pickle dataserialized_data
    pickled_data = pickle.dumps(serialized_model_coarse, 2) #second argument is the protocol and is NECESSARY (according to pybind11 docs)
    #overwrite the old serializer with the unpickled one
    serialized_model_coarse = pickle.loads(pickled_data)
    current_model = KratosMultiphysics.Model()
    serialized_model_coarse.Load("ModelSerialization",current_model)
    del(serialized_model_coarse)

    '''serialize model refined'''
    serialized_model_refined = KratosMultiphysics.StreamSerializer()
    serialized_model_refined.Save("ModelSerialization",current_model)

    # '''calculate NODAL_H'''
    # find_nodal_h = KratosMultiphysics.FindNodalHProcess(current_model.GetModelPart(name_model_part))
    # # find_nodal_h.Execute()

    # '''calculate the gradient of the distance variable'''
    # metric_param = KratosMultiphysics.Parameters(
    #     """{
    #         "hessian_strategy_parameters"              :{
    #                 "estimate_interpolation_error"     : false,
    #                 "interpolation_error"              : 0.004,
    #                 "mesh_dependent_constant"          : 0.28125
    #         },
    #         "enforce_current"                   : false,
    #         "anisotropy_remeshing"              : true,
    #         "anisotropy_parameters":{
    #             "reference_variable_name"          : "TEMPERATURE",
    #             "hmin_over_hmax_anisotropic_ratio" : 0.15,
    #             "boundary_layer_max_distance"      : 1.0,
    #             "interpolation"                    : "Linear"
    #         }
    #     }"""
    #     )
    # metric_param.AddEmptyValue("minimal_size")
    # metric_param["minimal_size"].SetDouble(minimal_size_value)
    # metric_param.AddEmptyValue("maximal_size")
    # metric_param["maximal_size"].SetDouble(maximal_size_value)

    # local_gradient = MeshingApplication.ComputeHessianSolMetricProcess2D(current_model.GetModelPart(name_model_part), KratosMultiphysics.TEMPERATURE, metric_param)
    # # local_gradient.Execute()

    # '''create the remeshing process'''
    # remesh_param = KratosMultiphysics.Parameters("""{ }""")
    # MmgProcess = MeshingApplication.MmgProcess2D(current_model.GetModelPart(name_model_part), remesh_param)
    # # MmgProcess.Execute()

    # # Finally we export to GiD
    # from gid_output_process import GiDOutputProcess
    # gid_output = GiDOutputProcess(current_model.GetModelPart(name_model_part),
    #                             "gid_output",
    #                             KratosMultiphysics.Parameters("""
    #                                 {
    #                                     "result_file_configuration" : {
    #                                         "gidpost_flags": {
    #                                             "GiDPostMode": "GiD_PostBinary",
    #                                             "WriteDeformedMeshFlag": "WriteUndeformed",
    #                                             "WriteConditionsFlag": "WriteConditions",
    #                                             "MultiFileFlag": "SingleFile"
    #                                         },        
    #                                         "nodal_results"       : ["TEMPERATURE","NODAL_H"]
    #                                     }
    #                                 }
    #                                 """)
    #                             )

    # gid_output.ExecuteInitialize()
    # gid_output.ExecuteBeforeSolutionLoop()
    # gid_output.ExecuteInitializeSolutionStep()
    # gid_output.PrintOutput()
    # gid_output.ExecuteFinalizeSolutionStep()
    # gid_output.ExecuteFinalize()

    # print(current_model.GetModelPart(name_model_part))

    # for prop in current_model.GetModelPart(name_model_part).Properties:
    #     print(prop.GetValue(KratosMultiphysics.CONDUCTIVITY))


    # '''the refinement process empties the coarse model part object and fill it with the refined model part
    # the solution on the refined grid is obtained from the interpolation of the coarse solution
    # there are not other operations, so to build the new model serialized we just need to change the model part of the coarse model serialized'''

    # '''try to reset the process info, GetPreviousSolutionStepInfo seems to allow to write again a gid file'''
    # # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0)
    # # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
    # # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)
    # # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.PRINTED_STEP, 0)
    # # current_model.GetModelPart(name_model_part).ProcessInfo.GetPreviousSolutionStepInfo()

    # '''serialize model refined'''
    # serialized_model_refined = KratosMultiphysics.StreamSerializer()
    # serialized_model_refined.Save("ModelSerialization",current_model)

    # '''build serializer refined parameters'''
    serialized_parameters_refined = serialized_parameters_coarse # check there is no reference (apparently not)

    # '''exploiting this second way I add manually the nodes, elements, etc, the problem is I am not able to load automatically the SubModelParts'''
    # # current_model_refined = KratosMultiphysics.Model()
    # # serialized_model_coarse.Load("ModelSerialization",current_model_refined)
    # # '''remove coarse model part and add refined model part'''
    # # current_model_refined.DeleteModelPart("MLMCLaplacianModelPart")
    # # current_model_refined.CreateModelPart("MLMCLaplacianModelPart")
    
    # # '''model part'''
    # # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetConditions(simulation_coarse._GetSolver().main_model_part.GetConditions())
    # # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetElements(simulation_coarse._GetSolver().main_model_part.GetElements())
    # # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetNodes(simulation_coarse._GetSolver().main_model_part.GetNodes())
    # # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetProperties(simulation_coarse._GetSolver().main_model_part.GetProperties())
    # # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetBufferSize(simulation_coarse._GetSolver().main_model_part.GetBufferSize())
    # # '''here I should do the same for all the sub model parts'''    
    # # ...

    return serialized_model_refined,serialized_parameters_refined


'''see /Examples/mmg_remeshing_examples/validation/hessian2D/source/test_hessian.py for details'''
def compute_refinement_from_analysisstage_object(simulation_coarse,minimal_size_value,maximal_size_value,serialized_parameters_coarse):

    simulation_coarse._GetSolver().print_on_rank_zero("::[compute_refinement]:: ", "refinement started")
    '''calculate NODAL_H'''
    find_nodal_h = KratosMultiphysics.FindNodalHProcess(simulation_coarse._GetSolver().main_model_part)
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
                "reference_variable_name"          : "DISTANCE",
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

    '''create the remeshing process'''
    remesh_param = KratosMultiphysics.Parameters("""{ }""")
    MmgProcess = MeshingApplication.MmgProcess2D(simulation_coarse._GetSolver().main_model_part, remesh_param)
    MmgProcess.Execute()

    '''the refinement process empties the coarse model part object and fill it with the refined model part
    the solution on the refined grid is obtained from the interpolation of the coarse solution
    there are not other operations, so to build the new model serialized we just need to change the model part of the coarse model serialized'''

    simulation_coarse._GetSolver().print_on_rank_zero("::[compute_refinement]:: ", "start saving refined model and parameters")

    '''build serializer refined model'''
    serialized_model_refined = KratosMultiphysics.StreamSerializer()
    current_model_refined = simulation_coarse.model
    
    '''try to reset the process info, GetPreviousSolutionStepInfo seems to allow to write again a gid file'''
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0)
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.PRINTED_STEP, 0)
    current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.GetPreviousSolutionStepInfo()
    '''serialize model refined'''
    serialized_model_refined.Save("ModelSerialization",current_model_refined)

    '''build serializer refined parameters'''
    serialized_parameters_refined = serialized_parameters_coarse # check there is no reference (apparently not)

    simulation_coarse._GetSolver().print_on_rank_zero("::[compute_refinement]:: ", "end saving refined model and parameters")
    simulation_coarse._GetSolver().print_on_rank_zero("::[compute_refinement]:: ", "refinement ended")
    del(simulation_coarse)

    '''exploiting this second way I add manually the nodes, elements, etc, the problem is I am not able to load automatically the SubModelParts'''
    # current_model_refined = KratosMultiphysics.Model()
    # serialized_model_coarse.Load("ModelSerialization",current_model_refined)
    # '''remove coarse model part and add refined model part'''
    # current_model_refined.DeleteModelPart("MLMCLaplacianModelPart")
    # current_model_refined.CreateModelPart("MLMCLaplacianModelPart")
    
    # '''model part'''
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetConditions(simulation_coarse._GetSolver().main_model_part.GetConditions())
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetElements(simulation_coarse._GetSolver().main_model_part.GetElements())
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetNodes(simulation_coarse._GetSolver().main_model_part.GetNodes())
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetProperties(simulation_coarse._GetSolver().main_model_part.GetProperties())
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetBufferSize(simulation_coarse._GetSolver().main_model_part.GetBufferSize())
    # '''here I should do the same for all the sub model parts'''    
    # ...

    return serialized_model_refined,serialized_parameters_refined

'''see /Examples/mmg_remeshing_examples/validation/hessian2D/source/test_hessian.py for details'''
def compute_refinement_run_inside(serialized_model_coarse,minimal_size_value,maximal_size_value,serialized_parameters_coarse,name_model_part):

    #pickle dataserialized_data
    pickled_data = pickle.dumps(serialized_model_coarse, 2) #second argument is the protocol and is NECESSARY (according to pybind11 docs)
    #overwrite the old serializer with the unpickled one
    serialized_model_coarse = pickle.loads(pickled_data)
    current_model = KratosMultiphysics.Model()
    serialized_model_coarse.Load("ModelSerialization",current_model)
    del(serialized_model_coarse)
    

    '''calculate NODAL_H'''
    find_nodal_h = KratosMultiphysics.FindNodalHProcess(current_model.GetModelPart(name_model_part))
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

    local_gradient = MeshingApplication.ComputeHessianSolMetricProcess2D(current_model.GetModelPart(name_model_part), KratosMultiphysics.TEMPERATURE, metric_param)
    local_gradient.Execute()

    '''create the remeshing process'''
    remesh_param = KratosMultiphysics.Parameters("""{ }""")
    MmgProcess = MeshingApplication.MmgProcess2D(current_model.GetModelPart(name_model_part), remesh_param)
    MmgProcess.Execute()

    # Finally we export to GiD
    from gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(current_model.GetModelPart(name_model_part),
                                "gid_output",
                                KratosMultiphysics.Parameters("""
                                    {
                                        "result_file_configuration" : {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostBinary",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "WriteConditionsFlag": "WriteConditions",
                                                "MultiFileFlag": "SingleFile"
                                            },        
                                            "nodal_results"       : ["TEMPERATURE","NODAL_H"]
                                        }
                                    }
                                    """)
                                )

    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

    print(current_model.GetModelPart(name_model_part))

    for prop in current_model.GetModelPart(name_model_part).Properties:
        print(prop.GetValue(KratosMultiphysics.CONDUCTIVITY))


    '''the refinement process empties the coarse model part object and fill it with the refined model part
    the solution on the refined grid is obtained from the interpolation of the coarse solution
    there are not other operations, so to build the new model serialized we just need to change the model part of the coarse model serialized'''

    '''try to reset the process info, GetPreviousSolutionStepInfo seems to allow to write again a gid file'''
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0)
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.TIME, 0)
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").ProcessInfo.SetValue(KratosMultiphysics.PRINTED_STEP, 0)
    # current_model.GetModelPart(name_model_part).ProcessInfo.GetPreviousSolutionStepInfo()

    '''serialize model refined'''
    serialized_model_refined = KratosMultiphysics.StreamSerializer()
    serialized_model_refined.Save("ModelSerialization",current_model)

    '''build serializer refined parameters'''
    serialized_parameters_refined = serialized_parameters_coarse # check there is no reference (apparently not)

    '''exploiting this second way I add manually the nodes, elements, etc, the problem is I am not able to load automatically the SubModelParts'''
    # current_model_refined = KratosMultiphysics.Model()
    # serialized_model_coarse.Load("ModelSerialization",current_model_refined)
    # '''remove coarse model part and add refined model part'''
    # current_model_refined.DeleteModelPart("MLMCLaplacianModelPart")
    # current_model_refined.CreateModelPart("MLMCLaplacianModelPart")
    
    # '''model part'''
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetConditions(simulation_coarse._GetSolver().main_model_part.GetConditions())
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetElements(simulation_coarse._GetSolver().main_model_part.GetElements())
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetNodes(simulation_coarse._GetSolver().main_model_part.GetNodes())
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetProperties(simulation_coarse._GetSolver().main_model_part.GetProperties())
    # current_model_refined.GetModelPart("MLMCLaplacianModelPart").SetBufferSize(simulation_coarse._GetSolver().main_model_part.GetBufferSize())
    # '''here I should do the same for all the sub model parts'''    
    # ...

    return serialized_model_refined,serialized_parameters_refined
