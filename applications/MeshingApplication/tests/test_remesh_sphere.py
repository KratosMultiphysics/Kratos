
# We import the libraies
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

# Import stuff
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

class TestRemeshMMG3D(KratosUnittest.TestCase):

    def test_remesh_sphere(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)
        #main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 1)

        # We add the variables needed
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

        # We import the model main_model_part
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO(file_path + "/mmg_eulerian_test/coarse_sphere_test").ReadModelPart(main_model_part)

        # We calculate the gradient of the distance variable
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, main_model_part.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero the metric
        ZeroVector = KratosMultiphysics.Vector(6)
        ZeroVector[0] = 0.0
        ZeroVector[1] = 0.0
        ZeroVector[2] = 0.0
        ZeroVector[3] = 0.0
        ZeroVector[4] = 0.0
        ZeroVector[5] = 0.0

        for node in main_model_part.Nodes:
            node.SetValue(MeshingApplication.METRIC_TENSOR_3D, ZeroVector)

        # We define a metric using the ComputeLevelSetSolMetricProcess
        MetricParameters = KratosMultiphysics.Parameters("""
        {
            "minimal_size"                      : 1.0e-1,
            "enforce_current"                   : false,
            "anisotropy_remeshing"              : false,
            "anisotropy_parameters"             :{
                "hmin_over_hmax_anisotropic_ratio"  : 0.15,
                "boundary_layer_max_distance"       : 1.0e-4,
                "interpolation"                     : "Linear"
            }
        }
        """)
        metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess3D(main_model_part, KratosMultiphysics.DISTANCE_GRADIENT, MetricParameters)
        metric_process.Execute()

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "filename"                         : "mmg_eulerian_test/coarse_sphere_test",
            "save_external_files"              : true,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        mmg_parameters["filename"].SetString(file_path + "/" + mmg_parameters["filename"].GetString())
        mmg_process = MeshingApplication.MmgProcess3D(main_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()

        # Finally we export to GiD
        gid_output = GiDOutputProcess(main_model_part,
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
                                                "nodal_results"       : ["DISTANCE"]
                                            }
                                        }
                                        """)
                                    )

        #gid_output.ExecuteInitialize()
        #gid_output.ExecuteBeforeSolutionLoop()
        #gid_output.ExecuteInitializeSolutionStep()
        #gid_output.PrintOutput()
        #gid_output.ExecuteFinalizeSolutionStep()
        #gid_output.ExecuteFinalize()

        check_parameters = KratosMultiphysics.Parameters("""
                            {
                                "reference_file_name"   : "mmg_eulerian_test/coarse_sphere_test_result.sol",
                                "output_file_name"      : "mmg_eulerian_test/coarse_sphere_test_step=0.sol",
                                "dimension"             : 3,
                                "comparison_type"       : "sol_file"
                            }
                            """)
        check_parameters["reference_file_name"].SetString(file_path + "/" + check_parameters["reference_file_name"].GetString())
        check_parameters["output_file_name"].SetString(file_path + "/" + check_parameters["output_file_name"].GetString())
        check_files = CompareTwoFilesCheckProcess(check_parameters)

        check_files.ExecuteInitialize()
        check_files.ExecuteBeforeSolutionLoop()
        check_files.ExecuteInitializeSolutionStep()
        check_files.ExecuteFinalizeSolutionStep()
        check_files.ExecuteFinalize()

        check_parameters = KratosMultiphysics.Parameters("""
        {
            "check_variables"      : ["DISTANCE"],
            "input_file_name"      : "mmg_eulerian_test/distante_extrapolation.json",
            "model_part_name"      : "MainModelPart",
            "time_frequency"       : 0.0
        }
        """)

        check_parameters["input_file_name"].SetString(file_path + "/" + check_parameters["input_file_name"].GetString())
        check = FromJsonCheckResultProcess(current_model, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()
        check.ExecuteFinalizeSolutionStep()

        ## The following is used to create the solution database
        #import json_output_process

        #out_parameters = KratosMultiphysics.Parameters("""
        #{
            #"output_variables"     : ["DISTANCE"],
            #"output_file_name"     : "mmg_eulerian_test/distante_extrapolation.json",
            #"model_part_name"      : "MainModelPart",
            #"time_frequency"       : 0.0
        #}
        #""")

        #out = json_output_process.JsonOutputProcess(current_model, out_parameters)
        #out.ExecuteInitialize()
        #out.ExecuteBeforeSolutionLoop()
        #out.ExecuteFinalizeSolutionStep()

    def test_remesh_sphere_skin(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)
        #main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 1)

        # We add the variables needed
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        # We import the model main_model_part
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO(file_path + "/mmg_eulerian_test/coarse_sphere_skin_test").ReadModelPart(main_model_part)

        for node in main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, abs(node.X))

        # We calculate the gradient of the distance variable
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, main_model_part.Nodes)
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()

        # We set to zero the metric
        metric_vector = KratosMultiphysics.Vector(6)
        metric_vector[0] = 1.0
        metric_vector[1] = 1.0
        metric_vector[2] = 1.0
        metric_vector[3] = 0.0
        metric_vector[4] = 0.0
        metric_vector[5] = 0.0

        for node in main_model_part.Nodes:
            node.SetValue(MeshingApplication.METRIC_TENSOR_3D, metric_vector)

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "filename"                         : "mmg_eulerian_test/coarse_sphere_skin_test",
            "save_external_files"              : true,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        mmg_parameters["filename"].SetString(file_path + "/" + mmg_parameters["filename"].GetString())
        mmg_process = MeshingApplication.MmgProcess3DSurfaces(main_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()

        # Finally we export to GiD
        gid_output = GiDOutputProcess(main_model_part,
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
                                                "nodal_results"       : ["DISTANCE"]
                                            }
                                        }
                                        """)
                                    )

        #gid_output.ExecuteInitialize()
        #gid_output.ExecuteBeforeSolutionLoop()
        #gid_output.ExecuteInitializeSolutionStep()
        #gid_output.PrintOutput()
        #gid_output.ExecuteFinalizeSolutionStep()
        #gid_output.ExecuteFinalize()

        check_parameters = KratosMultiphysics.Parameters("""
                            {
                                "reference_file_name"   : "mmg_eulerian_test/coarse_sphere_skin_test_result.sol",
                                "output_file_name"      : "mmg_eulerian_test/coarse_sphere_skin_test_step=0.sol",
                                "dimension"             : 3,
                                "comparison_type"       : "sol_file"
                            }
                            """)
        check_parameters["reference_file_name"].SetString(file_path + "/" + check_parameters["reference_file_name"].GetString())
        check_parameters["output_file_name"].SetString(file_path + "/" + check_parameters["output_file_name"].GetString())
        check_files = CompareTwoFilesCheckProcess(check_parameters)

        check_files.ExecuteInitialize()
        check_files.ExecuteBeforeSolutionLoop()
        check_files.ExecuteInitializeSolutionStep()
        check_files.ExecuteFinalizeSolutionStep()
        check_files.ExecuteFinalize()

        check_parameters = KratosMultiphysics.Parameters("""
        {
            "check_variables"      : ["DISTANCE"],
            "input_file_name"      : "mmg_eulerian_test/distante_extrapolation_skin.json",
            "model_part_name"      : "MainModelPart",
            "time_frequency"       : 0.0
        }
        """)

        check_parameters["input_file_name"].SetString(os.path.join(file_path, check_parameters["input_file_name"].GetString()))
        check = FromJsonCheckResultProcess(current_model, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()
        check.ExecuteFinalizeSolutionStep()

        ## The following is used to create the solution database
        #import json_output_process

        #out_parameters = KratosMultiphysics.Parameters("""
        #{
            #"output_variables"     : ["DISTANCE"],
            #"output_file_name"     : "mmg_eulerian_test/distante_extrapolation_skin.json",
            #"model_part_name"      : "MainModelPart",
            #"time_frequency"       : 0.0
        #}
        #""")

        #out = json_output_process.JsonOutputProcess(current_model, out_parameters)
        #out.ExecuteInitialize()
        #out.ExecuteBeforeSolutionLoop()
        #out.ExecuteFinalizeSolutionStep()

    def test_remesh_sphere_skin_prisms(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)
        #main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 1)

        # We add the variables needed
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        # We import the model main_model_part
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO(file_path + "/mmg_eulerian_test/coarse_sphere_skin_prisms_test").ReadModelPart(main_model_part)

        for node in main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, abs(node.X))

        # We calculate the gradient of the distance variable
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, main_model_part.Nodes)
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()

        # We set to zero the metric
        metric_vector = KratosMultiphysics.Vector(6)
        metric_vector[0] = 1.0
        metric_vector[1] = 1.0
        metric_vector[2] = 1.0
        metric_vector[3] = 0.0
        metric_vector[4] = 0.0
        metric_vector[5] = 0.0

        for node in main_model_part.Nodes:
            node.SetValue(MeshingApplication.METRIC_TENSOR_3D, metric_vector)

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "filename"                         : "mmg_eulerian_test/coarse_sphere_skin_prisms_test",
            "collapse_prisms_elements"         : true,
            "save_external_files"              : true,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        mmg_parameters["filename"].SetString(file_path + "/" + mmg_parameters["filename"].GetString())
        mmg_process = MeshingApplication.MmgProcess3DSurfaces(main_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()

        ## Finally we export to GiD
        #from gid_output_process import GiDOutputProcess
        #gid_output = GiDOutputProcess(main_model_part,
                                    #"gid_output",
                                    #KratosMultiphysics.Parameters("""
                                        #{
                                            #"result_file_configuration" : {
                                                #"gidpost_flags": {
                                                    #"GiDPostMode": "GiD_PostBinary",
                                                    #"WriteDeformedMeshFlag": "WriteUndeformed",
                                                    #"WriteConditionsFlag": "WriteConditions",
                                                    #"MultiFileFlag": "SingleFile"
                                                #},
                                                #"nodal_results"       : ["DISTANCE"]
                                            #}
                                        #}
                                        #""")
                                    #)

        #gid_output.ExecuteInitialize()
        #gid_output.ExecuteBeforeSolutionLoop()
        #gid_output.ExecuteInitializeSolutionStep()
        #gid_output.PrintOutput()
        #gid_output.ExecuteFinalizeSolutionStep()
        #gid_output.ExecuteFinalize()

        check_parameters = KratosMultiphysics.Parameters("""
                            {
                                "reference_file_name"   : "mmg_eulerian_test/coarse_sphere_skin_prisms_test_result.sol",
                                "output_file_name"      : "mmg_eulerian_test/coarse_sphere_skin_prisms_test_step=0.sol",
                                "dimension"             : 3,
                                "comparison_type"       : "sol_file"
                            }
                            """)
        check_parameters["reference_file_name"].SetString(file_path + "/" + check_parameters["reference_file_name"].GetString())
        check_parameters["output_file_name"].SetString(file_path + "/" + check_parameters["output_file_name"].GetString())
        check_files = CompareTwoFilesCheckProcess(check_parameters)

        check_files.ExecuteInitialize()
        check_files.ExecuteBeforeSolutionLoop()
        check_files.ExecuteInitializeSolutionStep()
        check_files.ExecuteFinalizeSolutionStep()
        check_files.ExecuteFinalize()

    def test_isosurface_remesh_sphere(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)
        #main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 1)

        # We add the variables needed
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

        # We import the model main_model_part
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO(file_path + "/mmg_eulerian_test/test_sphere_isosurface").ReadModelPart(main_model_part)

        # Set manually a distance function
        circle_radious = 0.25
        center_coordinates = [0.5, 0.5, 0.5]

        for node in main_model_part.Nodes:
            distance = ((node.X-center_coordinates[0])**2+(node.Y-center_coordinates[1])**2+(node.Z-center_coordinates[2])**2)**0.5 - circle_radious
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

        # We calculate the gradient of the distance variable
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "discretization_type"              : "Isosurface",
            "filename"                         : "mmg_eulerian_test/test_sphere_isosurface",
            "save_external_files"              : true,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        mmg_parameters["filename"].SetString(file_path + "/" + mmg_parameters["filename"].GetString())
        mmg_process = MeshingApplication.MmgProcess3D(main_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()

        # Finally we export to GiD
        gid_output = GiDOutputProcess(main_model_part,
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
                                                "nodal_results"       : ["DISTANCE"]
                                            }
                                        }
                                        """)
                                    )

        #gid_output.ExecuteInitialize()
        #gid_output.ExecuteBeforeSolutionLoop()
        #gid_output.ExecuteInitializeSolutionStep()
        #gid_output.PrintOutput()
        #gid_output.ExecuteFinalizeSolutionStep()
        #gid_output.ExecuteFinalize()

        check_parameters = KratosMultiphysics.Parameters("""
                            {
                                "reference_file_name"   : "mmg_eulerian_test/test_sphere_isosurface_result.sol",
                                "output_file_name"      : "mmg_eulerian_test/test_sphere_isosurface_step=0.o.sol",
                                "dimension"             : 3,
                                "comparison_type"       : "sol_file"
                            }
                            """)
        check_parameters["reference_file_name"].SetString(file_path + "/" + check_parameters["reference_file_name"].GetString())
        check_parameters["output_file_name"].SetString(file_path + "/" + check_parameters["output_file_name"].GetString())
        check_files = CompareTwoFilesCheckProcess(check_parameters)

        check_files.ExecuteInitialize()
        check_files.ExecuteBeforeSolutionLoop()
        check_files.ExecuteInitializeSolutionStep()
        check_files.ExecuteFinalizeSolutionStep()
        check_files.ExecuteFinalize()

if __name__ == '__main__':
    KratosUnittest.main()
