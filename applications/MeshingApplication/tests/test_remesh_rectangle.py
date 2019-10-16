
# We import the libraies
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

class TestRemeshMMG2D(KratosUnittest.TestCase):

    def test_remesh_rectangle_hessian(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart", 2)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)
        #main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 1)

        # We add the variables needed
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        # We import the model main_model_part
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO(file_path + "/mmg_lagrangian_test/remesh_rectangle").ReadModelPart(main_model_part)

        # We calculate the gradient of the distance variable
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, main_model_part.Nodes)

        main_model_part.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISTANCE,  9.86358e-08)
        main_model_part.Nodes[2].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 4.88117e-07)
        main_model_part.Nodes[3].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 8.21649e-07)
        main_model_part.Nodes[4].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 2.98426e-06)
        main_model_part.Nodes[5].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.9462e-06)
        main_model_part.Nodes[6].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 3.20072e-06)
        main_model_part.Nodes[7].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 8.13038e-07)
        main_model_part.Nodes[8].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.74975e-06)
        main_model_part.Nodes[9].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 8.66819e-06)
        main_model_part.Nodes[10].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 6.24329e-06)
        main_model_part.Nodes[11].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.10461e-05)
        main_model_part.Nodes[12].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 7.01043e-07)
        main_model_part.Nodes[13].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.90256e-05)
        main_model_part.Nodes[14].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 2.51694e-06)
        main_model_part.Nodes[15].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.51568e-05)
        main_model_part.Nodes[16].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 3.0177e-05)
        main_model_part.Nodes[17].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 8.48471e-06)
        main_model_part.Nodes[18].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 8.37749e-08)
        main_model_part.Nodes[19].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 4.64653e-07)
        main_model_part.Nodes[20].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 2.64854e-05)
        main_model_part.Nodes[21].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 3.71692e-05)
        main_model_part.Nodes[22].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 3.28506e-06)
        main_model_part.Nodes[23].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.97626e-05)
        main_model_part.Nodes[24].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 6.49145e-05)
        main_model_part.Nodes[25].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1,20e-02)
        main_model_part.Nodes[26].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 4.05303e-05)
        main_model_part.Nodes[27].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 6.49053e-05)
        main_model_part.Nodes[28].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 3.71567e-05)
        main_model_part.Nodes[29].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 3.01782e-05)
        main_model_part.Nodes[30].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000118602)
        main_model_part.Nodes[31].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 5.74756e-05)
        main_model_part.Nodes[32].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 6.26979e-05)
        main_model_part.Nodes[33].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 6.33597e-05)
        main_model_part.Nodes[34].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000101589)
        main_model_part.Nodes[35].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000192486)
        main_model_part.Nodes[36].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 9.93803e-05)
        main_model_part.Nodes[37].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000116324)
        main_model_part.Nodes[38].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 7.68763e-05)
        main_model_part.Nodes[39].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000149051)
        main_model_part.Nodes[40].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000289484)
        main_model_part.Nodes[41].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000192959)
        main_model_part.Nodes[42].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000145925)
        main_model_part.Nodes[43].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 9.79472e-05)
        main_model_part.Nodes[44].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000201518)
        main_model_part.Nodes[45].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000401176)
        main_model_part.Nodes[46].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000291243)
        main_model_part.Nodes[47].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000205528)
        main_model_part.Nodes[48].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000119595)
        main_model_part.Nodes[49].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000262135)
        main_model_part.Nodes[50].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000415149)
        main_model_part.Nodes[51].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.00053798)
        main_model_part.Nodes[52].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000268931)
        main_model_part.Nodes[53].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000139653)
        main_model_part.Nodes[54].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000331001)
        main_model_part.Nodes[55].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000546243)
        main_model_part.Nodes[56].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000675476)
        main_model_part.Nodes[57].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000327043)
        main_model_part.Nodes[58].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000158685)
        main_model_part.Nodes[59].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.00039571)
        main_model_part.Nodes[60].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000651727)
        main_model_part.Nodes[61].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000370579)
        main_model_part.Nodes[62].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000809699)
        main_model_part.Nodes[63].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000170989)
        main_model_part.Nodes[64].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000445764)
        main_model_part.Nodes[65].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000784852)
        main_model_part.Nodes[66].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000433683)
        main_model_part.Nodes[67].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000914023)
        main_model_part.Nodes[68].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000179197)
        main_model_part.Nodes[69].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.00047077)
        main_model_part.Nodes[70].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000904417)
        main_model_part.Nodes[71].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000464902)
        main_model_part.Nodes[72].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000994274)
        main_model_part.Nodes[73].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000173952)
        main_model_part.Nodes[74].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000466276)
        main_model_part.Nodes[75].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000971188)
        main_model_part.Nodes[76].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000474027)
        main_model_part.Nodes[77].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000946149)
        main_model_part.Nodes[78].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000160469)
        main_model_part.Nodes[79].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000397979)
        main_model_part.Nodes[80].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000996743)
        main_model_part.Nodes[81].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000455294)
        main_model_part.Nodes[82].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000925016)
        main_model_part.Nodes[83].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000132904)
        main_model_part.Nodes[84].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000373155)
        main_model_part.Nodes[85].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000963824)
        main_model_part.Nodes[86].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000386671)
        main_model_part.Nodes[87].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000858186)
        main_model_part.Nodes[88].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 9.47787e-05)
        main_model_part.Nodes[89].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000851174)
        main_model_part.Nodes[90].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000282759)
        main_model_part.Nodes[91].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000285314)
        main_model_part.Nodes[92].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000688848)
        main_model_part.Nodes[93].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 4.35228e-05)
        main_model_part.Nodes[94].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000686866)
        main_model_part.Nodes[95].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000179265)
        main_model_part.Nodes[96].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000175715)
        main_model_part.Nodes[97].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000386179)
        main_model_part.Nodes[98].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.59849e-05)
        main_model_part.Nodes[99].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.7304e-05)
        main_model_part.Nodes[100].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0.000382607)
        main_model_part.Nodes[101].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0)
        main_model_part.Nodes[102].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0)
        main_model_part.Nodes[103].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0)
        main_model_part.Nodes[104].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0)
        main_model_part.Nodes[105].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0)

        # We set to zero the metric
        ZeroVector = KratosMultiphysics.Vector(3)
        ZeroVector[0] = 0.0
        ZeroVector[1] = 0.0
        ZeroVector[2] = 0.0

        for node in main_model_part.Nodes:
            node.SetValue(MeshingApplication.METRIC_TENSOR_2D, ZeroVector)

        # We define a metric using the ComputeLevelSetSolMetricProcess
        MetricParameters = KratosMultiphysics.Parameters("""
        {
            "hessian_strategy_parameters"              :{
                "estimate_interpolation_error"     : false,
                "interpolation_error"              : 1.0e-6,
                "mesh_dependent_constant"          : 0.28125
            },
            "minimal_size"                     : 0.015,
            "maximal_size"                     : 0.5,
            "sizing_parameters":
            {
                "boundary_layer_max_distance"          : 1.0,
                "interpolation"                        : "constant"
            },
            "enforce_current"                      : false,
            "anisotropy_remeshing"                 : false,
            "enforce_anisotropy_relative_variable" : false
        }
        """)
        metric_process = MeshingApplication.ComputeHessianSolMetricProcess(main_model_part, KratosMultiphysics.DISTANCE, MetricParameters)
        metric_process.Execute()

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "filename"                         : "mmg_lagrangian_test/remesh_rectangle",
            "save_external_files"              : true,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        mmg_parameters["filename"].SetString(file_path + "/" + mmg_parameters["filename"].GetString())
        mmg_process = MeshingApplication.MmgProcess2D(main_model_part, mmg_parameters)

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
                                                #"nodal_results"       : ["DISTANCE"],
                                                #"nodal_nonhistorical_results": ["METRIC_TENSOR_2D","AUXILIAR_GRADIENT","AUXILIAR_HESSIAN"]
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

        # We check the solution
        check_parameters = KratosMultiphysics.Parameters("""
        {
            "reference_file_name"   : "mmg_lagrangian_test/remesh_rectangle_result.sol",
            "output_file_name"      : "mmg_lagrangian_test/remesh_rectangle_step=0.sol",
            "dimension"             : 2,
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
