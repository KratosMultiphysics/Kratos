from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import math
import time
import os
import sys

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities

class FluidDynamicsAnalysisSloshingTank(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super().__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def ModifyInitialGeometry(self):

        h0=0.5
        hmax=0.1
        for node in self._GetSolver().GetComputingModelPart().Nodes:            
            d0=node.Y-(h0+hmax*math.cos(math.pi*(1-node.X)))
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,d0)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

# Class derived from the UnitTest (KratosMultiphysics.KratosUnittest) class
class TwoFluidMassConservationTest(UnitTest.TestCase):

    def setUp(self):
        self.work_folder = "TwoFluidMassConservationSourceTest"
        self.check_absolute_tolerance = 1.0e-7
        self.check_relative_tolerance = 1.0e-5
        self.print_output = False
        self.print_reference_values = False

    # runs the two dimensinal test case
    def testTwoFluidMassConservationTest2D(self):
        self._AuxiliaryRunTest("ProjectParameters2D.json")

    # runs the three dimensional test case
    def testTwoFluidMassConservationTest3D(self):
        self._AuxiliaryRunTest("ProjectParameter3D.json")

    def tearDown(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('.time')

    def _AuxiliaryRunTest(self, project_parameters_name):
        with UnitTest.WorkFolderScope(self.work_folder,__file__):
            with open(project_parameters_name,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
                model = KratosMultiphysics.Model()

                if self.print_output:
                    self._AddOutput()

                if self.print_reference_values:
                    self._AddReferenceValuesOutput()
                else:
                    self._AddReferenceValuesCheck()

                # running
                self.simulation = FluidDynamicsAnalysisSloshingTank(model, self.parameters)
                self.simulation.Run()

    def _AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "FluidModelPart",
                "output_name"            : "TO_BE_DEFINED",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "step",
                        "output_control_type"         : "step",
                        "output_interval"             : 1.0,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["DISTANCE","VELOCITY","PRESSURE"]
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        output_name = "mass_conservation_{0}d".format(self.parameters["solver_settings"]["domain_size"].GetInt())
        gid_output_settings["Parameters"]["output_name"].SetString(output_name)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

    def _AddReferenceValuesOutput(self):
        dim = self.parameters["solver_settings"]["domain_size"].GetInt()
        json_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "json_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["DISTANCE"],
                "output_file_name" : "TO_BE_DEFINED",
                "model_part_name"  : "FluidModelPart.FluidParts_Fluid",
                "time_frequency"   : 0.01
            }
        }""")
        output_file_name = "mass_conservation_{0}d_results.json".format(dim)
        json_output_settings["Parameters"]["output_file_name"].SetString(output_file_name)
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def _AddReferenceValuesCheck(self):
        dim = self.parameters["solver_settings"]["domain_size"].GetInt()
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : ["DISTANCE"],
                "input_file_name"      : "TO_BE_DEFINED",
                "model_part_name"      : "FluidModelPart.FluidParts_Fluid",
                "tolerance"            : 0.0,
                "relative_tolerance"   : 0.0,
                "time_frequency"       : 0.01
            }
        }""")
        input_file_name = "mass_conservation_{0}d_results.json".format(dim)
        json_check_settings["Parameters"]["input_file_name"].SetString(input_file_name)
        json_check_settings["Parameters"]["tolerance"].SetDouble(self.check_absolute_tolerance)
        json_check_settings["Parameters"]["relative_tolerance"].SetDouble(self.check_relative_tolerance)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)


if __name__ == "__main__":
    UnitTest.main()
    test=TwoFluidMassConservationTest()
    test.testTwoFluidMassConservationTest2D()
    test.testTwoFluidMassConservationTest3D()