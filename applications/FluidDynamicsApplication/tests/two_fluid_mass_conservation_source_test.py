
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import math
import time
import sys

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities

class FluidDynamicsAnalysisMassConservation(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,test_inlet_case,flush_frequency=10.0):
        super().__init__(model,project_parameters)
        self.test_inlet_case = test_inlet_case
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def ModifyInitialGeometry(self):

        if self.test_inlet_case:

            h0=0.1
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                d0=node.X-h0
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,d0)
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,1,d0)
        else:
            h0=0.5
            hmax=0.1
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                d0=node.Y-(h0+hmax*math.cos(math.pi*(1-node.X)))
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,d0)
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,1,d0)

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

    # runs the two dimensional test case
    @UnitTest.expectedFailure #To be fixed in #9592
    def testTwoFluidMassConservationTest2D(self):
        self._has_inlet = False
        self._AuxiliaryRunTest("ProjectParameters2D.json")

    # runs the three dimensional test case
    @UnitTest.expectedFailure #To be fixed in #9592
    def testTwoFluidMassConservationTest3D(self):
        self._has_inlet = False
        self._AuxiliaryRunTest("ProjectParameter3D.json")

    # runs the three dimensional test case considering a inlet
    @UnitTest.expectedFailure #To be fixed in #9592
    def testTwoFluidMassConservationTest3DInlet(self):
        self._has_inlet = True
        self._AuxiliaryRunTest("ProjectParameters3DVariableInlet.json")

    def tearDown(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('.time')

    def _AuxiliaryRunTest(self, project_parameters_name):
        with UnitTest.WorkFolderScope(self.work_folder,__file__):
            with open(project_parameters_name,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
                model = KratosMultiphysics.Model()

                if self.parameters["problem_data"]["problem_name"].GetString()=="inlet_test":
                    self.test_inlet_case=True
                else:
                    self.test_inlet_case=False

                if self.print_output:
                    self._AddOutput()

                if self.print_reference_values:
                    self._AddReferenceValuesOutput()
                else:
                    self._AddReferenceValuesCheck()

                # running
                self.simulation = FluidDynamicsAnalysisMassConservation(model, self.parameters, self._has_inlet)
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
        output_name = "mass_conservation_{0}d_{1}".format(self.parameters["solver_settings"]["domain_size"].GetInt(),self.parameters["problem_data"]["problem_name"].GetString())
        gid_output_settings["Parameters"]["output_name"].SetString(output_name)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

    def _AddReferenceValuesOutput(self):
        dim = self.parameters["solver_settings"]["domain_size"].GetInt()
        name= self.parameters["problem_data"]["problem_name"].GetString()
        json_output_settings = KratosMultiphysics.Parameters("""{
            "name" : "Processes.KratosMultiphysics.JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["DISTANCE"],
                "output_file_name" : "TO_BE_DEFINED",
                "model_part_name"  : "FluidModelPart.FluidParts_Fluid",
                "time_frequency"   : 0.01
            }
        }""")
        output_file_name = "mass_conservation_{0}d_{1}_results.json".format(dim,name)
        json_output_settings["Parameters"]["output_file_name"].SetString(output_file_name)
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def _AddReferenceValuesCheck(self):
        dim = self.parameters["solver_settings"]["domain_size"].GetInt()
        name= self.parameters["problem_data"]["problem_name"].GetString()
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
        input_file_name = "mass_conservation_{0}d_{1}_results.json".format(dim,name)
        json_check_settings["Parameters"]["input_file_name"].SetString(input_file_name)
        json_check_settings["Parameters"]["tolerance"].SetDouble(self.check_absolute_tolerance)
        json_check_settings["Parameters"]["relative_tolerance"].SetDouble(self.check_relative_tolerance)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)


if __name__ == "__main__":
    UnitTest.main()