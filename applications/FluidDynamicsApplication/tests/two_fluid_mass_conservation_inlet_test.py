from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import math
import time
import os
import sys

import KratosMultiphysics.KratosUnittest as UnitTest

# Class to navigate through the folders
class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

# Class derived from the UnitTest (KratosMultiphysics.KratosUnittest) class
class TwoFluidMassConservationInletTest(UnitTest.TestCase):

    def __init__(self):
        self.work_folder = "TwoFluidMassInletTest"
        self.check_absolute_tolerance = 1.0e-7
        self.check_relative_tolerance = 1.0e-5
        self.print_output =False
        self.print_reference_values =False

    # runs the three dimensional test case
    def runTwoFluidMassConservationTestInlet3D(self):
        self._AuxiliaryRunTest("ProjectParameters3D.json")
    
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
                self.simulation = FluidDynamicsAnalysisChannel(model, self.parameters)
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
                "model_part_name"  : "FluidModelPart.FluidParts_Boundary-Fluid",
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
                "model_part_name"      : "FluidParts_Boundary-Fluid"",
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
class FluidDynamicsAnalysisChannel(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super().__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def ModifyInitialGeometry(self):

        h0=1
        
        for node in self._GetSolver().GetComputingModelPart().Nodes:            
            d0=node.X-h0
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,d0)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

if __name__ == "__main__":

    test = TwoFluidMassConservationInletTest()
    test.runTwoFluidMassConservationTestInlet3D()
