
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utils

try:
    import KratosMultiphysics.LinearSolversApplication
    have_external_solvers = True
except ImportError as e:
    have_external_solvers = False

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
class TwoFluidMassConservationTest(UnitTest.TestCase):

    def __init__(self):
        self.waterLevel = 1.0
        self.work_folder = "TwoFluidMassConservationProcTest"
        self.check_tolerance = 1e-7
        self.gravitationalAcceleration = 9.81
        # switch here for output
        self.print_output = False

    # runs the two dimensional test case
    def runTwoFluidMassConservationTest2D(self):
        with open("TwoFluidMassConservationProcTest/2Dtest.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()

            if self.print_output:
                parameters["output_processes"].AddValue("gid_output", KratosMultiphysics.Parameters(R'''[{
                    "python_module" : "gid_output_process",
                    "kratos_module" : "KratosMultiphysics",
                    "process_name"  : "GiDOutputProcess",
                    "help"          : "This process writes postprocessing files for GiD",
                    "Parameters"    : {
                        "model_part_name"        : "FluidModelPart",
                        "output_name"            : "FluidModelPart",
                        "postprocess_parameters" : {
                            "result_file_configuration" : {
                                "gidpost_flags"       : {
                                    "GiDPostMode"           : "GiD_PostBinary",
                                    "WriteDeformedMeshFlag" : "WriteDeformed",
                                    "WriteConditionsFlag"   : "WriteConditions",
                                    "MultiFileFlag"         : "SingleFile"
                                },
                                "file_label"          : "time",
                                "output_control_type" : "time",
                                "output_interval"     : 0.1,
                                "body_output"         : true,
                                "node_output"         : false,
                                "skin_output"         : false,
                                "plane_output"        : [],
                                "nodal_results"       : ["VELOCITY","PRESSURE","DISTANCE","DENSITY","DYNAMIC_VISCOSITY"],
                                "gauss_point_results" : []
                            },
                            "point_data_configuration"  : []
                        }
                    }
                }]''') )

            # running
            self.simulation = FluidDynamicsAnalysisWithFlush2D(model,parameters)
            self.simulation.Run()

            with open('TwoFluidMassConservationProcTest/mass_cons_2D.log', 'r') as file1, open('TwoFluidMassConservationProcTest/mass_cons_2D_referenceFile.ref', 'r') as file2:
                for line1, line2 in zip(file1, file2):
                    if line1 == line2:
                        # file are equal
                        self.assertAlmostEqual( 1.0, 1.0 )
                    else:
                        # deviation exist
                        self.assertAlmostEqual( 1.0, 0.0 )

            kratos_utils.DeleteFileIfExisting('TwoFluidMassConservationProcTest/mass_cons_2D.log')
            kratos_utils.DeleteFileIfExisting('TwoFluidMassConservationProcTest/2Dtest.time')


    # runs the three dimensional test case
    def runTwoFluidMassConservationTest3D(self):
        with open("TwoFluidMassConservationProcTest/3Dtest.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()

            if self.print_output:
                parameters["output_processes"].AddValue("gid_output", KratosMultiphysics.Parameters(R'''[{
                    "python_module" : "gid_output_process",
                    "kratos_module" : "KratosMultiphysics",
                    "process_name"  : "GiDOutputProcess",
                    "help"          : "This process writes postprocessing files for GiD",
                    "Parameters"    : {
                        "model_part_name"        : "FluidModelPart",
                        "output_name"            : "FluidModelPart",
                        "postprocess_parameters" : {
                            "result_file_configuration" : {
                                "gidpost_flags"       : {
                                    "GiDPostMode"           : "GiD_PostBinary",
                                    "WriteDeformedMeshFlag" : "WriteDeformed",
                                    "WriteConditionsFlag"   : "WriteConditions",
                                    "MultiFileFlag"         : "SingleFile"
                                },
                                "file_label"          : "time",
                                "output_control_type" : "time",
                                "output_interval"     : 0.1,
                                "body_output"         : true,
                                "node_output"         : false,
                                "skin_output"         : false,
                                "plane_output"        : [],
                                "nodal_results"       : ["VELOCITY","PRESSURE","DISTANCE","DENSITY","DYNAMIC_VISCOSITY"],
                                "gauss_point_results" : []
                            },
                            "point_data_configuration"  : []
                        }
                    }
                }]''') )

            # running
            self.simulation = FluidDynamicsAnalysisWithFlush3D(model,parameters)
            self.simulation.Run()

            with open('TwoFluidMassConservationProcTest/mass_cons_3D.log', 'r') as file1, open('TwoFluidMassConservationProcTest/mass_cons_3D_referenceFile.ref', 'r') as file2:
                for line1, line2 in zip(file1, file2):
                    if line1 == line2:
                        # file are equal
                        self.assertAlmostEqual( 1.0, 1.0 )
                    else:
                        # deviation exist
                        self.assertAlmostEqual( 1.0, 0.0 )

            kratos_utils.DeleteFileIfExisting('TwoFluidMassConservationProcTest/mass_cons_3D.log')
            kratos_utils.DeleteFileIfExisting('TwoFluidMassConservationProcTest/3Dtest.time')


class FluidDynamicsAnalysisWithFlush2D(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(FluidDynamicsAnalysisWithFlush2D,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def ModifyInitialGeometry(self):

        self.init_h = 1.5
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            distance = node.Y - self.init_h
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def ApplyBoundaryConditions(self):

        for node in self._GetSolver().GetComputingModelPart().Nodes:
            if ( node.Y > 1.8 ):
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.PRESSURE)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0.0)
            else:
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 1.0)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)

            node.Fix(KratosMultiphysics.DISTANCE)
            distance = node.Y - self.init_h
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def FinalizeSolutionStep(self):
        super(FluidDynamicsAnalysisWithFlush2D,self).FinalizeSolutionStep()
        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now


class FluidDynamicsAnalysisWithFlush3D(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(FluidDynamicsAnalysisWithFlush3D,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def ModifyInitialGeometry(self):

        self.init_h = 1.0
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            distance = node.Z - self.init_h
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def ApplyBoundaryConditions(self):

        for node in self._GetSolver().GetComputingModelPart().Nodes:
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 1.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, 0.0)
            node.Fix(KratosMultiphysics.DISTANCE)
            distance = node.Z - self.init_h
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def FinalizeSolutionStep(self):
        super(FluidDynamicsAnalysisWithFlush3D,self).FinalizeSolutionStep()
        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now


if __name__ == "__main__":

    test = TwoFluidMassConservationTest()

    test.runTwoFluidMassConservationTest2D()

    test.runTwoFluidMassConservationTest3D()