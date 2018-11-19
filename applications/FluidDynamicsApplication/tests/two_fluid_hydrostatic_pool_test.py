from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
from fluid_dynamics_analysis import FluidDynamicsAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utils

try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True

import sys
import time
import os

import KratosMultiphysics.KratosUnittest as UnitTest

# Class to navigate through the folders
class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

# Class derived from the UnitTest (KratosMultiphysics.KratosUnittest) class
class TwoFluidHydrostaticPoolTest(UnitTest.TestCase):

    def __init__(self):
        self.waterLevel = 0.5
        self.work_folder = "TwoFluidStaticPoolTest"
        self.settings = "TwoFluidStaticPoolTest2D.json"
        self.check_tolerance = 1e-10
        self.check_toleranceDistance = 0.03
        self.gravitationalAcceleration = 9.81
        self.domainHeight = 1.0
        self.rho1 = 1000.0
        self.rho2 = 1.0

    # runs the two dimensinal test case
    def runTwoFluidHydrostaticTest2D(self):
        with open("TwoFluidStaticPoolTest/TwoFluidStaticPoolTest2D.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            # running
            self.simulation = FluidDynamicsAnalysisWithFlush2D(model,parameters)
            self.simulation.Run()

            # testing
            for node in self.simulation._GetSolver().GetComputingModelPart().Nodes:
                velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                self.assertAlmostEqual(0.0, velocity[0], delta = self.check_tolerance)
                self.assertAlmostEqual(0.0, velocity[1], delta = self.check_tolerance)
                self.assertAlmostEqual(0.0, velocity[2], delta = self.check_tolerance)

                pressure = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
                if node.Y > self.waterLevel:
                    pressureAnalytic = (self.domainHeight-node.Y)*self.gravitationalAcceleration*self.rho2
                else:
                    pressureAnalytic = (self.domainHeight-self.waterLevel)*self.gravitationalAcceleration*self.rho2
                    pressureAnalytic += (self.waterLevel-node.Y)*self.gravitationalAcceleration*self.rho1
                self.assertAlmostEqual(pressureAnalytic, pressure, delta = self.check_tolerance)

                distance = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                distanceAnalytic = (node.Y - self.waterLevel)
                self.assertAlmostEqual(distanceAnalytic, distance, delta = self.check_toleranceDistance)

            kratos_utils.DeleteFileIfExisting('TwoFluidStaticPoolTest2D.post.bin')
            kratos_utils.DeleteFileIfExisting('tests.post.lst')

    # runs the three dimensional test case
    def runTwoFluidHydrostaticTest3D(self):
        with open("TwoFluidStaticPoolTest/TwoFluidStaticPoolTest3D.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            # running
            self.simulation = FluidDynamicsAnalysisWithFlush3D(model,parameters)
            self.simulation.Run()

            # testing
            for node in self.simulation._GetSolver().GetComputingModelPart().Nodes:
                velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                self.assertAlmostEqual(0.0, velocity[0], delta = self.check_tolerance)
                self.assertAlmostEqual(0.0, velocity[1], delta = self.check_tolerance)
                self.assertAlmostEqual(0.0, velocity[2], delta = self.check_tolerance)

                pressure = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
                if node.Z > self.waterLevel:
                    pressureAnalytic = (self.domainHeight-node.Z)*self.gravitationalAcceleration*self.rho2
                else:
                    pressureAnalytic = (self.domainHeight-self.waterLevel)*self.gravitationalAcceleration*self.rho2
                    pressureAnalytic += (self.waterLevel-node.Z)*self.gravitationalAcceleration*self.rho1
                self.assertAlmostEqual(pressureAnalytic, pressure, delta = self.check_tolerance)

                distance = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                distanceAnalytic = (node.Z - self.waterLevel)
                self.assertAlmostEqual(distanceAnalytic, distance, delta = self.check_toleranceDistance)
            
            kratos_utils.DeleteFileIfExisting('TwoFluidStaticPoolTest3D.post.bin')
            kratos_utils.DeleteFileIfExisting('tests.post.lst')


class FluidDynamicsAnalysisWithFlush2D(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(FluidDynamicsAnalysisWithFlush2D,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def ModifyInitialGeometry(self):
        
        init_h = 0.5
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            distance = node.Y - init_h
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def ApplyBoundaryConditions(self):

        v_zero = KratosMultiphysics.Vector(3,0.0)
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            if abs(node.X) > 0.499 and abs(node.X) < 0.501:
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)
                # node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, v_zero)
            if abs(node.Y) > -0.001 and abs(node.Y) < 0.001:
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
            	# node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, v_zero)
            if abs(node.Y) > 0.999 and abs(node.Y) < 1.001:
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
            	# node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, v_zero)

        nFix = 1
        v_zero = KratosMultiphysics.Vector(3,0.0)
        self._GetSolver().GetComputingModelPart().GetNode(nFix).Fix(KratosMultiphysics.VELOCITY_X)
        self._GetSolver().GetComputingModelPart().GetNode(nFix).Fix(KratosMultiphysics.VELOCITY_Y)
        self._GetSolver().GetComputingModelPart().GetNode(nFix).SetSolutionStepValue(KratosMultiphysics.VELOCITY, v_zero)

        self._GetSolver().GetComputingModelPart().GetNode(nFix).Fix(KratosMultiphysics.PRESSURE)
        self._GetSolver().GetComputingModelPart().GetNode(nFix).SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0.0)

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
        
        init_h = 0.5
        zero_vect = KratosMultiphysics.Vector(3,0.0)
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            distance = node.Z - init_h
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)

    def ApplyBoundaryConditions(self):

        # print("This function is executed")    
        v_zero = KratosMultiphysics.Vector(3,0.0)
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            if abs(node.X) > 0.499 and abs(node.X) < 0.501:
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)
                # node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, v_zero)
            if abs(node.Y) > 0.499 and abs(node.Y) < 0.501:
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
                # node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, v_zero)
            if abs(node.Z) < 0.001 or abs(node.Z) > 0.999:
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, 0.0)
                # node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, v_zero)

        nFix = 1726
        v_zero = KratosMultiphysics.Vector(3,0.0)
        self._GetSolver().GetComputingModelPart().GetNode(nFix).Fix(KratosMultiphysics.VELOCITY_X)
        self._GetSolver().GetComputingModelPart().GetNode(nFix).Fix(KratosMultiphysics.VELOCITY_Y)
        self._GetSolver().GetComputingModelPart().GetNode(nFix).Fix(KratosMultiphysics.VELOCITY_Z)
        self._GetSolver().GetComputingModelPart().GetNode(nFix).SetSolutionStepValue(KratosMultiphysics.VELOCITY, v_zero)

        self._GetSolver().GetComputingModelPart().GetNode(nFix).Fix(KratosMultiphysics.PRESSURE)
        self._GetSolver().GetComputingModelPart().GetNode(nFix).SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0.0)

    def FinalizeSolutionStep(self):
        super(FluidDynamicsAnalysisWithFlush3D,self).FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now


if __name__ == "__main__":

    test = TwoFluidHydrostaticPoolTest()
    test.runTwoFluidHydrostaticTest2D()
    test.runTwoFluidHydrostaticTest3D()




