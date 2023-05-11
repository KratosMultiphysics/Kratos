from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import math
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
# import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
# import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff
# import KratosMultiphysics.MeshingApplication as MeshApp
import KratosMultiphysics.PfemMeltingApplication as PfemM

import time as timer


# Importing the base class


import KratosMultiphysics.PfemMeltingApplication.coupled_fluid_thermal_solverwithoutmeshgeneration
BaseClass = KratosMultiphysics.PfemMeltingApplication.coupled_fluid_thermal_solverwithoutmeshgeneration.PfemCoupledFluidThermalSolver

def CreateSolver(main_model_part, custom_settings):

    return FemtolaserSolver(main_model_part, custom_settings)

class FemtolaserSolver(BaseClass):
   def SolveSolutionStep(self):

        for node in self.fluid_solver.main_model_part.Nodes:
            T = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            if False and node.Z > 0.0049:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 400)
                node.Fix(KratosMultiphysics.TEMPERATURE)

        t7=timer.time()

        # fluid_is_converged = self.fluid_solver.SolveSolutionStep()
        #self.Streamline.RungeKutta4ElementbasedSI(self.fluid_solver.main_model_part,100)

        for node in self.fluid_solver.main_model_part.Nodes:
            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY, velocity)

        self.HeatSource.Heat_Source(self.fluid_solver.main_model_part) #heat source for the thermal problem

        thermal_is_converged = self.thermal_solver.SolveSolutionStep()
        # self.CalculateViscosityaux()
        t8=timer.time()
        self.problemsolution=self.problemsolution + t8 - t7
        step=self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        self.outputfile6.write(str(step)+" "+ str(self.streamlineintegration) +" "+ str(self.meshingprocedure) +" "+ str(self.fillingsubmodelparts) +" "+ str(self.initializeSolutionStep)+" "+ str(self.problemsolution) +"\n")

        return thermal_is_converged