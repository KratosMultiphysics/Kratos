from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import sys
import time as timer
import os
import KratosMultiphysics as Kratos
from Kratos import Logger
import KratosMultiphysics.DEMApplication as Dem
sys.path.insert(0, '')
Logger.Print("Running under OpenMP........", label="DEM")
import DEM_procedures
import DEM_material_test_script
import KratosMultiphysics.StructuralMechanicsApplication as Structural

class Algorithm(object):

    def __init__(self):
        self.model = Kratos.Model()

        import dem_main_script_ready_for_coupling_with_fem
        self.dem_solution = dem_main_script_ready_for_coupling_with_fem.Solution()

        import structural_mechanics_analysis
        structural_parameters_file_name = "ProjectParameters.json"

        with open(structural_parameters_file_name,'r') as parameter_file:
            parameters = Kratos.Parameters(parameter_file.read())
        self.structural_solution = structural_mechanics_analysis.StructuralMechanicsAnalysis(self.model, parameters)


    def Run(self):
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def Initialize(self):
        self.dem_solution.Initialize()
        self.structural_solution.Initialize()

    def RunSolutionLoop(self):
        self.dem_solution.step = 0
        self.dem_solution.time = 0.0
        self.dem_solution.time_old_print = 0.0
        self.time_dem   = 0.0

        while self.structural_solution.time < self.structural_solution.end_time:
            self.structural_solution.time = self.structural_solution.solver.AdvanceInTime(self.structural_solution.time)
            self.structural_solution.InitializeSolutionStep()
            self.structural_solution.solver.Predict()
            self.structural_solution.solver.SolveSolutionStep()
            self.structural_solution.FinalizeSolutionStep()
            self.structural_solution.OutputSolutionStep()
            time_final_DEM_substepping = self.structural_solution.time

            self.Dt_DEM = self.dem_solution.spheres_model_part.ProcessInfo.GetValue(Kratos.DELTA_TIME)

            for self.dem_solution.time_dem in self.yield_DEM_time(self.time_dem, time_final_DEM_substepping, self.Dt_DEM):
                self.dem_solution.InitializeTimeStep()
                self.dem_solution.time = self.dem_solution.time + self.dem_solution.dt
                self.dem_solution.step += 1

                self.dem_solution.DEMFEMProcedures.UpdateTimeInModelParts(self.dem_solution.all_model_parts, self.dem_solution.time, self.dem_solution.dt, self.dem_solution.step)

                self.dem_solution.BeforeSolveOperations(self.dem_solution.time)

                self.dem_solution.SolverSolve()

                self.dem_solution.AfterSolveOperations()

                self.dem_solution.DEMFEMProcedures.MoveAllMeshes(self.dem_solution.all_model_parts, self.dem_solution.time, self.dem_solution.dt)
                #DEMFEMProcedures.MoveAllMeshesUsingATable(rigid_face_model_part, time, dt)

                ##### adding DEM elements by the inlet ######
                if self.dem_solution.DEM_parameters["dem_inlet_option"].GetBool():
                    self.dem_solution.DEM_inlet.CreateElementsFromInletMesh(self.dem_solution.spheres_model_part, self.dem_solution.cluster_model_part, self.dem_solution.creator_destructor)  # After solving, to make sure that neighbours are already set.

                stepinfo = self.dem_solution.report.StepiReport(timer, self.dem_solution.time, self.dem_solution.step)
                if stepinfo:
                    self.dem_solution.KRATOSprint(stepinfo)

                #### PRINTING GRAPHS ####
                os.chdir(self.dem_solution.graphs_path)
                self.dem_solution.post_utils.ComputeMeanVelocitiesInTrap("Average_Velocity.txt", self.dem_solution.time)

                self.dem_solution.materialTest.MeasureForcesAndPressure()
                self.dem_solution.materialTest.PrintGraph(self.dem_solution.time)

                self.dem_solution.DEMFEMProcedures.PrintGraph(self.dem_solution.time)
                self.dem_solution.DEMFEMProcedures.PrintBallsGraph(self.dem_solution.time)

                self.dem_solution.DEMEnergyCalculator.CalculateEnergyAndPlot(self.dem_solution.time)

                self.dem_solution.BeforePrintingOperations(self.dem_solution.time)

                #### GiD IO ##########################################
                if self.dem_solution.IsTimeToPrintPostProcess():
                    self.dem_solution.PrintResultsForGid(self.dem_solution.time)
                    self.dem_solution.time_old_print = self.dem_solution.time

                self.dem_solution.FinalizeTimeStep(self.dem_solution.time)


    def Finalize(self):
        self.dem_solution.Finalize()
        self.dem_solution.CleanUpOperations()
        self.structural_solution.Finalize()

    def yield_DEM_time(self, current_time, current_time_plus_increment, delta_time):
        current_time += delta_time

        tolerance = 0.0001
        while current_time < (current_time_plus_increment - tolerance * delta_time):
            yield current_time
            current_time += delta_time

        current_time = current_time_plus_increment
        yield current_time


if __name__ == "__main__":
    Algorithm().Run()