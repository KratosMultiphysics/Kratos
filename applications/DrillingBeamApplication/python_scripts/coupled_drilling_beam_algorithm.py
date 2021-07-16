from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import sys
import time as timer
import os
import weakref
import KratosMultiphysics as Kratos
from Kratos import Logger
import KratosMultiphysics.DEMApplication as Dem
sys.path.insert(0, '')
Logger.Print("Running under OpenMP........", label="DEM")
import KratosMultiphysics.DEMApplication.DEM_procedures as DEM_procedures
import KratosMultiphysics.DEMApplication.DEM_material_test_script as DEM_material_test_script
import KratosMultiphysics.StructuralMechanicsApplication as Structural
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem
import KratosMultiphysics.DemStructuresCouplingApplication.dem_structures_coupling_gid_output as dem_structures_coupling_gid_output
import KratosMultiphysics.DemStructuresCouplingApplication.dem_fem_coupling_algorithm as dem_fem_coupling_algorithm
import KratosMultiphysics.DrillingBeamApplication as DBA


class CoupledDrillingBeamAlgorithm(dem_fem_coupling_algorithm.Algorithm):

    def Initialize(self):
        self.structural_solution.Initialize()    # Reading mdpa
        self.structural_main_model_part = self.structural_solution._GetSolver().main_model_part

        self.dem_solution.Initialize()    # Adding DEM variables and reading

        self._DetectStructuresSkin()
        self._TransferStructuresSkinToDem()
        self.dem_solution._GetSolver().Initialize()

        self.sandwich_simulation = False

        self.beam_solids_utility = DBA.BeamSolidsUtility(self.structural_main_model_part)

        mixed_mp = self.model.CreateModelPart('MixedPart')
        filename = os.path.join(self.dem_solution.post_path, self.dem_solution.DEM_parameters["problem_name"].GetString())
        self.gid_output = dem_structures_coupling_gid_output.DemStructuresCouplingGiDOutput(
                            filename,
                            True,
                            "Binary",
                            "Multiples",
                            True,
                            True,
                            self.structural_solution._GetSolver().GetComputingModelPart(),
                            self.dem_solution.spheres_model_part,
                            self.dem_solution.cluster_model_part,
                            self.dem_solution.rigid_face_model_part,
                            self.dem_solution.contact_model_part,
                            mixed_mp
                            )

        structures_nodal_results = ["DEM_SURFACE_LOAD", "REACTION"]
        dem_nodal_results = ["TOTAL_FORCES"]
        clusters_nodal_results = []
        rigid_faces_nodal_results = []
        contact_model_part_results = []
        mixed_nodal_results = ["DISPLACEMENT"]
        gauss_points_results = ["CAUCHY_STRESS_TENSOR"]
        self.gid_output.initialize_dem_fem_results(structures_nodal_results,
                                                   dem_nodal_results,
                                                   clusters_nodal_results,
                                                   rigid_faces_nodal_results,
                                                   contact_model_part_results,
                                                   mixed_nodal_results,
                                                   gauss_points_results)




    def RunSolutionLoop(self):

        self.dem_solution.step = 0
        self.dem_solution.time = 0.0
        self.dem_solution.time_old_print = 0.0
        self.time_dem = 0.0
        self.Dt_structural = self.structural_solution._GetSolver().settings["time_stepping"]["time_step"].GetDouble()

        while self.structural_solution.time < self.structural_solution.end_time:

            portion_of_the_force_which_is_new = 1.0
            DemFem.DemStructuresCouplingUtilities().SmoothLoadTrasferredToFem(self.dem_solution.rigid_face_model_part, portion_of_the_force_which_is_new)

            self.structural_solution.time = self.structural_solution._GetSolver().AdvanceInTime(self.structural_solution.time)

            self.structural_solution.InitializeSolutionStep()
            self.structural_solution._GetSolver().Predict()
            self.structural_solution._GetSolver().SolveSolutionStep()
            self.structural_solution.FinalizeSolutionStep()

            self.structural_solution.OutputSolutionStep()

            time_final_DEM_substepping = self.structural_solution.time

            self.Dt_DEM = self.dem_solution.spheres_model_part.ProcessInfo.GetValue(Kratos.DELTA_TIME)

            DemFem.InterpolateStructuralSolutionForDEM().SaveStructuralSolution(self.structural_mp)

            DemFem.ComputeDEMFaceLoadUtility().ClearDEMFaceLoads(self.skin_mp)

            for self.dem_solution.time_dem in self.yield_DEM_time(self.dem_solution.time, time_final_DEM_substepping, self.Dt_DEM):
                self.dem_solution.InitializeSolutionStep()
                self.dem_solution.time = self.dem_solution.time + self.dem_solution._GetSolver().dt

                self.dem_solution.step += 1

                self.dem_solution.DEMFEMProcedures.UpdateTimeInModelParts(self.dem_solution.all_model_parts, self.dem_solution.time, self.dem_solution._GetSolver().dt, self.dem_solution.step)

                DemFem.InterpolateStructuralSolutionForDEM().InterpolateStructuralSolution(self.structural_mp, self.Dt_structural, self.structural_solution.time, self.dem_solution._GetSolver().dt, self.dem_solution.time)

                self.dem_solution.SolverSolve()

                DemFem.ComputeDEMFaceLoadUtility().CalculateDEMFaceLoads(self.skin_mp, self.dem_solution._GetSolver().dt, self.Dt_structural)

                ##### adding DEM elements by the inlet ######
                if self.dem_solution.DEM_parameters["dem_inlet_option"].GetBool():
                    self.dem_solution.DEM_inlet.CreateElementsFromInletMesh(self.dem_solution.spheres_model_part, self.dem_solution.cluster_model_part, self.dem_solution.creator_destructor)  # After solving, to make sure that neighbours are already set.

                stepinfo = self.dem_solution.report.StepiReport(timer, self.dem_solution.time, self.dem_solution.step)
                if stepinfo:
                    self.dem_solution.KratosPrintInfo(stepinfo)

                #### PRINTING GRAPHS ####
                os.chdir(self.dem_solution.graphs_path)
                self.dem_solution.post_utils.ComputeMeanVelocitiesInTrap("Average_Velocity.txt", self.dem_solution.time, self.dem_solution.graphs_path)

                self.dem_solution.materialTest.MeasureForcesAndPressure()
                self.dem_solution.materialTest.PrintGraph(self.dem_solution.time)

                self.dem_solution.DEMFEMProcedures.PrintGraph(self.dem_solution.time)
                self.dem_solution.DEMFEMProcedures.PrintBallsGraph(self.dem_solution.time)

                self.dem_solution.DEMEnergyCalculator.CalculateEnergyAndPlot(self.dem_solution.time)

                self.dem_solution.BeforePrintingOperations(self.dem_solution.time)

                #### GiD IO ##########################################
                if self.dem_solution.IsTimeToPrintPostProcess():
                    self.dem_solution._GetSolver().PrepareElementsForPrinting()
                    if self.dem_solution.DEM_parameters["ContactMeshOption"].GetBool():
                        self.dem_solution._GetSolver().PrepareContactElementsForPrinting()
                    self.dem_solution.PrintResultsForGid(self.dem_solution.time)
                    self.dem_solution.demio.PrintMultifileLists(self.dem_solution.time, self.dem_solution.post_path)
                    self.dem_solution.time_old_print = self.dem_solution.time

                self.dem_solution.FinalizeSolutionStep()

            DemFem.InterpolateStructuralSolutionForDEM().RestoreStructuralSolution(self.structural_mp)

            self.beam_solids_utility.ComputeEdgePointOfBeamSolids(self.structural_main_model_part)
            self.beam_solids_utility.ComputeDirectiveLineOfBeamSolids(self.structural_main_model_part)
            self.beam_solids_utility.ComputeReactionsOfBeamSolids(self.structural_main_model_part)

if __name__ == "__main__":
    CoupledDrillingBeamAlgorithm().Run()