from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import time as timer
import os
import KratosMultiphysics as Kratos
from Kratos import Logger
import KratosMultiphysics.DEMApplication as Dem
Logger.Print("Running under OpenMP........", label="DEM")
import KratosMultiphysics.StructuralMechanicsApplication as Structural
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem
from KratosMultiphysics.DemStructuresCouplingApplication.dem_fem_coupling_algorithm import Algorithm

class SPAlgorithm(Algorithm):

    def __init__(self):
        super(SPAlgorithm,self).__init__()

        sp_parameters_file_name = "SPParameters.json"

        with open(sp_parameters_file_name,'r') as parameter_file:
            self.sp_parameters = Kratos.Parameters(parameter_file.read())

        self.ValidateSettings()

        self.test_number = self.sp_parameters["problem_data"]["test_number"].GetInt()
        # Test types (4 different options):
        # Test number 0: no test simulation
        # Test number 1: CTW16 specimen
        # Test number 2: CTW10 specimen
        # Test number 3: Blind test specimen

        self.post_process_step_count = 0
        self.post_process_frequency = self.sp_parameters["post_process_tool"]["output_frequency"].GetInt()
        self.post_process_write_count = self.post_process_frequency

    @classmethod
    def GetDefaultSettings(cls):
        """This function returns the default-settings used by this class
        """
        return Kratos.Parameters("""{
            "problem_data"     : {
                "test_number"  : 1,
                "center" : [0.0,0.0,0.0],
                "axis"    : [0.0,0.0,1.0]
            },
            "post_process_tool":{
                "output_frequency": 0
            }
        }""")

    def ValidateSettings(self):
        """This function validates the settings of the solver
        """
        default_settings = self.GetDefaultSettings()
        self.sp_parameters.ValidateAndAssignDefaults(default_settings)

    def Initialize(self):
        super(SPAlgorithm,self).Initialize()

        self.InitializeAdditionalProcessInfoVars()

        if self.test_number:
            from KratosMultiphysics.DemStructuresCouplingApplication.control_module_fem_dem_utility import ControlModuleFemDemUtility
            self.control_module_fem_dem_utility = ControlModuleFemDemUtility(self.model, self.dem_solution.spheres_model_part, self.test_number)
            self.control_module_fem_dem_utility.ExecuteInitialize()

        # Create Postprocess tool for SP
        from KratosMultiphysics.DemStructuresCouplingApplication.sand_production_post_process_tool import SandProductionPostProcessTool
        self.sp_post_process_tool = SandProductionPostProcessTool(self.structural_solution._GetSolver().GetComputingModelPart(),
                                                                                                    self.dem_solution.spheres_model_part,
                                                                                                    self.test_number)

        from KratosMultiphysics.DemStructuresCouplingApplication import stress_failure_check_utility
        self.stress_failure_check_utility = stress_failure_check_utility.StressFailureCheckUtility(self.dem_solution.spheres_model_part, self.test_number)

    def InitializeAdditionalProcessInfoVars(self):
        self.dem_solution.spheres_model_part.ProcessInfo.SetValue(Dem.SIGMA_3_AVERAGE, 0.0)

    def RunSolutionLoop(self):

        self.dem_solution.step = 0
        self.dem_solution.time = 0.0
        self.dem_solution.time_old_print = 0.0
        self.time_dem = 0.0
        self.Dt_structural = self.structural_solution._GetSolver().settings["time_stepping"]["time_step"].GetDouble()

        while self.structural_solution.time < self.structural_solution.end_time:

            portion_of_the_force_which_is_new = 0.4
            DemFem.DemStructuresCouplingUtilities().SmoothLoadTrasferredToFem(self.dem_solution.rigid_face_model_part, portion_of_the_force_which_is_new)

            self.structural_solution.time = self.structural_solution._GetSolver().AdvanceInTime(self.structural_solution.time)

            self.structural_solution.InitializeSolutionStep()
            if self.test_number:
                self.control_module_fem_dem_utility.ExecuteInitializeSolutionStep()
            self.structural_solution._GetSolver().Predict()
            self.structural_solution._GetSolver().SolveSolutionStep()
            self.structural_solution.FinalizeSolutionStep()
            self.structural_solution.OutputSolutionStep()

            time_final_DEM_substepping = self.structural_solution.time

            self.Dt_DEM = self.dem_solution.spheres_model_part.ProcessInfo.GetValue(Kratos.DELTA_TIME)

            DemFem.InterpolateStructuralSolutionForDEM().SaveStructuralSolution(self.structural_mp)

            DemFem.ComputeDEMFaceLoadUtility().ClearDEMFaceLoads(self.skin_mp)

            if self.test_number == 1 or self.test_number == 2:
                self.outer_walls_model_part = self.model["Structure.SurfacePressure3D_lateral_pressure"]
                #DemFem.DemStructuresCouplingUtilities().ComputeSandProductionWithDepthFirstSearch(self.dem_solution.spheres_model_part, self.outer_walls_model_part, self.structural_solution.time)
                DemFem.DemStructuresCouplingUtilities().ComputeSandProductionWithDepthFirstSearchNonRecursiveImplementation(self.dem_solution.spheres_model_part, self.outer_walls_model_part, self.structural_solution.time)
                DemFem.DemStructuresCouplingUtilities().ComputeSandProduction(self.dem_solution.spheres_model_part, self.outer_walls_model_part, self.structural_solution.time)
            elif self.test_number == 3:
                self.outer_walls_model_part_1 = self.model["Structure.SurfacePressure3D_sigmaXpos"]
                self.outer_walls_model_part_2 = self.model["Structure.SurfacePressure3D_sigmaYpos"]
                DemFem.DemStructuresCouplingUtilities().ComputeTriaxialSandProduction(self.dem_solution.spheres_model_part, self.outer_walls_model_part_1, self.outer_walls_model_part_2, self.structural_solution.time)

            for self.dem_solution.time_dem in self.yield_DEM_time(self.dem_solution.time, time_final_DEM_substepping, self.Dt_DEM):
                self.dem_solution.time = self.dem_solution.time + self.dem_solution._GetSolver().dt

                self.dem_solution.step += 1

                self.dem_solution.DEMFEMProcedures.UpdateTimeInModelParts(self.dem_solution.all_model_parts, self.dem_solution.time, self.dem_solution._GetSolver().dt, self.dem_solution.step)

                self.dem_solution.InitializeSolutionStep()

                self.dem_solution._GetSolver().Predict()

                DemFem.InterpolateStructuralSolutionForDEM().InterpolateStructuralSolution(self.structural_mp, self.Dt_structural, self.structural_solution.time, self.dem_solution._GetSolver().dt, self.dem_solution.time)

                self.dem_solution.SolverSolve()

                DemFem.DemStructuresCouplingUtilities().MarkBrokenSpheres(self.dem_solution.spheres_model_part)

                center = Kratos.Array3()
                center[0] = self.sp_parameters["problem_data"]["center"][0].GetDouble()
                center[1] = self.sp_parameters["problem_data"]["center"][1].GetDouble()
                center[2] = self.sp_parameters["problem_data"]["center"][2].GetDouble()
                axis = Kratos.Array3()
                axis[0] = self.sp_parameters["problem_data"]["axis"][0].GetDouble()
                axis[1] = self.sp_parameters["problem_data"]["axis"][1].GetDouble()
                axis[2] = self.sp_parameters["problem_data"]["axis"][2].GetDouble()

                radius = 0
                if self.test_number == 1:
                    radius = 0.0036195; #95% of the real hole. CTW16 specimen
                elif self.test_number == 2:
                    radius = 0.012065; #95% of the real hole. CTW10 specimen
                elif self.test_number == 3:
                    radius = 0.036195; #95% of the real hole. Blind Test

                self.dem_solution.creator_destructor.MarkParticlesForErasingGivenCylinder(self.dem_solution.spheres_model_part, center, axis, radius)

                self.dem_solution.FinalizeSolutionStep()

                DemFem.ComputeDEMFaceLoadUtility().CalculateDEMFaceLoads(self.skin_mp, self.dem_solution._GetSolver().dt, self.Dt_structural)

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

                    if self.test_number:
                        self.stress_failure_check_utility.ExecuteFinalizeSolutionStep()

                self.dem_solution.FinalizeTimeStep(self.dem_solution.time)

            DemFem.InterpolateStructuralSolutionForDEM().RestoreStructuralSolution(self.structural_mp)

            if self.test_number:
                self.control_module_fem_dem_utility.ExecuteFinalizeSolutionStep()

            # Write SP data
            if self.IsPostProcessWriteStep():
                self.sp_post_process_tool.WriteData()

    def IsPostProcessWriteStep(self):
        self.post_process_step_count += 1
        if self.post_process_step_count == self.post_process_write_count:
            self.post_process_write_count += self.post_process_frequency
            return True
        else:
            return False

if __name__ == "__main__":
    SPAlgorithm().Run()
