
import KratosMultiphysics as KM
import KratosMultiphysics.FemToDemApplication as FEMDEM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem_for_PFEM_coupling as MainCouplingFemDem_for_PFEM_coupling
import KratosMultiphysics.FemToDemApplication.MainPFEM_for_coupling as MainPFEM_for_coupling
import os
import KratosMultiphysics.FemToDemApplication.fem_dem_coupled_gid_output as gid_output

def Wait():
    input("PFEM-FEMDEM -> Press Something")

def KratosPrintInfo(message):
    """This function prints info on screen
    """
    KM.Logger.Print("", message)
    KM.Logger.Flush()
#============================================================================================================================
class MainCouplingPfemFemDem_Solution:
#============================================================================================================================

    def __init__(self, Model, PFEMparameters, path = ""):
        # Initialize solutions of the FEMDEM and PFEM
        self.model = Model
        self.FEMDEM_Solution = MainCouplingFemDem_for_PFEM_coupling.MainCoupledFemDem_for_PFEM_coupling_Solution(Model, path)
        self.FEMDEM_Solution.is_slave = True
        self.FEMDEM_Solution.Initialize()

        self.PFEM_Solution = MainPFEM_for_coupling.MainPFEM_for_coupling_solution(Model, 
                                                                                  self.FEMDEM_Solution.FEM_Solution.main_model_part,
                                                                                  PFEMparameters)
        KratosPrintInfo("    ___                  _          _            _ _   _         ___  ___  __       "  + "\n" +
                       "    / __\___  _   _ _ __ | | ___  __| | __      _(_) |_| |__     / _ \/ __\/__\/\/\   " + "\n" +
                       "   / /  / _ \| | | | '_ \| |/ _ \/ _` | \ \ /\ / / | __| '_ \   / /_)/ _\ /_\ /    \  " + "\n" +
                       "  / /__| (_) | |_| | |_) | |  __/ (_| |  \ V  V /| | |_| | | | / ___/ /  //__/ /\/\ \ " + "\n" +
                       "  \____/\___/ \__,_| .__/|_|\___|\__,_|   \_/\_/ |_|\__|_| |_| \/   \/   \__/\/    \/ " + "\n")

#============================================================================================================================
    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

#============================================================================================================================
    def Initialize(self):
        # Now we storage the FEM modelpart in the PFEM_Solution
        self.PFEM_Solution.Initialize()

        if self.FEMDEM_Solution.CreateInitialSkin:
            self.FEMDEM_Solution.ComputeSkinSubModelPart()
            if self.FEMDEM_Solution.DEMFEM_contact:
                self.FEMDEM_Solution.TransferFEMSkinToDEM()
            FEMDEM.GenerateInitialSkinDEMProcess(self.FEMDEM_Solution.FEM_Solution.main_model_part, self.FEMDEM_Solution.SpheresModelPart).Execute()

        # We copy the output params in the PFEM
        self.PFEM_Solution.graphical_output = self.PFEM_Solution.SetCustomGraphicalOutput(self.FEMDEM_Solution.FEM_Solution.ProjectParameters)
        self.PFEM_Solution.GraphicalOutputExecuteInitialize()

        # Initialize the coupled post process
        self.InitializePostProcess()

#============================================================================================================================
    def RunMainTemporalLoop(self):
        self.RunBeforeSolutionLoopFEMDEM()
        while self.FEMDEM_Solution.FEM_Solution.time <= self.FEMDEM_Solution.FEM_Solution.end_time:
            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()

#============================================================================================================================
    def SolveSolutionStep(self):

        # It's necessary to Fix in order to maintain the FEMDEM Kinematics
        self.FixNodesModelPart(self.FEMDEM_Solution.FEM_Solution.main_model_part)

        KratosPrintInfo("==============================================" + "\n" +
                        " ==== SOLVING PFEM PART OF THE CALCULATION ====" + "\n" +
                        " ==============================================")
        self.SolveSolutionStepPFEM()

        # Now we Free the nodes to be calculated by the FEMDEM
        self.FreeNodesModelPart(self.FEMDEM_Solution.FEM_Solution.main_model_part)

        # Transfer pressure forces
        self.RegenerateAndUpdatePFEMPressureConditions()

        KratosPrintInfo("================================================" + "\n" +
                       " ==== SOLVING FEMDEM PART OF THE CALCULATION ====" + "\n" +
                       " ================================================")
        self.SolveSolutionStepFEMDEM()
        self.PFEM_Solution.main_model_part.RemoveNodes(KM.TO_ERASE)

        # We update the NO_MESH flag in the FEMDEM skin
        self.UpdateFEMDEMBoundary()


#============================================================================================================================
    def SolveSolutionStepPFEM(self):
        is_converged = self.PFEM_Solution._GetSolver().SolveSolutionStep()
        # self.PFEM_Solution.__CheckIfSolveSolutionStepReturnsAValue(is_converged)

#============================================================================================================================
    def SolveSolutionStepFEMDEM(self):
        self.FEMDEM_Solution.SolveSolutionStep()

#============================================================================================================================
    def RunBeforeSolutionLoopFEMDEM(self):
        # Solving the problem (time integration)
        self.FEMDEM_Solution.DEM_Solution.step           = 0
        self.FEMDEM_Solution.DEM_Solution.time           = 0.0
        self.FEMDEM_Solution.DEM_Solution.time_old_print = 0.0

        if self.FEMDEM_Solution.DoRemeshing:
            self.FEMDEM_Solution.RemeshingProcessMMG.ExecuteBeforeSolutionLoop()

#============================================================================================================================
    def FinalizeSolutionStep(self):
        self.PFEM_Solution.FinalizeSolutionStep()
        self.PFEM_Solution.OutputSolutionStep()
        self.FEMDEM_Solution.FinalizeSolutionStep()
        KM.PfemFluidDynamicsApplication.PostProcessUtilities().RebuildPostProcessModelPart(self.PFEM_Solution.post_process_model_part, self.PFEM_Solution.main_model_part)
        self.PrintResults()

#============================================================================================================================
    def Finalize(self):
        self.FEMDEM_Solution.Finalize()
        self.PFEM_Solution.Finalize()

#============================================================================================================================
    def InitializeSolutionStep(self):
        KratosPrintInfo("")
        KratosPrintInfo("_____________________________________________________________________________________")
        KratosPrintInfo("_____________________________________________________________________________________")
        KratosPrintInfo("")
        self.UpdateDeltaTimeInSolutions()
        self.PFEM_Solution.time = self.PFEM_Solution._GetSolver().AdvanceInTime(self.PFEM_Solution.time)
        self.PFEM_Solution.InitializeSolutionStep()
        self.PFEM_Solution._GetSolver().Predict()
        self.FEMDEM_Solution.InitializeSolutionStep()

        # We set the same delta time to all the solutions
        self.UpdateDeltaTimeInSolutions()

#============================================================================================================================
    def UpdateDeltaTimeInSolutions(self):
        FEM_delta_time = self.FEMDEM_Solution.FEM_Solution.delta_time
        self.PFEM_Solution.delta_time = FEM_delta_time
        self.PFEM_Solution.main_model_part.ProcessInfo[KM.DELTA_TIME] = FEM_delta_time

#============================================================================================================================
    def RegenerateAndUpdatePFEMPressureConditions(self):
        if self.FEMDEM_Solution.domain_size == 2:
            regenerate_cond_process = FEMDEM.RegeneratePfemPressureConditionsProcess2D(self.FEMDEM_Solution.FEM_Solution.main_model_part)
            update_cond_process     = FEMDEM.UpdatePressureValuePfemConditionsProcess2D(self.FEMDEM_Solution.FEM_Solution.main_model_part)
        else:
            regenerate_cond_process = FEMDEM.RegeneratePfemPressureConditionsProcess3D(self.FEMDEM_Solution.FEM_Solution.main_model_part)
            update_cond_process     = FEMDEM.UpdatePressureValuePfemConditionsProcess3D(self.FEMDEM_Solution.FEM_Solution.main_model_part)

        regenerate_cond_process.Execute()
        update_cond_process.Execute()

#============================================================================================================================
    def FixNodesModelPart(self, ModelPart):
        fix_nodes_model_part = FEMDEM.FixFreeVelocityOnNodesProcess(ModelPart, 0)
        fix_nodes_model_part.Execute()

    def FreeNodesModelPart(self, ModelPart):
        free_nodes_model_part = FEMDEM.FixFreeVelocityOnNodesProcess(ModelPart, 1)
        free_nodes_model_part.Execute()

#============================================================================================================================
    def ComputeSkinFEMDEMBoundary(self):
        if self.FEMDEM_Solution.domain_size == 2:
            skin_detection_process = KM.SkinDetectionProcess2D(self.FEMDEM_Solution.FEM_Solution.main_model_part,
                                                                               self.FEMDEM_Solution.SkinDetectionProcessParameters)
        else: # 3D
            skin_detection_process = KM.SkinDetectionProcess3D(self.FEMDEM_Solution.FEM_Solution.main_model_part,
                                                                               self.FEMDEM_Solution.SkinDetectionProcessParameters)    
        skin_detection_process.Execute()

#============================================================================================================================
    def UpdateFEMDEMBoundary(self):
        self.ComputeSkinFEMDEMBoundary()
        update_process = FEMDEM.UpdateFlagNoRemeshFemDemBoundaryProcess(self.FEMDEM_Solution.FEM_Solution.main_model_part)
        update_process.Execute()

#============================================================================================================================

    def InitializePostProcess(self):
        mixed_fluid_solid_mp       = self.model.CreateModelPart('mixed_fluid_solid_mp')
        mixed_fluid_solid_balls_mp = self.model.CreateModelPart('mixed_fluid_solid_balls_mp')
        mixed_solid_balls_mp       = self.model.CreateModelPart('mixed_solid_balls_mp')

        filename = os.path.join(self.FEMDEM_Solution.DEM_Solution.post_path, self.FEMDEM_Solution.DEM_Solution.DEM_parameters["problem_name"].GetString())
        self.gid_output = gid_output.FemDemCoupledGiDOutput(
                            filename,
                            True,
                            "Binary",
                            "Multiples",
                            True,
                            True,
                            self.FEMDEM_Solution.FEM_Solution.main_model_part,
                            self.PFEM_Solution.main_model_part,
                            self.FEMDEM_Solution.DEM_Solution.spheres_model_part,
                            self.FEMDEM_Solution.DEM_Solution.cluster_model_part,
                            self.FEMDEM_Solution.DEM_Solution.rigid_face_model_part,
                            mixed_fluid_solid_mp,
                            mixed_solid_balls_mp,
                            mixed_fluid_solid_balls_mp)

        solid_nodal_results = ["DISPLACEMENT", "ACCELERATION"]
        dem_nodal_results = ["TOTAL_FORCES", "RADIUS"]
        fluid_nodal_results = ["PRESSURE"]
        clusters_nodal_results = []
        rigid_faces_nodal_results = []
        mixed_solid_fluid_nodal_results = ["VELOCITY"]
        mixed_solid_balls_nodal_results = []
        mixed_solid_balls_fluid_nodal_results = []

        gauss_points_results = ["CAUCHY_STRESS_VECTOR", "DAMAGE_ELEMENT", "STRESS_VECTOR_INTEGRATED", "GREEN_LAGRANGE_STRAIN_VECTOR"]

        self.gid_output.initialize_dem_fem_results(solid_nodal_results,
                                                   fluid_nodal_results,
                                                   dem_nodal_results,
                                                   clusters_nodal_results,
                                                   rigid_faces_nodal_results,
                                                   mixed_solid_fluid_nodal_results,
                                                   mixed_solid_balls_nodal_results,
                                                   mixed_solid_balls_fluid_nodal_results,
                                                   gauss_points_results)


#PrintResults============================================================================================================================
    def PrintResults(self):
        print_parameters = self.FEMDEM_Solution.FEM_Solution.ProjectParameters["output_configuration"]["result_file_configuration"]
        if self.FEMDEM_Solution.FEM_Solution.step == 1: # always print the 1st step
            self.gid_output.Writeresults(self.FEMDEM_Solution.FEM_Solution.time)
            self.FEMDEM_Solution.FEM_Solution.time_old_print = self.FEMDEM_Solution.FEM_Solution.time
            self.FEMDEM_Solution.FEM_Solution.step_old_print = self.FEMDEM_Solution.FEM_Solution.step
        else:
            time_to_print = 0
            if print_parameters["output_control_type"].GetString() == "step":
                time_to_print = self.FEMDEM_Solution.FEM_Solution.step - self.FEMDEM_Solution.FEM_Solution.step_old_print
            else:
                time_to_print = self.FEMDEM_Solution.FEM_Solution.time - self.FEMDEM_Solution.FEM_Solution.time_old_print

            if print_parameters["output_control_type"].GetString() == "step":
                if print_parameters["output_interval"].GetInt() - time_to_print == 0:
                    self.gid_output.Writeresults(self.FEMDEM_Solution.FEM_Solution.time)
                    self.FEMDEM_Solution.FEM_Solution.step_old_print = self.FEMDEM_Solution.FEM_Solution.step
            else:
                if print_parameters["output_interval"].GetDouble() - time_to_print < 1e-2 * self.FEMDEM_Solution.FEM_Solution.delta_time:
                    self.gid_output.Writeresults(self.FEMDEM_Solution.FEM_Solution.time)
                    self.FEMDEM_Solution.FEM_Solution.time_old_print = self.FEMDEM_Solution.FEM_Solution.time