
import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem

def Wait():
    input("Press Something")

#============================================================================================================================
class MainCoupledFemDemSubstepping_Solution(MainCouplingFemDem.MainCoupledFemDem_Solution):
#============================================================================================================================

#============================================================================================================================
    def SolveSolutionStep(self):  # Method to perform the coupling FEM <-> DEM

        self.FEM_Solution.clock_time = self.FEM_Solution.StartTimeMeasuring()

        #### SOLVE FEM #########################################
        self.FEM_Solution.solver.Solve()
        ########################################################

        self.ExecuteBeforeGeneratingDEM()
        self.GenerateDEM() # we create the new DEM of this time step
        self.ExecuteAfterGeneratingDEM()

        FEMDEM_utilities = KratosFemDem.FEMDEMCouplingUtilities()
        FEMDEM_utilities.SaveStructuralSolution(self.FEM_Solution.main_model_part)
        # FEMDEM_utilities.IdentifyFreeParticles(self.FEM_Solution.main_model_part, self.DEM_Solution.spheres_model_part)

        # Perform substepping
        pseudo_substepping_time = 0
        if self.DEM_Solution.spheres_model_part.NumberOfElements() > 0:
            self.FEM_Solution.KratosPrintInfo("Performing DEM Substepping... Explicit time step: " + str(self.DEM_Solution._GetSolver().dt))
            self.FEM_Solution.KratosPrintInfo("")
            while pseudo_substepping_time <= self.FEM_Solution.delta_time:
                ### Begin Substepping

                self.BeforeSolveDEMOperations()
                FEMDEM_utilities.InterpolateStructuralSolution(self.FEM_Solution.main_model_part,
                                                               self.FEM_Solution.delta_time,
                                                               self.FEM_Solution.time,
                                                               self.DEM_Solution._GetSolver().dt,
                                                               self.DEM_Solution.time)
                self.UpdateDEMVariables()

                #### SOLVE DEM #########################################
                self.DEM_Solution.SolverSolve()
                ########################################################

                FEMDEM_utilities.AddExplicitImpulses(self.FEM_Solution.main_model_part, self.DEM_Solution._GetSolver().dt)

                self.DEM_Solution.FinalizeSolutionStep()
                # self.DEM_Solution.solver._MoveAllMeshes(self.DEM_Solution.time, self.DEM_Solution._GetSolver().dt)

                # We reset the position of the slave DEM
                self.UpdateDEMVariables()

                # Advancing in DEM explicit scheme
                pseudo_substepping_time += self.DEM_Solution._GetSolver().dt
            ### End Substepping

            # Reset the data base for the FEM
            FEMDEM_utilities.RestoreStructuralSolution(self.FEM_Solution.main_model_part)

        else: # In case there are no DEM yet
            self.OnlyUpdateTimeAndStepInDEM()

        # Transfer the contact forces of the DEM to the FEM nodes
        if self.TransferDEMContactForcesToFEM:
            FEMDEM_utilities.ComputeAndTranferAveragedContactTotalForces(self.FEM_Solution.main_model_part, self.FEM_Solution.delta_time)

#============================================================================================================================
    def FinalizeSolutionStep(self):

        # Transfer the contact forces of the DEM to the FEM nodes
        if self.TransferDEMContactForcesToFEM:
            self.TransferNodalForcesToFEM()

        self.FEM_Solution.StopTimeMeasuring(self.FEM_Solution.clock_time,"Solving", False)

        # Print required info
        # self.PrintPlotsFiles()

        # MODIFIED FOR THE REMESHING
        self.FEM_Solution.GraphicalOutputExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        self.FEM_Solution.model_processes.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output
        self.FEM_Solution.model_processes.ExecuteBeforeOutputStep()

        # processes to be executed after writting the output
        self.FEM_Solution.model_processes.ExecuteAfterOutputStep()

        if not self.is_slave:
            self.PrintResults()

#============================================================================================================================
    def BeforeSolveDEMOperations(self):
        self.DEM_Solution.time += self.DEM_Solution._GetSolver().dt
        self.DEM_Solution._GetSolver().AdvanceInTime(self.DEM_Solution.time)
        self.DEM_Solution.step += 1
        self.DEM_Solution._GetSolver().Predict()

#============================================================================================================================
    def OnlyUpdateTimeAndStepInDEM(self):
        self.DEM_Solution.time = self.FEM_Solution.time
        self.DEM_Solution.step = self.FEM_Solution.step
        self.DEM_Solution.DEMFEMProcedures.UpdateTimeInModelParts(self.DEM_Solution.all_model_parts,
                                                                   self.DEM_Solution.time,
                                                                   self.DEM_Solution.solver.dt,
                                                                   self.DEM_Solution.step,
                                                                   self.DEM_Solution.IsTimeToPrintPostProcess())
