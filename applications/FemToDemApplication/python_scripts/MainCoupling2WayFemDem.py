

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.StructuralMechanicsApplication as KratosSMA

import KratosMultiphysics.FemToDemApplication.MainDEM_for_coupling as DEM
import KratosMultiphysics.FemToDemApplication.MainFEM_for_coupling as FEM
import KratosMultiphysics.FemToDemApplication.FEMDEMParticleCreatorDestructor as PCD
import KratosMultiphysics.DEMApplication as KratosDEM
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem

def Wait():
    input("Press Something")

#============================================================================================================================
class MainCoupled2WayFemDem_Solution(MainCouplingFemDem.MainCoupledFemDem_Solution):
#============================================================================================================================
    def __init__(self, Model, path = ""):
        self.model = Model
        # Initialize solutions

        if path == "":
            DEMProjectParametersFile = open("ProjectParametersDEM.json", 'r')
        else:
            DEMProjectParametersFile = open(os.path.join(path, "ProjectParametersDEM.json"), 'r')
        DEM_project_parameters = KratosMultiphysics.Parameters(DEMProjectParametersFile.read())
        self.DEM_Solution = DEM.DEM_for_coupling_Solution(Model, DEM_project_parameters)

        self.DEM_Solution.Initialize()
        self.FEM_Solution = FEM.FEM_for_coupling_Solution(Model, path, self.DEM_Solution._GetSolver())

        self.domain_size = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # self.InitializePlotsFiles()
        self.echo_level = 0
        self.is_slave = False

#============================================================================================================================
    def Initialize(self):

        self.FEM_Solution.Initialize()
        self.InitializeProcessesAndVariables()

        self.FEM_Solution.KratosPrintInfo("")
        self.FEM_Solution.KratosPrintInfo("    ______                 ___    ____                 ")
        self.FEM_Solution.KratosPrintInfo("   / ____/___   ____ ___  |__ \  / __ \ ___   ____ ___ ")
        self.FEM_Solution.KratosPrintInfo("  / /_   / _ \ / __ `__ \ __/ / / / / // _ \ / __ `__ \ ")
        self.FEM_Solution.KratosPrintInfo(" / __/  /  __// / / / / // __/ / /_/ //  __// / / / / /")
        self.FEM_Solution.KratosPrintInfo("/_/     \___//_/ /_/ /_//____//_____/ \___//_/ /_/ /_/ 2 Way-Coupled Application")
        self.FEM_Solution.KratosPrintInfo("                           Developed by Alejandro Cornejo")
        self.FEM_Solution.KratosPrintInfo("")

#============================================================================================================================
    def InitializeSolutionStep(self):
        # Modified for the remeshing
        self.FEM_Solution.delta_time = self.ComputeDeltaTime()
        self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.FEM_Solution.delta_time
        self.FEM_Solution.time = self.FEM_Solution.time + self.FEM_Solution.delta_time
        self.FEM_Solution.main_model_part.CloneTimeStep(self.FEM_Solution.time)
        self.FEM_Solution.step = self.FEM_Solution.step + 1
        self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.FEM_Solution.step

        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: InitializeSolutionStep of the FEM part")

        self.FEM_Solution.InitializeSolutionStep()
        self.DEM_Solution._GetSolver().AdvanceInTime(self.FEM_Solution.time)
        self.DEM_Solution.time = self.FEM_Solution.time
        self.DEM_Solution.step = self.FEM_Solution.step

#============================================================================================================================
    def SolveSolutionStep(self):  # Method to perform the coupling FEM <-> DEM

        self.FEM_Solution.clock_time = self.FEM_Solution.StartTimeMeasuring()

        #### SOLVE FEM #########################################
        self.FEM_Solution.solver.InitializeSolutionStep()
        self.FEM_Solution.solver.Predict()
        self.FEM_Solution.solver.SolveSolutionStep()
        self.FEM_Solution.solver.FinalizeSolutionStep()
        ########################################################

        self.ExecuteBeforeGeneratingDEM()
        self.GenerateDEM() # we create the new DEM of this time step
        self.ExecuteAfterGeneratingDEM()
        self.BeforeSolveDEMOperations()


#============================================================================================================================
    def FinalizeSolutionStep(self):

        self.FEM_Solution.StopTimeMeasuring(self.FEM_Solution.clock_time,"Solving", False)

        # Print required info
        # self.PrintPlotsFiles()

        # MODIFIED FOR THE REMESHING
        self.FEM_Solution.GraphicalOutputExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        self.FEM_Solution.model_processes.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output
        self.FEM_Solution.model_processes.ExecuteBeforeOutputStep()


        # processes to be executed after writing the output
        self.FEM_Solution.model_processes.ExecuteAfterOutputStep()

        if not self.is_slave:
            self.PrintResults()

#InitializeDummyNodalForces============================================================================================================================
    def InitializeDummyNodalForces(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: InitializeDummyNodalForces")

        # we fill the submodel part with the nodes and dummy conditions
        max_id = self.GetMaximumConditionId()
        props = self.FEM_Solution.main_model_part.Properties[0]
        self.FEM_Solution.main_model_part.CreateSubModelPart("ContactForcesDEMConditions")
        self.FEM_Solution.main_model_part.GetSubModelPart("computing_domain").CreateSubModelPart("ContactForcesDEMConditions")
        for node in self.FEM_Solution.main_model_part.Nodes:
            self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").AddNode(node, 0)
            max_id += 1
            cond = self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").CreateNewCondition(
                                                                            "PointLoadCondition3D1N",
                                                                            max_id,
                                                                            [node.Id],
                                                                            props)
            self.FEM_Solution.main_model_part.GetSubModelPart("computing_domain").AddCondition(cond)
            self.FEM_Solution.main_model_part.GetSubModelPart("computing_domain").GetSubModelPart("ContactForcesDEMConditions").AddCondition(cond)
            self.FEM_Solution.main_model_part.GetCondition(max_id).SetValue(KratosSMA.POINT_LOAD, [0.0,0.0,0.0])

#TransferFEMSkinToDEM============================================================================================================================
    def TransferFEMSkinToDEM(self):
        fem_skin_mp = self.FEM_Solution.main_model_part.GetSubModelPart("SkinDEMModelPart")

        if self.DEM_Solution.rigid_face_model_part.HasSubModelPart("SkinTransferredFromStructure"):
            self.EraseConditionsAndNodesSubModelPart()
            dem_walls_mp = self.DEM_Solution.rigid_face_model_part.GetSubModelPart("SkinTransferredFromStructure")
            dem_walls_mp.SetValue(KratosDEM.RIGID_BODY_OPTION, False)
            props = self.DEM_Solution.spheres_model_part.GetProperties()[2]
            DemFem.DemStructuresCouplingUtilities().TransferStructuresSkinToDem(fem_skin_mp, dem_walls_mp, props)
        else: # have to create it
            # props = self.CreateFEMPropertiesForDEFEContact()
            props = self.DEM_Solution.spheres_model_part.GetProperties()[2]
            dem_walls_mp = self.DEM_Solution.rigid_face_model_part.CreateSubModelPart("SkinTransferredFromStructure")
            dem_walls_mp.SetValue(KratosDEM.RIGID_BODY_OPTION, False)
            # dem_walls_mp.AddProperties(props)
            DemFem.DemStructuresCouplingUtilities().TransferStructuresSkinToDem(fem_skin_mp, dem_walls_mp, props)