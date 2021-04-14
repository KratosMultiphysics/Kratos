

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.StructuralMechanicsApplication as KratosSMA

import KratosMultiphysics.FemToDemApplication.MainDEM_for_coupling as DEM
import KratosMultiphysics.FemToDemApplication.MainFEM_for_coupling as FEM
import KratosMultiphysics.FemToDemApplication.FEMDEMParticleCreatorDestructor as PCD
import KratosMultiphysics.MeshingApplication.mmg_process as MMG
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
        self.DEM_Solution = DEM.DEM_for_coupling_Solution(Model, path)
        self.DEM_Solution.Initialize()
        self.FEM_Solution = FEM.FEM_for_coupling_Solution(Model, path, self.DEM_Solution.solver)

        # Initialize Remeshing files
        self.DoRemeshing = self.FEM_Solution.ProjectParameters["AMR_data"]["activate_AMR"].GetBool()
        if self.DoRemeshing:
            self.mmg_parameter_file = open("MMGParameters.json",'r')
            self.mmg_parameters = KratosMultiphysics.Parameters(self.mmg_parameter_file.read())
            self.RemeshingProcessMMG = MMG.MmgProcess(Model, self.mmg_parameters)
        self.domain_size = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.InitializePlotsFiles()
        self.echo_level = 0
        self.is_slave = False

#============================================================================================================================
    def Initialize(self):
        if self.domain_size == 2:
            self.number_of_nodes_element = 3
        else: # 3D
            self.number_of_nodes_element = 4
            self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.ERASED_VOLUME] = 0.0 # Sand Production Calculations

        self.FEM_Solution.Initialize()


        nodes = self.FEM_Solution.main_model_part.Nodes
        utils = KratosMultiphysics.VariableUtils()
        # Initialize the "flag" IS_DEM in all the nodes
        utils.SetNonHistoricalVariable(KratosFemDem.IS_DEM, False, nodes)
        # Initialize the "flag" NODAL_FORCE_APPLIED in all the nodes
        utils.SetNonHistoricalVariable(KratosFemDem.NODAL_FORCE_APPLIED, False, nodes)
        # Initialize the "flag" RADIUS in all the nodes
        utils.SetNonHistoricalVariable(KratosMultiphysics.RADIUS, 0.0, nodes)

        # Initialize the var to track volume erased for each pressure
        utils.SetNonHistoricalVariable(KratosFemDem.PRESSURE_VOLUME, 0.0, nodes)
        utils.SetNonHistoricalVariable(KratosFemDem.PRESSURE_INITIAL_VOLUME, 0.0, nodes)

        self.SpheresModelPart = self.DEM_Solution.spheres_model_part
        self.DEMParameters = self.DEM_Solution.DEM_parameters
        self.DEMProperties = self.SpheresModelPart.GetProperties()[1]
        self.ParticleCreatorDestructor = PCD.FemDemParticleCreatorDestructor(self.SpheresModelPart,
                                                                           self.DEMProperties,
                                                                           self.DEMParameters)

        if self.domain_size == 3:
            self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part)

        if self.DoRemeshing:
            self.InitializeMMGvariables()
            self.RemeshingProcessMMG.ExecuteInitialize()

        if self.FEM_Solution.ProjectParameters.Has("transfer_dem_contact_forces") == False:
            self.TransferDEMContactForcesToFEM = True
        else:
            self.TransferDEMContactForcesToFEM = self.FEM_Solution.ProjectParameters["transfer_dem_contact_forces"].GetBool()

        if self.FEM_Solution.ProjectParameters.Has("pressure_load_extrapolation") == False:
            self.PressureLoad = False
        else:
            self.PressureLoad = self.FEM_Solution.ProjectParameters["pressure_load_extrapolation"].GetBool()

        if self.FEM_Solution.ProjectParameters.Has("DEM_FEM_contact") == False:
            self.DEMFEM_contact = False
        else:
            self.DEMFEM_contact = self.FEM_Solution.ProjectParameters["DEM_FEM_contact"].GetBool()
        self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.DEMFEM_CONTACT] = self.DEMFEM_contact


        # Initialize IP variables to zero
        self.InitializeIntegrationPointsVariables()

        if self.PressureLoad:
            KratosFemDem.AssignPressureIdProcess(self.FEM_Solution.main_model_part).Execute()
            KratosFemDem.ComputeInitialVolumeProcess(self.FEM_Solution.main_model_part).Execute()

        if self.FEM_Solution.ProjectParameters.Has("tangent_operator") == True:
            # 0 -> Elastic , 1 -> Secant , 2 -> Tangent , 3 -> Tangent 2nd Order
            tangent_type = self.FEM_Solution.ProjectParameters["tangent_operator"].GetInt()
            self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.TANGENT_CONSTITUTIVE_TENSOR] = tangent_type
        else:
            self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.TANGENT_CONSTITUTIVE_TENSOR] = 2

        self.SkinDetectionProcessParameters = KratosMultiphysics.Parameters("""
        {
            "name_auxiliar_model_part" : "SkinDEMModelPart",
            "name_auxiliar_condition"  : "Condition",
            "echo_level"               : 0
        }""")


        # for the dem contact forces coupling
        self.InitializeDummyNodalForces()

        # Just to find neighbours the 1st time
        self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM] = True
        if self.domain_size == 3:
            self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECOMPUTE_NEIGHBOURS] = True

        self.FEM_Solution.KratosPrintInfo("")
        self.FEM_Solution.KratosPrintInfo("    ______                 ___    ____                 ")
        self.FEM_Solution.KratosPrintInfo("   / ____/___   ____ ___  |__ \  / __ \ ___   ____ ___ ")
        self.FEM_Solution.KratosPrintInfo("  / /_   / _ \ / __ `__ \ __/ / / / / // _ \ / __ `__ \ ")
        self.FEM_Solution.KratosPrintInfo(" / __/  /  __// / / / / // __/ / /_/ //  __// / / / / /")
        self.FEM_Solution.KratosPrintInfo("/_/     \___//_/ /_/ /_//____//_____/ \___//_/ /_/ /_/ 2 Way-Coupled Application")
        self.FEM_Solution.KratosPrintInfo("                           Developed by Alejandro Cornejo")
        self.FEM_Solution.KratosPrintInfo("")

        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM Solution initialized")

        if self.domain_size == 3: # only in 3D
            # We assign the flag to recompute neighbours inside the 3D elements the 1st time
            utils = KratosMultiphysics.VariableUtils()
            utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)
            # We assign the flag to recompute neighbours inside the 3D elements the 1st time
            utils = KratosMultiphysics.VariableUtils()
            utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)

        if self.FEM_Solution.ProjectParameters.Has("create_initial_skin") == False:
            self.CreateInitialSkin = False
        else:
            self.CreateInitialSkin = self.FEM_Solution.ProjectParameters["create_initial_skin"].GetBool()

        if self.CreateInitialSkin:
            self.ComputeSkinSubModelPart()
            if self.DEMFEM_contact:
                self.TransferFEMSkinToDEM()
            KratosFemDem.GenerateInitialSkinDEMProcess(self.FEM_Solution.main_model_part, self.SpheresModelPart).Execute()

        # Initialize the coupled post process
        if not self.is_slave:
            self.InitializePostProcess()
        
        self.FindNeighboursIfNecessary()


#============================================================================================================================
    def InitializeSolutionStep(self):
        # Modified for the remeshing
        self.FEM_Solution.delta_time = self.ComputeDeltaTime()
        self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.FEM_Solution.delta_time
        self.FEM_Solution.time = self.FEM_Solution.time + self.FEM_Solution.delta_time
        self.FEM_Solution.main_model_part.CloneTimeStep(self.FEM_Solution.time)
        self.FEM_Solution.step = self.FEM_Solution.step + 1
        self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.FEM_Solution.step

        # self.FindNeighboursIfNecessary()
        self.PerformRemeshingIfNecessary()

        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: InitializeSolutionStep of the FEM part")

        self.FEM_Solution.InitializeSolutionStep()

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
        self.PrintPlotsFiles()
        
        # MODIFIED FOR THE REMESHING
        self.FEM_Solution.GraphicalOutputExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        self.FEM_Solution.model_processes.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output
        self.FEM_Solution.model_processes.ExecuteBeforeOutputStep()

        # write output results GiD: (frequency writing is controlled internally)
        # self.FEM_Solution.GraphicalOutputPrintOutput()

        # processes to be executed after writting the output
        self.FEM_Solution.model_processes.ExecuteAfterOutputStep()

        if self.DoRemeshing:
             self.RemeshingProcessMMG.ExecuteFinalizeSolutionStep()
        
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
            props = self.DEM_Solution.rigid_face_model_part.GetProperties(self.created_props_id,0)
            DemFem.DemStructuresCouplingUtilities().TransferStructuresSkinToDem(fem_skin_mp, dem_walls_mp, props)
        else: # have to create it
            props = self.CreateFEMPropertiesForDEFEContact()
            dem_walls_mp = self.DEM_Solution.rigid_face_model_part.CreateSubModelPart("SkinTransferredFromStructure")
            dem_walls_mp.SetValue(KratosDEM.RIGID_BODY_OPTION, False)
            dem_walls_mp.AddProperties(props)
            DemFem.DemStructuresCouplingUtilities().TransferStructuresSkinToDem(fem_skin_mp, dem_walls_mp, props)