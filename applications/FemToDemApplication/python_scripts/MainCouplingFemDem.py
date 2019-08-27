from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import FEMDEMParticleCreatorDestructor as PCD
import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import CouplingFemDem
import math
import KratosMultiphysics.MeshingApplication as MeshingApplication

def Wait():
    input("Press Something")

class MainCoupledFemDem_Solution:

#============================================================================================================================
    def __init__(self, Model):
        # Initialize solutions
        self.FEM_Solution = FEM.FEM_for_coupling_Solution(Model)
        self.DEM_Solution = DEM.DEM_for_coupling_Solution(Model)

        # Initialize Remeshing files
        self.DoRemeshing = self.FEM_Solution.ProjectParameters["AMR_data"]["activate_AMR"].GetBool()
        if self.DoRemeshing:
            self.mmg_parameter_file = open("MMGParameters.json",'r')
            self.mmg_parameters = KratosMultiphysics.Parameters(self.mmg_parameter_file.read())
            self.RemeshingProcessMMG = MMG.MmgProcess(Model, self.mmg_parameters)
        self.InitializePlotsFiles()
        self.echo_level = 0
        self.domain_size = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

#============================================================================================================================
    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

#============================================================================================================================
    def Initialize(self):
        if self.domain_size == 2:
            self.number_of_nodes_element = 3
        else: # 3D
            self.number_of_nodes_element = 4
            self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.ERASED_VOLUME] = 0.0 # Sand Production Calculations
        self.FEM_Solution.Initialize()
        self.DEM_Solution.Initialize()

        nodes = self.FEM_Solution.main_model_part.Nodes
        utils = KratosMultiphysics.VariableUtils()
        # Initialize the "flag" IS_DEM in all the nodes
        utils.SetNonHistoricalVariable(KratosFemDem.IS_DEM, False, nodes)
        # Initialize the "flag" NODAL_FORCE_APPLIED in all the nodes
        utils.SetNonHistoricalVariable(KratosFemDem.NODAL_FORCE_APPLIED, False, nodes)
        # Initialize the "flag" RADIUS in all the nodes
        utils.SetNonHistoricalVariable(KratosMultiphysics.RADIUS, 0.0, nodes)

        # Initialize IP variables to zero
        self.InitializeIntegrationPointsVariables()

        self.SpheresModelPart = self.DEM_Solution.spheres_model_part
        self.DEMParameters = self.DEM_Solution.DEM_parameters
        self.DEMProperties = self.SpheresModelPart.GetProperties()[1]
        self.ParticleCreatorDestructor = PCD.FemDemParticleCreatorDestructor(self.SpheresModelPart,
                                                                           self.DEMProperties,
                                                                           self.DEMParameters)

        if self.domain_size == 3:
            self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part, 4, 5)

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
        if self.PressureLoad:
            KratosFemDem.AssignPressureIdProcess(self.FEM_Solution.main_model_part).Execute()

        # if self.FEM_Solution.ProjectParameters.Has("displacement_perturbed_tangent") == False:
        #     self.DisplacementPerturbedTangent = False
        # else:
        #     self.DisplacementPerturbedTangent = self.FEM_Solution.ProjectParameters["displacement_perturbed_tangent"].GetBool()

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

        self.FEM_Solution.KratosPrintInfo(" /$$$$$$$$ /$$$$$$$$ /$$      /$$  /$$$$$$  /$$$$$$$  /$$$$$$$$ /$$      /$$")
        self.FEM_Solution.KratosPrintInfo("| $$_____/| $$_____/| $$$    /$$$ /$$__  $$| $$__  $$| $$_____/| $$$    /$$$")
        self.FEM_Solution.KratosPrintInfo("| $$      | $$      | $$$$  /$$$$|__/  \ $$| $$  \ $$| $$      | $$$$  /$$$$")
        self.FEM_Solution.KratosPrintInfo("| $$$$$   | $$$$$   | $$ $$/$$ $$  /$$$$$$/| $$  | $$| $$$$$   | $$ $$/$$ $$")
        self.FEM_Solution.KratosPrintInfo("| $$__/   | $$__/   | $$  $$$| $$ /$$____/ | $$  | $$| $$__/   | $$  $$$| $$")
        self.FEM_Solution.KratosPrintInfo("| $$      | $$      | $$\  $ | $$| $$      | $$  | $$| $$      | $$\  $ | $$")
        self.FEM_Solution.KratosPrintInfo("| $$      | $$$$$$$$| $$ \/  | $$| $$$$$$$$| $$$$$$$/| $$$$$$$$| $$ \/  | $$")
        self.FEM_Solution.KratosPrintInfo("|__/      |________/|__/     |__/|________/|_______/ |________/|__/     |__/ Application")
        self.FEM_Solution.KratosPrintInfo("                                                    Developed by Alejandro Cornejo")
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

#============================================================================================================================
    def RunMainTemporalLoop(self):
        # Solving the problem (time integration)
        self.DEM_Solution.step           = 0
        self.DEM_Solution.time           = 0.0
        self.DEM_Solution.time_old_print = 0.0

        if self.DoRemeshing:
            self.RemeshingProcessMMG.ExecuteBeforeSolutionLoop()

        # Temporal loop
        while self.FEM_Solution.time <= self.FEM_Solution.end_time:
            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()

#============================================================================================================================
    def InitializeSolutionStep(self):
        # Modified for the remeshing
		self.FEM_Solution.delta_time = self.ComputeDeltaTime()
		self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.FEM_Solution.delta_time
		self.FEM_Solution.time = self.FEM_Solution.time + self.FEM_Solution.delta_time
		self.FEM_Solution.main_model_part.CloneTimeStep(self.FEM_Solution.time)
		self.FEM_Solution.step = self.FEM_Solution.step + 1
		self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.FEM_Solution.step

		self.FindNeighboursIfNecessary()	
		self.PerformRemeshingIfNecessary()

		if self.echo_level > 0:
			self.FEM_Solution.KratosPrintInfo("FEM-DEM:: InitializeSolutionStep of the FEM part")

		self.FEM_Solution.InitializeSolutionStep()










#============================================================================================================================
    def InitializeIntegrationPointsVariables(self):
        utils = KratosMultiphysics.VariableUtils()
        elements = self.FEM_Solution.main_model_part.Elements
        if self.domain_size == 3:
            utils.SetNonHistoricalVariable(KratosFemDem.VOLUME_COUNTED, False, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.STRESS_THRESHOLD, 0.0, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.DAMAGE_ELEMENT, 0.0, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.PRESSURE_EXPANDED, 0, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.IS_SKIN, 0, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.SMOOTHING, 0, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.STRESS_VECTOR, [0.0,0.0,0.0,0.0,0.0,0.0], elements)
            utils.SetNonHistoricalVariable(KratosFemDem.STRAIN_VECTOR, [0.0,0.0,0.0,0.0,0.0,0.0], elements)
            utils.SetNonHistoricalVariable(KratosFemDem.STRESS_VECTOR_INTEGRATED, [0.0,0.0,0.0,0.0,0.0,0.0], elements)
        else: # 2D
            elements = self.FEM_Solution.main_model_part.Elements
            utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.STRESS_THRESHOLD, 0.0, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.DAMAGE_ELEMENT, 0.0, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.PRESSURE_EXPANDED, 0, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.IS_SKIN, 0, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.SMOOTHING, 0, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.STRESS_VECTOR, [0.0,0.0,0.0], elements)
            utils.SetNonHistoricalVariable(KratosFemDem.STRAIN_VECTOR, [0.0,0.0,0.0], elements)
            utils.SetNonHistoricalVariable(KratosFemDem.STRESS_VECTOR_INTEGRATED, [0.0,0.0,0.0], elements)

#============================================================================================================================
    def InitializeMMGvariables(self):

        ZeroVector3 = KratosMultiphysics.Vector(3)
        ZeroVector3[0] = 0.0
        ZeroVector3[1] = 0.0
        ZeroVector3[2] = 0.0

        utils = KratosMultiphysics.VariableUtils()
        nodes = self.FEM_Solution.main_model_part.Nodes
        utils.SetNonHistoricalVariable(MeshingApplication.AUXILIAR_GRADIENT, ZeroVector3, nodes)

#============================================================================================================================
    def InitializeDummyNodalForces(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: InitializeDummyNodalForces")

        # we fill the submodel part with the nodes and dummy conditions
        max_id = self.GetMaximumConditionId()
        props = self.FEM_Solution.main_model_part.Properties[0]
        self.FEM_Solution.main_model_part.CreateSubModelPart("ContactForcesDEMConditions")
        for node in self.FEM_Solution.main_model_part.Nodes:
            self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").AddNode(node, 0)
            max_id += 1
            cond = self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").CreateNewCondition(
                                                                            "PointLoadCondition2D1N",
                                                                            max_id,
                                                                            [node.Id],
                                                                            props)
            self.FEM_Solution.main_model_part.GetSubModelPart("computing_domain").AddCondition(cond)
            self.FEM_Solution.main_model_part.GetCondition(max_id).SetValue(Solid.FORCE_LOAD, [0.0,0.0,0.0])

#===================================================================================================================================
	def FindNeighboursIfNecessary(self):
		if self.echo_level > 0:
			self.FEM_Solution.KratosPrintInfo("FEM-DEM:: ComputeNeighboursIfNecessary")

        if self.domain_size == 3:
            if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM]: # The neighbours have changed
                self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part, 4, 5)
                self.nodal_neighbour_finder.Execute()
                # We reset the flag
                self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM] = False
        else: # 2D
            if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM]: # The neighbours have changed
                neighbour_elemental_finder =  KratosMultiphysics.FindElementalNeighboursProcess(self.FEM_Solution.main_model_part, 2, 5)
                neighbour_elemental_finder.Execute()
                # We reset the flag
                self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM] = False

#============================================================================================================================

    def PerformRemeshingIfNecessary(self):

        debug_metric = False
        if debug_metric:
            params = KratosMultiphysics.Parameters("""{}""")
            KratosFemDem.ComputeNormalizedFreeEnergyOnNodesProcess(self.FEM_Solution.main_model_part, self.FEM_Solution.ProjectParameters["AMR_data"]["hessian_variable_parameters"]).Execute()
            MeshingApplication.ComputeHessianSolMetricProcess(self.FEM_Solution.main_model_part, KratosFemDem.EQUIVALENT_NODAL_STRESS, params).Execute()

        if self.DoRemeshing:
            is_remeshing = self.CheckIfHasRemeshed()

            if is_remeshing:
                if self.echo_level > 0:
                    self.FEM_Solution.KratosPrintInfo("FEM-DEM:: ComputeNormalizedFreeEnergyOnNodesProcess")
                # Extrapolate the free energy as a remeshing criterion
                parameters = self.FEM_Solution.ProjectParameters["AMR_data"]["hessian_variable_parameters"]
                KratosFemDem.ComputeNormalizedFreeEnergyOnNodesProcess(self.FEM_Solution.main_model_part, parameters).Execute()

                # we eliminate the nodal DEM forces
                self.RemoveDummyNodalForces()

            # Perform remeshing
            self.RemeshingProcessMMG.ExecuteInitializeSolutionStep()

            if is_remeshing:
                if self.echo_level > 0:
                    self.FEM_Solution.KratosPrintInfo("FEM-DEM:: InitializeSolutionAfterRemeshing")

                if self.domain_size == 3:
                    self.RefineMappedVariables()
                    self.InitializeSolutionAfterRemeshing()
                    self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part, 4, 5)
                    self.nodal_neighbour_finder.Execute()
                    # We assign the flag to recompute neighbours inside the 3D elements
                    utils = KratosMultiphysics.VariableUtils()
                    utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)
                else: # 2D
                    self.InitializeSolutionAfterRemeshing()
                    neighbour_elemental_finder =  KratosMultiphysics.FindElementalNeighboursProcess(self.FEM_Solution.main_model_part, 2, 5)
                    neighbour_elemental_finder.ClearNeighbours()
                    neighbour_elemental_finder.Execute()

#============================================================================================================================
    def RefineMappedVariables(self):
        for elem in self.FEM_Solution.main_model_part.Elements:
            if elem.GetValue(KratosFemDem.DAMAGE_ELEMENT) < 0.0:
                elem.SetValue(KratosFemDem.DAMAGE_ELEMENT, 0.0)

#============================================================================================================================
    def InitializeSolutionAfterRemeshing(self):
        utils = KratosMultiphysics.VariableUtils()
        nodes = self.FEM_Solution.main_model_part.Nodes
        # Initialize the "flag" IS_DEM in all the nodes
        utils.SetNonHistoricalVariable(KratosFemDem.IS_DEM, False, nodes)
        # Initialize the "flag" NODAL_FORCE_APPLIED in all the nodes
        utils.SetNonHistoricalVariable(KratosFemDem.NODAL_FORCE_APPLIED, False, nodes)
        # Initialize the "flag" RADIUS in all the nodes
        utils.SetNonHistoricalVariable(KratosMultiphysics.RADIUS, 0.0, nodes)

        if self.FEM_Solution.ProjectParameters.Has("pressure_load_extrapolation") == False:
            self.PressureLoad = False
        else:
            self.PressureLoad = self.FEM_Solution.ProjectParameters["pressure_load_extrapolation"].GetBool()
        if self.PressureLoad:
            KratosFemDem.AssignPressureIdProcess(self.FEM_Solution.main_model_part).Execute()

        # Remove DEMS from previous mesh
        self.SpheresModelPart.Elements.clear()
        self.SpheresModelPart.Nodes.clear()

        self.InitializeDummyNodalForces()

        self.InitializeMMGvariables()
        self.FEM_Solution.model_processes = self.FEM_Solution.AddProcesses()
        self.FEM_Solution.model_processes.ExecuteInitialize()
        self.FEM_Solution.model_processes.ExecuteBeforeSolutionLoop()
        self.FEM_Solution.model_processes.ExecuteInitializeSolutionStep()

        # Search the skin nodes for the remeshing
        if self.domain_size == 2:
            skin_detection_process = KratosMultiphysics.SkinDetectionProcess2D(self.FEM_Solution.main_model_part,
                                                                                self.SkinDetectionProcessParameters)
        else: # 3D
            skin_detection_process = KratosMultiphysics.SkinDetectionProcess3D(self.FEM_Solution.main_model_part,
                                                                            skin_detection_process_param)    
        skin_detection_process.Execute()
        self.GenerateDemAfterRemeshing()