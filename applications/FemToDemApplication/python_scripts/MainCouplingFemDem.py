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
        # Initialize the "flag" IS_DEM in all the nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosFemDem.IS_DEM, False, nodes)
        # Initialize the "flag" NODAL_FORCE_APPLIED in all the nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosFemDem.NODAL_FORCE_APPLIED, False, nodes)
        # Initialize the "flag" RADIUS in all the nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.RADIUS, 0.0, nodes)

        # Initialize IP variables to zero
        self.InitializeIntegrationPointsVariables()

        self.SpheresModelPart = self.DEM_Solution.spheres_model_part
        self.DEMParameters = self.DEM_Solution.DEM_parameters
        self.DEMProperties = self.SpheresModelPart.GetProperties()[1]
        self.ParticleCreatorDestructor = PCD.FemDemParticleCreatorDestructor(self.SpheresModelPart,
                                                                           self.DEMProperties,
                                                                           self.DEMParameters)

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

        if self.FEM_Solution.ProjectParameters.Has("displacement_perturbed_tangent") == False:
            self.DisplacementPerturbedTangent = False
        else:
            self.DisplacementPerturbedTangent = self.FEM_Solution.ProjectParameters["displacement_perturbed_tangent"].GetBool()

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

        self.FEM_Solution.KratosPrintInfo(" /$$$$$$$$ /$$$$$$$$ /$$      /$$  /$$$$$$  /$$$$$$$  /$$$$$$$$ /$$      /$$")
        self.FEM_Solution.KratosPrintInfo("| $$_____/| $$_____/| $$$    /$$$ /$$__  $$| $$__  $$| $$_____/| $$$    /$$$")
        self.FEM_Solution.KratosPrintInfo("| $$      | $$      | $$$$  /$$$$|__/  \ $$| $$  \ $$| $$      | $$$$  /$$$$")
        self.FEM_Solution.KratosPrintInfo("| $$$$$   | $$$$$   | $$ $$/$$ $$  /$$$$$$/| $$  | $$| $$$$$   | $$ $$/$$ $$")
        self.FEM_Solution.KratosPrintInfo("| $$__/   | $$__/   | $$  $$$| $$ /$$____/ | $$  | $$| $$__/   | $$  $$$| $$")
        self.FEM_Solution.KratosPrintInfo("| $$      | $$      | $$\  $ | $$| $$      | $$  | $$| $$      | $$\  $ | $$")
        self.FEM_Solution.KratosPrintInfo("| $$      | $$$$$$$$| $$ \/  | $$| $$$$$$$$| $$$$$$$/| $$$$$$$$| $$ \/  | $$")
        self.FEM_Solution.KratosPrintInfo("|__/      |________/|__/     |__/|________/|_______/ |________/|__/     |__/ Application")

        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM Solution initialized")

        if self.echo_level > 0:
            KratosMultiphysics.Logger.PrintInfo("FEM-DEM Solution initialized")












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