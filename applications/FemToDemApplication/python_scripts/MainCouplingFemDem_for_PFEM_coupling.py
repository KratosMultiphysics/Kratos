

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem

import KratosMultiphysics.FemToDemApplication.MainFEM_for_PFEM_coupling as FEM
import KratosMultiphysics.FemToDemApplication.MainDEM_for_coupling as DEM
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem
import KratosMultiphysics.FemToDemApplication.FEMDEMParticleCreatorDestructor as PCD
import math
import os

def Wait():
    input("Press Something")

#============================================================================================================================
class MainCoupledFemDem_for_PFEM_coupling_Solution(MainCouplingFemDem.MainCoupledFemDem_Solution):
#============================================================================================================================
    def __init__(self, Model, path = ""):
        # Initialize solutions

        if path == "":
            DEMProjectParametersFile = open("ProjectParametersDEM.json", 'r')
        else:
            DEMProjectParametersFile = open(os.path.join(path, "ProjectParametersDEM.json"), 'r')
        DEM_project_parameters = KratosMultiphysics.Parameters(DEMProjectParametersFile.read())

        self.FEM_Solution = FEM.FEM_for_PFEM_coupling_Solution(Model, path)
        self.DEM_Solution = DEM.DEM_for_coupling_Solution(Model, DEM_project_parameters)

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
        self.DEM_Solution.Initialize()

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
            self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part, 4, 5)

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
        self.FEM_Solution.KratosPrintInfo("/_/     \___//_/ /_/ /_//____//_____/ \___//_/ /_/ /_/ Application")
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

        if self.FEM_Solution.ProjectParameters.Has("do_stabilization_solve") == False:
            self.do_stabilization_solve = False
        else:
            self.do_stabilization_solve = self.FEM_Solution.ProjectParameters["do_stabilization_solve"].GetBool()

        if self.CreateInitialSkin:
            self.ComputeSkinSubModelPart()
            if self.DEMFEM_contact:
                self.TransferFEMSkinToDEM()
            KratosFemDem.GenerateInitialSkinDEMProcess(self.FEM_Solution.main_model_part, self.SpheresModelPart).Execute()

        # Initialize the coupled post process
        if not self.is_slave:
            self.InitializePostProcess()

        self.FindNeighboursIfNecessary()