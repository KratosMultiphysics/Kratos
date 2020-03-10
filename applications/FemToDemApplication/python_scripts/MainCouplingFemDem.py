from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7


import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem

import KratosMultiphysics.FemToDemApplication.MainDEM_for_coupling as DEM
import KratosMultiphysics.FemToDemApplication.MainFEM_for_coupling as FEM
import KratosMultiphysics.FemToDemApplication.FEMDEMParticleCreatorDestructor as PCD
import math
import os
import KratosMultiphysics.MeshingApplication as MeshingApplication
import KratosMultiphysics.SolidMechanicsApplication as Solid
import KratosMultiphysics.MeshingApplication.mmg_process as MMG
import KratosMultiphysics.DEMApplication as KratosDEM
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

def Wait():
    input("Press Something")

#============================================================================================================================
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

        if self.CreateInitialSkin:
            self.ComputeSkinSubModelPart()
            if self.DEMFEM_contact:
                self.TransferFEMSkinToDEM()
            KratosFemDem.GenerateInitialSkinDEMProcess(self.FEM_Solution.main_model_part, self.SpheresModelPart).Execute()

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
    def SolveSolutionStep(self):  # Method to perform the coupling FEM <-> DEM

        self.FEM_Solution.clock_time = self.FEM_Solution.StartTimeMeasuring()

        #### SOLVE FEM #########################################
        self.FEM_Solution.solver.Solve()
        ########################################################

        self.ExecuteBeforeGeneratingDEM()
        self.GenerateDEM() # we create the new DEM of this time step
        self.ExecuteAfterGeneratingDEM()
        self.BeforeSolveDEMOperations()

        #### SOLVE DEM #########################################
        self.DEM_Solution.solver.Solve()
        ########################################################

        self.DEM_Solution.FinalizeSolutionStep()
        self.DEM_Solution.solver._MoveAllMeshes(self.DEM_Solution.time, self.DEM_Solution.solver.dt)

        # to print DEM with the FEM coordinates
        self.UpdateDEMVariables()

        # DEM GiD print output
        self.PrintDEMResults()

        self.DEM_Solution.FinalizeTimeStep(self.DEM_Solution.time)

        # Transfer the contact forces of the DEM to the FEM nodes
        if self.TransferDEMContactForcesToFEM:
            self.TransferNodalForcesToFEM()

        self.FEM_Solution.StopTimeMeasuring(self.FEM_Solution.clock_time,"Solving", False)

        # Update Coupled Postprocess file for Gid (post.lst)
        self.WritePostListFile()

        # Print required info
        self.PrintPlotsFiles()

#============================================================================================================================
    def FinalizeSolutionStep(self):
        # MODIFIED FOR THE REMESHING
        self.FEM_Solution.GraphicalOutputExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        self.FEM_Solution.model_processes.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output
        self.FEM_Solution.model_processes.ExecuteBeforeOutputStep()

        # write output results GiD: (frequency writing is controlled internally)
        self.FEM_Solution.GraphicalOutputPrintOutput()

        # processes to be executed after writting the output
        self.FEM_Solution.model_processes.ExecuteAfterOutputStep()

        if self.DoRemeshing:
             self.RemeshingProcessMMG.ExecuteFinalizeSolutionStep()

#============================================================================================================================
    def Finalize(self):
        self.FEM_Solution.Finalize()
        self.DEM_Solution.Finalize()
        self.DEM_Solution.CleanUpOperations()

        if self.DoRemeshing:
            self.RemeshingProcessMMG.ExecuteFinalize()

#InitializeIntegrationPointsVariables============================================================================================================================
    def InitializeIntegrationPointsVariables(self):
        utils = KratosMultiphysics.VariableUtils()
        elements = self.FEM_Solution.main_model_part.Elements
        nodes = self.FEM_Solution.main_model_part.Nodes
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
            utils.SetNonHistoricalVariable(KratosFemDem.STRESS_VECTOR_INTEGRATED, [0.0, 0.0, 0.0], elements)
        
        if self.PressureLoad:
            utils.SetNonHistoricalVariable(KratosFemDem.PRESSURE_ID, 0, nodes)

#InitializeMMGvariables============================================================================================================================
    def InitializeMMGvariables(self):

        ZeroVector3 = KratosMultiphysics.Vector(3)
        ZeroVector3[0] = 0.0
        ZeroVector3[1] = 0.0
        ZeroVector3[2] = 0.0

        utils = KratosMultiphysics.VariableUtils()
        nodes = self.FEM_Solution.main_model_part.Nodes
        utils.SetNonHistoricalVariable(MeshingApplication.AUXILIAR_GRADIENT, ZeroVector3, nodes)

#InitializeDummyNodalForces============================================================================================================================
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

#FindNeighboursIfNecessary===================================================================================================================================
    def FindNeighboursIfNecessary(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: ComputeNeighboursIfNecessary")

        if self.domain_size == 3:
            if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM]: # The neighbours have changed
                self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part)
                self.nodal_neighbour_finder.Execute()
                # We reset the flag
                self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM] = False
        else: # 2D
            if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM]: # The neighbours have changed
                neighbour_elemental_finder =  KratosMultiphysics.FindElementalNeighboursProcess(self.FEM_Solution.main_model_part, 2, 5)
                neighbour_elemental_finder.Execute()
                # We reset the flag
                self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM] = False

#PerformRemeshingIfNecessary============================================================================================================================
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
                    self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part)
                    self.nodal_neighbour_finder.Execute()
                    # We assign the flag to recompute neighbours inside the 3D elements
                    utils = KratosMultiphysics.VariableUtils()
                    utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)
                else: # 2D
                    self.InitializeSolutionAfterRemeshing()
                    neighbour_elemental_finder =  KratosMultiphysics.FindElementalNeighboursProcess(self.FEM_Solution.main_model_part, 2, 5)
                    neighbour_elemental_finder.ClearNeighbours()
                    neighbour_elemental_finder.Execute()

#RefineMappedVariables============================================================================================================================
    def RefineMappedVariables(self):
        for elem in self.FEM_Solution.main_model_part.Elements:
            if elem.GetValue(KratosFemDem.DAMAGE_ELEMENT) < 0.0:
                elem.SetValue(KratosFemDem.DAMAGE_ELEMENT, 0.0)

#InitializeSolutionAfterRemeshing============================================================================================================================
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
        self.ComputeSkinSubModelPart()
        if self.DEMFEM_contact:
            self.TransferFEMSkinToDEM()
        self.GenerateDemAfterRemeshing()

#ComputeSkinSubModelPart============================================================================================================================
    def ComputeSkinSubModelPart(self):
        # Search the skin nodes for the remeshing
        if self.domain_size == 2:
            skin_detection_process = KratosMultiphysics.SkinDetectionProcess2D(self.FEM_Solution.main_model_part,
                                                                               self.SkinDetectionProcessParameters)
        else: # 3D
            skin_detection_process = KratosMultiphysics.SkinDetectionProcess3D(self.FEM_Solution.main_model_part,
                                                                               self.SkinDetectionProcessParameters)
        skin_detection_process.Execute()

#ComputeDeltaTime============================================================================================================================
    def ComputeDeltaTime(self):
        if self.FEM_Solution.ProjectParameters["problem_data"].Has("time_step"):
            return self.FEM_Solution.ProjectParameters["problem_data"]["time_step"].GetDouble()

        elif self.FEM_Solution.ProjectParameters["problem_data"].Has("variable_time_steps"):

            current_time = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
            for key in self.FEM_Solution.ProjectParameters["problem_data"]["variable_time_steps"].keys():
                interval_settings = self.FEM_Solution.ProjectParameters["problem_data"]["variable_time_steps"][key]
                interval = KratosMultiphysics.IntervalUtility(interval_settings)

                 # Getting the time step of the interval
                if interval.IsInInterval(current_time):
                    return interval_settings["time_step"].GetDouble()
            # If we arrive here we raise an error because the intervals are not well defined
            raise Exception("::[MechanicalSolver]:: Time stepping not well defined!")
        else:
            raise Exception("::[MechanicalSolver]:: Time stepping not defined!")

#ExpandWetNodes============================================================================================================================
    def ExpandWetNodes(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: ExpandWetNodes")

        if self.PressureLoad:
            # This must be called before Generating DEM
            self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECONSTRUCT_PRESSURE_LOAD] = 0 # It is modified inside
            extend_wet_nodes_process = KratosFemDem.ExpandWetNodesProcess(self.FEM_Solution.main_model_part)
            extend_wet_nodes_process.Execute()

#GenerateDEM============================================================================================================================
    def GenerateDEM(self): # This method creates the DEM elements and remove the damaged FEM, Additionally remove the isolated elements
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: GenerateDEM")

        # If we want to compute sand production
        # self.CountErasedVolume()

        if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM]:
            dem_generator_process = KratosFemDem.GenerateDemProcess(self.FEM_Solution.main_model_part, self.SpheresModelPart)
            dem_generator_process.Execute()

            # We remove the inactive DEM associated to fem_nodes
            # self.RemoveAloneDEMElements()
            # self.RemoveIsolatedFiniteElements()
            element_eliminator = KratosMultiphysics.AuxiliarModelPartUtilities(self.FEM_Solution.main_model_part)
            element_eliminator.RemoveElementsAndBelongings(KratosMultiphysics.TO_ERASE)

            if self.domain_size == 3:
                # We assign the flag to recompute neighbours inside the 3D elements
                utils = KratosMultiphysics.VariableUtils()
                utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)
            
            # We update the skin for the DE-FE contact
            self.ComputeSkinSubModelPart()
            if self.DEMFEM_contact:
                self.TransferFEMSkinToDEM()


#RemoveIsolatedFiniteElements============================================================================================================================
    def RemoveIsolatedFiniteElements(self):
        FEM_Elements = self.FEM_Solution.main_model_part.Elements
        FEM_Nodes    = self.FEM_Solution.main_model_part.Nodes

        for node in FEM_Nodes:
            node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, 0)

        for Element in FEM_Elements:
            is_active = True
            if Element.IsDefined(KratosMultiphysics.ACTIVE):
                is_active = Element.Is(KratosMultiphysics.ACTIVE)

            if is_active == True:
                for i in range(0,3): # Loop over nodes of the element
                    node = Element.GetNodes()[i]
                    number_active_elements = node.GetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS)
                    number_active_elements += 1
                    node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, number_active_elements)

        for Element in FEM_Elements:
            total_elements_on_nodes = 0
            for i in range(0,3): # Loop over nodes of the element
                node = Element.GetNodes()[i]
                number_active_elements = node.GetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS)
                total_elements_on_nodes = total_elements_on_nodes + number_active_elements
            if total_elements_on_nodes == 3:
                Element.Set(KratosMultiphysics.TO_ERASE, True)

#ExtrapolatePressureLoad============================================================================================================================
    def ExtrapolatePressureLoad(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: ExtrapolatePressureLoad")

        if self.PressureLoad:
            # we reconstruct the pressure load if necessary
            if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECONSTRUCT_PRESSURE_LOAD] == 1:
                self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.INTERNAL_PRESSURE_ITERATION] = 1
                while self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.INTERNAL_PRESSURE_ITERATION] > 0:
                    if self.domain_size == 2:
                        KratosFemDem.ExtendPressureConditionProcess2D(self.FEM_Solution.main_model_part).Execute()
                    else:
                        KratosFemDem.ExtendPressureConditionProcess3D(self.FEM_Solution.main_model_part).Execute()

#UpdateDEMVariables============================================================================================================================
    def UpdateDEMVariables(self):
        update_de_kinematics_process = KratosFemDem.UpdateDemKinematicsProcess(self.FEM_Solution.main_model_part, self.SpheresModelPart)
        update_de_kinematics_process.Execute()

#TransferNodalForcesToFEM============================================================================================================================
    def TransferNodalForcesToFEM(self):
        tranfer_nodal_forces_process = KratosFemDem.TransferNodalForcesToFem(self.FEM_Solution.main_model_part, self.SpheresModelPart)
        tranfer_nodal_forces_process.Execute()

#WritePostListFile============================================================================================================================
    def WritePostListFile(self):
        post_file_name = self.FEM_Solution.problem_name + ".post.lst"
        time_label = round(self.FEM_Solution.step, 0)
        PostListFile = open(post_file_name, "w")
        PostListFile.write("Merge\n\n")
        PostListFile.write(self.FEM_Solution.problem_name + "_" + str(time_label) + ".post.res\n")
        PostListFile.write(self.FEM_Solution.problem_name + "_" + str(time_label) + ".post.msh\n")
        PostListFile.write(os.path.join(self.FEM_Solution.problem_name + "_Post_Files", self.FEM_Solution.problem_name + "_" + str(time_label) + ".post.bin"))
        PostListFile.close()

#InitializePlotsFiles============================================================================================================================
    def InitializePlotsFiles(self):
        # open general Displ/Reaction File
        self.PlotFile = open("PlotFile.txt","w")
        self.PlotFile.write("This File Plots the SUM of the displacement and reactions of the nodes selected in the lists!\n\n")
        self.PlotFile.write("       time           displ_x        displ_y      Reaction_x     Reaction_y    \n")
        self.PlotFile.close()
        self.TimePreviousPlotting = 0.0
        self.plot_files_nodes_list    = []
        self.plot_files_elements_list = []
        self.plot_files_nodes_id_list    = []
        self.plot_files_elements_id_list = []

        # open plots for nodes selected
        if self.FEM_Solution.ProjectParameters["watch_nodes_list"].size() != 0:
            number_nodes = self.FEM_Solution.ProjectParameters["watch_nodes_list"].size()
            for node in range(0, number_nodes):

                Id = self.FEM_Solution.ProjectParameters["watch_nodes_list"][node].GetInt()
                i_plot_file_node = open("PlotNode_" + str(Id) + ".txt","w")
                i_plot_file_node.write("\n")
                if self.domain_size == 2:
                    i_plot_file_node.write("       time          displ_x        displ_y         vel_x           vel_y         acc_x          acc_y        Reaction_x     Reaction_y    \n")
                else:
                    i_plot_file_node.write("       time          displ_x        displ_y        displ_z         vel_x           vel_y         vel_z           acc_x          acc_y          acc_z       Reaction_x     Reaction_y     Reaction_Z    \n")
                i_plot_file_node.close()
                self.plot_files_nodes_list.append(i_plot_file_node)
                self.plot_files_nodes_id_list.append(Id)

        # open plots for elements selected
        if self.FEM_Solution.ProjectParameters["watch_elements_list"].size() != 0:

            number_elems = self.FEM_Solution.ProjectParameters["watch_elements_list"].size()
            for elem in range(0, number_elems):
                Id = self.FEM_Solution.ProjectParameters["watch_elements_list"][elem].GetInt()
                i_plot_file_elem = open("PlotElement_" + str(Id) + ".txt","w")
                i_plot_file_elem.write("\n")
                if self.domain_size == 2:
                    i_plot_file_elem.write("          time                       Sxx                   Syy                      Sxy                    Exx                     Eyy                   Exy                Damage  \n")
                else:
                    i_plot_file_elem.write("       time             Sxx           Syy             Szz           Sxy            Syz            Sxz            Exx            Eyy            Ezz             Exy           Eyz            Exz          Damage  \n")
                i_plot_file_elem.close()
                self.plot_files_elements_list.append(i_plot_file_elem)
                self.plot_files_elements_id_list.append(Id)

#PrintPlotsFiles============================================================================================================================
    def PrintPlotsFiles(self):
        # Print the general file
        time = self.FEM_Solution.time
        total_reaction_x     = 0.0
        total_displacement_x = 0.0
        total_reaction_y     = 0.0
        total_displacement_y = 0.0
        total_displacement_z = 0.0
        total_reaction_z     = 0.0
        interval = self.FEM_Solution.ProjectParameters["interval_of_watching"].GetDouble()

        if self.FEM_Solution.time - self.TimePreviousPlotting >= interval:
            if self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size() > 0:
                if self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][0].IsInt():
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size()):
                        id_node = self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][index].GetInt()
                        node = self.FEM_Solution.main_model_part.GetNode(id_node)
                        total_displacement_x += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                        total_displacement_y += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                        total_displacement_z += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
                else:
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size()):
                        submodel_name = self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][index].GetString()
                        for node in self.FEM_Solution.main_model_part.GetSubModelPart(submodel_name).Nodes:
                            total_displacement_x += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                            total_displacement_y += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                            total_displacement_z += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)

                if self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][0].IsInt():
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"].size()):
                        id_node = self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][index].GetInt()
                        node = self.FEM_Solution.main_model_part.GetNode(id_node)
                        total_reaction_x += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
                        total_reaction_y += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y)
                        total_reaction_z += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Z)
                else:
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"].size()):
                        submodel_name = self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][index].GetString()
                        for node in self.FEM_Solution.main_model_part.GetSubModelPart(submodel_name).Nodes:
                            total_reaction_x += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
                            total_reaction_y += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y)
                            total_reaction_z += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Z)

                self.PlotFile = open("PlotFile.txt","a")
                if self.domain_size == 2:
                    self.PlotFile.write("    " + "{0:.4e}".format(time).rjust(11) + "    " + "{0:.4e}".format(total_displacement_x).rjust(11) +
                                        "    " + "{0:.4e}".format(total_displacement_y).rjust(11) + "    " + "{0:.4e}".format(total_reaction_x).rjust(11) +
                                        "    " + "{0:.4e}".format(total_reaction_y).rjust(11) + "\n")
                else:
                    self.PlotFile.write("    " + "{0:.4e}".format(time).rjust(11) + "    " + "{0:.4e}".format(total_displacement_x).rjust(11) +
                        "    " + "{0:.4e}".format(total_displacement_y).rjust(11) + "    " + "{0:.4e}".format(total_displacement_z).rjust(11) +
                        "    " + "{0:.4e}".format(total_reaction_x).rjust(11) + "    " + "{0:.4e}".format(total_reaction_y).rjust(11) + "    " +
                        "{0:.4e}".format(total_reaction_z).rjust(11) + "\n")
                self.PlotFile.close()


            # Print the selected nodes files
            if self.FEM_Solution.ProjectParameters["watch_nodes_list"].size() != 0:
                NumNodes = self.FEM_Solution.ProjectParameters["watch_nodes_list"].size()
                for inode in range(0, NumNodes):
                    id_node = self.PlotFilesNodesIdList[inode]
                    node = self.FEM_Solution.main_model_part.GetNode(id_node)
                    self.PlotFilesNodesList[inode] = open("PlotNode_" + str(id_node) + ".txt","a")

                    displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
                    velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                    reaction = node.GetSolutionStepValue(KratosMultiphysics.REACTION)
                    acceleration = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION)

                    dx = displacement[0]
                    dy = displacement[1]
                    dz = displacement[2]
                    Rx = reaction[0]
                    Ry = reaction[1]
                    Rz = reaction[2]
                    vx = velocity[0]
                    vy = velocity[1]
                    vz = velocity[2]
                    ax = acceleration[0]
                    ay = acceleration[1]
                    az = acceleration[2]

                    if self.domain_size == 2:
                        self.PlotFilesNodesList[inode].write("    " + "{0:.4e}".format(time).rjust(11) + "    " +
                            "{0:.4e}".format(dx).rjust(11) + "    " + "{0:.4e}".format(dy).rjust(11) + "    " +
                            "{0:.4e}".format(vx).rjust(11) + "    " + "{0:.4e}".format(vy).rjust(11) + "    " +
                            "{0:.4e}".format(ax).rjust(11) + "    " + "{0:.4e}".format(ay).rjust(11) + "    " +
                            "{0:.4e}".format(Rx).rjust(11) + "    " + "{0:.4e}".format(Ry).rjust(11) + "\n")
                    else:
                        self.PlotFilesNodesList[inode].write("    " + "{0:.4e}".format(time).rjust(11) + "    " +
                            "{0:.4e}".format(dx).rjust(11) + "    " + "{0:.4e}".format(dy).rjust(11) + "    " + "{0:.4e}".format(dz).rjust(11) + "    " +
                            "{0:.4e}".format(vx).rjust(11) + "    " + "{0:.4e}".format(vy).rjust(11) + "    " + "{0:.4e}".format(vz).rjust(11) + "    " +
                            "{0:.4e}".format(ax).rjust(11) + "    " + "{0:.4e}".format(ay).rjust(11) + "    " + "{0:.4e}".format(az).rjust(11) + "    " +
                            "{0:.4e}".format(Rx).rjust(11) + "    " + "{0:.4e}".format(Ry).rjust(11) + "    " + "{0:.4e}".format(Rz).rjust(11) + "\n")

                    self.PlotFilesNodesList[inode].close()

            # print the selected element files
            if self.FEM_Solution.ProjectParameters["watch_elements_list"].size() != 0:
                NumElem = self.FEM_Solution.ProjectParameters["watch_elements_list"].size()
                for iElem in range(0, NumElem):
                    Idelem = self.PlotFilesElementsIdList[iElem]
                    Elem = self.FEM_Solution.main_model_part.GetElement(Idelem)
                    self.PlotFilesElementsList[iElem] = open("PlotElement_" + str(Idelem) + ".txt","a")

                    stress_tensor = Elem.GetValuesOnIntegrationPoints(KratosFemDem.STRESS_VECTOR_INTEGRATED, self.FEM_Solution.main_model_part.ProcessInfo)
                    strain_vector = Elem.GetValue(KratosFemDem.STRAIN_VECTOR)

                    damage = Elem.GetValue(KratosFemDem.DAMAGE_ELEMENT)
                    
                    if self.domain_size == 2:
                        Sxx = stress_tensor[0][0]
                        Syy = stress_tensor[0][1]
                        Sxy = stress_tensor[0][2]
                        Exx = strain_vector[0]
                        Eyy = strain_vector[1]
                        Exy = strain_vector[2]
                        self.PlotFilesElementsList[iElem].write("    " + "{0:.4e}".format(time).rjust(11) + "    " +
                            "{0:.4e}".format(Sxx).rjust(11) + "    " + "{0:.4e}".format(Syy).rjust(11) + "    " +
                            "{0:.4e}".format(Sxy).rjust(11) + "    " + "{0:.4e}".format(Exx).rjust(11) +
                            "    " + "{0:.4e}".format(Eyy).rjust(11) + "    " + "{0:.4e}".format(Exy).rjust(11) +
                            "   " + "{0:.4e}".format(damage).rjust(11) + "\n")
                    else:
                        Sxx = stress_tensor[0][0]
                        Syy = stress_tensor[0][1]
                        Szz = stress_tensor[0][2]
                        Sxy = stress_tensor[0][3]
                        Syz = stress_tensor[0][4]
                        Sxz = stress_tensor[0][5]
                        Exx = strain_tensor[0]
                        Eyy = strain_tensor[1]
                        Ezz = strain_tensor[2]
                        Exy = strain_tensor[3]
                        Eyz = strain_tensor[4]
                        Exz = strain_tensor[5]
                        self.PlotFilesElementsList[iElem].write("    " + "{0:.4e}".format(time).rjust(11) + "    " +
                        "{0:.4e}".format(Sxx).rjust(11) + "    " + "{0:.4e}".format(Syy).rjust(11) + "    " +
                        "{0:.4e}".format(Szz).rjust(11) + "    " + "{0:.4e}".format(Sxy).rjust(11) + "    " +
                        "{0:.4e}".format(Syz).rjust(11) + "    " + "{0:.4e}".format(Sxz).rjust(11) + "    " +
                        "{0:.4e}".format(Exx).rjust(11) +
                        "    " + "{0:.4e}".format(Eyy).rjust(11) + "    " + "{0:.4e}".format(Ezz).rjust(11) +
                        "    " + "{0:.4e}".format(Exy).rjust(11) + "    " + "{0:.4e}".format(Eyz).rjust(11) +
                        "    " + "{0:.4e}".format(Exz).rjust(11) +
                        "   "  + "{0:.4e}".format(damage).rjust(11) + "\n")

                    self.PlotFilesElementsList[iElem].close()
            self.TimePreviousPlotting = time

#CountErasedVolume===================================================================================================================================
    def CountErasedVolume(self):
        count_erased_vol = True
        if count_erased_vol:
            erased_vol_process = KratosFemDem.ComputeSandProduction(self.FEM_Solution.main_model_part)
            erased_vol_process.Execute()

            self.ErasedVolume = open("ErasedVolume.txt","a")
            erased_vol = self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.ERASED_VOLUME]
            self.ErasedVolume.write("    " + "{0:.4e}".format(self.FEM_Solution.time).rjust(11) + "    " + "{0:.4e}".format(erased_vol).rjust(11) + "\n")
            self.ErasedVolume.close()

#GetMaximumConditionId============================================================================================================================
    def GetMaximumConditionId(self):
        max_id = 0
        for condition in self.FEM_Solution.main_model_part.Conditions:
            if condition.Id > max_id:
                max_id = condition.Id
        return max_id

#PrintDEMResults============================================================================================================================
    def PrintDEMResults(self):
        if self.DEM_Solution.step == 1: # always print the 1st step
            self.DEM_Solution.PrintResultsForGid(self.DEM_Solution.time)
            self.DEM_Solution.time_old_print = self.DEM_Solution.time

        else:
            time_to_print = self.DEM_Solution.time - self.DEM_Solution.time_old_print

            if self.DEM_Solution.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self.DEM_Solution.solver.dt:
                self.DEM_Solution.PrintResultsForGid(self.DEM_Solution.time)
                self.DEM_Solution.time_old_print = self.DEM_Solution.time

#CheckIfHasRemeshed============================================================================================================================
    def CheckIfHasRemeshed(self):
        execute_remesh = False
        step = self.RemeshingProcessMMG.step
        if not self.RemeshingProcessMMG.remesh_executed:
            if not self.RemeshingProcessMMG.initial_remeshing:
                # We need to check if the model part has been modified recently
                if self.RemeshingProcessMMG.main_model_part.Is(KratosMultiphysics.MODIFIED):
                    step = 0  # Reset (just to be sure)
                else:
                    step += 1
                    if self.RemeshingProcessMMG.step_frequency > 0:
                        if self.RemeshingProcessMMG.main_model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.RemeshingProcessMMG.initial_step:
                            if not self.RemeshingProcessMMG.initial_step_done:
                                    execute_remesh = True
                            else:
                                if step >= self.RemeshingProcessMMG.step_frequency:
                                    execute_remesh = True
        return execute_remesh

#RemoveDummyNodalForces============================================================================================================================
    def RemoveDummyNodalForces(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: RemoveDummyNodalForces")

        for condition in self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").Conditions:
            condition.Set(KratosMultiphysics.TO_ERASE, True)

        self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.FEM_Solution.main_model_part.RemoveSubModelPart("ContactForcesDEMConditions")

#GenerateDemAfterRemeshing============================================================================================================================
    def GenerateDemAfterRemeshing(self):
        # we extrapolate the damage to the nodes
        KratosFemDem.DamageToNodesProcess(self.FEM_Solution.main_model_part, 2).Execute()

        # we create a submodelpart containing the nodes and radius of the corresponding DEM
        KratosFemDem.DemAfterRemeshIdentificatorProcess(self.FEM_Solution.main_model_part, 0.95).Execute()

        # Loop over the elements of the Submodelpart to create the DEM
        for node in self.FEM_Solution.main_model_part.GetSubModelPart("DemAfterRemeshingNodes").Nodes:
            Id = node.Id
            R = node.GetValue(KratosFemDem.DEM_RADIUS)
            Coordinates = self.GetNodeCoordinates(node)
            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates, R, Id)
            node.SetValue(KratosFemDem.IS_DEM, True)

#RemoveAloneDEMElements============================================================================================================================
    def RemoveAloneDEMElements(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: RemoveAloneDEMElements")

        remove_alone_DEM_elements_process = KratosFemDem.RemoveAloneDEMElementsProcess(
                                                         self.FEM_Solution.main_model_part, 
                                                         self.SpheresModelPart)
        remove_alone_DEM_elements_process.Execute()

#ExecuteBeforeGeneratingDEM============================================================================================================================
    def ExecuteBeforeGeneratingDEM(self): 
        """Here the erased are labeled as INACTIVE so you can access to them. After calling
           GenerateDEM they are totally erased """
        if self.PressureLoad:
            self.ExpandWetNodes()
            KratosFemDem.UpdatePressureVolumeProcess(self.FEM_Solution.main_model_part).Execute()
            self.ExpandWetNodes()

#ExecuteAfterGeneratingDEM============================================================================================================================
    def ExecuteAfterGeneratingDEM(self):
        self.ExtrapolatePressureLoad()
        self.SpheresModelPart = self.ParticleCreatorDestructor.GetSpheresModelPart()
        # We update coordinates, displ and velocities of the DEM according to FEM
        self.UpdateDEMVariables()

#BeforeSolveDEMOperations============================================================================================================================
    def BeforeSolveDEMOperations(self):
        self.DEM_Solution.time = self.FEM_Solution.time
        self.DEM_Solution.step = self.FEM_Solution.step
        self.DEM_Solution.DEMFEMProcedures.UpdateTimeInModelParts(self.DEM_Solution.all_model_parts,
                                                                   self.DEM_Solution.time,
                                                                   self.DEM_Solution.solver.dt,
                                                                   self.DEM_Solution.step,
                                                                   self.DEM_Solution.IsTimeToPrintPostProcess())

#TransferFEMSkinToDEM============================================================================================================================
    def TransferFEMSkinToDEM(self):
        fem_skin_mp = self.FEM_Solution.main_model_part.GetSubModelPart("SkinDEMModelPart")

        if self.DEM_Solution.rigid_face_model_part.HasSubModelPart("SkinTransferredFromStructure"):
            self.EraseConditionsAndNodesSubModelPart()
            dem_walls_mp = self.DEM_Solution.rigid_face_model_part.GetSubModelPart("SkinTransferredFromStructure")
            props = self.DEM_Solution.rigid_face_model_part.GetProperties(self.created_props_id,0)
            DemFem.DemStructuresCouplingUtilities().TransferStructuresSkinToDem(fem_skin_mp, dem_walls_mp, props)
        else: # have to create it
            props = self.CreateFEMPropertiesForDEFEContact()
            dem_walls_mp = self.DEM_Solution.rigid_face_model_part.CreateSubModelPart("SkinTransferredFromStructure")
            dem_walls_mp.AddProperties(props)
            DemFem.DemStructuresCouplingUtilities().TransferStructuresSkinToDem(fem_skin_mp, dem_walls_mp, props)
    
    #-----------------------------------
    def EraseConditionsAndNodesSubModelPart(self):
        DEM_sub_model_part = self.DEM_Solution.rigid_face_model_part.GetSubModelPart("SkinTransferredFromStructure")
        self.DEM_Solution.rigid_face_model_part.Conditions.clear()
        self.DEM_Solution.rigid_face_model_part.Nodes.clear()
    #-----------------------------------
    def CreateFEMPropertiesForDEFEContact(self):
        max_id_properties = 0
        young = 0
        poisson = 0
        for prop in self.FEM_Solution.main_model_part.Properties:
            young = prop[KratosMultiphysics.YOUNG_MODULUS]
            poisson = prop[KratosMultiphysics.POISSON_RATIO]
            if max_id_properties < prop.Id:
                max_id_properties = prop.Id
        props = KratosMultiphysics.Properties(max_id_properties + 1)
        self.created_props_id = max_id_properties + 1
        props[KratosDEM.FRICTION] =  -0.5773502691896257
        props[KratosDEM.WALL_COHESION] = 0.0
        props[KratosDEM.COMPUTE_WEAR] = False
        props[KratosDEM.SEVERITY_OF_WEAR] = 0.001
        props[KratosDEM.IMPACT_WEAR_SEVERITY] = 0.001
        props[KratosDEM.BRINELL_HARDNESS] = 200.0
        props[KratosMultiphysics.YOUNG_MODULUS] = young # the PENALTY
        props[KratosMultiphysics.POISSON_RATIO] = poisson
        return props

#============================================================================================================================
