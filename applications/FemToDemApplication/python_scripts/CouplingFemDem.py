from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import MainDEM_for_coupling as DEM
import MainFEM_for_coupling as FEM
import FEMDEMParticleCreatorDestructor as PCD
import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.MeshingApplication as MeshingApplication
import KratosMultiphysics.SolidMechanicsApplication as Solid
import math
import os
import mmg_process as MMG

def Wait():
    input("Press Something")


# Main script of the coupled FEM-DEM Application 2D
class FEMDEM_Solution:
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

#============================================================================================================================
    def Run(self):

        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

#============================================================================================================================
    def Initialize(self):
        self.number_of_nodes_element = 3
        self.FEM_Solution.Initialize()
        self.DEM_Solution.Initialize()

        # Initialize the "flag" IS_DEM in all the nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosFemDem.IS_DEM, False, self.FEM_Solution.main_model_part.Nodes)
        # Initialize the "flag" NODAL_FORCE_APPLIED in all the nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosFemDem.NODAL_FORCE_APPLIED, False, self.FEM_Solution.main_model_part.Nodes)
        # Initialize the "flag" RADIUS in all the nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.RADIUS, False, self.FEM_Solution.main_model_part.Nodes)

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

        if self.FEM_Solution.ProjectParameters.Has("pressure_load_extrapolation") == False:
            self.PressureLoad = False
        else:
            self.PressureLoad = self.FEM_Solution.ProjectParameters["pressure_load_extrapolation"].GetBool()
        if self.PressureLoad:
            KratosFemDem.AssignPressureIdProcess(self.FEM_Solution.main_model_part).Execute()
        
        self.SkinDetectionProcessParameters = KratosMultiphysics.Parameters("""
        {
            "name_auxiliar_model_part" : "SkinDEMModelPart",
            "name_auxiliar_condition"  : "Condition",
            "echo_level"               : 0
        }""")

        # for the dem contact forces coupling
        self.InitializeDummyNodalForces()

        KratosMultiphysics.Logger.PrintInfo(" /$$$$$$$$ /$$$$$$$$ /$$      /$$  /$$$$$$  /$$$$$$$  /$$$$$$$$ /$$      /$$")
        KratosMultiphysics.Logger.PrintInfo("| $$_____/| $$_____/| $$$    /$$$ /$$__  $$| $$__  $$| $$_____/| $$$    /$$$")
        KratosMultiphysics.Logger.PrintInfo("| $$      | $$      | $$$$  /$$$$|__/  \ $$| $$  \ $$| $$      | $$$$  /$$$$")
        KratosMultiphysics.Logger.PrintInfo("| $$$$$   | $$$$$   | $$ $$/$$ $$  /$$$$$$/| $$  | $$| $$$$$   | $$ $$/$$ $$")
        KratosMultiphysics.Logger.PrintInfo("| $$__/   | $$__/   | $$  $$$| $$ /$$____/ | $$  | $$| $$__/   | $$  $$$| $$")
        KratosMultiphysics.Logger.PrintInfo("| $$      | $$      | $$\  $ | $$| $$      | $$  | $$| $$      | $$\  $ | $$")
        KratosMultiphysics.Logger.PrintInfo("| $$      | $$$$$$$$| $$ \/  | $$| $$$$$$$$| $$$$$$$/| $$$$$$$$| $$ \/  | $$")
        KratosMultiphysics.Logger.PrintInfo("|__/      |________/|__/     |__/|________/|_______/ |________/|__/     |__/ 2D Application")

#============================================================================================================================
    def RunMainTemporalLoop(self):

        # Solving the problem (time integration)
        self.DEM_Solution.step           = 0
        self.DEM_Solution.time           = 0.0
        self.DEM_Solution.time_old_print = 0.0

        if self.DoRemeshing:
            self.RemeshingProcessMMG.ExecuteBeforeSolutionLoop()

        while(self.FEM_Solution.time <= self.FEM_Solution.end_time):
            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()

#============================================================================================================================
    def InitializeSolutionStep(self):

        # modified for the remeshing
        self.FEM_Solution.delta_time = self.ComputeDeltaTime()
        self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.FEM_Solution.delta_time
        self.FEM_Solution.time = self.FEM_Solution.time + self.FEM_Solution.delta_time
        self.FEM_Solution.main_model_part.CloneTimeStep(self.FEM_Solution.time)
        self.FEM_Solution.step = self.FEM_Solution.step + 1
        self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.FEM_Solution.step

        neighbour_elemental_finder =  KratosMultiphysics.FindElementalNeighboursProcess(self.FEM_Solution.main_model_part, 2, 5)
        neighbour_elemental_finder.Execute()

        if self.DoRemeshing:
            is_remeshing = self.CheckIfHasRemeshed()

            if is_remeshing:
                # Extrapolate the free energy as a remeshing criterion
                KratosFemDem.ComputeNormalizedFreeEnergyOnNodesProcess(self.FEM_Solution.main_model_part, 2).Execute()

                # we eliminate the nodal DEM forces
                self.RemoveDummyNodalForces()

            # Perform remeshing
            self.RemeshingProcessMMG.ExecuteInitializeSolutionStep()

            if is_remeshing:
                self.InitializeSolutionAfterRemeshing()
                neighbour_elemental_finder =  KratosMultiphysics.FindElementalNeighboursProcess(self.FEM_Solution.main_model_part, 2, 5)
                neighbour_elemental_finder.ClearNeighbours()
                neighbour_elemental_finder.Execute()

        self.FEM_Solution.InitializeSolutionStep()

        # Create initial skin of DEM's
        self.create_initial_dem_skin = False  # Hard Coded TODO
        if self.create_initial_dem_skin and self.FEM_Solution.step == 1:
            self.CreateInitialSkinDEM()

#============================================================================================================================
    def SolveSolutionStep(self): # Function to perform the coupling FEM <-> DEM

        self.FEM_Solution.clock_time = self.FEM_Solution.StartTimeMeasuring()

        #### SOLVE FEM #########################################
        self.FEM_Solution.solver.Solve()
        ########################################################

        if self.PressureLoad:
            # This must be called before Generating DEM
            self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECONSTRUCT_PRESSURE_LOAD] = 0 # It is modified inside
            extend_wet_nodes_process = KratosFemDem.ExpandWetNodesProcess(self.FEM_Solution.main_model_part)
            extend_wet_nodes_process.Execute()

        # we create the new DEM of this time step
        self.GenerateDEM()

        if self.PressureLoad:
            # we reconstruct the pressure load if necessary
            if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECONSTRUCT_PRESSURE_LOAD] == 1:
                self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.INTERNAL_PRESSURE_ITERATION] = 1
                while self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.INTERNAL_PRESSURE_ITERATION] > 0:
                    KratosFemDem.ExtendPressureConditionProcess2D(self.FEM_Solution.main_model_part).Execute()
            
        self.SpheresModelPart = self.ParticleCreatorDestructor.GetSpheresModelPart()
        self.CheckForPossibleIndentations()

        # We update coordinates, displ and velocities of the DEM according to FEM
        self.UpdateDEMVariables()

        self.DEM_Solution.InitializeTimeStep()
        self.DEM_Solution.time = self.FEM_Solution.time
        self.DEM_Solution.step = self.FEM_Solution.step
        self.DEM_Solution.DEMFEMProcedures.UpdateTimeInModelParts(self.DEM_Solution.all_model_parts,
                                                                   self.DEM_Solution.time,
                                                                   self.DEM_Solution.solver.dt,
                                                                   self.DEM_Solution.step,
                                                                   self.DEM_Solution.IsTimeToPrintPostProcess())
        self.DEM_Solution._BeforeSolveOperations(self.DEM_Solution.time)

        #### SOLVE DEM #########################################
        self.DEM_Solution.solver.Solve()
        ########################################################

        self.DEM_Solution.AfterSolveOperations()
        self.DEM_Solution.solver._MoveAllMeshes(self.DEM_Solution.time, self.DEM_Solution.solver.dt)
        
        # to print DEM with the FEM coordinates
        self.UpdateDEMVariables()

        # DEM GiD print output
        self.PrintDEMResults()

        self.DEM_Solution.FinalizeTimeStep(self.DEM_Solution.time)

        # Transfer the contact forces of the DEM to the FEM nodes
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

        # Remove the submodel to be recomputed at each dt
        if self.FEM_Solution.main_model_part.HasSubModelPart("SkinDEMModelPart"):
            for cond in self.FEM_Solution.main_model_part.GetSubModelPart("SkinDEMModelPart").Conditions:
                cond.Set(KratosMultiphysics.TO_ERASE)

        self.FEM_Solution.main_model_part.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.FEM_Solution.main_model_part.RemoveSubModelPart("SkinDEMModelPart")


#============================================================================================================================
    def Finalize(self):

        self.FEM_Solution.Finalize()
        self.DEM_Solution.Finalize()
        self.DEM_Solution.CleanUpOperations()

        if self.DoRemeshing:
            self.RemeshingProcessMMG.ExecuteFinalize()


#============================================================================================================================
    def GenerateDEM(self): # This method creates the DEM elements and remove the damaged FEM, Additionally remove the isolated elements

        FEM_elements = self.FEM_Solution.main_model_part.Elements

        # Loop Over Elements to find the INACTIVE ones and generate the DEM only once
        for Element in FEM_elements:
            is_active     = True
            DEM_Generated = Element.GetValue(KratosFemDem.DEM_GENERATED)

            if Element.IsDefined(KratosMultiphysics.ACTIVE):
                is_active = Element.Is(KratosMultiphysics.ACTIVE)

                NumberOfDEM = 0         # Number of nodes with DEM Associated
                for node in range(0, self.number_of_nodes_element): # Loop over nodes of the FE
                    Node = Element.GetNodes()[node]
                    if Node.GetValue(KratosFemDem.IS_DEM) == True:
                        NumberOfDEM += 1

                if is_active == False and DEM_Generated == False: # Let's generate the remaining DEM of this FE
                    # print("elemento borrado: ", Element.Id)

                    # --------------------- 1ST SCENARIO -----------------------------
                    if NumberOfDEM == 0: # we must create 3 DEM

                        dist01  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[1])
                        dist02  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[2])
                        dist12  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[1], Element.GetNodes()[2])

                        # look to the node 1 --------------
                        Radius1 = self.GetMinimumValue(dist01, dist02) * 0.5
                        Coordinates1 = self.GetNodeCoordinates(Element.GetNodes()[0])
                        Id1 = Element.GetNodes()[0].Id

                        self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, Radius1, Id1)
                        Element.GetNodes()[0].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                        Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, Radius1)

                        # look to the node 2 --------------
                        # Radius2 = self.GetMinimumValue(dist01-Radius1, dist12*0.5)
                        Radius2 = dist01-Radius1
                        Coordinates2 = self.GetNodeCoordinates(Element.GetNodes()[1])
                        Id2 = Element.GetNodes()[1].Id

                        self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, Radius2, Id2)
                        Element.GetNodes()[1].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                        Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, Radius2)

                        # look to the node 3 --------------
                        Radius3 = self.GetMinimumValue(dist02-Radius1, dist12-Radius2)
                        Coordinates3 = self.GetNodeCoordinates(Element.GetNodes()[2])
                        Id3 = Element.GetNodes()[2].Id

                        self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates3, Radius3, Id3)
                        Element.GetNodes()[2].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                        Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, Radius3)

                        # DEM generated for this Element
                        Element.SetValue(KratosFemDem.DEM_GENERATED, True)
                        Element.Set(KratosMultiphysics.TO_ERASE, True)

                    # --------------------- 2ND SCENARIO -----------------------------
                    elif NumberOfDEM == 2: # we must create 1 DEM

                        dist01  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[1])
                        dist02  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[2])
                        dist12  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[1], Element.GetNodes()[2])

                        localId = 0 # Local Id of the node without DEM
                        for index in range(0,3):
                            if Element.GetNodes()[index].GetValue(KratosFemDem.IS_DEM) == False:
                                localId = index
                                break

                        Coordinates = self.GetNodeCoordinates(Element.GetNodes()[localId])
                        Id = Element.GetNodes()[localId].Id

                        if localId == 0:
                            R1 = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)
                            R2 = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)
                            R0 = self.GetMinimumValue(dist01-R1, dist02-R2)
                            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates, R0, Id)
                            Element.GetNodes()[index].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                            Element.GetNodes()[index].SetValue(KratosMultiphysics.RADIUS, R0)

                        elif localId == 1:
                            R0 = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
                            R2 = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)
                            R1 = self.GetMinimumValue(dist01-R0, dist12-R2)
                            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates, R1, Id)
                            Element.GetNodes()[index].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                            Element.GetNodes()[index].SetValue(KratosMultiphysics.RADIUS, R1)

                        elif localId == 2:
                            R1 = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)
                            R0 = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
                            R2 = self.GetMinimumValue(dist02-R0, dist12-R1)
                            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates, R2, Id)
                            Element.GetNodes()[index].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                            Element.GetNodes()[index].SetValue(KratosMultiphysics.RADIUS, R2)

                        # DEM generated for this Element
                        Element.SetValue(KratosFemDem.DEM_GENERATED, True)
                        Element.Set(KratosMultiphysics.TO_ERASE, True)

                    # --------------------- 3RD SCENARIO -----------------------------
                    elif NumberOfDEM == 1: # we must create 2 DEM

                        dist01  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[1])
                        dist02  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[2])
                        dist12  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[1], Element.GetNodes()[2])

                        localId = 0 # Local Id of the node with DEM
                        for index in range(0,3):
                            if Element.GetNodes()[index].GetValue(KratosFemDem.IS_DEM) == True:
                                localId = index
                                break

                        RadiusOfDem = Element.GetNodes()[localId].GetValue(KratosMultiphysics.RADIUS)

                        # --------------
                        if localId == 0:

                            # R1 = self.GetMinimumValue(dist01-RadiusOfDem, dist12*0.5)
                            R1 = dist01-RadiusOfDem
                            R2 = self.GetMinimumValue(dist02-RadiusOfDem, dist12-R1)

                            Coordinates1 = self.GetNodeCoordinates(Element.GetNodes()[1])
                            Id1 = Element.GetNodes()[1].Id
                            Coordinates2 = self.GetNodeCoordinates(Element.GetNodes()[2])
                            Id2 = Element.GetNodes()[2].Id

                            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, R1, Id1)
                            Element.GetNodes()[1].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                            Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

                            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, R2, Id2)
                            Element.GetNodes()[2].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                            Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

                        # --------------
                        elif localId == 1:

                            # R0 = self.GetMinimumValue(dist01-RadiusOfDem, dist02*0.5)
                            R0 = dist01-RadiusOfDem
                            R2 = self.GetMinimumValue(dist12-RadiusOfDem, dist02-R0)

                            Coordinates0 = self.GetNodeCoordinates(Element.GetNodes()[0])
                            Id0 = Element.GetNodes()[0].Id
                            Coordinates2 = self.GetNodeCoordinates(Element.GetNodes()[2])
                            Id2 = Element.GetNodes()[2].Id

                            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates0, R0, Id0)
                            Element.GetNodes()[0].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                            Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)

                            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, R2, Id2)
                            Element.GetNodes()[2].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                            Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

                        # --------------
                        elif localId == 2:

                            # R0 = self.GetMinimumValue(dist02-RadiusOfDem, dist01*0.5)
                            R0 = dist02-RadiusOfDem
                            R1 = self.GetMinimumValue(dist12-RadiusOfDem, dist01-R0)

                            Coordinates0 = self.GetNodeCoordinates(Element.GetNodes()[0])
                            Id0 = Element.GetNodes()[0].Id
                            Coordinates1 = self.GetNodeCoordinates(Element.GetNodes()[1])
                            Id1 = Element.GetNodes()[1].Id

                            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates0, R0, Id0)
                            Element.GetNodes()[0].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                            Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)

                            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, R1, Id1)
                            Element.GetNodes()[1].SetValue(KratosFemDem.IS_DEM, True)        # Has an asociated DEM now
                            Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

                        # DEM generated for this Element
                        Element.SetValue(KratosFemDem.DEM_GENERATED, True)
                        Element.Set(KratosMultiphysics.TO_ERASE, True)

                    # --------------------- 4TH SCENARIO -----------------------------
                    elif NumberOfDEM == 3: # We must avoid possible indentations

                        dist01  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[1])
                        dist02  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[2])
                        dist12  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[1], Element.GetNodes()[2])

                        R0  = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
                        R1  = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)
                        R2  = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)
                        Id0 = Element.GetNodes()[0].Id
                        Id1 = Element.GetNodes()[1].Id
                        Id2 = Element.GetNodes()[2].Id

                        # Check the 3 edges of the element
                        if R0 + R1 > dist01:
                            # R0 = self.GetMinimumValue(R0, dist01*0.5)
                            R1 = dist01 - R0

                            # assign the new radius to the DEM nodes
                            self.DEM_Solution.spheres_model_part.GetNode(Id0).SetSolutionStepValue(KratosMultiphysics.RADIUS, R0)
                            self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
                            Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)
                            Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

                        if R0 + R2 > dist02:
                            # R0 = self.GetMinimumValue(R0, 0.5*dist02)
                            R2 = dist02 - R0

                            # assign the new radius to the DEM nodes
                            self.DEM_Solution.spheres_model_part.GetNode(Id0).SetSolutionStepValue(KratosMultiphysics.RADIUS, R0)
                            self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
                            Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)
                            Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

                        if R1 + R2 > dist12:
                            # R1 = self.GetMinimumValue(R1, 0.5*dist12)
                            R2 = dist12 - R1

                            # assign the new radius to the DEM nodes
                            self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
                            self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
                            Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)
                            Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

                        # DEM generated for this Element
                        Element.SetValue(KratosFemDem.DEM_GENERATED, True)
                        Element.Set(KratosMultiphysics.TO_ERASE, True)

                    else:
                        raise Exception("Not possible")

                elif is_active == False and DEM_Generated == True:
                    Element.Set(KratosMultiphysics.TO_ERASE, True)

                # Remove the isolated Elements
                self.RemoveIsolatedFiniteElements()

        # We remove the inactive DEM associated to fem_nodes
        self.RemoveAloneDEMElements()
        element_eliminator = KratosMultiphysics.AuxiliarModelPartUtilities(self.FEM_Solution.main_model_part)
        element_eliminator.RemoveElementsAndBelongings(KratosMultiphysics.TO_ERASE)

#============================================================================================================================
    def CheckForPossibleIndentations(self): # Verifies if an element has indentations between its DEM

        FEM_elements = self.FEM_Solution.main_model_part.Elements

        for Element in FEM_elements:
            is_active     = True
            DEM_Generated = Element.GetValue(KratosFemDem.DEM_GENERATED)
            if Element.IsDefined(KratosMultiphysics.ACTIVE):
                is_active = Element.Is(KratosMultiphysics.ACTIVE)

            NumberOfDEM = 0         # Number of nodes with DEM Associated
            for node in range(0, 3): # Loop over nodes of the FE
                Node = Element.GetNodes()[node]

                if Node.GetValue(KratosFemDem.IS_DEM) == True:
                    NumberOfDEM += 1

            if NumberOfDEM == 3 and is_active == True and DEM_Generated == False:  # Case in which the DEM have been generated by its neighbours
                # Just avoid the initial indentations
                dist01  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[1])
                dist02  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[0], Element.GetNodes()[2])
                dist12  = self.CalculateDistanceBetweenNodes(Element.GetNodes()[1], Element.GetNodes()[2])
                R0  = Element.GetNodes()[0].GetValue(KratosMultiphysics.RADIUS)
                R1  = Element.GetNodes()[1].GetValue(KratosMultiphysics.RADIUS)
                R2  = Element.GetNodes()[2].GetValue(KratosMultiphysics.RADIUS)
                Id0 = Element.GetNodes()[0].Id
                Id1 = Element.GetNodes()[1].Id
                Id2 = Element.GetNodes()[2].Id

                # Check the 3 edges of the element
                if R0 + R1 > dist01:
                    R0 = self.GetMinimumValue(R0, dist01*0.5)
                    R1 = dist01 - R0
                    # assign the new radius to the DEM nodes
                    self.DEM_Solution.spheres_model_part.GetNode(Id0).SetSolutionStepValue(KratosMultiphysics.RADIUS, R0)
                    self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
                    Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)
                    Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

                if R0 + R2 > dist02:
                    R0 = self.GetMinimumValue(R0, 0.5*dist02)
                    R2 = dist02 - R0
                    # assign the new radius to the DEM nodes
                    self.DEM_Solution.spheres_model_part.GetNode(Id0).SetSolutionStepValue(KratosMultiphysics.RADIUS, R0)
                    self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
                    Element.GetNodes()[0].SetValue(KratosMultiphysics.RADIUS, R0)
                    Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)

                if R1 + R2 > dist12:
                    R1 = self.GetMinimumValue(R1, 0.5*dist12)
                    R2 = dist12 - R1
                    # assign the new radius to the DEM nodes
                    self.DEM_Solution.spheres_model_part.GetNode(Id1).SetSolutionStepValue(KratosMultiphysics.RADIUS, R1)
                    self.DEM_Solution.spheres_model_part.GetNode(Id2).SetSolutionStepValue(KratosMultiphysics.RADIUS, R2)
                    Element.GetNodes()[2].SetValue(KratosMultiphysics.RADIUS, R2)
                    Element.GetNodes()[1].SetValue(KratosMultiphysics.RADIUS, R1)

                # DEM generated for this Element
                Element.SetValue(KratosFemDem.DEM_GENERATED, True)

#============================================================================================================================
    def CalculateDistanceBetweenNodes(self, Node1, Node2):
        # only in 2D
        X1 = Node1.X
        X2 = Node2.X
        Y1 = Node1.Y
        Y2 = Node2.Y
        return math.sqrt((X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2))

    def GetMinimumValue(self, val1, val2):
        res = val1
        if val2 < val1:
            res = val2
        return res

    def GetNodeCoordinates(self, Node):
        X = Node.X
        Y = Node.Y
        Z = Node.Z
        coord = KratosMultiphysics.Array3()
        coord[0] = X
        coord[1] = Y
        coord[2] = Z
        return coord

#============================================================================================================================
    def UpdateDEMVariables(self):

        FEM_Nodes = self.FEM_Solution.main_model_part.Nodes
        for fem_node in FEM_Nodes:
            if fem_node.GetValue(KratosFemDem.IS_DEM):
                id_node = fem_node.Id
                associated_dem = self.SpheresModelPart.GetNode(id_node)

                Coordinates    = self.GetNodeCoordinates(fem_node)
                velocity = fem_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                displacement = fem_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)

                # Update Coordinates
                associated_dem.X = Coordinates[0]
                associated_dem.Y = Coordinates[1]

                # Update Displacements
                associated_dem.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, displacement[0])
                associated_dem.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, displacement[1])

                # Update Velocities
                associated_dem.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, velocity[0])
                associated_dem.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, velocity[1])

#============================================================================================================================
    def CheckInactiveNodes(self):

        FEM_Elements = self.FEM_Solution.main_model_part.Elements
        FEM_Nodes    = self.FEM_Solution.main_model_part.Nodes
        erased_nodes_id = []
        conditions_to_erase_id = []

        for node in FEM_Nodes:
            node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, 0)

        for Element in FEM_Elements:
            is_active = True
            if Element.IsDefined(KratosMultiphysics.ACTIVE):
                is_active = Element.Is(KratosMultiphysics.ACTIVE)

            if is_active == True:
                for i in range(0,3): # Loop over nodes of the element
                    node = Element.GetNodes()[i]
                    NumberOfActiveElements = node.GetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS)
                    NumberOfActiveElements += 1
                    node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, NumberOfActiveElements)

        NumberOfActiveElements = 0

        for node in FEM_Nodes:
            NumberOfActiveElements = node.GetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS)
            if NumberOfActiveElements == 0 and node.GetValue(KratosFemDem.INACTIVE_NODE) == False:
                Id = node.Id
                node.SetValue(KratosFemDem.INACTIVE_NODE, True)
                node.Set(KratosMultiphysics.TO_ERASE, True) # added
                DEMnode = self.SpheresModelPart.GetNode(Id)
                DEMnode.SetValue(KratosFemDem.INACTIVE_NODE, True)
                DEMnode.Set(KratosMultiphysics.TO_ERASE, True)
                erased_nodes_id.append(Id)

                for condition in self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").Conditions:
                    if condition.GetNodes()[0].Id == Id:
                        conditions_to_erase_id.append(condition.Id)

        # let's remove the nodal dem conditions according to inactive nodes
        for Id in conditions_to_erase_id:
            self.FEM_Solution.main_model_part.RemoveCondition(Id)

        # Remove inactive nodes
        self.SpheresModelPart.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.FEM_Solution.main_model_part.GetRootModelPart().RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE) # added

#============================================================================================================================
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

#============================================================================================================================

    def TransferNodalForcesToFEM(self):

        for condition in self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").Conditions:
            id_node = condition.GetNodes()[0].Id

            if self.FEM_Solution.main_model_part.GetNode(id_node).GetValue(KratosFemDem.IS_DEM):
                dem_forces = self.SpheresModelPart.GetNode(id_node).GetSolutionStepValue(KratosMultiphysics.TOTAL_FORCES)
                condition.SetValue(Solid.FORCE_LOAD, dem_forces)

#============================================================================================================================
    def WritePostListFile(self):

        post_file_name = self.FEM_Solution.problem_name + ".post.lst"
        time_label = round(self.FEM_Solution.step, 0)
        PostListFile = open(post_file_name, "w")
        PostListFile.write("Merge\n\n")
        PostListFile.write(self.FEM_Solution.problem_name + "_" + str(time_label) + ".post.res\n")
        PostListFile.write(self.FEM_Solution.problem_name + "_" + str(time_label) + ".post.msh\n")
        PostListFile.write(os.path.join(self.FEM_Solution.problem_name + "_Post_Files", self.FEM_Solution.problem_name + "_" + str(time_label) + ".post.bin"))
        PostListFile.close()

#============================================================================================================================
    def PrintPlotsFiles(self):

        # Print the general file
        time = self.FEM_Solution.time
        total_reaction_x     = 0.0
        total_displacement_x = 0.0
        total_reaction_y     = 0.0
        total_displacement_y = 0.0
        interval = self.FEM_Solution.ProjectParameters["interval_of_watching"].GetDouble()
        
        if self.FEM_Solution.time - self.TimePreviousPlotting >= interval:
            if self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size() > 0:
                if self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][0].IsInt():
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size()):
                        IdNode = self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][index].GetInt()
                        node = self.FEM_Solution.main_model_part.GetNode(IdNode)
                        total_displacement_x += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                        total_displacement_y += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                else:
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size()):
                        submodel_name = self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][index].GetString()
                        for node in self.FEM_Solution.main_model_part.GetSubModelPart(submodel_name).Nodes:
                            total_displacement_x += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                            total_displacement_y += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y) 

                if self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][0].IsInt():
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"].size()):
                        IdNode = self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][index].GetInt()
                        node = self.FEM_Solution.main_model_part.GetNode(IdNode)
                        total_reaction_x += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
                        total_reaction_y += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y)
                else:
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"].size()):
                        submodel_name = self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][index].GetString()
                        for node in self.FEM_Solution.main_model_part.GetSubModelPart(submodel_name).Nodes:
                            total_reaction_x += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
                            total_reaction_y += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y) 

                self.PlotFile = open("PlotFile.txt","a")
                self.PlotFile.write("    " + "{0:.4e}".format(time).rjust(11) + "    " + "{0:.4e}".format(total_displacement_x).rjust(11) +
                                    "    " + "{0:.4e}".format(total_displacement_y).rjust(11) + "    " + "{0:.4e}".format(total_reaction_x).rjust(11) +
                                    "    " + "{0:.4e}".format(total_reaction_y).rjust(11) + "\n")
                self.PlotFile.close()


            # Print the selected nodes files
            if self.FEM_Solution.ProjectParameters["watch_nodes_list"].size() != 0:
                NumNodes = self.FEM_Solution.ProjectParameters["watch_nodes_list"].size()
                for inode in range(0, NumNodes):
                    IdNode = self.PlotFilesNodesIdList[inode]
                    node = self.FEM_Solution.main_model_part.GetNode(IdNode)
                    self.PlotFilesNodesList[inode] = open("PlotNode_" + str(IdNode) + ".txt","a")

                    displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
                    velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                    reaction = node.GetSolutionStepValue(KratosMultiphysics.REACTION)
                    acceleration = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION)

                    dx = displacement[0]
                    dy = displacement[1]
                    Rx = reaction[0]
                    Ry = reaction[1]
                    vx = velocity[0]
                    vy = velocity[1]
                    ax = acceleration[0]
                    ay = acceleration[1]

                    self.PlotFilesNodesList[inode].write("    " + "{0:.4e}".format(time).rjust(11) + "    " +
                        "{0:.4e}".format(dx).rjust(11) + "    " + "{0:.4e}".format(dy).rjust(11) + "    " +
                        "{0:.4e}".format(vx).rjust(11) + "    " + "{0:.4e}".format(vy).rjust(11) + "    " +
                        "{0:.4e}".format(ax).rjust(11) + "    " + "{0:.4e}".format(ay).rjust(11) + "    " +
                        "{0:.4e}".format(Rx).rjust(11) + "    " + "{0:.4e}".format(Ry).rjust(11) + "\n")

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

                    Sxx = stress_tensor[0][0]
                    Syy = stress_tensor[0][1]
                    Sxy = stress_tensor[0][2]

                    Exx = strain_vector[0]
                    Eyy = strain_vector[1]
                    Exy = strain_vector[2]

                    damage = Elem.GetValue(KratosFemDem.DAMAGE_ELEMENT)

                    self.PlotFilesElementsList[iElem].write("    " + "{0:.4e}".format(time).rjust(11) + "    " +
                        "{0:.4e}".format(Sxx).rjust(11) + "    " + "{0:.4e}".format(Syy).rjust(11) + "    " +
                        "{0:.4e}".format(Sxy).rjust(11) + "    " + "{0:.4e}".format(Exx).rjust(11) +
                        "    " + "{0:.4e}".format(Eyy).rjust(11) + "    " + "{0:.4e}".format(Exy).rjust(11) +
                        "   " + "{0:.4e}".format(damage).rjust(11) + "\n")

                    self.PlotFilesElementsList[iElem].close()
            self.TimePreviousPlotting = time

#============================================================================================================================
    def InitializePlotsFiles(self):

        # open general Displ/Reaction File
        self.PlotFile = open("PlotFile.txt","w")
        self.PlotFile.write("This File Plots the SUM of the displacement and reactions of the nodes selected in the lists!\n\n")
        self.PlotFile.write("       time           displ_x        displ_y      Reaction_x     Reaction_y    \n")
        self.PlotFile.close()
        self.TimePreviousPlotting = 0.0
        self.PlotFilesNodesList    = []
        self.PlotFilesElementsList = []
        self.PlotFilesNodesIdList    = []
        self.PlotFilesElementsIdList = []

        # open plots for nodes selected
        if self.FEM_Solution.ProjectParameters["watch_nodes_list"].size() != 0:
            NumNodes = self.FEM_Solution.ProjectParameters["watch_nodes_list"].size()
            for node in range(0, NumNodes):

                Id = self.FEM_Solution.ProjectParameters["watch_nodes_list"][node].GetInt()
                iPlotFileNode = open("PlotNode_" + str(Id) + ".txt","w")
                iPlotFileNode.write("\n")
                iPlotFileNode.write("       time          displ_x        displ_y         vel_x           vel_y         acc_x          acc_y        Reaction_x     Reaction_y    \n")
                iPlotFileNode.close()
                self.PlotFilesNodesList.append(iPlotFileNode)
                self.PlotFilesNodesIdList.append(Id)

        # open plots for elements selected
        if self.FEM_Solution.ProjectParameters["watch_elements_list"].size() != 0:

            NumNElements = self.FEM_Solution.ProjectParameters["watch_elements_list"].size()
            for elem in range(0, NumNElements):
                Id = self.FEM_Solution.ProjectParameters["watch_elements_list"][elem].GetInt()
                iPlotFileElem = open("PlotElement_" + str(Id) + ".txt","w")
                iPlotFileElem.write("\n")
                iPlotFileElem.write("          time                       Sxx                   Syy                      Sxy                    Exx                     Eyy                   Exy                Damage  \n")
                iPlotFileElem.close()
                self.PlotFilesElementsList.append(iPlotFileElem)
                self.PlotFilesElementsIdList.append(Id)

#============================================================================================================================
    def InitializeMMGvariables(self):

        ZeroVector3 = KratosMultiphysics.Vector(3)
        ZeroVector3[0] = 0.0
        ZeroVector3[1] = 0.0
        ZeroVector3[2] = 0.0

        for node in self.FEM_Solution.main_model_part.Nodes:
            node.SetValue(MeshingApplication.AUXILIAR_GRADIENT, ZeroVector3)

#============================================================================================================================
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

#============================================================================================================================
    def CheckIfHasRemeshed(self):

        is_remeshed = False

        if (self.RemeshingProcessMMG.initial_remeshing == False):
            step = self.RemeshingProcessMMG.step + 1
            # We need to check if the model part has been modified recently
            if self.RemeshingProcessMMG.step_frequency > 0:
                if step >= self.RemeshingProcessMMG.step_frequency:
                        if self.RemeshingProcessMMG.main_model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.RemeshingProcessMMG.initial_step:
                            # Has remeshed
                            is_remeshed = True
        return is_remeshed

#============================================================================================================================
    def PrintDEMResults(self):

        if self.DEM_Solution.step == 1: # always print the 1st step
            self.DEM_Solution.PrintResultsForGid(self.DEM_Solution.time)
            self.DEM_Solution.time_old_print = self.DEM_Solution.time

        else:
            time_to_print = self.DEM_Solution.time - self.DEM_Solution.time_old_print

            if (self.DEM_Solution.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self.DEM_Solution.solver.dt):
                self.DEM_Solution.PrintResultsForGid(self.DEM_Solution.time)
                self.DEM_Solution.time_old_print = self.DEM_Solution.time


#============================================================================================================================
    def CreateInitialSkinDEM(self):

        initial_dem_skin_process = KratosFemDem.InitialDemSkinProcess(self.FEM_Solution.main_model_part)
        initial_dem_skin_process.Execute()

        # Loop over the elements of the Submodelpart to create the DEM
        for node in self.FEM_Solution.main_model_part.GetSubModelPart("InitialDemSkin").Nodes:

            Id = node.Id
            R = node.GetValue(KratosFemDem.DEM_RADIUS)
            Coordinates = self.GetNodeCoordinates(node)
            self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates, R, Id)
            node.SetValue(KratosFemDem.IS_DEM, True)

#============================================================================================================================

    def InitializeIntegrationPointsVariables(self):

        for elem in self.FEM_Solution.main_model_part.Elements:
            elem.SetValue(KratosFemDem.STRESS_THRESHOLD, 0.0)
            elem.SetValue(KratosFemDem.DAMAGE_ELEMENT, 0.0)
            elem.SetValue(KratosFemDem.PRESSURE_EXPANDED, 0)
            elem.SetValue(KratosFemDem.IS_SKIN, 0)
            elem.SetValue(KratosFemDem.SMOOTHING, 0)
            elem.SetValue(KratosFemDem.STRESS_VECTOR, [0.0,0.0,0.0])
            elem.SetValue(KratosFemDem.STRAIN_VECTOR, [0.0,0.0,0.0])

#============================================================================================================================

    def GetMaximumConditionId(self):
        max_id = 0
        for condition in self.FEM_Solution.main_model_part.Conditions:
            if condition.Id > max_id:
                max_id = condition.Id
        return max_id

#============================================================================================================================

    def InitializeDummyNodalForces(self):
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

#============================================================================================================================
    def RemoveDummyNodalForces(self):

        for condition in self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").Conditions:
            condition.Set(KratosMultiphysics.TO_ERASE, True)

        self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.FEM_Solution.main_model_part.RemoveSubModelPart("ContactForcesDEMConditions")

#============================================================================================================================
    def RemoveAloneDEMElements(self):
        # method to remove the dem corresponding to inactive nodes
        FEM_Nodes = self.FEM_Solution.main_model_part.Nodes
        FEM_Elements = self.FEM_Solution.main_model_part.Elements

        for node in FEM_Nodes:
            node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, 0)

        for Element in FEM_Elements:
            for i in range(0, self.number_of_nodes_element): # Loop over nodes of the element
                if Element.IsNot(KratosMultiphysics.TO_ERASE):
                    node = Element.GetNodes()[i]
                    NumberOfActiveElements = node.GetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS)
                    NumberOfActiveElements += 1
                    node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, NumberOfActiveElements)

        NumberOfActiveElements = 0
        for node in FEM_Nodes:
            NumberOfActiveElements = node.GetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS)
            if NumberOfActiveElements == 0:
                self.SpheresModelPart.GetNode(node.Id).Set(KratosMultiphysics.TO_ERASE, True)

        self.SpheresModelPart.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)

#============================================================================================================================

    def InitializeSolutionAfterRemeshing(self):
        # Initialize the "flag" IS_DEM in all the nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosFemDem.IS_DEM, False, self.FEM_Solution.main_model_part.Nodes)
        # Initialize the "flag" NODAL_FORCE_APPLIED in all the nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosFemDem.NODAL_FORCE_APPLIED, False, self.FEM_Solution.main_model_part.Nodes)
        # Initialize the "flag" RADIUS in all the nodes
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.RADIUS, False, self.FEM_Solution.main_model_part.Nodes)

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
        skin_detection_process = KratosMultiphysics.SkinDetectionProcess2D(self.FEM_Solution.main_model_part,
                                                                            self.SkinDetectionProcessParameters)
        skin_detection_process.Execute()
        self.GenerateDemAfterRemeshing()

#============================================================================================================================

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


            
