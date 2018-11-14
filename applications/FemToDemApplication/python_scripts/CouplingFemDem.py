from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import MainDEM_for_coupling as DEM
import MainFEM_for_coupling as FEM
import FEMDEMParticleCreatorDestructor as PCD
import KratosMultiphysics
import KratosMultiphysics.DEMApplication as DEMApplication
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.MeshingApplication as MeshingApplication
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
            Model = {self.mmg_parameters["model_part_name"].GetString(): self.FEM_Solution.main_model_part}
            self.RemeshingProcessMMG = MMG.MmgProcess(Model, self.mmg_parameters)

        self.InitializePlotsFiles()


#============================================================================================================================
    def Run(self):
        
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

#============================================================================================================================
    def Initialize(self):

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

        self.pressure_load = True #hard coded
        if self.pressure_load:
            KratosFemDem.AssignPressureIdProcess(self.FEM_Solution.main_model_part).Execute()



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
        self.FEM_Solution.delta_time = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self.FEM_Solution.time = self.FEM_Solution.time + self.FEM_Solution.delta_time
        self.FEM_Solution.main_model_part.CloneTimeStep(self.FEM_Solution.time)
        self.FEM_Solution.step = self.FEM_Solution.step + 1
        self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.FEM_Solution.step

        if self.DoRemeshing:
            is_remeshing = self.CheckIfHasRemeshed()
            
            if is_remeshing:
                # Extrapolate the VonMises normalized stress to nodes (remeshing)
                KratosFemDem.StressToNodesProcess(self.FEM_Solution.main_model_part, 2).Execute()
            # Perform remeshing
            self.RemeshingProcessMMG.ExecuteInitializeSolutionStep()

            if is_remeshing:
                # Initialize the "flag" IS_DEM in all the nodes
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosFemDem.IS_DEM, False, self.FEM_Solution.main_model_part.Nodes)
                # Initialize the "flag" NODAL_FORCE_APPLIED in all the nodes
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosFemDem.NODAL_FORCE_APPLIED, False, self.FEM_Solution.main_model_part.Nodes)
                # Initialize the "flag" RADIUS in all the nodes
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.RADIUS, False, self.FEM_Solution.main_model_part.Nodes)

                self.pressure_load = True #hard coded
                if self.pressure_load:
                    KratosFemDem.AssignPressureIdProcess(self.FEM_Solution.main_model_part).Execute()

                # Remove DEMS from previous mesh
                self.SpheresModelPart.Elements.clear()
                self.SpheresModelPart.Nodes.clear()

                self.InitializeMMGvariables()
                self.FEM_Solution.model_processes = self.FEM_Solution.AddProcesses()
                self.FEM_Solution.model_processes.ExecuteInitialize()
                self.FEM_Solution.model_processes.ExecuteBeforeSolutionLoop()
                self.FEM_Solution.model_processes.ExecuteInitializeSolutionStep()

		# Search the skin nodes for the remeshing
        if self.DoRemeshing:
            skin_detection_process_param = KratosMultiphysics.Parameters("""
            {
                "name_auxiliar_model_part" : "SkinDEMModelPart",
                "name_auxiliar_condition"  : "Condition",
                "echo_level"               : 0
            }""")
            skin_detection_process = KratosMultiphysics.SkinDetectionProcess2D(self.FEM_Solution.main_model_part,
                                                                                skin_detection_process_param)
            skin_detection_process.Execute()

        self.FEM_Solution.InitializeSolutionStep()

        # Create initial skin of DEM's
        self.create_initial_dem_skin = False  # Hard Coded TODO
        if self.create_initial_dem_skin and self.FEM_Solution.step == 1:
            self.CreateInitialSkinDEM()

        # Create the DEM after the remeshing
        if self.DoRemeshing and is_remeshing:
            self.GenerateDemAfterRemeshing()



#============================================================================================================================
    def SolveSolutionStep(self): # Function to perform the coupling FEM <-> DEM

        self.FEM_Solution.clock_time = self.FEM_Solution.StartTimeMeasuring()

        #### SOLVE FEM #########################################
        self.FEM_Solution.solver.Solve()
        ########################################################

        if self.pressure_load:
            # we reconstruct the pressure load
            self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.ITER] = 1
            itera = 0
            while self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.ITER] > 0: 
                print("--------------------iteracion while--------------------------", itera)
                KratosFemDem.ExtendPressureConditionProcess2D(self.FEM_Solution.main_model_part,).Execute()
                itera += 1

        # we create the new DEM of this time step
        self.GenerateDEM()
        self.SpheresModelPart = self.ParticleCreatorDestructor.GetSpheresModelPart()
        self.CheckForPossibleIndentations()
        self.CheckInactiveNodes()

        # We update coordinates, displ and velocities of the DEM according to FEM
        self.UpdateDEMVariables()

        # Extrapolate the VonMises normalized stress to nodes (remeshing)
        KratosFemDem.StressToNodesProcess(self.FEM_Solution.main_model_part, 2).Execute()

        self.DEM_Solution.InitializeTimeStep()
        self.DEM_Solution.time = self.FEM_Solution.time
        self.DEM_Solution.step = self.FEM_Solution.step
        self.DEM_Solution.DEMFEMProcedures.UpdateTimeInModelParts(self.DEM_Solution.all_model_parts,
                                                                   self.DEM_Solution.time,
                                                                   self.DEM_Solution.dt,
                                                                   self.DEM_Solution.step,
                                                                   self.DEM_Solution.IsTimeToPrintPostProcess(self.DEM_Solution.time))
        self.DEM_Solution.BeforeSolveOperations(self.DEM_Solution.time)

        #### SOLVE DEM #########################################
        self.DEM_Solution.solver.Solve()
        ########################################################

        self.DEM_Solution.AfterSolveOperations()
        self.DEM_Solution.DEMFEMProcedures.MoveAllMeshes(self.DEM_Solution.all_model_parts, self.DEM_Solution.time, self.DEM_Solution.dt)
        
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
    def GenerateDEM(self):

        FEM_elements = self.FEM_Solution.main_model_part.Elements

        # Loop Over Elements to find the INACTIVE ones and generate the DEM only once
        for Element in FEM_elements:

            is_active     = True
            DEM_Generated = Element.GetValue(KratosFemDem.DEM_GENERATED)

            if Element.IsDefined(KratosMultiphysics.ACTIVE):
                is_active = Element.Is(KratosMultiphysics.ACTIVE)

                NumberOfDEM = 0         # Number of nodes with DEM Associated
                for node in range(0,3): # Loop over nodes of the FE
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

        self.FEM_Solution.main_model_part.GetRootModelPart().RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)

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
                # Element.Set(KratosMultiphysics.TO_ERASE, True)
        
        # self.FEM_Solution.main_model_part.GetRootModelPart().RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)


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
        DEM_Nodes = self.SpheresModelPart.Nodes

        for DEM_Node in DEM_Nodes:  # Loop over DEM nodes
            if (DEM_Node.GetValue(KratosFemDem.INACTIVE_NODE) == False):

                Id = DEM_Node.Id
                Corresponding_FEM_Node = self.FEM_Solution.main_model_part.GetNode(Id)
                Coordinates    = self.GetNodeCoordinates(Corresponding_FEM_Node)
                Velocity_x     = Corresponding_FEM_Node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                Velocity_y     = Corresponding_FEM_Node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                Displacement_x = Corresponding_FEM_Node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                Displacement_y = Corresponding_FEM_Node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)

                # Update Coordinates
                DEM_Node.X = Coordinates[0]
                DEM_Node.Y = Coordinates[1]

                # Update Displacements
                DEM_Node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, Displacement_x)
                DEM_Node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, Displacement_y)

                # Update Velocities
                DEM_Node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, Velocity_x)
                DEM_Node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, Velocity_y)


#============================================================================================================================
    def CheckInactiveNodes(self):

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

        # Remove inactive nodes
        self.SpheresModelPart.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.FEM_Solution.main_model_part.GetRootModelPart().RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE) # added

#============================================================================================================================
    def TransferNodalForcesToFEM(self):

        DEM_Nodes = self.SpheresModelPart.Nodes

        for DEM_node in DEM_Nodes:

            # Get the contact forces from the DEM
            Id = DEM_node.Id
            Force_X = DEM_node.GetSolutionStepValue(KratosMultiphysics.TOTAL_FORCES_X)
            Force_Y = DEM_node.GetSolutionStepValue(KratosMultiphysics.TOTAL_FORCES_Y)

            # Tranfer the Force information to the FEM nodes
            self.FEM_Solution.main_model_part.GetNode(Id).SetValue(KratosFemDem.NODAL_FORCE_X, Force_X)
            self.FEM_Solution.main_model_part.GetNode(Id).SetValue(KratosFemDem.NODAL_FORCE_Y, Force_Y)


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
        TotalReaction_x     = 0.0
        TotalDisplacement_x = 0.0
        TotalReaction_y     = 0.0
        TotalDisplacement_y = 0.0
        interval = self.FEM_Solution.ProjectParameters["interval_of_watching"].GetDouble()


        if self.FEM_Solution.time - self.TimePreviousPlotting >= interval:

            for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size()):

                IdNode = self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][index].GetInt()
                node = self.FEM_Solution.main_model_part.GetNode(IdNode)
                TotalDisplacement_x += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                TotalDisplacement_y += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)

            for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"].size()):

                IdNode = self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][index].GetInt()
                node = self.FEM_Solution.main_model_part.GetNode(IdNode)
                TotalReaction_x += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
                TotalReaction_y += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y)

            self.PlotFile = open("PlotFile.txt","a")
            self.PlotFile.write("    " + "{0:.4e}".format(time).rjust(11) + "    " + "{0:.4e}".format(TotalDisplacement_x).rjust(11) +
                                "    " + "{0:.4e}".format(TotalDisplacement_y).rjust(11) + "    " + "{0:.4e}".format(TotalReaction_x).rjust(11) +
                                "    " + "{0:.4e}".format(TotalReaction_y).rjust(11) + "\n")
            self.PlotFile.close()

            # Print the selected nodes files
            if self.FEM_Solution.ProjectParameters["watch_nodes_list"].size() != 0:

                NumNodes = self.FEM_Solution.ProjectParameters["watch_nodes_list"].size()

                for inode in range(0, NumNodes):

                    IdNode = self.PlotFilesNodesIdList[inode]
                    node = self.FEM_Solution.main_model_part.GetNode(IdNode)

                    self.PlotFilesNodesList[inode] = open("PlotNode_" + str(IdNode) + ".txt","a")

                    dx = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                    dy = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                    Rx = node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
                    Ry = node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y)
                    vx = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                    vy = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    ax = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION_X)
                    ay = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION_Y)

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

                    Sxx = Elem.GetValuesOnIntegrationPoints(KratosFemDem.STRESS_VECTOR_INTEGRATED, self.FEM_Solution.main_model_part.ProcessInfo)[0][0]
                    Syy = Elem.GetValuesOnIntegrationPoints(KratosFemDem.STRESS_VECTOR_INTEGRATED, self.FEM_Solution.main_model_part.ProcessInfo)[0][1]
                    Sxy = Elem.GetValuesOnIntegrationPoints(KratosFemDem.STRESS_VECTOR_INTEGRATED, self.FEM_Solution.main_model_part.ProcessInfo)[0][2]

                    Exx = Elem.GetValue(KratosFemDem.STRAIN_VECTOR)[0]
                    Eyy = Elem.GetValue(KratosFemDem.STRAIN_VECTOR)[1]
                    Exy = Elem.GetValue(KratosFemDem.STRAIN_VECTOR)[2]

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
        self.PlotFile.write("       time          displ_x        displ_y       Reaction_x     Reaction_y    \n")
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
                iPlotFileElem.write("       time             Sxx           Syy             Sxy           Exx             Eyy          Exy           Damage  \n")
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
                        if self.RemeshingProcessMMG.model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.RemeshingProcessMMG.initial_step:
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

            if (self.DEM_Solution.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self.DEM_Solution.dt):
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
            elem.SetValue(KratosFemDem.STRESS_VECTOR, [0.0,0.0,0.0])
            elem.SetValue(KratosFemDem.STRAIN_VECTOR, [0.0,0.0,0.0])