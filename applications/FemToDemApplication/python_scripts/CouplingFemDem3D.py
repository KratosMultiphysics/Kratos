from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import FEMDEMParticleCreatorDestructor as PCD
import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import CouplingFemDem
import math
import KratosMultiphysics.MeshingApplication as MeshingApplication

def Wait():
	input("Press Something")

# Main script of the coupled FEM-DEM Application 3D
class FEMDEM3D_Solution(CouplingFemDem.FEMDEM_Solution):

#============================================================================================================================
	def Info(self):
		print("Coupling of the 3D FEMDEM App")

#============================================================================================================================
	def Initialize(self):
		self.number_of_nodes_element = 4
		self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.ERASED_VOLUME] = 0.0 #Sand Production Calculations
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

		self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part, 4, 5)

		if self.DoRemeshing:
			self.InitializeMMGvariables()
			self.RemeshingProcessMMG.ExecuteInitialize()

		if self.FEM_Solution.ProjectParameters.Has("pressure_load_extrapolation") == False:
			self.PressureLoad = False
		else:
			self.PressureLoad = self.FEM_Solution.ProjectParameters["pressure_load_extrapolation"].GetBool()
		if self.PressureLoad:
			KratosFemDem.AssignPressureIdProcess(self.FEM_Solution.main_model_part).Execute()

        # for the dem contact forces coupling
		self.InitializeDummyNodalForces()

		# Just to find neighbours the 1st time
		self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM] = True
		self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECOMPUTE_NEIGHBOURS] = True

		self.FEM_Solution.KratosPrintInfo(" /$$$$$$$$ /$$$$$$$$ /$$      /$$  /$$$$$$  /$$$$$$$  /$$$$$$$$ /$$      /$$")
		self.FEM_Solution.KratosPrintInfo("| $$_____/| $$_____/| $$$    /$$$ /$$__  $$| $$__  $$| $$_____/| $$$    /$$$")
		self.FEM_Solution.KratosPrintInfo("| $$      | $$      | $$$$  /$$$$|__/  \ $$| $$  \ $$| $$      | $$$$  /$$$$")
		self.FEM_Solution.KratosPrintInfo("| $$$$$   | $$$$$   | $$ $$/$$ $$  /$$$$$$/| $$  | $$| $$$$$   | $$ $$/$$ $$")
		self.FEM_Solution.KratosPrintInfo("| $$__/   | $$__/   | $$  $$$| $$ /$$____/ | $$  | $$| $$__/   | $$  $$$| $$")
		self.FEM_Solution.KratosPrintInfo("| $$      | $$      | $$\  $ | $$| $$      | $$  | $$| $$      | $$\  $ | $$")
		self.FEM_Solution.KratosPrintInfo("| $$      | $$$$$$$$| $$ \/  | $$| $$$$$$$$| $$$$$$$/| $$$$$$$$| $$ \/  | $$")
		self.FEM_Solution.KratosPrintInfo("|__/      |________/|__/     |__/|________/|_______/ |________/|__/     |__/ 3D Application")
		self.FEM_Solution.KratosPrintInfo("")

		if self.echo_level > 0:
			self.FEM_Solution.KratosPrintInfo("FEM-DEM Solution initialized")

		# We assign the flag to recompute neighbours inside the 3D elements the 1st time
		utils = KratosMultiphysics.VariableUtils()
		utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)

		if self.echo_level > 0:
			KratosMultiphysics.Logger.PrintInfo("FEM-DEM Solution initialized")

		# We assign the flag to recompute neighbours inside the 3D elements the 1st time
		utils = KratosMultiphysics.VariableUtils()
		utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)

#============================================================================================================================
	def InitializeSolutionStep(self):

        # modified for the remeshing
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
	def SolveSolutionStep(self):

		# Function to perform the coupling FEM <-> DEM
		self.FEM_Solution.clock_time = self.FEM_Solution.StartTimeMeasuring()

		#### SOLVE FEM #########################################
		self.FEM_Solution.solver.Solve()
		########################################################

		self.ExpandWetNodes()
		self.GenerateDEM() # we create the new DEM of this time step
		self.ExtrapolatePressure()

		self.SpheresModelPart = self.ParticleCreatorDestructor.GetSpheresModelPart()

		self.UpdateDEMVariables()     # We update coordinates, displ and velocities of the DEM according to FEM

		self.DEM_Solution.InitializeTimeStep()

		self.DEM_Solution.time = self.FEM_Solution.time
		self.DEM_Solution.step = self.FEM_Solution.step

		self.DEM_Solution.DEMFEMProcedures.UpdateTimeInModelParts(self.DEM_Solution.all_model_parts, self.DEM_Solution.time,self.DEM_Solution.solver.dt,self.DEM_Solution.step, self.DEM_Solution.IsTimeToPrintPostProcess())
		self.DEM_Solution._BeforeSolveOperations(self.DEM_Solution.time)

		#### SOLVE DEM #########################################
		self.DEM_Solution.solver.Solve()
		########################################################
		self.DEM_Solution.AfterSolveOperations()

		self.DEM_Solution.solver._MoveAllMeshes(self.DEM_Solution.time, self.DEM_Solution.solver.dt)
		self.UpdateDEMVariables() # to print DEM with the FEM coordinates

		self.PrintDEMResultsForGid()

		self.DEM_Solution.FinalizeTimeStep(self.DEM_Solution.time)

		# Transfer the contact forces of the DEM to the FEM nodes
		self.TransferNodalForcesToFEM()

		self.FEM_Solution.StopTimeMeasuring(self.FEM_Solution.clock_time,"Solving", False)

		# Update Coupled Postprocess file for Gid (post.lst)
		self.WritePostListFile()

		# Print required info
		self.PrintPlotsFiles()

#============================================================================================================================
	def GenerateDEM(self): # 3D version
		if self.echo_level > 0:
			self.FEM_Solution.KratosPrintInfo("FEM-DEM:: GenerateDEM")
		self.CountErasedVolume()

		if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM]:
			dem_generator_process = KratosFemDem.GenerateDemProcess(self.FEM_Solution.main_model_part, self.SpheresModelPart)
			dem_generator_process.Execute()

			self.RemoveAloneDEMElements()
			element_eliminator = KratosMultiphysics.AuxiliarModelPartUtilities(self.FEM_Solution.main_model_part)
			element_eliminator.RemoveElementsAndBelongings(KratosMultiphysics.TO_ERASE)

			# We assign the flag to recompute neighbours inside the 3D elements
			utils = KratosMultiphysics.VariableUtils()
			utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)

#============================================================================================================================
	def CheckInactiveNodes(self):

		FEM_Elements = self.FEM_Solution.main_model_part.Elements
		FEM_Nodes    = self.FEM_Solution.main_model_part.Nodes
		erased_nodes_id = []
		conditions_to_erase_id = []

		for node in FEM_Nodes:
			node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, 0)

		for Element in FEM_Elements:
			if Element.IsNot(KratosMultiphysics.TO_ERASE):
				for i in range(0, 4): # Loop over nodes of the element
					node = Element.GetNodes()[i]
					NumberOfActiveElements = node.GetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS)
					NumberOfActiveElements += 1
					node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, NumberOfActiveElements)

		NumberOfActiveElements = 0
		for node in FEM_Nodes:
			NumberOfActiveElements = node.GetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS)
			if NumberOfActiveElements == 0 and node.GetValue(KratosFemDem.INACTIVE_NODE) == False:
				Id = node.Id
				# print("nodo eliminado: ", Id)
				DEMnode = self.SpheresModelPart.GetNode(Id)
				node.SetValue(KratosFemDem.INACTIVE_NODE, True)
				node.Set(KratosMultiphysics.TO_ERASE, True) # added
				DEMnode.SetValue(KratosFemDem.INACTIVE_NODE, True)
				DEMnode.Set(KratosMultiphysics.TO_ERASE, True)
				erased_nodes_id.append(Id)

				for condition in self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").Conditions:
					if condition.GetNodes()[0].Id == Id:
						conditions_to_erase_id.append(condition.Id)

			# Reset the value to the next step
			node.SetValue(KratosFemDem.NUMBER_OF_ACTIVE_ELEMENTS, 0)

        # let's remove the nodal dem conditions according to inactive nodes
		for Id in conditions_to_erase_id:
			self.FEM_Solution.main_model_part.RemoveCondition(Id)

		# Remove inactive nodes
		self.SpheresModelPart.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)
		self.FEM_Solution.main_model_part.GetRootModelPart().RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE) # added

#============================================================================================================================
	def UpdateDEMVariables(self):
		update_de_kinematics_process = KratosFemDem.UpdateDemKinematicsProcess(self.FEM_Solution.main_model_part, self.SpheresModelPart)
		update_de_kinematics_process.Execute()
#============================================================================================================================
	def PrintPlotsFiles(self):

		# Print the general file
		time = self.FEM_Solution.time
		total_reaction_x     = 0.0
		total_displacement_x = 0.0
		total_reaction_y     = 0.0
		total_displacement_y = 0.0
		total_reaction_z     = 0.0
		total_displacement_z = 0.0
		interval = self.FEM_Solution.ProjectParameters["interval_of_watching"].GetDouble()

		if self.FEM_Solution.time - self.TimePreviousPlotting >= interval:
			if self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size() > 0:
				if self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][0].IsInt():
					for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size()):
						IdNode = self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][index].GetInt()
						node = self.FEM_Solution.main_model_part.GetNode(IdNode)
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
						IdNode = self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][index].GetInt()
						node = self.FEM_Solution.main_model_part.GetNode(IdNode)
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
				self.PlotFile.write("    " + "{0:.4e}".format(time).rjust(11) + "    " + "{0:.4e}".format(total_displacement_x).rjust(11) +
					"    " + "{0:.4e}".format(total_displacement_y).rjust(11) + "    " + "{0:.4e}".format(total_displacement_z).rjust(11) +
					"    " + "{0:.4e}".format(total_reaction_x).rjust(11) + "    " + "{0:.4e}".format(total_reaction_y).rjust(11) + "    " +
					"{0:.4e}".format(total_reaction_z).rjust(11) + "\n")
				self.PlotFile.close()

			# Print the selected nodes files
			if self.FEM_Solution.ProjectParameters["watch_nodes_list"].size() != 0:
				NumNodes = self.FEM_Solution.ProjectParameters["watch_nodes_list"].size()
				for inode in range(0, NumNodes):
					IdNode = self.PlotFilesNodesIdList[inode]
					node = self.FEM_Solution.main_model_part.GetNode(IdNode)
					self.PlotFilesNodesList[inode] = open("PlotNode_" + str(IdNode) + ".txt", "a")

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
					strain_tensor = Elem.GetValue(KratosFemDem.STRAIN_VECTOR)

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

					damage = Elem.GetValue(KratosFemDem.DAMAGE_ELEMENT)

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

#============================================================================================================================
	def InitializePlotsFiles(self):

		# open general Displ/Reaction File
		self.PlotFile = open("PlotFile.txt","w")
		self.PlotFile.write("This File Plots the SUM of the displacement and reactions of the nodes selected in the lists!\n\n")
		self.PlotFile.write("       time          displ_x        displ_y        displ_z       Reaction_x     Reaction_y     Reaction_z    \n")
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
				iPlotFileNode.write("       time          displ_x        displ_y        displ_z         vel_x           vel_y         vel_z           acc_x          acc_y          acc_z       Reaction_x     Reaction_y     Reaction_Z    \n")
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
				iPlotFileElem.write("       time             Sxx           Syy             Szz           Sxy            Syz            Sxz            Exx            Eyy            Ezz             Exy           Eyz            Exz          Damage  \n")
				iPlotFileElem.close()
				self.PlotFilesElementsList.append(iPlotFileElem)
				self.PlotFilesElementsIdList.append(Id)

#============================================================================================================================
	def RefineMappedVariables(self):
		for elem in self.FEM_Solution.main_model_part.Elements:
			if elem.GetValue(KratosFemDem.DAMAGE_ELEMENT) < 0.0:
				elem.SetValue(KratosFemDem.DAMAGE_ELEMENT, 0.0)

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
		skin_detection_process_param = KratosMultiphysics.Parameters("""
		{
			"name_auxiliar_model_part" : "SkinDEMModelPart",
			"name_auxiliar_condition"  : "Condition",
			"echo_level"               : 0
		}""")
		skin_detection_process = KratosMultiphysics.SkinDetectionProcess3D(self.FEM_Solution.main_model_part,
																		skin_detection_process_param)
		skin_detection_process.Execute()
		self.GenerateDemAfterRemeshing()

#============================================================================================================================

	def PrintDEMResultsForGid(self):

		# DEM GiD print output
		if self.DEM_Solution.step == 1: # always print the 1st step
			self.DEM_Solution.PrintResultsForGid(self.DEM_Solution.time)
			self.DEM_Solution.time_old_print = self.DEM_Solution.time
		else:
			time_to_print = self.DEM_Solution.time - self.DEM_Solution.time_old_print

			if (self.DEM_Solution.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self.DEM_Solution.solver.dt):

				self.DEM_Solution.PrintResultsForGid(self.DEM_Solution.time)
				self.DEM_Solution.time_old_print = self.DEM_Solution.time

#============================================================================================================================

	def InitializeIntegrationPointsVariables(self):

		utils = KratosMultiphysics.VariableUtils()
		utils.SetNonHistoricalVariable(KratosFemDem.VOLUME_COUNTED, False, self.FEM_Solution.main_model_part.Elements)
		utils.SetNonHistoricalVariable(KratosFemDem.STRESS_THRESHOLD, 0.0, self.FEM_Solution.main_model_part.Elements)
		utils.SetNonHistoricalVariable(KratosFemDem.DAMAGE_ELEMENT, 0.0, self.FEM_Solution.main_model_part.Elements)
		utils.SetNonHistoricalVariable(KratosFemDem.PRESSURE_EXPANDED, 0, self.FEM_Solution.main_model_part.Elements)
		utils.SetNonHistoricalVariable(KratosFemDem.IS_SKIN, 0, self.FEM_Solution.main_model_part.Elements)
		utils.SetNonHistoricalVariable(KratosFemDem.SMOOTHING, 0, self.FEM_Solution.main_model_part.Elements)
		utils.SetNonHistoricalVariable(KratosFemDem.STRESS_VECTOR, [0.0,0.0,0.0,0.0,0.0,0.0], self.FEM_Solution.main_model_part.Elements)
		utils.SetNonHistoricalVariable(KratosFemDem.STRAIN_VECTOR, [0.0,0.0,0.0,0.0,0.0,0.0], self.FEM_Solution.main_model_part.Elements)
		utils.SetNonHistoricalVariable(KratosFemDem.STRESS_VECTOR_INTEGRATED, [0.0,0.0,0.0,0.0,0.0,0.0], self.FEM_Solution.main_model_part.Elements)

#===================================================================================================================================

	def ExpandWetNodes(self):
		if self.PressureLoad:
			# This must be called before Generating DEM
			self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECONSTRUCT_PRESSURE_LOAD] = 0 # It is modified inside
			extend_wet_nodes_process = KratosFemDem.ExpandWetNodesProcess(self.FEM_Solution.main_model_part)
			extend_wet_nodes_process.Execute()

#===================================================================================================================================

	def ExtrapolatePressure(self):
		if self.echo_level > 0:
			self.FEM_Solution.KratosPrintInfo("FEM-DEM:: ExtrapolatePressureLoad")
		if self.PressureLoad:
			# we reconstruct the pressure load if necessary
			if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECONSTRUCT_PRESSURE_LOAD] == 1:
				self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.INTERNAL_PRESSURE_ITERATION] = 1
				while self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.INTERNAL_PRESSURE_ITERATION] > 0:
					KratosFemDem.ExtendPressureConditionProcess3D(self.FEM_Solution.main_model_part).Execute()

#===================================================================================================================================

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

				# Extrapolate the VonMises normalized stress to nodes (remeshing)
				parameters = self.FEM_Solution.ProjectParameters["AMR_data"]["hessian_variable_parameters"]
				KratosFemDem.ComputeNormalizedFreeEnergyOnNodesProcess(self.FEM_Solution.main_model_part, parameters).Execute()

				# we eliminate the nodal DEM forces
				self.RemoveDummyNodalForces()

			# Perform remeshing
			self.RemeshingProcessMMG.ExecuteInitializeSolutionStep()

			if is_remeshing:
				if self.echo_level > 0:
					self.FEM_Solution.KratosPrintInfo("FEM-DEM:: InitializeSolutionAfterRemeshing")

				self.RefineMappedVariables()
				self.InitializeSolutionAfterRemeshing()
				self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part, 4, 5)
				self.nodal_neighbour_finder.Execute()
				# We assign the flag to recompute neighbours inside the 3D elements
				utils = KratosMultiphysics.VariableUtils()
				utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)

#===================================================================================================================================

	def FindNeighboursIfNecessary(self):
		if self.echo_level > 0:
			self.FEM_Solution.KratosPrintInfo("FEM-DEM:: ComputeNeighboursIfNecessary")

		if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM]: # The neighbours have changed
			self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part, 4, 5)
			self.nodal_neighbour_finder.Execute()
			# We reset the flag
			self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM] = False

	def CountErasedVolume(self):
		count_erased_vol = True
		if count_erased_vol:
			erased_vol_process = KratosFemDem.ComputeSandProduction(self.FEM_Solution.main_model_part)
			erased_vol_process.Execute()

			self.ErasedVolume = open("ErasedVolume.txt","a")
			erased_vol = self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.ERASED_VOLUME]
			self.ErasedVolume.write("    " + "{0:.4e}".format(self.FEM_Solution.time).rjust(11) + "    " + "{0:.4e}".format(erased_vol).rjust(11) + "\n")
			self.ErasedVolume.close()
		
