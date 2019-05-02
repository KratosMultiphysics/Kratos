from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import CouplingFemDem3D
import KratosMultiphysics


def Wait():
	input("Press Something")

# Main script of the coupled FEM-DEM Application 3D for hexahedrons
class FEMDEM3DHexahedrons_Solution(CouplingFemDem3D.FEMDEM3D_Solution):

#============================================================================================================================
	def GenerateDEM(self): # 3D version for hexahedrons
		pass
#============================================================================================================================
	def CheckForPossibleIndentations(self): # Verifies if an element has indentations between its DEM
		pass
#============================================================================================================================

	def CheckInactiveNodes(self):
		pass
#============================================================================================================================

	def InitializeSolutionStep(self):

        # modified for the remeshing
		self.FEM_Solution.delta_time = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
		self.FEM_Solution.time = self.FEM_Solution.time + self.FEM_Solution.delta_time
		self.FEM_Solution.main_model_part.CloneTimeStep(self.FEM_Solution.time)
		self.FEM_Solution.step = self.FEM_Solution.step + 1
		self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.FEM_Solution.step

		self.FEM_Solution.InitializeSolutionStep()

		# Create initial skin of DEM's
		self.create_initial_dem_skin = False  # Hard Coded TODO
		if self.create_initial_dem_skin and self.FEM_Solution.step == 1:
			self.CreateInitialSkinDEM()


#============================================================================================================================
