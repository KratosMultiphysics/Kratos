from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
CheckForPreviousImport()

# ======================================================================================
# ======================================================================================
# ======================================================================================

class ResultsIO:
	
	# ==================================================================================
	
	def __init__(self, 
				 ModelPart, 
				 OutputFileName,
				 ResultsOnNodes = [], 
				 ResultsOnGaussPoints = []):
		
		self.ModelPart = ModelPart
		self.OutputFileName = OutputFileName
		self.MeshName = 0.0
		
		self.ResultsOnNodes = ResultsOnNodes
		self.ResultsOnGaussPoints = ResultsOnGaussPoints
		
		self.myIO = GidIO(
			OutputFileName,
			GiDPostMode.GiD_PostBinary,
			MultiFileFlag.SingleFile,
			WriteDeformedMeshFlag.WriteUndeformed, 
			WriteConditionsFlag.WriteElementsOnly)
		
		self.IsInitialized = False
		
	# ==================================================================================
	
	def Initialize(self):
		
		if(self.IsInitialized == False):
			
			self.myIO.InitializeMesh(self.MeshName);
			self.myIO.WriteMesh(self.ModelPart.GetMesh());
			self.myIO.FinalizeMesh()
			self.myIO.InitializeResults(self.MeshName, self.ModelPart.GetMesh())
			
			self.IsInitialized = True
			
	# ==================================================================================
	
	def Finalize(self):
		
		if(self.IsInitialized):
			
			self.myIO.FinalizeResults()
			self.IsInitialized = False
			
	# ==================================================================================
	
	def Write(self, theCurrentTime):
		
		if(self.IsInitialized):
			
			for result in self.ResultsOnNodes:
				self.myIO.WriteNodalResults(result, 
											self.ModelPart.Nodes, 
											theCurrentTime, 0)
				
			for result in self.ResultsOnGaussPoints:
				self.myIO.PrintOnGaussPoints(result, 
											 self.ModelPart, 
											 theCurrentTime)