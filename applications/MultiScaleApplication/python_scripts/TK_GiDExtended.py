from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.MultiScaleApplication import *
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
				 ResultsOnGaussPoints = [],
				 Frequency = None):
		
		self.ModelPart = ModelPart
		self.OutputFileName = OutputFileName
		self.MeshName = 0.0
		
		self.ResultsOnNodes = ResultsOnNodes
		self.ResultsOnGaussPoints = ResultsOnGaussPoints
		
		self.myIO = GidIO_Extended(
			OutputFileName,
			GiDPostMode.GiD_PostBinary,
			MultiFileFlag.SingleFile,
			WriteDeformedMeshFlag.WriteUndeformed, 
			WriteConditionsFlag.WriteElementsOnly)
		
		self.IsInitialized = False
		
		self.Frequency = Frequency
		self.Tn = None
		
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
	
	def __check_write_freq(self,t):
		r = True
		f = self.Frequency
		if(f is not None):
			if(self.Tn is None):
				self.Tn = t
			else:
				dt = t-self.Tn
				if(dt >= f):
					self.Tn = t
				else:
					r = False
		return r
	
	def Write(self, theCurrentTime):
		
		if(self.IsInitialized and self.__check_write_freq(theCurrentTime)):
			
			for result in self.ResultsOnNodes:
				self.myIO.WriteNodalResults(result, 
											self.ModelPart.Nodes, 
											theCurrentTime, 0)
			
			for result in self.ResultsOnGaussPoints:
				self.myIO.PrintOnGaussPoints(result, 
											 self.ModelPart, 
											 theCurrentTime)
			
			self.myIO.Flush();