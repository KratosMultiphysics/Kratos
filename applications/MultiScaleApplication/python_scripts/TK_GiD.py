from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *

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

		self.myIO = GidIO(
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

class ResultsIO_2Physics:

	# ==================================================================================

	def __init__(self,
				 ModelPartA,
				 ModelPartB,
				 OutputFileName,
				 ResultsOnNodes = [],
				 ResultsOnGaussPoints_ModA = [],
				 ResultsOnGaussPoints_ModB = [],
				 Frequency = None):

		self.ModelPartA = ModelPartA
		self.ModelPartB = ModelPartB
		self.OutputFileName = OutputFileName
		self.MeshName = 0.0

		self.ResultsOnNodes = ResultsOnNodes
		self.ResultsOnGaussPoints_ModA = ResultsOnGaussPoints_ModA
		self.ResultsOnGaussPoints_ModB = ResultsOnGaussPoints_ModB

		self.myIO = GidIO(
			OutputFileName,
			GiDPostMode.GiD_PostBinary, #GiD_PostAscii, #
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
			self.myIO.WriteMesh(self.ModelPartA.GetMesh());
			self.myIO.FinalizeMesh()
			self.myIO.InitializeResults(self.MeshName, self.ModelPartA.GetMesh())
			self.myIO.InitializeResults(self.MeshName, self.ModelPartB.GetMesh())

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

	# ==================================================================================

	def Write(self, theCurrentTime):

		if(self.IsInitialized and self.__check_write_freq(theCurrentTime)):

			for result in self.ResultsOnNodes:
				self.myIO.WriteNodalResults(result,
											self.ModelPartA.Nodes,
											theCurrentTime, 0)
			# print('Print on Nodes Done')
			# print('self.ModelPartA',self.ModelPartA)
			# print('self.ModelPartB',self.ModelPartB)

			for result in self.ResultsOnGaussPoints_ModA:
				self.myIO.PrintOnGaussPoints(result,
											 self.ModelPartA,
											 theCurrentTime)
			# print('Print on OnGaussPoints_ModA Done')

			for result in self.ResultsOnGaussPoints_ModB:
				self.myIO.PrintOnGaussPoints(result,
											 self.ModelPartB,
											 theCurrentTime)
			# print('Print on OnGaussPoints_ModB Done')

			self.myIO.Flush();