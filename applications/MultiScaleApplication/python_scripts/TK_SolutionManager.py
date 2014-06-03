from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# ======================================================================================
# ======================================================================================
# ======================================================================================

class Study:
	
	# ==================================================================================
	
	def __init__(
				self, 
				InitialTime):
		
		self.InitialTime = InitialTime
		self.SolutionStages = []
	
	# ==================================================================================
	
	def Append(self, SolutionStage):
		
		self.SolutionStages.append(SolutionStage)
	
	# ==================================================================================
	
	def Run(self):
		
		# set intial times for each solution stage
		initial_time = self.InitialTime
		for istage in self.SolutionStages:
			istage.SetInitialTime(initial_time)
			initial_time = istage.GetEndTime()
		
		# initialize ResultsIO
		for istage in self.SolutionStages:
			if(istage.ResultsIO != None):
				istage.ResultsIO.Initialize()
		
		# initialize stages
		for istage in self.SolutionStages:
			istage.Initialize()
		
		# run stages
		for istage in self.SolutionStages:
			istage.Solve()
			if(istage.IsConverged == False):
				break;
		
		# finalize stages
		for istage in self.SolutionStages:
			istage.Finalize()
		
		# finalize ResultsIO
		for istage in self.SolutionStages:
			if(istage.ResultsIO != None):
				istage.ResultsIO.Finalize()
		