from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import datetime
import time
import TK_TimeLines
from KratosMultiphysics import *
from KratosMultiphysics.MultiScaleApplication import *
CheckForPreviousImport()

# ======================================================================================
#
# Returns a list of all Variables required to use the SolutionStage object in
# this module.
#
# ======================================================================================

def Variables():
	return([
			DISPLACEMENT, 
			ROTATION,
			DISPLACEMENT_LAGRANGE,
			ROTATION_LAGRANGE,
			REACTION,
			MOMENTUM,
			FORCE,
			MOMENT,
			])

# ======================================================================================
#
# Returns a list of all (Physical) D.O.F.s required to use the SolutionStage object in
# this module.
#
# ======================================================================================

def Dofs():
	return ([
			DISPLACEMENT_LAGRANGE_X, 
			DISPLACEMENT_LAGRANGE_Y, 
			DISPLACEMENT_LAGRANGE_Z,
			ROTATION_LAGRANGE_X, 
			ROTATION_LAGRANGE_Y, 
			ROTATION_LAGRANGE_Z,
			])

# ======================================================================================
#
# Returns a list of tuples.
# Each tuple in the list is a Pair with the Lagrangian D.O.F. and the corresponding
# Variable used as reaction.
#
# ======================================================================================

def DofsWithReactions():
	return ([
			(DISPLACEMENT_X, REACTION_X),
			(DISPLACEMENT_Y, REACTION_Y),
			(DISPLACEMENT_Z, REACTION_Z),
			(ROTATION_X, MOMENTUM_X),
			(ROTATION_Y, MOMENTUM_Y),
			(ROTATION_Z, MOMENTUM_Z),
			])

# ======================================================================================
#
# The SolutionStage object.
# This class solves a stage of the analysis history using one or more time steps, 
# and it takes care of saving the results if required.
#
# ======================================================================================

class SolutionStage:

	# ==================================================================================
	
	def __init__(
				self, 
				ModelPart,
				TimeLine = TK_TimeLines.FixedTimeLine(
					Duration = 1.0,
					Increment = 1.0
					),
				Convergence = ResidualNormCriteria(1.0E-4, 1.0E-4),
				MaxIterations = 30,
				CalculateReactions = True,
				ReformDofSetAtEachStep = False,
				MoveMesh = True,
				ResultsIO = None,
				LinearSolver = SuperLUSolver(),
				CalculateTangent = True,
				Parallel = True,
				CustomOp = None):
		
		# Model Part
		self.ModelPart = ModelPart
		
		# Time line
		self.TimeLine = TimeLine
		
		# Time Scheme
		self.TimeScheme = StaticGeneralScheme()
		self.TimeScheme.Check(self.ModelPart)
		
		# Linear Solver
		self.LinearSolver = LinearSolver
		
		# Convergence Criteria
		self.Convergence = Convergence
		self.Convergence.Check(self.ModelPart)
		
		# Builder and Solver
		self.Parallel = Parallel
		if(self.Parallel):
			self.BuilderAndSolver = StaticGeneralBuilderAndSolver(self.LinearSolver)
		else:
			self.BuilderAndSolver = StaticGeneralBuilderAndSolverSequential(self.LinearSolver)
		
		# Misc. Solver parameters
		self.MaxIterations = MaxIterations
		self.CalculateReactions = CalculateReactions
		self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
		self.MoveMesh = MoveMesh
		
		# Results IO
		self.ResultsIO = ResultsIO
		
		# Create and initialize the solver
		self.Solver = StaticGeneralStrategy(
			self.ModelPart,
			self.TimeScheme,
			self.LinearSolver,
			self.Convergence,
			self.BuilderAndSolver,
			self.MaxIterations,
			self.CalculateReactions,
			self.ReformDofSetAtEachStep,
			self.MoveMesh)
		
		self.Solver.SetKeepSystemConstantDuringIterations(not CalculateTangent)
		
		self.Solver.Check();
		self.Solver.SetEchoLevel(1);
		
		self.IsConverged = False
		
		self.CustomOp = CustomOp
	
	# ==================================================================================
	
	def SetInitialTime(self, InitialTime):
		
		self.TimeLine.SetInitialTime(InitialTime)
	
	# ==================================================================================
	
	def GetEndTime(self):
		
		return self.TimeLine.EndTime
	
	# ==================================================================================
	
	def Initialize(self):
		if(self.CustomOp is not None):
			for cop in self.CustomOp:
				cop.Initialize(self.ModelPart)
	
	# ==================================================================================
	
	def Finalize(self):
		if(self.CustomOp is not None):
			for cop in self.CustomOp:
				cop.Finalize(self.ModelPart)
	
	# ==================================================================================
	
	def Solve(self):
		
		self.PrintHeader()
		
		self.IsConverged = True
		
		timer_beg = time.time()
		
		last_time_step = self.TimeLine.InitialTime
		time_step_counter = 0
		
		while True:
			
			hasNextTimeStep, NextTimeStep = self.TimeLine.NextTimeStep(
				LastIterationConverged = self.IsConverged
				)
			
			if(hasNextTimeStep == False):
				self.PrintNonConvergence()
				self.IsConverged = False
				break
			
			self.ModelPart.CloneTimeStep(NextTimeStep)
			
			time_step_counter += 1
			self.ModelPart.ProcessInfo[TIME_STEPS] = time_step_counter
			
			# custom operation : OnBeforeSolutionStep
			
			for cop in self.CustomOp:
				cop.OnBeforeSolutionStep(self.ModelPart)
			
			# solve
			
			self.Solver.Solve()
			
			if(self.Solver.IsConverged() == False):
				self.PrintNonConvergence()
				self.IsConverged = False
				# continue
			
			last_time_step = NextTimeStep
			self.IsConverged = True
			
			if(self.ResultsIO != None):
				self.ResultsIO.Write(NextTimeStep)
			
			# end solve
			
			# custom operation : OnSolutionStepCompleted
			
			for cop in self.CustomOp:
				cop.OnSolutionStepCompleted(self.ModelPart)
			
			if(self.TimeLine.Finished):
				self.IsConverged = True
				break
		
		timer_end = time.time()
		
		self.PrintFooter(timer_end - timer_beg)
	
	# ==================================================================================
	
	def PrintHeader(self):
		print ("")
		print ("")
		print ("========================================================================")
		print (" Solution Stage: STATIC GENERAL")
		print ("")
		print (" Initial Time:\t", self.TimeLine.InitialTime, "s.")
		print (" End Time:\t", self.TimeLine.EndTime, "s.")
		print ("")
		print (" Solution started at:")
		print (datetime.datetime.now().strftime(" %Y-%m-%d %H:%M:%S"))
		print ("========================================================================")
		print ("")
	
	# ==================================================================================
	
	def PrintFooter(self, elapsedTime):
		print ("========================================================================")
		print (" Solution terminated at:")
		print (datetime.datetime.now().strftime(" %Y-%m-%d %H:%M:%S"))
		print ("")
		print (" Elapsed: ", elapsedTime, " seconds")
		print ("========================================================================")
		print ("")
		
	# ==================================================================================
	
	def PrintNonConvergence(self):
		print ("")
		print ("========================================================================")
		print (" ERROR:")
		print ("------------------------------------------------------------------------")
		print (" The current Step did not converge.")
		print (" The results of the current step will not be calculated.")
		print ("========================================================================")
		print ("")