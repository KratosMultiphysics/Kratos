from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import datetime
import time
import TK_TimeLines
import TK_LoadFunctions
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.MultiScaleApplication import *

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
			# New...
			REACTION_DISPLACEMENT_LAGRANGE,
			REACTION_ROTATION_LAGRANGE,
			])

# ======================================================================================
#
# Returns a list of all (Physical) D.O.F.s required to use the SolutionStage object in
# this module.
#
# ======================================================================================

def Dofs():
	return ([
			# DISPLACEMENT_LAGRANGE_X,
			# DISPLACEMENT_LAGRANGE_Y,
			# DISPLACEMENT_LAGRANGE_Z,
			# ROTATION_LAGRANGE_X,
			# ROTATION_LAGRANGE_Y,
			# ROTATION_LAGRANGE_Z,
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
			# NEW...
			(DISPLACEMENT_LAGRANGE_X, REACTION_DISPLACEMENT_LAGRANGE_X),
			(DISPLACEMENT_LAGRANGE_Y, REACTION_DISPLACEMENT_LAGRANGE_Y),
			(DISPLACEMENT_LAGRANGE_Z, REACTION_DISPLACEMENT_LAGRANGE_Z),
			(ROTATION_LAGRANGE_X,     REACTION_ROTATION_LAGRANGE_X    ),
			(ROTATION_LAGRANGE_Y,     REACTION_ROTATION_LAGRANGE_Y    ),
			(ROTATION_LAGRANGE_Z,     REACTION_ROTATION_LAGRANGE_Z    ),
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
				DesiredIterations = 10,
				CalculateReactions = True,
				ReformDofSetAtEachStep = False,
				MoveMesh = True,
				ResultsIO = None,
				LinearSolver = SkylineLUFactorizationSolver(),
				ArcLength = 1.0,
				PsiU = 1.0,
				PsiF = 1.0,
				URef = 1.0,
				Parallel = True,
				CustomOp = None):

		# Parallelism
		self.Parallel = Parallel

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
		if(self.Parallel):
			self.BuilderAndSolver = StaticGeneralBuilderAndSolver(self.LinearSolver)
		else:
			self.BuilderAndSolver = StaticGeneralBuilderAndSolverSequential(self.LinearSolver)

		# Misc. Solver parameters
		self.MaxIterations = MaxIterations
		self.CalculateReactions = CalculateReactions
		self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
		self.MoveMesh = MoveMesh
		self.ArcLength = ArcLength
		self.PsiU      = PsiU
		self.PsiF      = PsiF
		self.URef      = URef

		# Results IO
		self.ResultsIO = ResultsIO

		# Create and initialize the solver
		self.Solver = ArcLengthDissNrgStrategy(
			self.ModelPart,
			self.TimeScheme,
			self.LinearSolver,
			self.Convergence,
			self.MaxIterations,
			self.CalculateReactions,
			self.ReformDofSetAtEachStep,
			self.MoveMesh)

		self.Solver.Check();
		self.Solver.SetEchoLevel(1);

		self.IsConverged = False

		self.CustomOp = CustomOp

		# @todo: put it in the c-tor
		self.AdaptiveIncrementation = True
		self.DesiredIterations = DesiredIterations
		if(self.DesiredIterations >= self.MaxIterations):
			self.DesiredIterations = self.MaxIterations/2;

	# ==================================================================================

	def SetInitialTime(self, InitialTime):
		self.TimeLine.SetInitialTime(InitialTime)

	# ==================================================================================

	def GetEndTime(self):
		return self.TimeLine.EndTime

	# ==================================================================================

	def SetTimeBoundsOnProcessInfo(self, initial_time, end_time):
		self.ModelPart.ProcessInfo[START_TIME] = initial_time
		self.ModelPart.ProcessInfo[END_TIME]   = end_time

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

		# stage initializations
		self.IsConverged = True
		timer_start = time.time()
		self.PrintHeader()

		# time incrementation parameters
		tstart = self.TimeLine.InitialTime
		tend = tstart + self.TimeLine.Duration
		delta_time = self.TimeLine.Increment
		current_time = tstart
		time_tolerance = abs(delta_time) * 1.0E-6
		increment_id = 0
		increment_mult = 1.0
		increment_min = self.TimeLine.MinIncrement
		increment_max = self.TimeLine.MaxIncrement

		self.ModelPart.ProcessInfo[TIME] = current_time

		# call custom operations on stage initialization
		# NOTE: before calling it, set LAMBDA to 1.
		# the arc length requires a starting value of LAMBDA equal to 1
		# to calculate the reference load vector.
		self.ModelPart.ProcessInfo[LAMBDA] = 1.0 # to calculate reference load vector
		if(self.CustomOp is not None):
			for cop in self.CustomOp:
				cop.OnBeforeSolutionStage(self.ModelPart) # here BC are set on nodes with LAMBDA=1

		# PREPARE FIRST SOVE==============================================================
		# set some data to the process info
		self.ModelPart.CloneTimeStep(current_time)
		self.ModelPart.ProcessInfo[TIME_STEPS] = increment_id
		self.ModelPart.ProcessInfo[TIME] = current_time
		self.ModelPart.ProcessInfo[DELTA_TIME] = delta_time
		self.Solver.PrepareFirstSolve()
		# PREPARE FIRST SOVE==============================================================

		# nodal values should be set to 0. arc length strategy
		# will use the reference load vector and an internal LAMBDA
		self.ModelPart.ProcessInfo[LAMBDA] = 0.0

		# begin time incrementation loop
		while True:

			# define the new time step size
			old_time = current_time
			new_time = current_time + delta_time # * increment_mult
			if(new_time > tend):
				new_time = tend
			current_delta_time = new_time - current_time
			current_time = new_time
			increment_id += 1

			# set some data to the process info
			self.ModelPart.CloneTimeStep(current_time)
			self.ModelPart.ProcessInfo[TIME_STEPS] = increment_id
			self.ModelPart.ProcessInfo[TIME] = current_time
			self.ModelPart.ProcessInfo[DELTA_TIME] = current_delta_time

			# save current LAMBDA, and then set it to 0 before calling OnBeforeSolutionStep on custom operations
			# lambda_saved = self.ModelPart.ProcessInfo[LAMBDA]
			self.ModelPart.ProcessInfo[LAMBDA] = 0.0

			# custom operations
			for cop in self.CustomOp:
				cop.OnBeforeSolutionStep(self.ModelPart)

			# reset the correct LAMBDA
			# self.ModelPart.ProcessInfo[LAMBDA] = lambda_saved

			# solve current step
			self.Solver.Solve()

			if(self.Solver.IsConverged()):

				# finalize the current step
				# write results
				if(self.ResultsIO != None):
					self.ResultsIO.Write(current_time)
				# custom operations
				for cop in self.CustomOp:
					cop.OnSolutionStepCompleted(self.ModelPart)

				# exit the incrementation loop if the end time has been reached
				self.IsConverged = True
				if(abs(current_time - tend) < time_tolerance):
					break

				# adapt the next time increment size
				if(self.AdaptiveIncrementation):
					target_iter = float(self.DesiredIterations)
					needed_iter = float(self.ModelPart.ProcessInfo[NL_ITERATION_NUMBER])
					if(needed_iter == 0): needed_iter = 1
					# factor = pow(0.5,0.25*(needed_iter-target_iter))
					factor = target_iter/needed_iter
					increment_mult_old = increment_mult
					increment_mult *= factor
					new_delta_time = delta_time * increment_mult
					if(new_delta_time > increment_max):
						#ndt = increment_mult_old*factor*dt = dtmax
						# factor = dtmax/(dt*mult_old)
						factor = increment_max/(delta_time*increment_mult_old)
						increment_mult = increment_mult_old*factor
					self.Solver.SetLoadFactors(factor, increment_mult)

			else:

				self.IsConverged = False

				# adapt the next time increment size
				if(self.AdaptiveIncrementation):

					target_iter = float(self.DesiredIterations)
					needed_iter = float(self.ModelPart.ProcessInfo[NL_ITERATION_NUMBER])
					if(needed_iter == 0): needed_iter = 1
					# factor = pow(0.5,0.25*(needed_iter-target_iter))
					factor = target_iter/needed_iter
					increment_mult *= factor
					new_delta_time = delta_time * increment_mult
					self.Solver.SetLoadFactors(factor, increment_mult)

					if(new_delta_time < increment_min):

						# exit the time incrementation loop
						print("")
						print("TIME", current_time, " - ERROR: ")
						print(" The iteration process did not converge.")
						print(" The current increment (", new_delta_time,
							") is less than the minimum (", increment_min, ")")
						print(" The solution cannot continue")
						print("")
						break;

					else:

						print("")
						print("TIME", current_time, " - WARNING: ")
						print(" The iteration process did not converge.")
						print(" The current increment (", current_delta_time,
							") will be reduced to (", new_delta_time, ")")
						print(" The solution cannot continue")
						print("")
						current_time = old_time

				else:

					# exit the time incrementation loop
					self.PrintNonConvergence()
					break

		# stage finalizations
		timer_end = time.time()
		self.PrintFooter(timer_end - timer_start)

		# call custom operations on stage finalization
		if(self.CustomOp is not None):
			for cop in self.CustomOp:
				cop.OnSolutionStageCompleted(self.ModelPart)

	# ==================================================================================

	def PrintHeader(self):
		print ("")
		print ("")
		print ("========================================================================")
		print (" Solution Stage: ARC LENGTH")
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