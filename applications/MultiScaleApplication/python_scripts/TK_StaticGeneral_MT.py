from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import datetime
import time
import TK_TimeLines
from KratosMultiphysics import *
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
			# Thermic
			TEMPERATURE,
			TEMPERATURE_REACTION,
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
			# Thermic
			(TEMPERATURE, TEMPERATURE_REACTION),
			])

# ======================================================================================
#
# The SolutionStage object.
# This class solves a stage of the analysis history using one or more time steps,
# and it takes care of saving the results if required.
#
# ======================================================================================

class ALGO_TYPE:
	FULL_NEWTON=1
	KRYLOV_NEWTON=2

class SolutionStage:

	# ==================================================================================

	def __init__(
				self,
				ModelPartMechanical,
				ModelPartThermal,
				TimeLine = TK_TimeLines.FixedTimeLine(
					Duration = 1.0,
					Increment = 1.0
					),
				ConvergenceCriteriaMechanical = ResidualNormCriteria(1.0E-4, 1.0E-4),
				ConvergenceCriteriaThermal = ResidualNormCriteria(1.0E-4, 1.0E-4),
				MaxIterations = 30,
				DesiredIterations = 10,
				CalculateReactions = True,
				ReformDofSetAtEachStep = False,
				MoveMesh = True,
				ResultsIO = None,
				LinearSolverMechanical = SkylineLUFactorizationSolver(),
				LinearSolverThermal = SkylineLUFactorizationSolver(),
				CalculateTangent = True,
				Parallel = True,
				CustomOpMechanical = None,
				CustomOpThermal = None,
				AlgorithmMechanical = ALGO_TYPE.FULL_NEWTON,
				AlgorithmThermal = ALGO_TYPE.FULL_NEWTON):

		# Parallelism
		self.Parallel = Parallel

		# Model Parts
		self.ModelPartMechanical = ModelPartMechanical
		self.ModelPartThermal = ModelPartThermal

		# Time line
		self.TimeLine = TimeLine

		# Time Scheme
		self.TimeSchemeMechanical = StaticGeneralScheme()
		self.TimeSchemeMechanical.Check(self.ModelPartMechanical)
		self.TimeSchemeThermal = StaticGeneralScheme() # per adesso va bene -> (i)->SolutionStepValue += Dx(i)
		self.TimeSchemeThermal.Check(self.ModelPartThermal)

		# LinearSolver
		self.LinearSolverMechanical = LinearSolverMechanical
		self.LinearSolverThermal = LinearSolverThermal

		# ConvergenceCriteria
		self.ConvergenceCriteriaMechanical = ConvergenceCriteriaMechanical
		self.ConvergenceCriteriaMechanical.Check(self.ModelPartMechanical)
		self.ConvergenceCriteriaThermal = ConvergenceCriteriaThermal
		self.ConvergenceCriteriaThermal.Check(self.ModelPartThermal)

		# BuilderAndSolver
		if(self.Parallel):
			self.BuilderAndSolverMechanical = StaticGeneralBuilderAndSolver(self.LinearSolverMechanical)
			self.BuilderAndSolverThermal = StaticGeneralBuilderAndSolver(self.LinearSolverThermal)# per adesso va bene
		else:
			self.BuilderAndSolverMechanical = StaticGeneralBuilderAndSolverSequential(self.LinearSolverMechanical)
			self.BuilderAndSolverThermal = StaticGeneralBuilderAndSolverSequential(self.LinearSolverThermal)# per adesso va bene

		# Misc. SolverParameters
		self.MaxIterations = MaxIterations
		self.CalculateReactions = CalculateReactions
		self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
		self.MoveMesh = MoveMesh

		# Results IO
		self.ResultsIO = ResultsIO

		# Create and initialize the solver
		if(AlgorithmMechanical == ALGO_TYPE.FULL_NEWTON):
			self.SolverMechanical = StaticGeneralStrategy(
				self.ModelPartMechanical,
				self.TimeSchemeMechanical,
				self.LinearSolverMechanical,
				self.ConvergenceCriteriaMechanical,
				self.BuilderAndSolverMechanical,
				self.MaxIterations,
				self.CalculateReactions,
				self.ReformDofSetAtEachStep,
				self.MoveMesh)
		elif(AlgorithmMechanical == ALGO_TYPE.KRYLOV_NEWTON):
			self.SolverMechanical = StaticGeneralStrategyKrylovNewton(
				self.ModelPartMechanical,
				self.TimeSchemeMechanical,
				self.LinearSolverMechanical,
				self.ConvergenceCriteriaMechanical,
				self.BuilderAndSolverMechanical,
				self.MaxIterations,
				self.CalculateReactions,
				self.ReformDofSetAtEachStep,
				self.MoveMesh)
		self.SolverMechanical.SetKeepSystemConstantDuringIterations(not CalculateTangent)
		self.SolverMechanical.Check();
		self.SolverMechanical.SetEchoLevel(1);

		# Create and initialize the solver
		if(AlgorithmThermal == ALGO_TYPE.FULL_NEWTON):
			self.SolverThermal = StaticGeneralStrategy(
				self.ModelPartThermal,
				self.TimeSchemeThermal,
				self.LinearSolverThermal,
				self.ConvergenceCriteriaThermal,
				self.BuilderAndSolverThermal,
				self.MaxIterations,
				self.CalculateReactions,
				self.ReformDofSetAtEachStep,
				self.MoveMesh)
		elif(AlgorithmThermal == ALGO_TYPE.KRYLOV_NEWTON):
			self.SolverThermal = StaticGeneralStrategyKrylovNewton(
				self.ModelPartThermal,
				self.TimeSchemeThermal,
				self.LinearSolverThermal,
				self.ConvergenceCriteriaThermal,
				self.BuilderAndSolverThermal,
				self.MaxIterations,
				self.CalculateReactions,
				self.ReformDofSetAtEachStep,
				self.MoveMesh)
		self.SolverThermal.SetKeepSystemConstantDuringIterations(not CalculateTangent)
		self.SolverThermal.Check();
		self.SolverThermal.SetEchoLevel(1);

		# misc
		self.IsConverged = False
		self.CustomOpMechanical = CustomOpMechanical
		self.CustomOpThermal = CustomOpThermal

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
		self.ModelPartMechanical.ProcessInfo[START_TIME] = initial_time
		self.ModelPartMechanical.ProcessInfo[END_TIME]   = end_time

	# ==================================================================================

	def Initialize(self):
		print(" ")
		print(" #------------------------ CustomOpMechanical ------------------------# ")
		if(self.CustomOpMechanical is not None):
			for cop in self.CustomOpMechanical:
				cop.Initialize(self.ModelPartMechanical)
		print(" ")
		print(" #------------------------ CustomOpThermal ------------------------# ")
		if(self.CustomOpThermal is not None):
			for cop in self.CustomOpThermal:
				cop.Initialize(self.ModelPartThermal)

	# ==================================================================================

	def Finalize(self):
		if(self.CustomOpMechanical is not None):
			for cop in self.CustomOpMechanical:
				cop.Finalize(self.ModelPartMechanical)
		if(self.CustomOpThermal is not None):
			for cop in self.CustomOpThermal:
				cop.Finalize(self.ModelPartThermal)

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

		# begin time incrementation loop
		while True:
			print(' Start WHILE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
			# define the new time step size
			old_time = current_time
			new_time = current_time + delta_time * increment_mult
			if(new_time > tend):
				new_time = tend
			current_delta_time = new_time - current_time
			current_time = new_time
			increment_id += 1

			# set some data to the process info
			self.ModelPartMechanical.CloneTimeStep(current_time)
			self.ModelPartMechanical.ProcessInfo[TIME_STEPS] = increment_id
			self.ModelPartMechanical.ProcessInfo[TIME] = current_time
			self.ModelPartMechanical.ProcessInfo[DELTA_TIME] = current_delta_time
			self.ModelPartThermal.ProcessInfo = self.ModelPartMechanical.ProcessInfo

			# custom operations
			print('START CUSTOM OPERATION OnBeforeSolutionStep')
			for cop in self.CustomOpMechanical:
				cop.OnBeforeSolutionStep(self.ModelPartMechanical)
			# for cop in self.CustomOpThermal:
				# cop.OnBeforeSolutionStep(self.ModelPartThermal)

			# solve the current time step
			temp_converged_flag = False
			print(" #------------------------ SOLVE THERMAL PART ------------------------# ")
			self.SolverThermal.Solve()
			print(" #----------------------- THERMAL PART SOLVED ------------------------# ")
			temp_converged_flag = self.SolverThermal.IsConverged()
			temp_nl_iter_number = float(self.ModelPartThermal.ProcessInfo[NL_ITERATION_NUMBER])
			if(temp_converged_flag):
				print(" #----------------------- SOLVE MECHANICAL PART ------------------------# ")
				self.SolverMechanical.Solve()
				print(" #---------------------- MECHANICAL PART SOLVED ------------------------# ")
				temp_converged_flag = self.SolverMechanical.IsConverged()
				# print('temp_converged_flag',temp_converged_flag)
				temp_nl_iter_number = max(temp_nl_iter_number,float(self.ModelPartMechanical.ProcessInfo[NL_ITERATION_NUMBER]))

			if(temp_converged_flag):

				# finalize the current step
				# write results
				if(self.ResultsIO != None):
					self.ResultsIO.Write(current_time)
				# custom operations
				for cop in self.CustomOpMechanical:
					# print('in self.CustomOpMechanical')
					cop.OnSolutionStepCompleted(self.ModelPartMechanical)
				# for cop in self.CustomOpThermal:
					# cop.OnSolutionStepCompleted(self.ModelPartThermal)

				# exit the incrementation loop if the end time has been reached
				self.IsConverged = True
				if(abs(current_time - tend) < time_tolerance):
					break

				# adapt the next time increment size
				if(self.AdaptiveIncrementation):
					target_iter = float(self.DesiredIterations)
					needed_iter = temp_nl_iter_number
					if(needed_iter == 0):
						needed_iter = 1
					if(needed_iter > 0):
						increment_mult *= target_iter/needed_iter
					new_delta_time = delta_time * increment_mult
					if(new_delta_time > increment_max):
						increment_mult = increment_max / delta_time

			else:

				self.IsConverged = False

				# adapt the next time increment size
				if(self.AdaptiveIncrementation):

					target_iter = float(self.DesiredIterations)
					needed_iter = temp_nl_iter_number
					increment_mult *= target_iter/needed_iter
					new_delta_time = delta_time * increment_mult

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
		print(' End Step <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
		# stage finalizations
		timer_end = time.time()
		self.PrintFooter(timer_end - timer_start)

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

class SolutionStageWithDamageInfo:

	# ==================================================================================

	def __init__(
				self,
				ModelPartMechanical,
				ModelPartThermal,
				TimeLine = TK_TimeLines.FixedTimeLine(
					Duration = 1.0,
					Increment = 1.0
					),
				ConvergenceCriteriaMechanical = ResidualNormCriteria(1.0E-4, 1.0E-4),
				ConvergenceCriteriaThermal = ResidualNormCriteria(1.0E-4, 1.0E-4),
				MaxIterations = 30,
				DesiredIterations = 10,
				CalculateReactions = True,
				ReformDofSetAtEachStep = False,
				MoveMesh = True,
				ResultsIO = None,
				LinearSolverMechanical = SkylineLUFactorizationSolver(),
				LinearSolverThermal = SkylineLUFactorizationSolver(),
				CalculateTangent = True,
				Parallel = True,
				CustomOpMechanical = None,
				CustomOpThermal = None,
				AlgorithmMechanical = ALGO_TYPE.FULL_NEWTON,
				AlgorithmThermal = ALGO_TYPE.FULL_NEWTON):

		# Parallelism
		self.Parallel = Parallel

		# Model Parts
		self.ModelPartMechanical = ModelPartMechanical
		self.ModelPartThermal = ModelPartThermal

		# Time line
		self.TimeLine = TimeLine

		# Time Scheme
		self.TimeSchemeMechanical = StaticGeneralScheme()
		self.TimeSchemeMechanical.Check(self.ModelPartMechanical)
		self.TimeSchemeThermal = StaticGeneralScheme() # per adesso va bene -> (i)->SolutionStepValue += Dx(i)
		self.TimeSchemeThermal.Check(self.ModelPartThermal)

		# LinearSolver
		self.LinearSolverMechanical = LinearSolverMechanical
		self.LinearSolverThermal = LinearSolverThermal

		# ConvergenceCriteria
		self.ConvergenceCriteriaMechanical = ConvergenceCriteriaMechanical
		self.ConvergenceCriteriaMechanical.Check(self.ModelPartMechanical)
		self.ConvergenceCriteriaThermal = ConvergenceCriteriaThermal
		self.ConvergenceCriteriaThermal.Check(self.ModelPartThermal)

		# BuilderAndSolver
		if(self.Parallel):
			self.BuilderAndSolverMechanical = StaticGeneralBuilderAndSolver(self.LinearSolverMechanical)
			self.BuilderAndSolverThermal = StaticGeneralBuilderAndSolver(self.LinearSolverThermal)# per adesso va bene
		else:
			self.BuilderAndSolverMechanical = StaticGeneralBuilderAndSolverSequential(self.LinearSolverMechanical)
			self.BuilderAndSolverThermal = StaticGeneralBuilderAndSolverSequential(self.LinearSolverThermal)# per adesso va bene

		# Misc. SolverParameters
		self.MaxIterations = MaxIterations
		self.CalculateReactions = CalculateReactions
		self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
		self.MoveMesh = MoveMesh

		# Results IO
		self.ResultsIO = ResultsIO

		# Create and initialize the solver
		if(AlgorithmMechanical == ALGO_TYPE.FULL_NEWTON):
			self.SolverMechanical = StaticGeneralStrategy(
				self.ModelPartMechanical,
				self.TimeSchemeMechanical,
				self.LinearSolverMechanical,
				self.ConvergenceCriteriaMechanical,
				self.BuilderAndSolverMechanical,
				self.MaxIterations,
				self.CalculateReactions,
				self.ReformDofSetAtEachStep,
				self.MoveMesh)
		elif(AlgorithmMechanical == ALGO_TYPE.KRYLOV_NEWTON):
			self.SolverMechanical = StaticGeneralStrategyKrylovNewton(
				self.ModelPartMechanical,
				self.TimeSchemeMechanical,
				self.LinearSolverMechanical,
				self.ConvergenceCriteriaMechanical,
				self.BuilderAndSolverMechanical,
				self.MaxIterations,
				self.CalculateReactions,
				self.ReformDofSetAtEachStep,
				self.MoveMesh)
		self.SolverMechanical.SetKeepSystemConstantDuringIterations(not CalculateTangent)
		self.SolverMechanical.Check();
		self.SolverMechanical.SetEchoLevel(1);

		# Create and initialize the solver
		if(AlgorithmThermal == ALGO_TYPE.FULL_NEWTON):
			self.SolverThermal = StaticGeneralStrategy(
				self.ModelPartThermal,
				self.TimeSchemeThermal,
				self.LinearSolverThermal,
				self.ConvergenceCriteriaThermal,
				self.BuilderAndSolverThermal,
				self.MaxIterations,
				self.CalculateReactions,
				self.ReformDofSetAtEachStep,
				self.MoveMesh)
		elif(AlgorithmThermal == ALGO_TYPE.KRYLOV_NEWTON):
			self.SolverThermal = StaticGeneralStrategyKrylovNewton(
				self.ModelPartThermal,
				self.TimeSchemeThermal,
				self.LinearSolverThermal,
				self.ConvergenceCriteriaThermal,
				self.BuilderAndSolverThermal,
				self.MaxIterations,
				self.CalculateReactions,
				self.ReformDofSetAtEachStep,
				self.MoveMesh)
		self.SolverThermal.SetKeepSystemConstantDuringIterations(not CalculateTangent)
		self.SolverThermal.Check();
		self.SolverThermal.SetEchoLevel(1);

		# misc
		self.IsConverged = False
		self.CustomOpMechanical = CustomOpMechanical
		self.CustomOpThermal = CustomOpThermal

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
		self.ModelPartMechanical.ProcessInfo[START_TIME] = initial_time
		self.ModelPartMechanical.ProcessInfo[END_TIME]   = end_time

	# ==================================================================================

	def Initialize(self):
		print(" ")
		print(" #------------------------ CustomOpMechanical ------------------------# ")
		if(self.CustomOpMechanical is not None):
			for cop in self.CustomOpMechanical:
				cop.Initialize(self.ModelPartMechanical)
		print(" ")
		print(" #------------------------ CustomOpThermal ------------------------# ")
		if(self.CustomOpThermal is not None):
			for cop in self.CustomOpThermal:
				cop.Initialize(self.ModelPartThermal)

	# ==================================================================================

	def Finalize(self):
		if(self.CustomOpMechanical is not None):
			for cop in self.CustomOpMechanical:
				cop.Finalize(self.ModelPartMechanical)
		if(self.CustomOpThermal is not None):
			for cop in self.CustomOpThermal:
				cop.Finalize(self.ModelPartThermal)

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

		counter = 0 #Counter for the damaged surfaces
		damage_limit_surf = [self.ModelPartMechanical.ProcessInfo[RVE_DAMAGE_SURFACE_LIMIT]];
		damage_reached_flag = 0
		num_surf = len(damage_limit_surf)
		# begin time incrementation loop
		while True:

			# define the new time step size
			old_time = current_time
			new_time = current_time + delta_time * increment_mult
			if(new_time > tend):
				new_time = tend
			current_delta_time = new_time - current_time
			current_time = new_time
			increment_id += 1

			# set some data to the process info
			self.ModelPartMechanical.CloneTimeStep(current_time)
			self.ModelPartMechanical.ProcessInfo[TIME_STEPS] = increment_id
			self.ModelPartMechanical.ProcessInfo[TIME] = current_time
			self.ModelPartMechanical.ProcessInfo[DELTA_TIME] = current_delta_time
			self.ModelPartThermal.ProcessInfo = self.ModelPartMechanical.ProcessInfo

			#Stefano mod for multiple damage surface
			self.ModelPartMechanical.ProcessInfo[RVE_DAMAGE_SURFACE_LIMIT] = damage_limit_surf[counter]
			filename = "d" + str(damage_limit_surf[counter]) + "_15.txt"
			print(filename)
			# filename = "d" + str(damage_limit_surf) + "_8.txt"
			########

			# custom operations
			for cop in self.CustomOpMechanical:
				cop.OnBeforeSolutionStep(self.ModelPartMechanical)
			# for cop in self.CustomOpThermal:
				# cop.OnBeforeSolutionStep(self.ModelPartThermal)

			# solve the current time step
			temp_converged_flag = False
			print(" #------------------------ SOLVE THERMAL PART ------------------------# ")
			self.SolverThermal.Solve()
			temp_converged_flag = self.SolverThermal.IsConverged()
			temp_nl_iter_number = float(self.ModelPartThermal.ProcessInfo[NL_ITERATION_NUMBER])
			if(temp_converged_flag):
				print(" #----------------------- SOLVE MECHANICAL PART ------------------------# ")
				self.SolverMechanical.Solve()
				temp_converged_flag = self.SolverMechanical.IsConverged()
				temp_nl_iter_number = max(temp_nl_iter_number,float(self.ModelPartMechanical.ProcessInfo[NL_ITERATION_NUMBER]))

			if(temp_converged_flag):

				# finalize the current step
				# write results
				if(self.ResultsIO != None):
					self.ResultsIO.Write(current_time)
				# custom operations
				for cop in self.CustomOpMechanical:
					cop.OnSolutionStepCompleted(self.ModelPartMechanical)
				# for cop in self.CustomOpThermal:
					# cop.OnSolutionStepCompleted(self.ModelPartThermal)

				# exit the incrementation loop if the end time has been reached
				self.IsConverged = True
				if(abs(current_time - tend) < time_tolerance):
					break

				# adapt the next time increment size
				if(self.AdaptiveIncrementation):
					target_iter = float(self.DesiredIterations)
					needed_iter = temp_nl_iter_number
					if(needed_iter == 0):
						needed_iter = 1
					if(needed_iter > 0):
						increment_mult *= target_iter/needed_iter
					new_delta_time = delta_time * increment_mult
					if(new_delta_time > increment_max):
						increment_mult = increment_max / delta_time

			else:

				self.IsConverged = False

				# adapt the next time increment size
				if(self.AdaptiveIncrementation):

					target_iter = float(self.DesiredIterations)
					needed_iter = temp_nl_iter_number
					increment_mult *= target_iter/needed_iter
					new_delta_time = delta_time * increment_mult

					if(new_delta_time < increment_min):
						ofile = open(filename, "a")
						ofile2 = open("FailedAnalysis_" + filename, "a")
						print("num_surf: ",num_surf,"	counter: ",counter)
						if(counter == num_surf-1):
						# exit the time incrementation loop
							print("")
							print("TIME", current_time, " - ERROR: ")
							print(" The iteration process did not converge.")
							print(" The current increment (", new_delta_time,
								") is less than the minimum (", increment_min, ")")
							print(" The solution cannot continue")
							print("")
							print ("------------------------------------------------------------------------")
							print("")
							print("")
							print("")
							# print("             Final Damage: ", damage_limit_surf[counter], " Reached!")
							print("             Final Damage: ", damage_limit_surf, " Reached!")
							print("")
							print("")
							print("")
							print ("------------------------------------------------------------------------")
							# Print Results
							strain = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_TENSOR, self.ModelPartMechanical.ProcessInfo)
							stress = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(RVE_GENERAL_STRESS_TENSOR, self.ModelPartMechanical.ProcessInfo)
							C 	   = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(HOMOGENIZED_CONST_TENS, self.ModelPartMechanical.ProcessInfo)
							tag = self.ModelPartMechanical.ProcessInfo[ACTUAL_TAG]
							reached_damage = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(EQUIVALENT_DAMAGE, self.ModelPartMechanical.ProcessInfo)
							print(reached_damage)
							if (reached_damage[0][0] >= 0.99*damage_limit_surf[counter]) & (reached_damage[0][0] <= 1.01*damage_limit_surf[counter]):
							# if (reached_damage[0][0] >= 0.99*damage_limit_surf) & (reached_damage[0][0] <= 1.01*damage_limit_surf):
								if (len(strain[0]) == 4):
									ofile.write(str(reached_damage[0]) + "   " + str(tag[0]) + "," + str(tag[1]) + "   " + str(strain[0][0]) + " " + str(strain[0][3]) + " " + str(strain[0][1]) + "   " + str(stress[0][0]) + " " + str(stress[0][3]) + " " + str(stress[0][1]) + "   " + str(C[0])+ '\n')
								else:
									ofile.write(str(reached_damage[0]) + "   " + str(tag[0]) + "," + str(tag[1]) + "," + str(tag[2]) + "," + str(tag[3]) + "," + str(tag[4]) + "   " + str(strain[0][0]) + " " + str(strain[0][4]) + " " + str(strain[0][8]) + " " + str(strain[0][1]) + " " + str(strain[0][5]) + " " + str(strain[0][2]) + "   " + str(stress[0][0]) + " " + str(stress[0][4]) + " " + str(stress[0][8]) + " " + str(stress[0][1])  + " " + str(stress[0][5])  + " " + str(stress[0][2])  + "   " + str(C[0])+ '\n')
								ofile.flush()
								ofile.close()
								damage_reached_flag = 1
							else:
								if (len(strain[0]) == 4):
									ofile2.write(str(reached_damage[0]) + "   " + str(tag[0]) + "," + str(tag[1]) + "   " + str(strain[0][0]) + " " + str(strain[0][3]) + " " + str(strain[0][1]) + "   " + str(stress[0][0]) + " " + str(stress[0][3]) + " " + str(stress[0][1]) + "   " + str(C[0])+ '\n')
								else:
									ofile2.write(str(reached_damage[0]) + "   " + str(tag[0]) + "," + str(tag[1]) + "," + str(tag[2]) + "," + str(tag[3]) + "," + str(tag[4]) + "   " + str(strain[0][0]) + " " + str(strain[0][4]) + " " + str(strain[0][8]) + " " + str(strain[0][1]) + " " + str(strain[0][5]) + " " + str(strain[0][2]) + "   " + str(stress[0][0]) + " " + str(stress[0][4]) + " " + str(stress[0][8]) + " " + str(stress[0][1])  + " " + str(stress[0][5])  + " " + str(stress[0][2])  + "   " + str(C[0])+ '\n')
								ofile2.flush()
								ofile2.close()
								damage_reached_flag = 1
							break;
						else:
							print ("------------------------------------------------------------------------")
							print("")
							print("")
							print("")
							print("             Damage: ", damage_limit_surf[counter], " Reached!")
							# print("             Damage: ", damage_limit_surf, " Reached!")
							print("")
							print("")
							print("")
							print ("------------------------------------------------------------------------")
							# Print Results
							strain = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_TENSOR, self.ModelPartMechanical.ProcessInfo)
							stress = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(RVE_GENERAL_STRESS_TENSOR, self.ModelPartMechanical.ProcessInfo)
							C 	   = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(HOMOGENIZED_CONST_TENS, self.ModelPartMechanical.ProcessInfo)
							tag = self.ModelPartMechanical.ProcessInfo[ACTUAL_TAG]
							reached_damage = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(EQUIVALENT_DAMAGE, self.ModelPartMechanical.ProcessInfo)
							print(reached_damage)
							if (reached_damage[0][0] >= 0.99*damage_limit_surf[counter]) & (reached_damage[0][0] <= 1.01*damage_limit_surf[counter]):
							# if (reached_damage[0][0] >= 0.99*damage_limit_surf) & (reached_damage[0][0] <= 1.01*damage_limit_surf):
								if (len(strain[0]) == 4):
									ofile.write(str(reached_damage[0]) + "   " + str(tag[0]) + "," + str(tag[1]) + "   " + str(strain[0][0]) + " " + str(strain[0][3]) + " " + str(strain[0][1]) + "   " + str(stress[0][0]) + " " + str(stress[0][3]) + " " + str(stress[0][1]) + "   " + str(C[0])+ '\n')
								else:
									ofile.write(str(reached_damage[0]) + "   " + str(tag[0]) + "," + str(tag[1]) + "," + str(tag[2]) + "," + str(tag[3]) + "," + str(tag[4]) + "   " + str(strain[0][0]) + " " + str(strain[0][4]) + " " + str(strain[0][8]) + " " + str(strain[0][1]) + " " + str(strain[0][5]) + " " + str(strain[0][2]) + "   " + str(stress[0][0]) + " " + str(stress[0][4]) + " " + str(stress[0][8]) + " " + str(stress[0][1])  + " " + str(stress[0][5])  + " " + str(stress[0][2])  + "   " + str(C[0])+ '\n')
								ofile.flush()
								ofile.close()
								damage_reached_flag = 1
							else:
								if (len(strain[0]) == 4):
									ofile2.write(str(reached_damage[0]) + "   " + str(tag[0]) + "," + str(tag[1]) + "   " + str(strain[0][0]) + " " + str(strain[0][3]) + " " + str(strain[0][1]) + "   " + str(stress[0][0]) + " " + str(stress[0][3]) + " " + str(stress[0][1]) + "   " + str(C[0])+ '\n')
								else:
									ofile2.write(str(reached_damage[0]) + "   " + str(tag[0]) + "," + str(tag[1]) + "," + str(tag[2]) + "," + str(tag[3]) + "," + str(tag[4]) + "   " + str(strain[0][0]) + " " + str(strain[0][4]) + " " + str(strain[0][8]) + " " + str(strain[0][1]) + " " + str(strain[0][5]) + " " + str(strain[0][2]) + "   " + str(stress[0][0]) + " " + str(stress[0][4]) + " " + str(stress[0][8]) + " " + str(stress[0][1])  + " " + str(stress[0][5])  + " " + str(stress[0][2])  + "   " + str(C[0])+ '\n')
								ofile2.flush()
								ofile2.close()
								damage_reached_flag = 1
							counter = counter + 1
							current_time = delta_time

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
		# Print Results
		if (damage_reached_flag == 0):
			fname = "NotEnoughStrain_" + str(damage_limit_surf[counter]) + ".txt"
			ofile = open(fname, "a")
			strain = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_TENSOR, self.ModelPartMechanical.ProcessInfo)
			stress = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(RVE_GENERAL_STRESS_TENSOR, self.ModelPartMechanical.ProcessInfo)
			C 	   = self.ModelPartMechanical.Elements[1].GetValuesOnIntegrationPoints(HOMOGENIZED_CONST_TENS, self.ModelPartMechanical.ProcessInfo)
			tag = self.ModelPartMechanical.ProcessInfo[ACTUAL_TAG]
			ofile.write("[" + str(damage_limit_surf[counter]) + "]" + "   " + str(tag[0]) + "," + str(tag[1]) + "   " + str(strain[0][0]) + " " + str(strain[0][3]) + " " + str(strain[0][1]) + "   " + str(stress[0][0]) + " " + str(stress[0][3]) + " " + str(stress[0][1]) + "   " + str(C[0])+ '\n')
			ofile.flush()
			ofile.close()


		# stage finalizations
		timer_end = time.time()
		self.PrintFooter(timer_end - timer_start)

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