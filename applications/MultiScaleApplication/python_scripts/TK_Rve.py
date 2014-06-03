## @package Rve
#  This module contains classes (Rve Modelers) that are used
#  to handle the generation, assignment and tracking of 
#  RveConstitutiveLaws.
#
#  More details

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import datetime
import time
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MultiScaleApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
import TK_Props
CheckForPreviousImport()

## RVEPropertyMap
#
# Detailed description...
class RVEPropertyMap:
	
	## Constructor
	def __init__(self,
				 PropertyID,
				 Values):
		self.PropertyID = PropertyID
		self.Values = Values


## RVEModelPartPrototype
#
# Detailed description...
class RVEModelPartPrototype:
	
	## Constructor
	def __init__(self,
				 ModelName,
				 NodalVariables = ([
					DISPLACEMENT,
					DISPLACEMENT_LAGRANGE,
					]),
				 DOFs = ([
					DISPLACEMENT_X,
					DISPLACEMENT_Y,
					DISPLACEMENT_Z,
					DISPLACEMENT_LAGRANGE_X,
					DISPLACEMENT_LAGRANGE_Y,
					DISPLACEMENT_LAGRANGE_Z,
					]),
				 BufferSize = 2,
				 RVEPropertyMapList = None):
		
		# create the model part
		self.Model = ModelPart(ModelName)
		
		# add nodal variables
		for ivar in NodalVariables:
			self.Model.AddNodalSolutionStepVariable(ivar)
		
		# read the model part
		model_part_io = ModelPartIO(ModelName)
		model_part_io.ReadModelPart(self.Model)
		
		# add all degrees of freedom
		for inode in self.Model.Nodes:
			for idof in DOFs:
				inode.AddDof(idof)
		
		# set buffer size
		self.Model.SetBufferSize(BufferSize)
		
		# set up all the properties
		for ipmap in RVEPropertyMapList:
			TK_Props.Property(
				Pro = self.Model.Properties[ipmap.PropertyID],
				Values = ipmap.Values
				)


## RVEStrainSize
#
# Detailed description...
class RVEStrainSize:
	RVE_PLANE_STRESS = 0
	RVE_PLANE_STRAIN = 1
	RVE_3D = 2

## The RveModeler for continuum elements.
#
#  This class is a specialized RveModeler for continuum elements.
class RVEModelerSolid:

	## Constructor.
	def __init__(self, 
				 MicroModelPartPrototype, 
				 StrainSize,
				 ResultsIOClass,
				 ResultsOnNodes = [], 
				 ResultsOnGaussPoints = [],
				 LinearSolverClass = SuperLUSolver,
				 MaxIterations = 20,
				 CalculateReactions = False,
				 ReformDofSetAtEachIteration = False,
				 MoveMesh = False,
				 ConvergenceCriteriaClass = ResidualNormCriteria,
				 ConvergenceRelativeTolerance = 1.0E-6,
				 ConvergenceAbsoluteTolerance = 1.0E-9,
				 ConvergenceIsVerbose = False,
				 TargetElementList = [],
				 OutputElementList = []):
		
		self.MicroModelPartPrototype = MicroModelPartPrototype
		self.StrainSize = StrainSize
		
		if(self.StrainSize == RVEStrainSize.RVE_PLANE_STRESS):
			self.RveBoundaryManagerClass = RveBoundary2D
			self.RveAdapterClass = RvePlaneStressAdapter
			self.RveMaterialClass = RveConstitutiveLawPlaneStress
		elif(self.StrainSize == RVEStrainSize.RVE_PLANE_STRAIN):
			raise Exception("Rve Plane Strain Not Yet Implemented")
		else: # RVEStrainSize.RVE_3D):
			self.RveBoundaryManagerClass = RveBoundary3D
			self.RveAdapterClass = Rve3DAdapter
			self.RveMaterialClass = RveConstitutiveLaw3D
		
		self.RveBoundaryManager = None
		
		self.ResultsIOClass = ResultsIOClass
		self.ResultsOnNodes = ResultsOnNodes
		self.ResultsOnGaussPoints = ResultsOnGaussPoints
		
		self.LinearSolverClass = LinearSolverClass
		self.MaxIterations = MaxIterations
		self.CalculateReactions = CalculateReactions
		self.ReformDofSetAtEachIteration = ReformDofSetAtEachIteration
		self.MoveMesh = MoveMesh
		
		self.ConvergenceCriteriaClass = ConvergenceCriteriaClass
		self.ConvergenceRelativeTolerance = ConvergenceRelativeTolerance
		self.ConvergenceAbsoluteTolerance = ConvergenceAbsoluteTolerance
		self.ConvergenceIsVerbose = ConvergenceIsVerbose
		
		self.TargetElementList = TargetElementList
		self.OutputElementList = OutputElementList
		
		self.TrackList = {}
		
		self.Initialized = False
	
	## Initialize
	#
	# called at the very beginning of the analysis history to
	# perform all initializations. This method should be called
	# only once.
	def Initialize(self, Model):
		if(self.Initialized == False):
			# initialize the boundary manager.
			self.RveBoundaryManager = self.RveBoundaryManagerClass(self.MicroModelPartPrototype.Model)
			
			# generate,assign and track all required rve's
			for elem_id in self.TargetElementList:
				elem = Model.Elements[elem_id]
				self.__assign_rve_constitutive_law(elem)
			
			# initialize the output
			self.__initialize_output()
			
			# set initialization flag
			self.Initialized = True
	
	## OnBeforeSolutionStep
	#
	# called before each solutions steps is solved
	def OnBeforeSolutionStep(self, Model):
		pass
	
	## OnSolutionStepCompleted
	#
	# called after each solutions steps is solved
	def OnSolutionStepCompleted(self, Model):
		# write the output for this time step
		self.__write_output(Model.ProcessInfo[TIME])
	
	## Finalize
	#
	# called at the end of the analysis history to
	# perform all finalizations. This method should be called
	# only once.
	def Finalize(self, Model):
		if(self.Initialized == True):
			
			# finalize the output
			self.__finalize_output()
	
	# private methods *******************************************************************************************
	
	## __generate_rve_constitutive_law
	#
	# This method generates a new rve constitutive law
	# by cloning the rve model part prototype and creating
	# a new rve constitutive law out of it.
	# This method is meant to be private, do NOT call it explicitly
	def __generate_rve_constitutive_law(self):
		
		modelPartClone = ModelPart(self.MicroModelPartPrototype.Model.Name + "_RVE")
		RveCloneModelPart(self.MicroModelPartPrototype.Model, modelPartClone) # clone the model part prototype
		
		msStatus = RveMacroscaleStatus() # the macro-scale status
		self.RveBoundaryManager.AddConditions(modelPartClone, msStatus) # create the rve boundary conditions
		
		linSolver = self.LinearSolverClass()
		timeScheme = StaticGeneralScheme() # TODO: NOT IN SETTINGS...
		timeScheme.Check(modelPartClone)
		convCriteria = self.ConvergenceCriteriaClass(
			self.ConvergenceRelativeTolerance,
			self.ConvergenceAbsoluteTolerance,
			self.ConvergenceIsVerbose,
			)
		
		rveStrategy = StaticGeneralStrategy( # TODO: NOT IN SETTINGS...
			modelPartClone,
			timeScheme,
			linSolver,
			convCriteria,
			StaticGeneralBuilderAndSolver(linSolver), # TODO: NOT IN SETTINGS...
			self.MaxIterations,
			self.CalculateReactions,
			self.ReformDofSetAtEachIteration,
			self.MoveMesh,
			)
		rveStrategy.SetEchoLevel(0) # TODO: NOT IN SETTINGS...
		
		adapter = self.RveAdapterClass() # generate the rve adapter
		adapter.SetRveData(modelPartClone, rveStrategy, msStatus) # set all data (just for testing...)
		
		rveLaw = self.RveMaterialClass(adapter) # finally generate the constitutive law adapter
		for i in range(modelPartClone.GetBufferSize()):
			modelPartClone.CloneTimeStep(0.0)
		return rveLaw
	
	## __track_rve_constitutive_law
	#
	# This method tracks a rve constitutive law
	# at a given element in a given gauss point.
	# This method is meant to be private, do NOT call it explicitly
	def __track_rve_constitutive_law(self, rveLaw, elemID, gpID):
		elInfo = SolidElementInfo(elemID, gpID)
		outputFileName = self.MicroModelPartPrototype.Model.Name + "__" + elInfo.GetStringExtension()
		rveLawIO = self.ResultsIOClass(rveLaw.GetModelPart(), outputFileName, self.ResultsOnNodes, self.ResultsOnGaussPoints)
		self.TrackList[elInfo] = (rveLaw, rveLawIO)
	
	## __assign_rve_constitutive_law
	#
	# This method assignes a rve constitutive law
	# at a given element.
	# This method is meant to be private, do NOT call it explicitly
	def __assign_rve_constitutive_law(self, Element):
		
		# get the number of integration points
		elemIntPoints = Element.GetIntegrationPoints()
		num_gp = len(elemIntPoints)
		elem_id = Element.Id
		
		# get a reference to the process into
		pinfo = self.MicroModelPartPrototype.Model.ProcessInfo
		
		# prepare the list of constitutive laws for the element
		constitutiveLaws = []
		
		# for each element integration point ...
		for gp_id in range(num_gp):
			
			# generate a new rve constitutive law
			aRveLaw = self.__generate_rve_constitutive_law()
			constitutiveLaws.append(aRveLaw)
			
			# TODO: check what rve law to track...
			# for the moment let's track them all
			self.__track_rve_constitutive_law(aRveLaw, elem_id, gp_id)
		
		# assign the list of constitutive laws
		Element.SetValuesOnIntegrationPoints(CONSTITUTIVE_LAW_POINTER, constitutiveLaws, pinfo)
	
	## Initializes the output for the tracked rves (only if required)
	def __initialize_output(self):
		for key, value in self.TrackList.items():
			rveIO = value[1]
			if(rveIO is not None):
				rveIO.Initialize()
	
	## Writes the output for the tracked rves (only if required)
	def __write_output(self, currentTime):
		for key, value in self.TrackList.items():
			rveIO = value[1]
			if(rveIO is not None):
				rveIO.Write(currentTime)
	
	## Finalizes the output for the tracked rves (only if required)
	def __finalize_output(self):
		for key, value in self.TrackList.items():
			rveIO = value[1]
			if(rveIO is not None):
				rveIO.Finalize()
	
	## Prints an extensive description of this object
	def __print_info(self):
		
		print ("")
		print ("====================================================")
		print ("RveModelerShell - Info:")
		print ("====================================================")
		print ("MODEL PART - PROTOTYPE:")
		print (self.MicroModelPartPrototype.Model)
		print ("====================================================")
		print ("TRACK LIST:")
		ii = 0
		print ("+--------------------------------------------------------+")
		for key, value in self.TrackList.items():
			print ("AT[", ii, "]")
			print ("Info:")
			print (key)
			print ("(RveMaterial, IO)")
			print (value)
			print ("Micro Model Clone:")
			micro = value[0].GetModelPart()
			print (hex(id(micro)))
			print (micro)
			print ("+--------------------------------------------------------+")
			ii+=1