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
from KratosMultiphysics.StructuralMechanicsApplication import *
import TK_Props

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
					REACTION,
					]),
				 DOFs = ([
					(DISPLACEMENT_X, REACTION_X),
					(DISPLACEMENT_Y, REACTION_Y),
					(DISPLACEMENT_Z, REACTION_Z),
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
				inode.AddDof(idof[0], idof[1])

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
	RVE_THICK_SHELL = 3
	RVE_THIN_SHELL = 4

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
				 RveConstraintHandlerClass = RveConstraintHandler_ZBF_SD,
				 RveHomogenizerClass = RveHomogenizer,
				 SchemeClass = RveStaticScheme,
				 LinearSolverClass = SuperLUSolver,
				 MaxIterations = 10,
				 CalculateReactions = False,
				 ReformDofSetAtEachIteration = False,
				 MoveMesh = False,
				 ConvergenceCriteriaClass = ResidualNormCriteria,
				 ConvergenceRelativeTolerance = 1.0E-6,
				 ConvergenceAbsoluteTolerance = 1.0E-9,
				 ConvergenceIsVerbose = False,
				 TargetElementList = [],
				 OutputElementList = [],
				 BoundingPolygonNodesID = None):

		self.MicroModelPartPrototype = MicroModelPartPrototype
		self.StrainSize = StrainSize

		self.BoundingPolygonNodesID = BoundingPolygonNodesID

		if(self.StrainSize == RVEStrainSize.RVE_PLANE_STRESS):
			self.RveAdapterClass = RvePlaneStressAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV2PlaneStress
		elif(self.StrainSize == RVEStrainSize.RVE_3D):
			self.RveAdapterClass = Rve3DAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV23D
		elif(self.StrainSize == RVEStrainSize.RVE_THICK_SHELL):
			self.RveAdapterClass = RveThickShellAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV2ThickShell
		else:
			raise Exception("Rve Type Not Yet Implemented")

		self.RveGeometryDescr = None

		self.ResultsIOClass = ResultsIOClass
		self.ResultsOnNodes = ResultsOnNodes
		self.ResultsOnGaussPoints = ResultsOnGaussPoints

		self.RveConstraintHandlerClass = RveConstraintHandlerClass
		self.RveHomogenizerClass       = RveHomogenizerClass
		self.SchemeClass               = SchemeClass

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

			# initialize the geometry descriptor
			self.RveGeometryDescr = RveGeometryDescriptor()
			if(self.BoundingPolygonNodesID is not None):
				self.RveGeometryDescr.SetUserCornerNodes(self.BoundingPolygonNodesID)
			self.RveGeometryDescr.Build(self.MicroModelPartPrototype.Model)
			print(self.RveGeometryDescr)

			# generate,assign and track all required rve's
			for elem_id in self.TargetElementList:
				elem = Model.Elements[elem_id]
				self.__assign_rve_constitutive_law(elem)

			# initialize the output
			self.__initialize_output()

			# set initialization flag
			self.Initialized = True

	## OnBeforeSolutionStage
	#
	# called before each solutions stage
	def OnBeforeSolutionStage(self, Model):
		pass

	## OnSolutionStageCompleted
	#
	# called after each solutions stage
	def OnSolutionStageCompleted(self, Model):
		pass

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

		msData = RveMacroscaleData()

		linSolver = self.LinearSolverClass()

		timeScheme = self.SchemeClass()
		timeScheme.Check(modelPartClone)

		convCriteria = self.ConvergenceCriteriaClass(
			self.ConvergenceRelativeTolerance,
			self.ConvergenceAbsoluteTolerance,
			self.ConvergenceIsVerbose,
			)

		constraint_handler = self.RveConstraintHandlerClass()

		homogenizer = self.RveHomogenizerClass()

		adapter = self.RveAdapterClass() # generate the rve adapter
		adapter.SetRveData(
			modelPartClone,
			msData,
			self.RveGeometryDescr,
			constraint_handler,
			RveLinearSystemOfEquations(linSolver),
			homogenizer,
			timeScheme,
			convCriteria
		) # set all data (just for testing...)

		rveLaw = self.RveMaterialClass(adapter) # finally generate the constitutive law adapter
		for i in range(modelPartClone.GetBufferSize()):
			modelPartClone.CloneTimeStep(0.0)
		return rveLaw

	## __generate_rve_constitutive_law
	#
	# This method generates a new rve constitutive law
	# by cloning the rve model part prototype and creating
	# a new rve constitutive law out of it.
	# This method is meant to be private, do NOT call it explicitly
	def __generate_rve_constitutive_law_2(self,elid):

		modelPartClone = ModelPart(self.MicroModelPartPrototype.Model.Name + "_RVE")
		RveCloneModelPart(self.MicroModelPartPrototype.Model, modelPartClone) # clone the model part prototype

		msData = RveMacroscaleData()

		linSolver = self.LinearSolverClass()

		timeScheme = self.SchemeClass()
		timeScheme.Check(modelPartClone)

		convCriteria = self.ConvergenceCriteriaClass(
			self.ConvergenceRelativeTolerance,
			self.ConvergenceAbsoluteTolerance,
			self.ConvergenceIsVerbose,
			)

		constraint_handler = None
		bel = [1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,16,17,24,25,32,33,40,41,48,49,56]
		if(elid in bel):
			constraint_handler = RveConstraintHandler_PBFZR_SD_ThickShell()
		else:
			constraint_handler = self.RveConstraintHandlerClass()

		homogenizer = self.RveHomogenizerClass()

		adapter = self.RveAdapterClass() # generate the rve adapter
		adapter.SetRveData(
			modelPartClone,
			msData,
			self.RveGeometryDescr,
			constraint_handler,
			RveLinearSystemOfEquations(linSolver),
			homogenizer,
			timeScheme,
			convCriteria
		) # set all data (just for testing...)

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
		if( next((x for x in self.OutputElementList if x == elemID), None) is not None ):
			outputFileName = self.MicroModelPartPrototype.Model.Name + "__" + elInfo.GetStringExtension()
			rveLawIO = self.ResultsIOClass(rveLaw.GetModelPart(), outputFileName, self.ResultsOnNodes, self.ResultsOnGaussPoints)
			self.TrackList[elInfo] = (rveLaw, rveLawIO)
		else:
			self.TrackList[elInfo] = (rveLaw, None)

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
			# aRveLaw = self.__generate_rve_constitutive_law_2(elem_id)
			constitutiveLaws.append(aRveLaw)

			# TODO: check what rve law to track...
			# for the moment let's track them all
			self.__track_rve_constitutive_law(aRveLaw, elem_id, gp_id)

		# assign the list of constitutive laws
		Element.SetValuesOnIntegrationPoints(CONSTITUTIVE_LAW, constitutiveLaws, pinfo)

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

	def GENERATE_RVE_LAW(self):
		return self.__generate_rve_constitutive_law()

	def TRACK_RVE_LAW(self,rveLaw,elemID,gpID):
		return self.__track_rve_constitutive_law(rveLaw,elemID,gpID)

	def INIT_OUTPUT(self):
		return self.__initialize_output()

	def WRITE_OUTPUT(self,time):
		return self.__write_output(time)

	def FIN_OUTPUT(self):
		return self.__finalize_output()

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


# class RVE_MPI_Utils ( 1) copy prototype; 2) partition target/output elements; 3) generate a modeler for each node)