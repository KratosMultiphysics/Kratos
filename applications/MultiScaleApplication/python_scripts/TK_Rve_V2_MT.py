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

class RVEModelPartPrototype_2Phisics:

	## Constructor
	def __init__(self,
				 ModelNameA,
				 ModelNameB,
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
				 RVEPropertyMapListA = None,
				 RVEPropertyMapListB = None):

		# create the model part
		self.ModelA = ModelPart(ModelNameA)
		self.ModelB = ModelPart(ModelNameB)

		# add nodal variables
		for ivar in NodalVariables:
			self.ModelA.AddNodalSolutionStepVariable(ivar)

		# read the model part
		model_part_io = ModelPartIO(ModelNameA)
		model_part_io.ReadModelPart(self.ModelA)

		# add all degrees of freedom
		for inode in self.ModelA.Nodes:
			for idof in DOFs:
				inode.AddDof(idof[0], idof[1])

		# set buffer size
		self.ModelA.SetBufferSize(BufferSize)
		self.ModelB.SetBufferSize(BufferSize)

		# preserve connectivities
		model_part_io_conn_preserv = ModelPartIOConnPreserver(ModelNameB, self.ModelA)
		model_part_io_conn_preserv.ReadModelPart(self.ModelB)

		# set up all the properties
		for ipmap in RVEPropertyMapListA:
			TK_Props.Property(
				Pro = self.ModelA.Properties[ipmap.PropertyID],
				Values = ipmap.Values
				)
		# set up all the properties
		for ipmap in RVEPropertyMapListB:
			TK_Props.Property(
				Pro = self.ModelB.Properties[ipmap.PropertyID],
				Values = ipmap.Values
				)

## RVEStrainSize
#
# Detailed description...
class RVEStrainSize:
	RVE_PLANE_STRESS = 0
	RVE_PLANE_STRAIN = 1
	RVE_3D = 2
	RVE_THERMAL_PLANE_STRESS = 3
	RVE_THERMAL_3D = 4

## The RveModeler for continuum elements.
#
#  This class is a specialized RveModeler for continuum elements in 1 Physics.
class RVEModelerSolid:

	## Constructor.
	def __init__(self,
				 MicroModelPart, # Actual Physics
				 MicroModelPartPrototype, # Database Of 2 Physics (In 1 Physics MicroModelPartPrototype.Model = MicroModelPartPrototype)
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
				 BoundingPolygonNodesID = None,
				 # NEW
				 SecondaryRveModeler = None,
				 IsSecondary = False):

		self.MicroModelPart = MicroModelPart
		self.MicroModelPartPrototype = MicroModelPartPrototype
		self.StrainSize = StrainSize

		self.BoundingPolygonNodesID = BoundingPolygonNodesID

		self.IsSecondary = IsSecondary

		if(self.StrainSize == RVEStrainSize.RVE_THERMAL_PLANE_STRESS):
			self.RveAdapterClass = RveThermal2DAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV2Thermal2D
		elif(self.StrainSize == RVEStrainSize.RVE_THERMAL_3D):
			self.RveAdapterClass = RveThermal3DAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV2Thermal3D
		elif(self.StrainSize == RVEStrainSize.RVE_PLANE_STRESS):
			self.RveAdapterClass = RvePlaneStressAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV2PlaneStress
		elif(self.StrainSize == RVEStrainSize.RVE_PLANE_STRAIN):
			raise Exception("Rve Plane Strain Not Yet Implemented")
		else: # RVEStrainSize.RVE_3D):
			self.RveAdapterClass = Rve3DAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV23D

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

		# NEW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		self.SecondaryRveModeler = SecondaryRveModeler
		self.IsSecondary         = IsSecondary
		# if(self.IsSecondary == False):
			# if(self.SecondaryRveModeler is None):
				# raise exeption(" -- RVEModelerSolid is the first physic and need SecondaryRveModeler for add RVE_Clone_ModelPart_List -- ")
		# NEW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
			self.RveGeometryDescr.Build(self.MicroModelPart.Model)
			#print(self.RveGeometryDescr)

			# generate,assign and track all required rve's
			if(self.IsSecondary == True):
				# il primario ha generato la lista di cloni
				# ho bisogno di sapere in quale id mi trovo della lista
				self.clone_list_counter = 0
				for elem_id in self.TargetElementList:
					elem = Model.Elements[elem_id]
					dummy = self.__assign_rve_constitutive_law(elem)
				self.clone_list_counter = 0 # non necessario ma per sicurazzo lo riazzeriamo
			else:
				# se sono il primario genero una lista di [nelem*ngauss] di rve clones...
				self.stored_rvemdpa_clones=[]
				for elem_id in self.TargetElementList:
					elem = Model.Elements[elem_id]
					elem_rvemdpa_clone_list = self.__assign_rve_constitutive_law(elem)
					for iclone in elem_rvemdpa_clone_list:
						self.stored_rvemdpa_clones.append(iclone)
				# ... e la copio nel modeler secondario (che non dovra generarla!!!!!)
				if(self.SecondaryRveModeler is not None):
					self.SecondaryRveModeler.stored_rvemdpa_clones = self.stored_rvemdpa_clones

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

		if(self.IsSecondary == True):
			current_rve_primary_clone = self.stored_rvemdpa_clones[ self.clone_list_counter ]
			modelPartClone = ModelPart(self.MicroModelPart.Model.Name + "_RVE")
			RveCloneModelPart2Physics(self.MicroModelPart.Model, current_rve_primary_clone, modelPartClone) # clone the model part prototype

			self.clone_list_counter = self.clone_list_counter + 1
		else:
			modelPartClone = ModelPart(self.MicroModelPart.Model.Name + "_RVE")
			RveCloneModelPart(self.MicroModelPart.Model, modelPartClone) # clone the model part prototype

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

		if (self.IsSecondary == False):
			for i in range(modelPartClone.GetBufferSize()):
				modelPartClone.CloneTimeStep(0.0)
		return  (rveLaw,modelPartClone) # return a tuple

	## __track_rve_constitutive_law
	#
	# This method tracks a rve constitutive law
	# at a given element in a given gauss point.
	# This method is meant to be private, do NOT call it explicitly
	def __track_rve_constitutive_law(self, rveLaw, elemID, gpID):
		elInfo = SolidElementInfo(elemID, gpID)

		if( next((x for x in self.OutputElementList if x == elemID), None) is not None ):
			outputFileName = self.MicroModelPart.Model.Name + "__" + elInfo.GetStringExtension()
			rveLawIO = self.ResultsIOClass(rveLaw.GetModelPart(), outputFileName, self.ResultsOnNodes, self.ResultsOnGaussPoints)
			# if (self.IsSecondary == False):
				# print ("ResultsIOClass Mechanical Mdpa")
				# rveLawIO = self.ResultsIOClass(rveLaw.GetModelPart(), self.MicroModelPartB.Model, outputFileName, self.ResultsOnNodes, self.ResultsOnGaussPoints_ModA, self.ResultsOnGaussPoints_ModB)
			self.TrackList[elInfo] = (rveLaw, rveLawIO)
		else:
			self.TrackList[elInfo] = (rveLaw, None)

	## __assign_rve_constitutive_law
	#
	# This method assignes a rve constitutive law
	# at a given element.
	# This method is meant to be private, do NOT call it explicitly
	def __assign_rve_constitutive_law(self, Element):

		# list of generated rve mdpa clones
		rve_mdpa_clones = []

		# get the number of integration points
		elemIntPoints = Element.GetIntegrationPoints()
		num_gp = len(elemIntPoints)
		elem_id = Element.Id

		# get a reference to the process into
		pinfo = self.MicroModelPart.Model.ProcessInfo

		# prepare the list of constitutive laws for the element
		constitutiveLaws = []

		# for each element integration point ...
		for gp_id in range(num_gp):

			# generate a new rve constitutive law
			rve_law__rve_mdpa__tuple = self.__generate_rve_constitutive_law()
			aRveLaw = rve_law__rve_mdpa__tuple[0]
			constitutiveLaws.append(aRveLaw)

			# TODO: check what rve law to track...
			# for the moment let's track them all
			self.__track_rve_constitutive_law(aRveLaw, elem_id, gp_id)

			# store the rve mdpa clone
			rve_mdpa_clones.append(rve_law__rve_mdpa__tuple[1])

		# assign the list of constitutive laws
		Element.SetValuesOnIntegrationPoints(CONSTITUTIVE_LAW_POINTER, constitutiveLaws, pinfo)

		return rve_mdpa_clones

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
		print (self.MicroModelPart.Model)
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

#
class RVEModelerSolidPredictor:

	## Constructor.
	def __init__(self,
				 MicroModelPart, # Actual Physics
				 MicroModelPartPrototype, # Database Of 2 Physics (In 1 Physics MicroModelPartPrototype.Model = MicroModelPartPrototype)
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
				 BoundingPolygonNodesID = None,
				 # NEW
				 C0FileName = [],
				 HashTagFileName = [],
				 JsonFileName = [],
				 SecondaryRveModeler = None,
				 IsSecondary = False):

		self.MicroModelPart = MicroModelPart
		self.MicroModelPartPrototype = MicroModelPartPrototype
		self.StrainSize = StrainSize

		self.BoundingPolygonNodesID = BoundingPolygonNodesID

		self.IsSecondary = IsSecondary

		self.C0FileName = C0FileName
		self.HashTagFileName = HashTagFileName
		self.JsonFileName = JsonFileName

		if(self.StrainSize == RVEStrainSize.RVE_THERMAL_PLANE_STRESS):
			self.RveAdapterClass = RveThermal2DAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV2Thermal2D
		elif(self.StrainSize == RVEStrainSize.RVE_THERMAL_3D):
			self.RveAdapterClass = RveThermal3DAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV2Thermal3D
		elif(self.StrainSize == RVEStrainSize.RVE_PLANE_STRESS):
			self.RveAdapterClass = RvePlaneStressAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV2PlaneStress
		elif(self.StrainSize == RVEStrainSize.RVE_PLANE_STRAIN):
			raise Exception("Rve Plane Strain Not Yet Implemented")
		else: # RVEStrainSize.RVE_3D):
			self.RveAdapterClass = Rve3DAdapterV2
			self.RveMaterialClass = RveConstitutiveLawV23D

		self.RveGeometryDescr = None

		self.RvePredictorCalc = None

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
		self.ActiveElementList = []
		self.AdapterDict = {}
		self.RveGenReq = False

		self.TrackList = {}

		self.Initialized = False

		# NEW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		self.SecondaryRveModeler = SecondaryRveModeler
		self.IsSecondary         = IsSecondary
		# if(self.IsSecondary == False):
			# if(self.SecondaryRveModeler is None):
				# raise exeption(" -- RVEModelerSolid is the first physic and need SecondaryRveModeler for add RVE_Clone_ModelPart_List -- ")
		# NEW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
			self.RveGeometryDescr.Build(self.MicroModelPart.Model)
			#print(self.RveGeometryDescr)

			# initialize the Predictor with HashTag in orter to create the maps
			if(self.IsSecondary == False):
				# print("self.HashTagFileName: ",self.HashTagFileName[0])
				self.RvePredictorCalc = RvePredictorCalculator(self.C0FileName[0],self.HashTagFileName[0],self.JsonFileName[0])
			else:
				self.RvePredictorCalc = RvePredictorCalculator("ThermalAnalysis"," do not need"," Prediction")
				#print(self.RvePredictorCalc)

			# generate,assign and track all required rve's
			if(self.IsSecondary == True):
				# il primario ha generato la lista di cloni
				# ho bisogno di sapere in quale id mi trovo della lista
				self.clone_list_counter = 0
				for elem_id in self.TargetElementList:
					elem = Model.Elements[elem_id]
					dummy = self.__assign_rve_constitutive_law(elem)
				self.clone_list_counter = 0 # non necessario ma per sicurezza lo riazzeriamo
			else:
				# se sono il primario genero una lista di [nelem*ngauss] di rve clones...
				self.stored_rvemdpa_clones=[]
				for elem_id in self.TargetElementList:
					elem = Model.Elements[elem_id]
					elem_rvemdpa_clone_list = self.__assign_rve_constitutive_law(elem)
					for iclone in elem_rvemdpa_clone_list:
						self.stored_rvemdpa_clones.append(iclone)
				# ... e la copio nel modeler secondario (che non dovra generarla!!!!!)
				self.SecondaryRveModeler.stored_rvemdpa_clones = self.stored_rvemdpa_clones

			# initialize the output
			# self.__initialize_output()

			# set initialization flag
			self.Initialized = True

	## OnBeforeSolutionStep
	#
	# called before each solutions steps is solved
	def OnBeforeSolutionStep(self, Model):

		print(' OnBeforeSolutionStep...')

		# generate,assign and track all required rve's
		if(self.IsSecondary == False):
			# print('1. Mechanical')
			# se sono il meccanico genero una lista di [nelem*ngauss] di rve clones...
			for elem_id in self.TargetElementList:
				#Take information of Elements that needs Rve
				non_linear_flags = []
				non_linear_flags = Model.Elements[elem_id].GetValuesOnIntegrationPoints(RVE_NON_LINEAR_FLAG,Model.ProcessInfo)
				for i_non_linear_flag in non_linear_flags:
					if Model.ProcessInfo[RVE_PREDICTION_FLAG] == -1:
						i_non_linear_flag[0] = 1.0
					# print('i_non_linear_flag: ',i_non_linear_flag)
					if i_non_linear_flag[0] == 1.0:
						# print(' RVE_NON_LINEAR_FLAG != 0.0 for the element ', elem_id)
						if self.ActiveElementList.count(elem_id) == 0:
							self.RveGenReq = True
							# set initialization flag
							self.ActiveElementList.append(elem_id)
							elem = Model.Elements[elem_id]
							self.__assign_rve_constitutive_law_from_adapter(elem)
		if (len(self.OutputElementList) <= 0):
			self.OutputElementList = self.ActiveElementList
		# self.OutputElementList = self.ActiveElementList
		# initialize the output
		self.__initialize_output()
		print(' Number of Active Elements: ',len(self.ActiveElementList), ' - Over ', len(self.TargetElementList), ' total elements')

		# pass

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

		if(self.IsSecondary == True):
			current_rve_primary_clone = self.stored_rvemdpa_clones[ self.clone_list_counter ]
			modelPartClone = ModelPart(self.MicroModelPart.Model.Name + "_RVE")
			RveCloneModelPart2Physics(self.MicroModelPart.Model, current_rve_primary_clone, modelPartClone) # clone the model part prototype

			self.clone_list_counter = self.clone_list_counter + 1
		else:
			modelPartClone = ModelPart(self.MicroModelPart.Model.Name + "_RVE")
			RveCloneModelPart(self.MicroModelPart.Model, modelPartClone) # clone the model part prototype

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

		# RveGenReq = adapter.RveGenerationRequested()
		# print('RveGenReq: ',self.RveGenReq)

		if (self.IsSecondary == False):
			if (self.RveGenReq == False):
				# print('SetPredictorData...')
				adapter.SetPredictorData(
					modelPartClone,
					self.RvePredictorCalc
				)
			else:
				adapter.SetRveDataAfterPredictor(
					# modelPartClone,
					msData,
					self.RveGeometryDescr,
					constraint_handler,
					RveLinearSystemOfEquations(linSolver),
					homogenizer,
					timeScheme,
					convCriteria
				) # set all data (just for testing...)
		else:
			adapter.SetPredictorData(
						modelPartClone,
						self.RvePredictorCalc
					)
			adapter.SetRveDataAfterPredictor(
				# modelPartClone,
				msData,
				self.RveGeometryDescr,
				constraint_handler,
				RveLinearSystemOfEquations(linSolver),
				homogenizer,
				timeScheme,
				convCriteria
			) # set all data (just for testing...)

		rveLaw = self.RveMaterialClass(adapter) # finally generate the constitutive law adapter

		if (self.IsSecondary == False):
			for i in range(modelPartClone.GetBufferSize()):
				modelPartClone.CloneTimeStep(0.0)

		return  (rveLaw,modelPartClone,adapter) # occhio return a tuple <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	## __generate_rve_constitutive_law_from_adapter
	#
	# This method generates a new rve constitutive law
	# by cloning the rve model part prototype and creating
	# a new rve constitutive law out of it.
	# This method is meant to be private, do NOT call it explicitly
	# ONLY MECHANICAL MODELPART USE IT
	def __generate_rve_constitutive_law_from_adapter(self, Adapter):
		#                        [0]         [1]        [2]
		# Adapter is a tuple: {aRveLaw,modelPartClone,adapter}

		modelPartClone = Adapter[1]

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

		adapt = Adapter[2]

		adapt.SetRveDataAfterPredictor(
			# modelPartClone,
			msData,
			self.RveGeometryDescr,
			constraint_handler,
			RveLinearSystemOfEquations(linSolver),
			homogenizer,
			timeScheme,
			convCriteria
		) # set all data (just for testing...)

		return

	## __track_rve_constitutive_law
	#
	# This method tracks a rve constitutive law
	# at a given element in a given gauss point.
	# This method is meant to be private, do NOT call it explicitly
	def __track_rve_constitutive_law(self, rveLaw, elemID, gpID):
		elInfo = SolidElementInfo(elemID, gpID)

		if( next((x for x in self.OutputElementList if x == elemID), None) is not None ):
			outputFileName = self.MicroModelPart.Model.Name + "__" + elInfo.GetStringExtension()
			rveLawIO = self.ResultsIOClass(rveLaw.GetModelPart(), outputFileName, self.ResultsOnNodes, self.ResultsOnGaussPoints)
			# if (self.IsSecondary == False):
				# print ("ResultsIOClass Mechanical Mdpa")
				# rveLawIO = self.ResultsIOClass(rveLaw.GetModelPart(), self.MicroModelPartB.Model, outputFileName, self.ResultsOnNodes, self.ResultsOnGaussPoints_ModA, self.ResultsOnGaussPoints_ModB)
			self.TrackList[elInfo] = (rveLaw, rveLawIO)
		else:
			self.TrackList[elInfo] = (rveLaw, None)

	## __assign_rve_constitutive_law
	#
	# This method assignes a rve constitutive law
	# at a given element.
	# This method is meant to be private, do NOT call it explicitly
	def __assign_rve_constitutive_law(self, Element):

		# list of generated rve mdpa clones
		rve_mdpa_clones = []

		# get the number of integration points
		elemIntPoints = Element.GetIntegrationPoints()
		num_gp = len(elemIntPoints)
		elem_id = Element.Id

		# get a reference to the process into
		pinfo = self.MicroModelPart.Model.ProcessInfo

		# prepare the list of constitutive laws for the element
		constitutiveLaws = []

		# for each element integration point ...
		for gp_id in range(num_gp):

			# generate a new rve constitutive law
			rve_law__rve_mdpa__tuple = self.__generate_rve_constitutive_law()
			aRveLaw = rve_law__rve_mdpa__tuple[0]
			constitutiveLaws.append(aRveLaw)

			# TODO: check what rve law to track...
			# for the moment let's track them all
			self.__track_rve_constitutive_law(aRveLaw, elem_id, gp_id)

			# store the rve mdpa clone
			rve_mdpa_clones.append(rve_law__rve_mdpa__tuple[1])

			#Fill Adapter Dict
			adapt = rve_law__rve_mdpa__tuple[2]
			self.AdapterDict[elem_id,gp_id] = (rve_law__rve_mdpa__tuple)

		# assign the list of constitutive laws
		Element.SetValuesOnIntegrationPoints(CONSTITUTIVE_LAW_POINTER, constitutiveLaws, pinfo)

		return rve_mdpa_clones

	## __assign_rve_constitutive_law_from_adapter
	#
	# This method assignes a rve constitutive law
	# at a given element.
	# This method is meant to be private, do NOT call it explicitly
	def __assign_rve_constitutive_law_from_adapter(self, Element):

		# get the number of integration points
		elemIntPoints = Element.GetIntegrationPoints()
		num_gp = len(elemIntPoints)
		elem_id = Element.Id

		# for each element integration point ...
		for gp_id in range(num_gp):
			# generate a new rve constitutive law
			self.__generate_rve_constitutive_law_from_adapter(self.AdapterDict[elem_id,gp_id])
			aRveLaw = self.AdapterDict[elem_id,gp_id][0]
			# TODO: check what rve law to track...
			# for the moment let's track them all
			self.__track_rve_constitutive_law(aRveLaw, elem_id, gp_id)

		return

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
		print (self.MicroModelPart.Model)
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