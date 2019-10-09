from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.MultiScaleApplication import *

# =====================================================================================
#
# Assigns an Finite Element Formulation to a set of elements
#
# =====================================================================================

class Builder:

	# =================================================================================
	#
	# Creates a new Builder
	#
	# =================================================================================

	def __init__(
				self,
				FileName,
				Model = None,
				Variables = [],
				Dofs = [],
				DofsWithReactions = [],
				BufferSize = 2
				):

		self.FileName = FileName
		self.Model = Model
		self.Variables = Variables
		self.Dofs = Dofs
		self.DofsWithReactions = DofsWithReactions
		self.BufferSize = BufferSize

		# create the model part object
		self.ModelPart = ModelPart(self.FileName)

		# add variables
		for iVar in self.Variables:
			self.ModelPart.AddNodalSolutionStepVariable(iVar)

	# =================================================================================
	#
	# Builds the ModelPart
	#
	# =================================================================================

	def Build(self):

		# read the model part
		if(self.Model is None):
			model_part_io = ModelPartIO(self.FileName)
			model_part_io.ReadModelPart(self.ModelPart)
		else:
			self.Model.BuildModelPart(self.ModelPart)

		# add degrees of freedom
		for idof in self.Dofs:
			for inode in self.ModelPart.Nodes:
				inode.AddDof(idof)

		# add degrees of freedom with reactions
		for idof in self.DofsWithReactions:
			for inode in self.ModelPart.Nodes:
				inode.AddDof(idof[0], idof[1])

		# set buffer size here
		self.ModelPart.SetBufferSize(self.BufferSize)