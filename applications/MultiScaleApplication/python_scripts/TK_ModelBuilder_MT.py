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
				FileName_Mechanical,
				FileName_Thermal,
				Variables = [],
				Dofs = [],
				DofsWithReactions = [],
				BufferSize = 2
				):

		self.FileName_Mechanical = FileName_Mechanical
		self.FileName_Thermal = FileName_Thermal
		self.Variables = Variables
		self.Dofs = Dofs
		self.DofsWithReactions = DofsWithReactions
		self.BufferSize = BufferSize

		# create the model part object
		self.ModelPart_Mechanical = ModelPart(self.FileName_Mechanical)
		self.ModelPart_Thermal = ModelPart(self.FileName_Thermal)

		# add variables
		for iVar in self.Variables:
			self.ModelPart_Mechanical.AddNodalSolutionStepVariable(iVar)

	# =================================================================================
	#
	# Builds the ModelPart
	#
	# =================================================================================

	def Build(self):

		# read the model part
		model_part_io = ModelPartIO(self.FileName_Mechanical)
		model_part_io.ReadModelPart(self.ModelPart_Mechanical)

		# add degrees of freedom
		for idof in self.Dofs:
			for inode in self.ModelPart_Mechanical.Nodes:
				inode.AddDof(idof)

		# add degrees of freedom with reactions
		for idof in self.DofsWithReactions:
			for inode in self.ModelPart_Mechanical.Nodes:
				inode.AddDof(idof[0], idof[1])

		# set buffer size here
		self.ModelPart_Mechanical.SetBufferSize(self.BufferSize)
		self.ModelPart_Thermal.SetBufferSize(self.BufferSize)

		model_part_io_conn_preserv = ModelPartIOConnPreserver(self.FileName_Thermal, self.ModelPart_Mechanical)
		model_part_io_conn_preserv.ReadModelPart(self.ModelPart_Thermal)