## @package BC
#  This module contains classes to handle
#  the definition, application and update of boundary conditions
#
#  More details

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.MultiScaleApplication import *
import TK_LoadFunctions
CheckForPreviousImport()

## PrescribedValueBC.
#
#  A class that manages the nodal variables for a group of nodes in a model
#  using a user-defines (or a default) load function
class PrescribedValueBC:
	
	def __init__(
				 self,
				 Model=None,
				 Nodes=None,
				 Values=None,
				 LoadFunction=None,
				 Fix = True):
		self.Model = Model
		self.Nodes = Nodes
		self.Values = Values
		self.LoadFunction = LoadFunction
		if(self.LoadFunction is None):
			self.LoadFunction = TK_LoadFunctions.ConstantLoadFunction()
		self.Fix = Fix
	
	def Initialize(self, Model):
		info = self.Model.ProcessInfo
		for ipair in self.Values:
			mult = self.LoadFunction.GetMultiplier(info)
			var = ipair[0]
			val = ipair[1] * mult
			for i in self.Nodes:
				inode = self.Model.Nodes[i]
				if(self.Fix):
					inode.Fix(var)
				inode.SetSolutionStepValue(var, val)
	
	def OnBeforeSolutionStep(self, Model):
		info = self.Model.ProcessInfo
		for ipair in self.Values:
			mult = self.LoadFunction.GetMultiplier(info)
			var = ipair[0]
			val = ipair[1] * mult
			for i in self.Nodes:
				inode = self.Model.Nodes[i]
				inode.SetSolutionStepValue(var, val)
	
	def OnSolutionStepCompleted(self, Model):
		pass
	
	def Finalize(self, Model):
		pass
	