## @package BC
#  This module contains classes to handle
#  the definition, application and update of boundary conditions
#
#  More details

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.MultiScaleApplication import *
import TK_LoadFunctions

def dummy_spatial_function(x,y,z):
	return 1.0

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
				 SpatialFunction=None,
				 Fix = True,
				 TimeBegin = 0.0,
				 TimeEnd = 1.0):
		self.Model = Model
		self.Nodes = Nodes
		self.Values = Values
		self.LoadFunction = LoadFunction
		if(self.LoadFunction is None):
			self.LoadFunction = TK_LoadFunctions.ConstantLoadFunction()
		self.SpatialFunction = SpatialFunction
		if(self.SpatialFunction is None):
			self.SpatialFunction = dummy_spatial_function
		self.Fix = Fix
		self.TimeBegin = TimeBegin
		self.TimeEnd = TimeEnd
		self.Initialized = False
		self.InitialValues = [[0 for i in range(0,len(self.Nodes))] for i in range(0,len(self.Values))]

	def Initialize(self, Model):
		pass

	def OnBeforeSolutionStage(self, Model):
		info = self.Model.ProcessInfo
		time = info[TIME]
		if(time >= self.TimeBegin):
			if(time <= self.TimeEnd):
				if(not self.Initialized):
					id_1=0
					for ipair in self.Values:
						id_2=0
						mult = self.LoadFunction.GetMultiplier(info)
						var = ipair[0]
						val = ipair[1] * mult
						for i in self.Nodes:
							inode = self.Model.Nodes[i]
							spatial_mult = self.SpatialFunction(inode.X0, inode.Y0, inode.Z0)
							initial_value = inode.GetSolutionStepValue(var)
							self.InitialValues[id_1][id_2]=initial_value
							inode.SetSolutionStepValue(var,0, val*spatial_mult+initial_value)
							if(self.Fix):
								inode.Fix(var)
							id_2 += 1
						id_1 += 1
					self.Initialized = True

	def OnSolutionStageCompleted(self, Model):
		pass

	def OnBeforeSolutionStep(self, Model):
		info = self.Model.ProcessInfo
		time = info[TIME]
		if(time > self.TimeBegin):
			if(time <= self.TimeEnd):
				if(not self.Initialized):
					self.OnBeforeSolutionStage(Model)
				if(self.Initialized):
					id_1=0
					for ipair in self.Values:
						mult = self.LoadFunction.GetMultiplier(info)
						id_2=0
						var = ipair[0]
						val = ipair[1] * mult
						for i in self.Nodes:
							inode = self.Model.Nodes[i]
							spatial_mult = self.SpatialFunction(inode.X0, inode.Y0, inode.Z0)
							initial_value = self.InitialValues[id_1][id_2]
							inode.SetSolutionStepValue(var,0, val*spatial_mult+initial_value)
							if(self.Fix):
								inode.Fix(var)
							id_2 += 1
						id_1 += 1
			else:
				for ipair in self.Values:
					var = ipair[0]
					val = 0.0
					for i in self.Nodes:
						inode = self.Model.Nodes[i]
						inode.SetSolutionStepValue(var,0, val)
						if(self.Fix):
							inode.Free(var)

	def OnSolutionStepCompleted(self, Model):
		self.OnBeforeSolutionStep(Model) # call it once again, just to make sure! it's useful with arc length strategy

	def Finalize(self, Model):
		pass
