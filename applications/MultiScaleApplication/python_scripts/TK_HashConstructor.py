## @package BC
#  This module contains classes to handle
#  the definition, application and update of boundary conditions
#
#  More details

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.MultiScaleApplication import *
import TK_LoadFunctions

from math import *

## PrescribedValueBC.
#
#  A class that manages the nodal variables for a group of nodes in a model
#  using a user-defines (or a default) load function
class HashDatabaseConstructor:

	def __init__(self,
				 NDivision=None,
				 Tags=None,
				 Matrix=None,
				 Vector=None,
				 Double=None,
				 ):
		self.Tags = Tags
		self.Matrix = Matrix
		self.Vector = Vector
		self.Double = Double

	def DefaultConstructor(self,
						   NDivision=None):
		self.NDivision = NDivision

		Tags = []; #zeros((m+1)^2,2);
		for i in range(-m,m+1):
			for j in range(-m,m+1):
				Tags.append([i,j]);
		return Tags

	def AddMatrixToHash(self,
							Tags=None,
							Matrix=None):
		self.Tags = Tags
		self.Matrix = Matrix

		HashMap = {}
		HashMap[Tags] = self.Matrix

		return HashMap

	def AddVectorToHash(self,
							Tags=None,
							Vector=None):
		self.Tags = Tags
		self.Vector = Vector

		HashMap = {}
		HashMap[Tags] = self.Vector

		return HashMap

	def AddDoubleToHash(self,
							Tags=None,
							Double=None):
		self.Tags = Tags
		self.Double = Double

		HashMap = {}
		HashMap[Tags] = self.Double

		return HashMap

	def GetHash(self,
				):

		pass

