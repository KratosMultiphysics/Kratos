from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.MultiScaleApplication import *

def Property(
			Pro,
			Values = []
			):

	for iValue in Values:
		try:
			Pro.SetValue(iValue[0], iValue[1])
		except:
			iValue[1].AddToProperties(iValue[0], iValue[1], Pro)

def ShellHomogeneousSection(
							Thickness,
							Material,
							Pro,
							Offset = 0.0,
							NumberOfIntegrationPoints = 5
							):

	sec = ShellCrossSection()
	sec.SetOffset(Offset)
	sec.BeginStack()
	sec.AddPly(Thickness, 0.0, NumberOfIntegrationPoints, Material)
	sec.EndStack()
	sec.AddToProperties(SHELL_CROSS_SECTION, sec, Pro)

class Ply:
	def __init__(
				self,
				Thickness,
				Material,
				NumberOfIntegrationPoints = 5,
				Orientation = 0.0):

		self.Thickness = Thickness
		self.Material = Material
		self.NumberOfIntegrationPoints = NumberOfIntegrationPoints
		self.Orientation = Orientation

def ShellCompositeSection(
						  Pro,
						  Offset = 0.0,
						  Layup = []
						  ):

	sec = ShellCrossSection()
	sec.SetOffset(Offset)
	sec.BeginStack()
	for iPly in Layup:
		sec.AddPly(
			iPly.Thickness,
			iPly.Orientation,
			iPly.NumberOfIntegrationPoints,
			iPly.Material
			)
	sec.EndStack()
	sec.AddToProperties(SHELL_CROSS_SECTION, sec, Pro)