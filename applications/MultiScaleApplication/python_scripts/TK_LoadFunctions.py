from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.MultiScaleApplication import *

class ProlongationType:
	CONSTANT = 0
	LINEAR = 2
	ZERO = 3


class ConstantLoadFunction:
	def GetMultiplier(self,
					  Info=None,
					  ):
		return 1.0

class LambdaLoadFunction:
	def GetMultiplier(self,
					  Info=None,
					  ):
		return Info[LAMBDA]

class RampLoadFunction:
	def GetMultiplier(self,
					  Info=None,
					  ):
		start_time = Info[START_TIME]
		end_time = Info[END_TIME]
		current_time = Info[TIME]
		if(current_time <= start_time):
			return 0.0
		elif(current_time >= end_time):
			return 1.0
		else:
			return (current_time-start_time)/(end_time-start_time)


class PieceWiseLoadFunction:

	def __init__(self,
				 X=None,
				 Y=None,
				 ProlongationLeft=ProlongationType.CONSTANT,
				 ProlongationRight=ProlongationType.CONSTANT
				 ):
		self.X = X
		self.Y = Y
		self.ProlongationLeft = ProlongationLeft
		self.ProlongationRight = ProlongationRight
		self.__check_input()

	def GetMultiplier(self,
					  Info=None,
					  ):
		n = len(self.X)
		currentTime = Info[TIME]
		if(currentTime < self.X[0]):
			if(self.ProlongationLeft == ProlongationType.ZERO):
				return 0.0
			elif(self.ProlongationLeft == ProlongationType.CONSTANT):
				return self.Y[0]
			elif(self.ProlongationLeft == ProlongationType.LINEAR):
				return self.__get_linear_factor_left(currentTime)
		elif(currentTime > (self.X[n-1] + self.Tolerance)):
			if(self.ProlongationRight == ProlongationType.ZERO):
				return 0.0
			elif(self.ProlongationRight == ProlongationType.CONSTANT):
				return self.Y[n-1]
			elif(self.ProlongationRight == ProlongationType.LINEAR):
				return self.__get_linear_factor_right(currentTime)
		else :
			return self.__get_linear_factor(currentTime)

	def __check_input(self):
		n = len(self.X)
		if(n != len(self.Y)):
			raise Exception("The size of X and Y must be the same")
		minStep = self.X[n-1] - self.X[0]
		for i in range(1,n):
			x0 = self.X[i-1]
			x1 = self.X[i]
			if(x1 <= x0):
				raise Exception("The X values should increase monotonically")
			istep = x1-x0
			if(minStep > istep):
				minStep = istep
		self.Tolerance = minStep*1.0E-6

	def __get_linear_factor(self, currentTime):
		for i in range(1,len(self.X)):
			if(currentTime <= self.X[i] + self.Tolerance):
				x0 = self.X[i-1]
				x1 = self.X[i]
				y0 = self.Y[i-1]
				y1 = self.Y[i]
				deltaTime = x1 - x0
				if(deltaTime == 0.0):
					return y1
				deltaFactor = y1 - y0
				relativeTime = currentTime - x0
				timeRatio = relativeTime / deltaTime
				return y0 + timeRatio * deltaFactor
		return 0.0;

	def __get_linear_factor_left(self, currentTime):
		x0 = self.X[0]
		x1 = self.X[1]
		y0 = self.Y[0]
		y1 = self.Y[1]
		deltaTime = x1 - x0
		if(deltaTime == 0.0):
			return 0.0;
		deltaFactor = y1 - y0
		relativeTime = currentTime - x0
		timeRatio = relativeTime / deltaTime
		return y0 + timeRatio * deltaFactor

	def __get_linear_factor_right(self, currentTime):
		i = len(self.X)
		x0 = self.X[i-1]
		x1 = self.X[i]
		y0 = self.Y[i-1]
		y1 = self.Y[i]
		deltaTime = x1 - x0;
		if(deltaTime == 0.0):
			return 0.0;
		deltaFactor = y1 - y0
		relativeTime = currentTime - x0
		timeRatio = relativeTime / deltaTime
		return y0 + timeRatio * deltaFactor

