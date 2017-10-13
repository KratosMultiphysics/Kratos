## @package TimeLines
#  This module contains classes (Rve Modelers) that are used
#  to handle the generation, assignment and tracking of 
#  RveConstitutiveLaws.
#
#  More details

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7


class FixedTimeLine:

	def __init__(
				self,
				Duration = 1.0,
				Increment = 1.0):
		
		self.InitialTime = 0.0
		self.Duration = Duration
		self.Increment = Increment
		
		self.EndTime = self.InitialTime + self.Duration
		self.CurrentTime = self.InitialTime
		
		self.Finished = False

	def SetInitialTime(self, InitialTime):
	
		self.InitialTime = InitialTime
		self.EndTime = self.InitialTime + self.Duration
		self.CurrentTime = self.InitialTime
	
	def NextTimeStep(self, LastIterationConverged = True):
		
		if(self.Finished):
			return (False, self.CurrentTime)
		else:
			next_time = self.CurrentTime + self.Increment
			if(next_time > self.EndTime):
				next_time = self.EndTime
				self.Finished = True
			self.CurrentTime = next_time
			return (True, next_time)


class NonUniformTimeLine:

	def __init__(
				self,
				DurationAndIncrements = []):
		
		self.InitialTime = 0.0
		self.DurationAndIncrements = DurationAndIncrements
		
		self.Duration = self.__get_duration()
		
		self.EndTime = self.InitialTime + self.Duration
		self.CurrentTime = self.InitialTime
		
		self.Finished = False

	def SetInitialTime(self, InitialTime):
	
		self.InitialTime = InitialTime
		self.EndTime = self.InitialTime + self.Duration
		self.CurrentTime = self.InitialTime
	
	def NextTimeStep(self, LastIterationConverged = True):
		
		if(self.Finished):
			return (False, self.CurrentTime)
		else:
			increment = 0.0
			next_stop = self.InitialTime
			for i in self.DurationAndIncrements:
				next_increment = i[1]
				next_stop += i[0]
				if(self.CurrentTime < next_stop):
					increment = next_increment
					break
			next_time = self.CurrentTime + increment
			if(next_time >= self.EndTime-1.0e-14):
				next_time = self.EndTime
				self.Finished = True
			self.CurrentTime = next_time
			return (True, next_time)
	
	def __get_duration(self):
		duration=0.0
		for i in self.DurationAndIncrements:
			duration += i[0]
		return duration

class AdaptiveTimeLine:
	
	def __init__(
				self,
				Duration = 1.0, 
				Increment = 1.0,  
				MinIncrement = 1.0, 
				MaxIncrement = 1.0):
		
		self.InitialTime = 0.0
		self.Duration = Duration
		
		self.Increment = Increment
		
		self.EndTime = self.InitialTime + self.Duration
		
		self.CurrentIncrement = self.Increment
		self.CurrentTime = self.InitialTime
		self.LastIterationConverged = True
		
		self.MinIncrement = MinIncrement
		self.MaxIncrement = MaxIncrement
		
		if(self.MinIncrement > self.MaxIncrement):
			temp = self.MinIncrement
			self.MinIncrement = self.MaxIncrement
			self.MaxIncrement = temp
		
		self.ASmallTolerance = self.Duration * 1.0E-10
		if(self.ASmallTolerance > self.MinIncrement):
			self.ASmallTolerance = self.MinIncrement * 1.0E-2
		
		self.Finished = False
	
	def SetInitialTime(self, InitialTime):
		
		self.InitialTime = InitialTime
		self.EndTime = self.InitialTime + self.Duration
		self.CurrentTime = self.InitialTime
	
	def NextTimeStep(self, LastIterationConverged = True):
		
		if(LastIterationConverged == True):
			
			self.LastIterationConverged = LastIterationConverged
			self.CurrentTime += self.CurrentIncrement
			self.CheckFinishedState()
			return (True, self.CurrentTime)
			
		else:
			
			self.LastIterationConverged = LastIterationConverged
			print ("")
			print (" Reducing Time Step due to NON CONVERGENCE: ")
			print ("   Previous Increment : ", self.CurrentIncrement)
			print ("   Current Increment  : ", self.CurrentIncrement * 0.5)
			print ("")
			self.CurrentTime -= self.CurrentIncrement
			self.CurrentIncrement *= 0.5
			if(self.CurrentIncrement < self.MinIncrement):
				print (" WARNING: The required increment is smaller than the Mininum increment!")
				print (" Current increment: ", self.CurrentIncrement, " < Min.Increment: ", self.MinIncrement)
				return (False, self.CurrentTime)
			else:
				self.CurrentTime += self.CurrentIncrement
				self.CheckFinishedState()
				return (True, self.CurrentTime)
		
		return (False, self.CurrentTime)
	
	def CheckFinishedState(self):
		
		self.Finished = False
		if(self.CurrentTime >= (self.EndTime - self.ASmallTolerance)):
			self.CurrentTime = self.EndTime
			self.Finished = True

# class CompositeTimeLine:
	
	# def __init__(
				# self,
				# TimeLines = []):
		
		# if(len(TimeLines) < 1):
			# raise Exception('At least one timeline is required')
		# self.TimeLines = TimeLines
		
		# self.InitialTime = 0.0
		# self.Duration = 0.0
		# for tl in self.TimeLines:
			# self.Duration = self.Duration + tl.Duration
		
		# self.TLID = 0
		
		# self.Increment = self.TimeLines[self.TLID].Increment
		
		# self.EndTime = self.InitialTime + self.Duration
		
		# self.CurrentIncrement = self.Increment
		# self.CurrentTime = self.InitialTime
		# self.LastIterationConverged = True
		
		# self.MinIncrement = self.TimeLines[self.TLID].MinIncrement
		# self.MaxIncrement = self.TimeLines[self.TLID].MaxIncrement
		
		# self.ASmallTolerance = self.Duration * 1.0E-10
		# if(self.ASmallTolerance > self.MinIncrement):
			# self.ASmallTolerance = self.MinIncrement * 1.0E-2
		
		# self.Finished = False
	
	# def SetInitialTime(self, InitialTime):
		
		# self.InitialTime = InitialTime
		# last_initial_time = self.InitialTime
		# for tl in self.TimeLines:
			# tl.SetInitialTime(last_initial_time)
			# last_initial_time = tl.EndTime
		# self.EndTime = self.InitialTime + self.Duration
		# self.CurrentTime = self.InitialTime
	
	# def NextTimeStep(self, LastIterationConverged = True):
		
		
		# if(LastIterationConverged == True):
			
			# self.LastIterationConverged = LastIterationConverged
			# self.CurrentTime += self.CurrentIncrement
			# self.CheckFinishedState()
			# return (True, self.CurrentTime)
			
		# else:
			
			# self.LastIterationConverged = LastIterationConverged
			# print ("")
			# print (" Reducing Time Step due to NON CONVERGENCE: ")
			# print ("   Previous Increment : ", self.CurrentIncrement)
			# print ("   Current Increment  : ", self.CurrentIncrement * 0.5)
			# print ("")
			# self.CurrentTime -= self.CurrentIncrement
			# self.CurrentIncrement *= 0.5
			# if(self.CurrentIncrement < self.MinIncrement):
				# print (" WARNING: The required increment is smaller than the Mininum increment!")
				# print (" Current increment: ", self.CurrentIncrement, " < Min.Increment: ", self.MinIncrement)
				# return (False, self.CurrentTime)
			# else:
				# self.CurrentTime += self.CurrentIncrement
				# self.CheckFinishedState()
				# return (True, self.CurrentTime)
		
		# return (False, self.CurrentTime)
	
	# def CheckFinishedState(self):
		
		# self.Finished = False
		# if(self.CurrentTime >= (self.EndTime - self.ASmallTolerance)):
			# self.CurrentTime = self.EndTime
			# self.Finished = True