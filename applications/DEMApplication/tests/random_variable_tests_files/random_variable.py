import numpy as np
import random
import KratosMultiphysics


def GaussianPDF(x, mu=0, sigma=1):
	return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))

class RandomVariable:
	def __init__(self, parameters):
		pass

	def Sample(self):
		return 0.0

class PiecewiseLinearRV(RandomVariable):
	def __init__(self, parameters):
		if not isinstance(parameters, KratosMultiphysics.Parameters):
			raise Exception("expected input shall be a Parameters object, encapsulating a json string")
		self.boundaries = np.array(parameters['breakpoints'].GetVector(), dtype=np.float32)
		self.heights = np.array(parameters['values'].GetVector(), dtype=np.float32)
		self.ranges = np.array([(self.boundaries[i+1] - self.boundaries[i] for i in range(len(self.boundaries)))])
		self.areas = np.zeros(len(self.heights) - 1)
		self.trapezoid_indices = range(len(self.areas))
		self.Normalize()
		self.Sample()
		super().__init__(parameters)

	def InterpolateHeight(self, x):
		if x < self.boundaries[0] or x > self.boundaries[-1]:
			raise Exception('x =' + str(x) + 'is out of bounds (' + str([self.boundaries[0], self.boundaries[-1]]) + ')')

		for i, x_boundary in enumerate(self.boundaries):
			if x < x_boundary:
				break

		b1 = self.heights[i - 1]
		b2 = self.heights[i]
		h = self.boundaries[i] - self.boundaries[i - 1]
		x1 = self.boundaries[i - 1]
		x2 = self.boundaries[i]
		return (b1 * (x2 - x) + b2 * (x - x1)) / h

	def Normalize(self):
		self.total_area = 0.0
		for i in range(len(self.heights[:-1])):
			area = 0.5 * (self.boundaries[i+1] - self.boundaries[i]) * (self.heights[i+1] + self.heights[i])
			self.areas[i] = area
			self.total_area += area
		self.heights /= self.total_area
		self.discrete_probabilities = self.areas / self.total_area

	def Sample(self):
		i_trapezoid = self.SampleTrapezoidChoice()[0]
		x0 = self.boundaries[i_trapezoid]
		H = self.boundaries[i_trapezoid + 1] - x0
		B1 = self.heights[i_trapezoid]
		B2 = self.heights[i_trapezoid + 1]
		x_within = self.__class__.SampleWithinTrapezoid(H, B1, B2)
		return x0 + x_within

	def SampleTrapezoidChoice(self, n_to_pick=1):
		draw = np.random.choice(len(self.areas), n_to_pick, p=self.discrete_probabilities)
		return draw

	@staticmethod
	def SampleWithinTrapezoid(H, B1, B2):
		if B1 == 0:
			x = PiecewiseLinearRV.SamplePositiveSlopingStandardTriangle()
		else:
			beta = B2/B1
			b = 2.0 / (1 + beta)
			x = PiecewiseLinearRV.SampleXWithinStandardTrapezoid(b)
		return H * x

	@staticmethod
	def SampleXWithinStandardTrapezoid(b):
		alpha = random.uniform(0, 1)
		if alpha < b/2:
			y = PiecewiseLinearRV.SampleNegativeSlopingStandardTriangle()
		else:
			y = PiecewiseLinearRV.SamplePositiveSlopingStandardTriangle()
		return y

	@staticmethod
	def SamplePositiveSlopingStandardTriangle():
		return np.sqrt(random.uniform(0, 1))

	@staticmethod
	def SampleNegativeSlopingStandardTriangle():
		return 1.0 - PiecewiseLinearRV.SamplePositiveSlopingStandardTriangle()

