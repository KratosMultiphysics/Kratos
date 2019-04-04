"""simple helper class to generate meaningful gradient fields"""
import math

class SumOfNodalDistancesResponse():
    def __init__(self, value=1.0):
        self.value = value

    def GetGradient(self, model_part, variable):
        for node in model_part.Nodes:
            _gradient = [0.0,0.0,0.0]
            _gradient[2] = math.copysign(1, node.Z - self.value)
            node.SetSolutionStepValue(variable, _gradient)

class SumOfNodalDistancesLinearResponse():

    def GetGradient(self, model_part, variable):
        for node in model_part.Nodes:
            _gradient = [0.0,0.0,0.0]
            _gradient[2] = node.Z - node.X
            node.SetSolutionStepValue(variable, _gradient)

class DiscreteVolumeResponse():
    def __init__(self, value=1.0):
        self.value = value

    def GetGradient(self, model_part, variable):
        for node in model_part.Nodes:
            node.SetSolutionStepValue(variable, [0.0, 0.0, 0.0])
        for condition in model_part.Conditions:
            for node in condition.GetNodes():
                area = condition.GetGeometry().Area()
                _gradient = [0.0,0.0,0.0]
                _gradient[2] = math.copysign(1, node.Z - self.value) * area
                orig = node.GetSolutionStepValue(variable)
                _gradient[0] += orig[0]
                _gradient[1] += orig[1]
                _gradient[2] += orig[2]
                node.SetSolutionStepValue(variable, _gradient)


class DiscreteLinearResponse():

    def GetGradient(self, model_part, variable):
        for node in model_part.Nodes:
            node.SetSolutionStepValue(variable, [0.0, 0.0, 0.0])
        for condition in model_part.Conditions:
            for node in condition.GetNodes():
                area = condition.GetGeometry().Area()
                _gradient = [0.0,0.0,0.0]
                _gradient[2] = (node.Z - node.X) * area
                orig = node.GetSolutionStepValue(variable)
                _gradient[0] += orig[0]
                _gradient[1] += orig[1]
                _gradient[2] += orig[2]
                node.SetSolutionStepValue(variable, _gradient)


