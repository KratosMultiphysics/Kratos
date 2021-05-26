# This test example is from M. Hojjat, E. Stavropoulou,
# K.-U. Bletzinger, The Vertex Morphing method for node-based
# shape optimization, Comput. Methods Appl. Mech. Engrg. 268
# (2014) 494-513.
#
# The target curve is the tent function illustrated below.
#
#                    z=1
#                     /\
#                    /  \
#                   /    \
#  |--> x          /      \           z=0
#  _______________/        \_______________
#  |----- 15 -----|-- 10 --|----- 15 -----|

from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_base import AnalyzerBaseClass

class CustomAnalyzer(AnalyzerBaseClass):

    # --------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        if communicator.isRequestingValueOf("targetDeviation"):
            communicator.reportValue("targetDeviation", self.__ObjectiveFunction(current_design))

        if communicator.isRequestingGradientOf("targetDeviation"):
            communicator.reportGradient("targetDeviation", self.__ObjectiveGradient(current_design))

    # --------------------------------------------------------------------------------------------------
    def __ObjectiveFunction(self, current_design):
        """ Returns the objective function to be minimized """
        objective = 0.0
        for node in current_design.Nodes:
            objective = objective + abs(self.__TentFunction(node.X) - node.Z)
        return objective

    # --------------------------------------------------------------------------------------------------
    def __ObjectiveGradient(self, current_design):
        """ Returns the gradient of the objective function """
        sensitivity = dict()
        for node in current_design.Nodes:
            delta = node.Z - self.__TentFunction(node.X)
            if abs(delta) == 0.0:
                sz = 0.0
            else:
                sz = delta / abs(delta)
            sensitivity[node.Id] = [0.0, 0.0, sz]
        return sensitivity

    # --------------------------------------------------------------------------------------------------
    @staticmethod
    def __TentFunction(x):
        """ Defines the target curve z=__TentFunction(x) """
        if x <= 15.0:
            return 0.0
        elif x<= 20.0:
            return (x - 15.0) / 5.0
        elif x <= 25.0:
            return 1.0 - (x - 20.0) / 5.0
        else:
            return 0.0