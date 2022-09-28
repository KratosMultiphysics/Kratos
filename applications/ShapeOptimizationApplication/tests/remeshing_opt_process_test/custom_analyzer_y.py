from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_base import AnalyzerBaseClass

class CustomAnalyzer(AnalyzerBaseClass):
    """
    Response is square sum of all y-coordinates.
    If minimized, the geometry is flatteened in y direction.
    """

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
            objective = objective + node.Y * node.Y
        return objective

    # --------------------------------------------------------------------------------------------------
    def __ObjectiveGradient(self, current_design):
        """ Returns the gradient of the objective function """
        sensitivity = dict()
        for node in current_design.Nodes:
            sensitivity[node.Id] = [0.0, 2*node.Y, 0.0]
        return sensitivity
