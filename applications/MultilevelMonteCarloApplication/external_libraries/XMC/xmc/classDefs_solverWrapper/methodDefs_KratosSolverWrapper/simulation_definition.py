# Import generic Python libraries
# Type hinting
from typing import Union, List, Optional

# Import Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage

# Define type hints
simulation_output_type = Optional[List[Union[float, List[Union[float, List[float]]]]]]
event_type = Union[List[float], List[List[float]]]


class SimulationScenario(AnalysisStage):
    """
    This class is an example and cannot be used to run any problem.

    Class defining the analysis stage, namely the problem, XMC is solving.
    The user should derive the customized analysis stage from the base class AnalysisStage
    or from other analyses of Kratos.
    """

    def __init__(self, model, parameters, random_variable: event_type):
        """
        The constructor of the AnalysisStage-Object.

        Inputs:

        model: KratosMultiphysics.Model object.
            The Model to be used
        parameters: KratosMultiphysics.Parameters object.
            The ProjectParameters to be used
        random_variable: list.
            The unknown parameter of the problem
        """

        raise RuntimeError(
            'Example "SimulationScenario" is being called. '
            'Instead, a specific "SimulationScenario" must be defined and '
            'the key "analysisStage" with the path to it added to the '
            '"KratosSolverWrapper" dictionary.'
        )

    def EvaluateQuantityOfInterest(self) -> simulation_output_type:
        """
        Method evaluating the output quantity/quantities of interest of the problem.
        The output list structure should agree with KratosSolverWrapper and moment estimators
        formats.
        This method is called once per every sample, at the end of the Kratos simulation.

        Inputs:

        Outputs:

        qoi_list: list.
            It contains scalars or list of scalars
        """

    def MappingAndEvaluateQuantityOfInterest(self) -> simulation_output_type:
        """
        Method mapping output quantities of interest onto reference model parts
        and evaluating and returning the output quantity/quantities of interest of the problem.
        The output list structure should agree with KratosSolverWrapper and moment estimators
        formats.
        This method is called once per every sample, at the end of the Kratos simulation.
        This method is called if mappingOutputQuantities is set to True,
        and internally can call the EvaluateQuantityOfInterest method.

        Inputs:

        Outputs:

        qoi_list: list.
            It contains scalars or list of scalars
        """
