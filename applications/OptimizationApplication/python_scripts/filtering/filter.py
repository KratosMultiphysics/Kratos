import typing
from abc import ABC, abstractmethod
from importlib import import_module

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class Filter(ABC):
    """Base filter class

    This class unifies the filters used in the optimization application.
    """
    def __init__(self) -> None:
        self.__component_data_view: 'typing.Optional[ComponentDataView]' = None

    def SetComponentDataView(self, component_data_view: ComponentDataView) -> None:
        self.__component_data_view = component_data_view

    def GetComponentDataView(self) -> ComponentDataView:
        if self.__component_data_view is None:
            raise RuntimeError("Please use SetComponentDataView first.")
        return self.__component_data_view

    @abstractmethod
    def Initialize(self) -> None:
        """Initializes the filter.

        This method initializes the filter. This is only called once in the whole optimization process.

        """
        pass

    @abstractmethod
    def Check(self) -> None:
        """Checks the filter.

        This method checks the filter. This is only called once in the whole optimization process after calling the initialize.

        """
        pass

    @abstractmethod
    def Finalize(self) -> None:
        """Finalizes the filter.

        This method checks the filter. This is only called once in the whole optimization process at the end.

        """
        pass

    @abstractmethod
    def FilterField(self, unfiltered_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        """Filter the input unfiltered_field.

        This method filters the passed input field. Unfiltered field is assumed to be a
        non-integrated field over entity domains.

        Args:
            unfiltered_field (ContainerExpressionTypes): Input field

        Returns:
            ContainerExpressionTypes: Filtered output field.
        """
        pass

    @abstractmethod
    def FilterIntegratedField(self, unfiltered_integrated_field: ContainerExpressionTypes) -> ContainerExpressionTypes:
        """Filter the input integrated field.

        This method filters the passed integrated input field. The unfiltered integrated field is
        assumed to be already integrated over entity domains.

        Args:
            unfiltered_integrated_field (ContainerExpressionTypes): Unfiltered integrated field.

        Returns:
            ContainerExpressionTypes: Filtered output field.
        """
        pass

    @abstractmethod
    def Update(self) -> None:
        """Updates the filter.

        This method is used to update the filter state. This will be called by the control, when control changes
        the optimization domain.

        """
        pass

def Factory(model: Kratos.Model, filtering_model_part_name: str, variable: SupportedSensitivityFieldVariableTypes, data_location: Kratos.Globals.DataLocation, settings: Kratos.Parameters) -> Filter:
    if not settings.Has("filter_type"):
        raise RuntimeError(f"\"filter_type\" not provided in the following filter settings:\n{settings}")

    filter_type = settings["filter_type"].GetString()
    filter_module_name = f"KratosMultiphysics.OptimizationApplication.filtering.{filter_type}"

    module = import_module(filter_module_name)
    if not hasattr(module, "Factory"):
        raise RuntimeError(f"Python module {filter_module_name} does not have a Factory method.")
    return getattr(module, "Factory")(model, filtering_model_part_name, variable, data_location, settings)
