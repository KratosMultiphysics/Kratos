import typing
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll

class MasterControl:
    """Master control class.

    This class is used to simplify working with many controls at once. Following responsibilities are assumed:
        1. Maps physical gradients from different domains to one @ref Kratos::CombinedTensorAdaptor (using MapGradient).
        2. Updates each respective domain from updates given by one @ref Kratos::CombinedTensorAdaptor (using Update).

    There should be only one master control class per optimization problem.

    """
    def __init__(self) -> None:
        self.__list_of_controls: 'list[Control]' = []

    def AddControl(self, control: Control) -> None:
        """Adds a given control to the master control.

        Args:
            control (Control): Control to be added
        """
        self.__list_of_controls.append(control)

    def GetListOfControls(self) -> 'list[Control]':
        """Returns the list of controls in the master control.

        Returns:
            list[Control]: List of controls.
        """
        return self.__list_of_controls

    def GetPhysicalKratosVariableCombinedTensorAdaptorsMap(self) -> 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]':
        """Returns map of physical variables and @ref Kratos::CombinedTensorAdaptor from each control.

        This returns a map of physical control variables and a @ref Kratos::CombinedTensorAdaptor. The @ref Kratos::CombinedTensorAdaptor will contain
        all the @ref Kratos::TensorAdaptor for respective control's control domains for each physical control variable. The repeated @ref Kratos::TensorAdaptor s
        for the same physical control variable is omitted to avoid double calculations.

        Returns:
            dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]: Physical control variable and @ref Kratos::CombinedTensorAdaptor map.
        """

        physical_variable_containers_dict: 'dict[SupportedSensitivityFieldVariableTypes, list[typing.Union[Kratos.NodesArray, Kratos.ElementsArray, Kratos.ConditionsArray]]]' = {}
        for control in self.__list_of_controls:
            for physical_variable in control.GetPhysicalKratosVariables():
                # check whether the physical variable is already there.
                if not physical_variable in physical_variable_containers_dict.keys():
                    physical_variable_containers_dict[physical_variable] = []

                # check whether the container for that physical variable is already there.
                control_ta = control.GetEmptyField()
                if control_ta.GetContainer() not in physical_variable_containers_dict[physical_variable]:
                    physical_variable_containers_dict[physical_variable].append(control_ta.GetContainer())

        physical_variable_combined_ta_dict: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]' = {}
        for variable, containers_list in physical_variable_containers_dict.items():
            physical_variable_combined_ta_dict[variable] = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([Kratos.TensorAdaptors.DoubleTensorAdaptor(Kratos.TensorAdaptors.VariableTensorAdaptor(container, variable), copy=False) for container in containers_list], perform_collect_data_recursively=False, perform_store_data_recursively=False)

        return physical_variable_combined_ta_dict

    def GetEmptyField(self) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        """Returns empty CollectiveExpression containing empty ContainerExpressions for each control.

        Returns:
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor: Empty CollectiveExpression
        """
        ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([control.GetEmptyField() for control in self.__list_of_controls], perform_collect_data_recursively=False, perform_store_data_recursively=False)
        ta.CollectData() # reading all the data from each ta adaptor to the combined tensor adaptor
        return ta

    def GetControlField(self) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        """Returns CollectiveExpression containing control field ContainerExpressions for each control.

        Returns:
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor: Control field CollectiveExpression
        """
        ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([control.GetControlField() for control in self.__list_of_controls], perform_collect_data_recursively=False, perform_store_data_recursively=False)
        ta.CollectData() # reading all the data from each ta adaptor to the combined tensor adaptor
        return ta

    @time_decorator()
    def MapGradient(self, physical_space_gradient_variable_and_combined_tensor_adaptor_map: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        """Maps physical space gradients to a combined tensor adaptor.

        This method maps sensitivities w.r.t. physical space variables to control space by using each control. It is done by converting input
        map with @ref Kratos::CombinedTensorAdaptor to a map with @ref Kratos::TensorAdaptor. The following is a pseudo example:

        input:

        physical_space_gradient_variable_and_combined_tensor_adaptor_map = {
            "YOUND_MODULUS": [ControlDomain1, ControlDomain2],
            "DENSITY"      : [ControlDomain1],
            "VISCOSITY"    : [ControlDomain2],
        }

        control info:
        control1: domain = ControlDomain1
                  physical_vars = YOUND_MODULUS, DENSITY

        control2: domain = ControlDomain2
                  physical_vars = VISCOSITY, DENSITY


        from above information, following two maps are created and passed to each control to obtain one @ref Kratos::TensorAdaptor from each control.
        control_specific_maps:
        for control1: {YOUND_MODULUS: ControlDomain1, DENSITY: ControlDomain1}
        for control2: {VISCOSITY: ControlDomain2, DENSITY: A zero valued ControlDomain}

        then returned mapped gradients are added to one @ref Kratos::CombinedTensorAdaptor.

        In here, the missing gradients for required physical variables will be assumed to be zero.

        This converts given map of sensitivities w.r.t. different physical space variables to one sensitivities in control space and aggregated
        to one @ref Kratos::CombinedTensorAdaptor.

        Args:
            physical_space_gradient_variable_and_combined_tensor_adaptor_map (dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]): Map of physical space sensitivities.

        Returns:
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor: Control space sensitivities.
        """


        mapped_gradients: 'list[Kratos.TensorAdaptors.DoubleTensorAdaptor]' = []

        for control in self.__list_of_controls:
            # iterate through each control to create its own container map from the combined tensor adaptor map given as input
            control_physical_sensitivities_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleTensorAdaptor]' = {}
            for physical_control_variable in control.GetPhysicalKratosVariables():
                # first assume the gradients for this physical_control_variable is zero, hence get the zero valued tensor adaptor.
                control_ta = control.GetEmptyField()

                # get the required physical variables from control.
                if physical_control_variable in physical_space_gradient_variable_and_combined_tensor_adaptor_map.keys():
                    # if the sensitivities for the given physical control variable exists, then try to find whether
                    # for this specific control the physical control variable sensitivities exists.

                    sensitivity_cta = physical_space_gradient_variable_and_combined_tensor_adaptor_map[physical_control_variable]
                    for ta in sensitivity_cta.GetTensorAdaptors():
                        if ta.GetContainer() == control_ta.GetContainer():
                            # there exists for this control's physical variables sensitivities.
                            control_physical_sensitivities_container_expression_map[physical_control_variable] = ta
                            break

                if physical_control_variable not in control_physical_sensitivities_container_expression_map.keys():
                    # If it is found from input gradients, the control_expression will have those values,
                    # otherwise it will have representative zero control_expression.
                    control_physical_sensitivities_container_expression_map[physical_control_variable] = control_ta

            # map the physical control variable sensitivities to one control space
            mapped_gradients.append(control.MapGradient(control_physical_sensitivities_container_expression_map))

        result = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([mapped_gradient for mapped_gradient in mapped_gradients], perform_collect_data_recursively=False, perform_store_data_recursively=False)
        result.CollectData()
        return result

    @time_decorator()
    def Update(self, update_combined_tensor_adaptor: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor) -> 'dict[Control, bool]':
        """Update each control with given combined tensor adaptor's respective tensor adaptor.

        Args:
            update_combined_tensor_adaptor (Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor): Update

        Raises:
            RuntimeError: If number of controls and number of tensor adaptors mismatch.

        Returns:
            dict[Control, bool]: A map with control and a boolean whether the update changed anything in that control.
        """
        if len(self.__list_of_controls) != len(update_combined_tensor_adaptor.GetTensorAdaptors()):
            raise RuntimeError(f"Controls size and update size mismatch [ number of controls: {len(self.__list_of_controls)}, number of tensor adaptors: {len(update_combined_tensor_adaptor.GetTensorAdaptors())} ].")

        update_map: 'dict[Control, bool]' = {}
        for control, ta in zip(self.__list_of_controls, update_combined_tensor_adaptor.GetTensorAdaptors()):
            update_map[control] = control.Update(ta)
            if update_map[control]:
                Kratos.Logger.PrintInfo(f"Control {control.GetName()} updated.")
            else:
                Kratos.Logger.PrintInfo(f"Control {control.GetName()} not updated")

        return update_map

    def Check(self) -> None:
        CallOnAll(self.__list_of_controls, Control.Check)

    def Initialize(self) -> None:
        CallOnAll(self.__list_of_controls, Control.Initialize)

    def Finalize(self) -> None:
        CallOnAll(self.__list_of_controls, Control.Finalize)

    def GetName(self) -> str:
        return "master_control"
