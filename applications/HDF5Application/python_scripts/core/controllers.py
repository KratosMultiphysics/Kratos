"""!@package HDF5Application

HDF5 controllers.

This module contains the controllers which control the frequency of HDF5 IO
operations.

license: HDF5Application/license.txt
"""


# Kratos imports
import KratosMultiphysics

##!@addtogroup HDF5Application
##!@{
##!@name Kratos classes
##!@{

class DefaultController(KratosMultiphysics.Controller):
    def __init__(self) -> None:
        KratosMultiphysics.Controller.__init__(self)

    def Evaluate(self) -> bool:
        return True

class SingleTimeController(KratosMultiphysics.Controller):
    def __init__(self, temporal_controller: KratosMultiphysics.OutputController) -> None:
        KratosMultiphysics.Controller.__init__(self)
        self.__temporal_controller = temporal_controller
        self.__is_evaluated = temporal_controller.GetNextPossibleEvaluateControlValue() > temporal_controller.GetInterval()

    def Evaluate(self) -> bool:
        return not self.__is_evaluated and self.__temporal_controller.Evaluate()

    def Update(self) -> None:
        if self.__temporal_controller.Evaluate():
            self.__is_evaluated = True
        self.__temporal_controller.Update()

def Factory(model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> KratosMultiphysics.Controller:
    """!@brief Return the controller specified by the setting 'controller_type'.
    @detail Empty settings will contain default values after returning from the
    function call.
    """
    parameters.AddMissingParameters(KratosMultiphysics.Parameters("""{
        "controller_type" : "default_controller"
    }"""))
    controller_type = parameters['controller_type'].GetString()
    if controller_type == 'default_controller':
        return DefaultController()
    elif controller_type == 'temporal_controller':
        return KratosMultiphysics.OutputController(model, parameters)
    elif controller_type == "single_time_controller":
        return SingleTimeController(KratosMultiphysics.OutputController(model, parameters))
    else:
        raise ValueError(f"Unsupported controller_type = \"{controller_type}\". Followings are supported controller types:\n\ttemporal_controller\n\tdefault_controller\n\tsingle_time_controller")
##!@}
