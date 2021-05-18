import KratosMultiphysics
from KratosMultiphysics.check_scalar_base_process import CheckScalarBaseProcess

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if not isinstance(settings,KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckScalarFromProcessInfoProcess(Model, settings["Parameters"])


class CheckScalarFromProcessInfoProcess(CheckScalarBaseProcess, KratosUnittest.TestCase):
    """This process checks analytically from a function the solution (scalar) in the process info belonging a certain submodelpart

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        # Copy from input
        self.settings = settings

        # Construct the base process.
        super(CheckScalarFromProcessInfoProcess, self).__init__(Model, self.settings)

        # Raise error in case of not time function
        if not self.is_time_function and not self.value_is_numeric:
            raise Exception("Must depend on time alone!: " + self.function_string)

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        process_info = self.model_part.ProcessInfo
        current_time = process_info[KratosMultiphysics.TIME] - process_info[KratosMultiphysics.DELTA_TIME]

        if current_time >= self.interval[0] and current_time < self.interval[1]:
            if self.value_is_numeric:
                value = process_info.GetValue(self.variable)
                self.assertAlmostEqual(self.value, value, self.tolerance_rank)
            else:
                self.value = self.aux_function.f(0.0, 0.0, 0.0, current_time)
                value = process_info.GetValue(self.variable)
                self.assertAlmostEqual(self.value, value, self.tolerance_rank)
