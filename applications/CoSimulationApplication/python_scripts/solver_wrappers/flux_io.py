from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

def Create(*args):
    """
    Create a FluxIO instance.

    Returns:
        FluxIO: An instance of the FluxIO class.
    """
    return FluxIO(*args)

class FluxIO(CoSimulationIO):
    """
    This is the IO wrapper for the Flux solver.
    All Import/Export functions are methods of the FluxWrapper class.

    Args:
        settings (Kratos.Parameters): The Kratos parameters containing the settings for this IO.
        model (Kratos.Model): The Kratos model containing the problem's data.
        solver_name (str): The name of the solver associated with this IO.
        data_communicator (Kratos.DataCommunicator): The data communicator for parallel communication.

    Methods:
        __init__(self, settings, model, solver_name, data_communicator): Initialize FluxIO.
        PrintInfo(self): Print information about the IO.
        Check(self): Check the state of the IO.

    Note:
        This class inherits from CoSimulationIO and provides methods for Import/Export operations.

    """

    def __init__(self, settings, model, solver_name, data_communicator):
        """
        Initialize the FluxIO object.

        Args:
            settings (Kratos.Parameters): The Kratos parameters containing the settings for this IO.
            model (Kratos.Model): The Kratos model containing the problem's data.
            solver_name (str): The name of the solver associated with this IO.
            data_communicator (Kratos.DataCommunicator): The data communicator for parallel communication.
        """
        super().__init__(settings, model, solver_name, data_communicator)

    def PrintInfo(self):
        """
        Print information about the FluxIO.
        """
        print("This is the IO dedicated to Flux")

    def Check(self):
        """
        Check the state of the FluxIO.
        """
        pass