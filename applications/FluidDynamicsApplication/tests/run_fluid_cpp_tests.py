import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing import utilities as testing_utils

def run():
    commander = testing_utils.Commander()
    commander.RunCppTests('FluidDynamicsApplication',int(900),{})

if __name__ == '__main__':
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    run()