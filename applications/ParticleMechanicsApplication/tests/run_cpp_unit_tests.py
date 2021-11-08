from KratosMultiphysics import *
from KratosMultiphysics.ParticleMechanicsApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosParticleMechanicsFastSuite")

if __name__ == '__main__':
    run()
