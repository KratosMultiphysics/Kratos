import os

# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis
import KratosMultiphysics.DEMApplication.DEM_procedures as DEM_procedures
# This utility will control the execution scope
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class PorosityConservationTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Setting parameters

            with open(self.file_parameters,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Create Model
            model = KratosMultiphysics.Model()

            self.test = PorosityConservationAnalysis(model, parameters)

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Run()

    def tearDown(self):
        pass


# General test factory
class TestFactory(KratosUnittest.TestCase):

    def setUp(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Setting parameters

            with open(self.file_parameters,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Create Model
            model = KratosMultiphysics.Model()

            self.test = SwimmingDEMAnalysis(model, parameters)

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Run()

    def tearDown(self):
        pass

class PorosityConservationAnalysis(SwimmingDEMAnalysis, DEM_procedures.Procedures):

    def __init__(self, model, parameters):
        super().__init__(model, parameters)

    def Initialize(self):
        super().Initialize()
        solid_domain_volume, fluid_domain_volume = self.ComputePhysicalVolumes()
        self.solid_domain_volume = solid_domain_volume
        self.fluid_domain_volume = fluid_domain_volume

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.CheckConservation()

    def ComputePhysicalVolumes(self):
        import numpy
        solid_volume = 0.0
        for node in self.spheres_model_part.Nodes:
            sphere_volume = node.GetSolutionStepValue(KratosMultiphysics.RADIUS) ** 3 * 4/3 * numpy.pi
            solid_volume += sphere_volume

        fluid_volume = 0.0
        for node in self.fluid_model_part.Nodes:
            fluid_volume += node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)

        return solid_volume, fluid_volume

    def CheckConservation(self):
        import numpy
        numerical_fluid_volume = 0.0
        for node in self.fluid_model_part.Nodes:
            numerical_fluid_volume +=  node.GetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION) * node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)

        assert numpy.abs(self.fluid_domain_volume - self.solid_domain_volume - numerical_fluid_volume) < 1e-17
