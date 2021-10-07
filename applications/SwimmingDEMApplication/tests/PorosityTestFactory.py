import os

# Importing the Kratos Library
import KratosMultiphysics
import numpy as np
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
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'porosity_tests/porosity_conservation/')

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.solid_domain_volume = self.ComputeSolidVolume()
        self.CheckConservation()

    def ComputeSolidVolume(self):
        solid_volume = 0.0
        spheres_model_part = self.model.GetModelPart('SpheresPart')
        for node in spheres_model_part.Nodes:
            sphere_volume = node.GetSolutionStepValue(KratosMultiphysics.RADIUS) ** 3 * 4/3 * np.pi
            solid_volume += sphere_volume
        return solid_volume

    def CheckConservation(self):
        numerical_fluid_volume = 0.0
        fluid_volume = 0.0
        for node in self.fluid_model_part.Nodes:
            fluid_volume += node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
            numerical_fluid_volume +=  node.GetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION) * node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)

        assert fluid_volume != numerical_fluid_volume
        assert np.abs(fluid_volume - self.solid_domain_volume - numerical_fluid_volume) < 1e-17
