from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_prebuckling_analysis import StructuralMechanicsPrebucklingAnalysis
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import kratos_utilities as kratos_utils
eigensolvers_application_available = kratos_utils.CheckIfApplicationsAvailable("EigenSolversApplication")

'''
Test description:
This test does an eigenvalue analysis on a simple cantilever beam
It compares the results between the scipy solver and the standard eigen solver.
'''

@KratosUnittest.skipUnless(eigensolvers_application_available,"Missing required application: EigenSolversApplication")
class TestPrebucklingAnalysis(KratosUnittest.TestCase):
    # muting the output
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    def test_prebuckling_analsis(self):
        reference = 92.80
        with open("prebuckling/ProjectParameters_Symmetric.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        model_symmetric = KratosMultiphysics.Model()
        simulation_symmetric = StructuralMechanicsPrebucklingAnalysis(model_symmetric,parameters)
        simulation_symmetric.Run()
        model_part_symmetric = model_symmetric.GetModelPart("Structure")

        with open("prebuckling/ProjectParameters_Full.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        model_full = KratosMultiphysics.Model()
        simulation_full = StructuralMechanicsPrebucklingAnalysis(model_full,parameters)
        simulation_full.Run()
        model_part_full = model_full.GetModelPart("Structure")

        self._check_load_multiplier(model_part_symmetric, model_part_full, reference)

    def _check_load_multiplier(self,mp_sym, mp_full, reference):
        #Check if value stays the same in the first and last loadstep
        #self.assertLess( abs(1-load_multiplier1[0]/load_multiplier1[7]), 1.0e-4)
        load_multiplier1 = mp_sym.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR][0]
        load_multiplier2 = mp_full.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR][0]
        print("sym: ", load_multiplier1)
        print("full: ", load_multiplier2)
        self.assertAlmostEqual(load_multiplier1, load_multiplier2, 5)

        #self.assertAlmostEqual(load_multiplier1[0], load_multiplier1[1], 1)
        #Check if both models give same values
        #self.assertLess( abs(1-load_multiplier1[0]/reference), 1.0e-2)
        #Compare value against reference from abaqus



if __name__ == '__main__':
    KratosUnittest.main()
