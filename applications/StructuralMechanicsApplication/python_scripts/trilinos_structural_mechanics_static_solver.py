from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.mpi as mpi
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.MetisApplication as MetisApplication
import trilinos_structural_mechanics_solver


# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return TrilinosStaticMechanicalSolver(main_model_part, custom_settings)


class TrilinosStaticMechanicalSolver(trilinos_structural_mechanics_solver.TrilinosMechanicalSolver):
    """The trilinos structural mechanics static solver.

    For more information see:
    structural_mechanics_solver.py
    trilinos_structural_mechanics_solver.py
    """
    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver.
        super(TrilinosStaticMechanicalSolver, self).__init__(main_model_part, custom_settings)
        print("::[TrilinosStaticMechanicalSolver]:: Construction finished")

    def _create_solution_scheme(self):
        return TrilinosApplication.TrilinosResidualBasedIncrementalUpdateStaticScheme()
