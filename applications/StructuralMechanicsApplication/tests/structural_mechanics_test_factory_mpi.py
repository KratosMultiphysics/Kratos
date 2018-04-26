from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

from  structural_mechanics_test_factory import StructuralMechanicsTestFactory

class StructuralMechanicsTestFactoryMPI(StructuralMechanicsTestFactory):

    def modify_parameters(self, project_parameters):
        """Setting MPI for the tests execution
        """
        project_parameters["problem_data"]["parallel_type"].SetString("MPI")

class ShellT3AndQ4LinearStaticStructPinchedCylinderTests(StructuralMechanicsTestFactoryMPI):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_pinched_cylinder"
