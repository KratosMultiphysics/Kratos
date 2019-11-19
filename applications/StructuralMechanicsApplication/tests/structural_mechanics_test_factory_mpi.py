from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

from  structural_mechanics_test_factory import StructuralMechanicsTestFactory
import KratosMultiphysics.kratos_utilities as kratos_utils

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class StructuralMechanicsTestFactoryMPI(StructuralMechanicsTestFactory):

    def modify_parameters(self, project_parameters):
        """Setting MPI for the tests execution
        The existing parameters are modified on the fly to execute in MPI
        """
        # Change parallel type to MPI
        project_parameters["problem_data"]["parallel_type"].SetString("MPI")
        model_import_settings = project_parameters["solver_settings"]["model_import_settings"]
        if not model_import_settings.Has("partition_in_memory"):
            model_import_settings.AddEmptyValue("partition_in_memory")
        model_import_settings["partition_in_memory"].SetBool(True)


class ShellT3AndQ4LinearStaticStructPinchedCylinderTests(StructuralMechanicsTestFactoryMPI):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_pinched_cylinder"
