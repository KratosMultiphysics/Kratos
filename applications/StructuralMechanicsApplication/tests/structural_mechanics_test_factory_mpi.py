from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

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

        # Here we append "_mpi" to the name of the results file
        # e.g. "test_results.json" => "test_results_mpi.json"
        # This is necessary becs the node numbering changes with the partitioning
        results_file_param = project_parameters["json_check_process"][0]["Parameters"]["input_file_name"]
        results_file_name = results_file_param.GetString()
        raw_path, file_name = os.path.split(results_file_name)
        raw_file_name, file_ext = os.path.splitext(file_name)
        mpi_file_name = raw_file_name + "_mpi" + file_ext
        full_new_path = os.path.join(raw_path, mpi_file_name)
        results_file_param.SetString(full_new_path)

        # saving the name of the input file such that the partitioned files can be deleted later
        self.input_filename = project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()

    def tearDown(self):
        super(StructuralMechanicsTestFactoryMPI, self).tearDown()

        # Now delete the partitioned mdpa files
        for i in range(KratosMPI.mpi.size):
            partitioned_mdpa_file_name = self.input_filename + "_" + str(i) + ".mdpa"
            kratos_utils.DeleteFileIfExisting(GetFilePath(partitioned_mdpa_file_name))


class ShellT3AndQ4LinearStaticStructPinchedCylinderTests(StructuralMechanicsTestFactoryMPI):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_pinched_cylinder"
