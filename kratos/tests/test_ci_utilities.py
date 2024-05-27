from pathlib import Path
from typing import List
from unittest.mock import patch

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing import ci_utilities

PATCH_NAME: str = "KratosMultiphysics.testing.ci_utilities.changed_files"

FILES_CHANGED_1: List[Path] = [Path("kratos/linear_solvers/amgcl_solver.h")]

FILES_CHANGED_2: List[Path] = [Path("kratos/mpi/python_scripts/distributed_import_model_part_utility.py")]

FILES_CHANGED_3: List[Path] = [
    Path("applications/FluidDynamicsApplication/custom_elements/vms.cpp"),
    Path("applications/LinearSolversApplication/custom_solvers/eigen_direct_solver.h"),
    Path("applications/StructuralMechanicsApplication/CMakeLists.txt"),
]


class TestCIUtilities(KratosUnittest.TestCase):
    @patch(PATCH_NAME, return_value=FILES_CHANGED_1)
    def test_core_changes(self, _):
        self.assertCountEqual(ci_utilities.get_changed_applications(), set())
        self.assertTrue(ci_utilities.is_core_changed())
        self.assertFalse(ci_utilities.is_mpi_core_changed())
        self.assertCountEqual(ci_utilities.get_changed_files_extensions(), {".h"})
        self.assertFalse(ci_utilities.are_only_python_files_changed())

    @patch(PATCH_NAME, return_value=FILES_CHANGED_2)
    def test_mpi_core_changes(self, _):
        self.assertCountEqual(ci_utilities.get_changed_applications(), set())
        self.assertTrue(ci_utilities.is_core_changed())
        self.assertTrue(ci_utilities.is_mpi_core_changed())
        self.assertCountEqual(ci_utilities.get_changed_files_extensions(), {".py"})
        self.assertTrue(ci_utilities.are_only_python_files_changed())

    @patch(PATCH_NAME, return_value=FILES_CHANGED_3)
    def test_application_changes(self, _):
        self.assertCountEqual(
            ci_utilities.get_changed_applications(),
            {"FluidDynamicsApplication", "LinearSolversApplication", "StructuralMechanicsApplication"},
        )
        self.assertFalse(ci_utilities.is_core_changed())
        self.assertFalse(ci_utilities.is_mpi_core_changed())
        self.assertCountEqual(ci_utilities.get_changed_files_extensions(), {".cpp", ".h", ".txt"})
        self.assertFalse(ci_utilities.are_only_python_files_changed())


if __name__ == "__main__":
    KratosUnittest.main()
