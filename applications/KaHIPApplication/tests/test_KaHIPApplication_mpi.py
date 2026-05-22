"""
MPI test launcher for KaHIPApplication.

Run with:
    mpirun -n 4 python test_KaHIPApplication_mpi.py
or via the Kratos testing infrastructure:
    python kratos_run_tests.py --level mpi_small
"""

import KratosMultiphysics
import KratosMultiphysics.KaHIPApplication  # register the application

if not KratosMultiphysics.IsDistributedRun():
    raise Exception(
        "test_KaHIPApplication_mpi.py must be run with mpirun (distributed mode).")

import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_kahip_mpi_partitioner import TestKaHIPMPIPartitioner


def AssembleTestSuites():
    """Return the test suites for the MPI tests."""
    suites = KratosUnittest.KratosSuites

    small_mpi = suites['mpi_small']
    small_mpi.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TestKaHIPMPIPartitioner,
        ]))

    nightly_mpi = suites['mpi_nightly']
    nightly_mpi.addTests(small_mpi)

    all_mpi = suites['mpi_all']
    all_mpi.addTests(nightly_mpi)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
