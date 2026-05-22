"""
Tests for KaHIPDivideHeterogeneousInputProcess.

These tests mirror test_quad_partition.py from MetisApplication, verifying that:
  - The partitioner produces valid partition assignments
  - Condition synchronisation places each boundary condition in the same
    partition as its parent element
  - The output files are written correctly and can be read back
"""

import pathlib

import KratosMultiphysics as KM
from KratosMultiphysics import (
    KratosUnittest,
    Model,
    ModelPartIO,
    PARTITION_INDEX,
    FRACTIONAL_STEP,
    Logger,
    Parameters,
)
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.KaHIPApplication import KaHIPDivideHeterogeneousInputProcess


WORK_DIR = pathlib.Path(__file__).parent.absolute()


class TestKaHIPPartitioner(KratosUnittest.TestCase):
    """Integration tests for KaHIPDivideHeterogeneousInputProcess."""

    # ── Condition synchronisation ─────────────────────────────────────────────

    def test_condition_partitioning_synchronized(self):
        """Every condition must reside in the same partition as its parent element."""
        test_file = "quads"
        partitions = 3
        sync_conditions = True
        domain_size = 2
        verbosity = 0

        self.addCleanup(
            DeleteDirectoryIfExisting,
            WORK_DIR / f'{test_file}_partitioned')

        with KratosUnittest.WorkFolderScope(WORK_DIR, ''):
            io_flags = ModelPartIO.READ | ModelPartIO.SKIP_TIMER | ModelPartIO.IGNORE_VARIABLES_ERROR
            io = ModelPartIO(test_file, io_flags)

            partitioner = KaHIPDivideHeterogeneousInputProcess(
                io, partitions, domain_size, verbosity, sync_conditions)
            partitioner.Execute()

            model = Model()
            for i in range(partitions):
                mp = model.CreateModelPart(f'part_{i}')
                mp.AddNodalSolutionStepVariable(PARTITION_INDEX)

                part_io = ModelPartIO(
                    str(WORK_DIR / f'{test_file}_partitioned/{test_file}_{i}'),
                    io_flags)
                part_io.ReadModelPart(mp)

                for condition in mp.Conditions:
                    parent_id = condition.GetValue(FRACTIONAL_STEP)
                    self.assertIn(
                        parent_id, mp.Elements,
                        msg=f"Partition {i}: condition parent element {parent_id} not in partition")

    # ── Partition index range ─────────────────────────────────────────────────

    def test_partition_index_range(self):
        """All nodes must receive a PARTITION_INDEX in [0, nparts)."""
        test_file = "quads"
        partitions = 4
        domain_size = 2
        verbosity = 0

        self.addCleanup(
            DeleteDirectoryIfExisting,
            WORK_DIR / f'{test_file}_partitioned')

        with KratosUnittest.WorkFolderScope(WORK_DIR, ''):
            io_flags = ModelPartIO.READ | ModelPartIO.SKIP_TIMER | ModelPartIO.IGNORE_VARIABLES_ERROR
            io = ModelPartIO(test_file, io_flags)

            partitioner = KaHIPDivideHeterogeneousInputProcess(
                io, partitions, domain_size, verbosity, False)
            partitioner.Execute()

            model = Model()
            node_counts = [0] * partitions
            for i in range(partitions):
                mp = model.CreateModelPart(f'range_test_{i}')
                mp.AddNodalSolutionStepVariable(PARTITION_INDEX)

                part_io = ModelPartIO(
                    str(WORK_DIR / f'{test_file}_partitioned/{test_file}_{i}'),
                    io_flags)
                part_io.ReadModelPart(mp)

                for node in mp.Nodes:
                    part_idx = node.GetSolutionStepValue(PARTITION_INDEX)
                    self.assertGreaterEqual(int(part_idx), 0)
                    self.assertLess(int(part_idx), partitions)
                    if int(part_idx) == i:
                        node_counts[i] += 1

            # Each partition must own at least one local node
            for i, count in enumerate(node_counts):
                self.assertGreater(count, 0,
                    msg=f"Partition {i} has no locally-owned nodes")

    # ── Parameters-based constructor ──────────────────────────────────────────

    def test_parameters_constructor(self):
        """Parameters-based constructor with ECO mode should produce a valid partition."""
        test_file = "quads"
        partitions = 3

        self.addCleanup(
            DeleteDirectoryIfExisting,
            WORK_DIR / f'{test_file}_partitioned')

        with KratosUnittest.WorkFolderScope(WORK_DIR, ''):
            io_flags = ModelPartIO.READ | ModelPartIO.SKIP_TIMER | ModelPartIO.IGNORE_VARIABLES_ERROR
            io = ModelPartIO(test_file, io_flags)

            settings = Parameters("""{
                "preconfiguration": "eco",
                "imbalance": 0.05,
                "seed": 1,
                "suppress_output": true,
                "num_trials": 3,
                "verbosity": 0,
                "synchronize_conditions": true
            }""")

            partitioner = KaHIPDivideHeterogeneousInputProcess(io, partitions, settings)
            partitioner.Execute()

            # Verify output files were created
            for i in range(partitions):
                part_file = WORK_DIR / f'{test_file}_partitioned/{test_file}_{i}.mdpa'
                self.assertTrue(
                    part_file.exists(),
                    msg=f"Expected partition file not found: {part_file}")

    # ── Multi-trial mode ──────────────────────────────────────────────────────

    def test_multi_trial(self):
        """Partitioner with num_trials=5 should still produce a valid result."""
        test_file = "quads"
        partitions = 2

        self.addCleanup(
            DeleteDirectoryIfExisting,
            WORK_DIR / f'{test_file}_partitioned')

        with KratosUnittest.WorkFolderScope(WORK_DIR, ''):
            io_flags = ModelPartIO.READ | ModelPartIO.SKIP_TIMER | ModelPartIO.IGNORE_VARIABLES_ERROR
            io = ModelPartIO(test_file, io_flags)

            settings = Parameters("""{
                "preconfiguration": "fast",
                "imbalance": 0.03,
                "seed": 0,
                "suppress_output": true,
                "num_trials": 5
            }""")

            partitioner = KaHIPDivideHeterogeneousInputProcess(io, partitions, settings)
            partitioner.Execute()

            # Verify that each partition file exists and has elements
            model = Model()
            for i in range(partitions):
                mp = model.CreateModelPart(f'multi_trial_{i}')
                mp.AddNodalSolutionStepVariable(PARTITION_INDEX)
                part_io = ModelPartIO(
                    str(WORK_DIR / f'{test_file}_partitioned/{test_file}_{i}'),
                    io_flags)
                part_io.ReadModelPart(mp)
                self.assertGreater(mp.NumberOfNodes(), 0)
                self.assertGreater(mp.NumberOfElements(), 0)

    # ── Preconfiguration modes ────────────────────────────────────────────────

    def _run_with_preconfiguration(self, mode):
        test_file = "quads"
        partitions = 3

        self.addCleanup(
            DeleteDirectoryIfExisting,
            WORK_DIR / f'{test_file}_partitioned')

        with KratosUnittest.WorkFolderScope(WORK_DIR, ''):
            io_flags = ModelPartIO.READ | ModelPartIO.SKIP_TIMER | ModelPartIO.IGNORE_VARIABLES_ERROR
            io = ModelPartIO(test_file, io_flags)

            settings = Parameters(f"""{{
                "preconfiguration": "{mode}",
                "imbalance": 0.03,
                "seed": 0,
                "suppress_output": true,
                "num_trials": 1
            }}""")

            partitioner = KaHIPDivideHeterogeneousInputProcess(io, partitions, settings)
            partitioner.Execute()

            for i in range(partitions):
                self.assertTrue(
                    (WORK_DIR / f'{test_file}_partitioned/{test_file}_{i}.mdpa').exists())

    def test_fast_mode(self):
        self._run_with_preconfiguration("fast")

    def test_eco_mode(self):
        self._run_with_preconfiguration("eco")

    def test_strong_mode(self):
        self._run_with_preconfiguration("strong")


if __name__ == '__main__':
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
