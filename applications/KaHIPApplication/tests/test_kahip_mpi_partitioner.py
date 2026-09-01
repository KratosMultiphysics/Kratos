"""
MPI tests for KaHIPDivideHeterogeneousInputInMemoryProcess.

These tests verify that the in-memory partitioner distributes mesh data correctly
across MPI ranks. They mirror the MetisApplication MPI test pattern.

Run with:  mpirun -n 4 python test_KaHIPApplication_mpi.py
"""

import pathlib
import os

import KratosMultiphysics as KM
from KratosMultiphysics import (
    KratosUnittest,
    Model,
    ModelPartIO,
    ReorderConsecutiveModelPartIO,
    Logger,
    Parameters,
    PARTITION_INDEX,
    TEMPERATURE,
)
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.KaHIPApplication  # ensure application is registered
from KratosMultiphysics.KaHIPApplication import (
    KaHIPDivideHeterogeneousInputProcess,
    KaHIPDivideHeterogeneousInputInMemoryProcess,
    KaHIPPartitioningModeler,
)

WORK_DIR = pathlib.Path(__file__).parent.parent.absolute() / "test_examples"


class TestKaHIPMPIPartitioner(KratosUnittest.TestCase):
    """Integration tests for MPI-parallel KaHIP partitioning."""

    def setUp(self):
        self.comm = KM.Testing.GetDefaultDataCommunicator()
        self.rank = self.comm.Rank()
        self.size = self.comm.Size()

    # ── In-memory process: basic correctness ──────────────────────────────────

    def test_in_memory_process_partitions_cube(self):
        """Each rank must receive local nodes and elements after in-memory partition."""
        file_name = str(WORK_DIR / "cube")
        io_flags = ModelPartIO.READ | ModelPartIO.SKIP_TIMER | ModelPartIO.IGNORE_VARIABLES_ERROR

        reorder_io = ReorderConsecutiveModelPartIO(file_name, io_flags)
        serial_io = ModelPartIO(file_name, io_flags)

        partitioner = KaHIPDivideHeterogeneousInputInMemoryProcess(
            reorder_io, serial_io, self.comm,
            dimension=3, verbosity=0, synchronize_conditions=True)
        partitioner.Execute()

        # After Execute, serial_io reads from the local in-memory buffer
        model = Model()
        model_part = model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        serial_io.ReadModelPart(model_part)

        # Each rank must have at least one local node and element
        local_nodes = model_part.GetCommunicator().LocalMesh().NumberOfNodes()
        local_elements = model_part.GetCommunicator().LocalMesh().NumberOfElements()
        if self.size <= 10:
            self.assertGreater(local_nodes, 0,
                msg=f"Rank {self.rank}: no local nodes after in-memory partitioning")
            self.assertGreater(local_elements, 0,
                msg=f"Rank {self.rank}: no local elements after in-memory partitioning")

    def test_in_memory_process_global_totals(self):
        """Global node/element counts must match the original mesh."""
        file_name = str(WORK_DIR / "cube")
        io_flags = ModelPartIO.READ | ModelPartIO.SKIP_TIMER | ModelPartIO.IGNORE_VARIABLES_ERROR

        reorder_io = ReorderConsecutiveModelPartIO(file_name, io_flags)
        serial_io = ModelPartIO(file_name, io_flags)

        settings = Parameters("""{
            "preconfiguration": "eco",
            "imbalance": 0.03,
            "seed": 0,
            "echo_level": 0,
            "num_trials": 1
        }""")
        partitioner = KaHIPDivideHeterogeneousInputInMemoryProcess(
            reorder_io, serial_io, self.comm, settings, True)
        partitioner.Execute()

        model = Model()
        model_part = model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        serial_io.ReadModelPart(model_part)
        KratosMPI.ParallelFillCommunicator(model_part.GetRootModelPart(), self.comm).Execute()

        # cube.mdpa has 413 nodes, 1191 elements, 780 conditions
        total_nodes = model_part.GetCommunicator().GlobalNumberOfNodes()
        total_elements = model_part.GetCommunicator().GlobalNumberOfElements()
        total_conditions = model_part.GetCommunicator().GlobalNumberOfConditions()

        self.assertEqual(total_nodes, 413,
            msg="Global node count mismatch after in-memory partitioning")
        self.assertEqual(total_elements, 1191,
            msg="Global element count mismatch after in-memory partitioning")
        self.assertEqual(total_conditions, 780,
            msg="Global condition count mismatch after in-memory partitioning")

    # ── Modeler: file-based MPI partitioning ──────────────────────────────────

    def test_modeler_file_based_partitioning(self):
        """KaHIPPartitioningModeler (file-based) must distribute the cube mesh."""
        file_name = str(WORK_DIR / "cube")
        model = Model()
        model_part = model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3

        settings = Parameters("""{
            "model_part_name": "Main",
            "input_filename": \"""" + file_name + """\",
            "number_of_partitions": 0,
            "dimension": 3,
            "verbosity": 0,
            "synchronize_conditions": true,
            "partition_in_memory": false,
            "perform_partitioning": true,
            "kahip_settings": {
                "preconfiguration": "eco",
                "imbalance": 0.03,
                "seed": 0,
                "echo_level": 0,
                "num_trials": 1
            }
        }""")

        modeler = KaHIPPartitioningModeler(model, settings)
        modeler.SetupModelPart()
        KratosMPI.ParallelFillCommunicator(model_part.GetRootModelPart(), self.comm).Execute()

        local_nodes = model_part.GetCommunicator().LocalMesh().NumberOfNodes()
        if self.size <= 10:
            self.assertGreater(local_nodes, 0,
                msg=f"Rank {self.rank}: no local nodes after file-based modeler partitioning")

        total_nodes = model_part.GetCommunicator().GlobalNumberOfNodes()
        self.assertEqual(total_nodes, 413)

        # Clean up partitioned directory
        self.comm.Barrier()
        if self.rank == 0:
            import shutil
            partitioned_dir = str(WORK_DIR / "cube_partitioned")
            if pathlib.Path(partitioned_dir).exists():
                shutil.rmtree(partitioned_dir)
        self.comm.Barrier()

    def test_modeler_in_memory_partitioning(self):
        """KaHIPPartitioningModeler (in-memory) must distribute the cube mesh."""
        file_name = str(WORK_DIR / "cube")
        model = Model()
        model_part = model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3

        settings = Parameters("""{
            "model_part_name": "Main",
            "input_filename": \"""" + file_name + """\",
            "number_of_partitions": 0,
            "dimension": 3,
            "verbosity": 0,
            "synchronize_conditions": true,
            "partition_in_memory": true,
            "perform_partitioning": true,
            "kahip_settings": {
                "preconfiguration": "fast",
                "imbalance": 0.05,
                "seed": 42,
                "echo_level": 0,
                "num_trials": 1
            }
        }""")

        modeler = KaHIPPartitioningModeler(model, settings)
        modeler.SetupModelPart()
        KratosMPI.ParallelFillCommunicator(model_part.GetRootModelPart(), self.comm).Execute()

        local_nodes = model_part.GetCommunicator().LocalMesh().NumberOfNodes()
        if self.size <= 10:
            self.assertGreater(local_nodes, 0,
                msg=f"Rank {self.rank}: no local nodes after in-memory modeler partitioning")

        total_nodes = model_part.GetCommunicator().GlobalNumberOfNodes()
        self.assertEqual(total_nodes, 413)


if __name__ == '__main__':
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
