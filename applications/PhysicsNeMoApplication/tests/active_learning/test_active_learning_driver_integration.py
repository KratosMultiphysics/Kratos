import queue
import sys
from pathlib import Path

import numpy

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.active_learning.sample_io import KratosALSample
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.in_process_backend import InProcessBackend
from KratosMultiphysics.PhysicsNeMoApplication.active_learning import kratos_label_strategy
from KratosMultiphysics.PhysicsNeMoApplication.active_learning import metrology
from KratosMultiphysics.PhysicsNeMoApplication.active_learning import query_strategies

try:
    import torch
    import physicsnemo
    from physicsnemo.active_learning.driver import Driver
    from physicsnemo.active_learning.config import DriverConfig, StrategiesConfig, TrainingConfig
    from physicsnemo.active_learning.protocols import QueryStrategy
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False

_AUX_FILES = Path(__file__).parent / "aux_files"
sys.path.insert(0, str(_AUX_FILES))

if have_physicsnemo:
    class _StubQueryStrategy(QueryStrategy):
        """Enqueues two fixed design points."""

        __protocol_name__ = "StubQueryStrategy"
        max_samples = 2

        def __init__(self):
            self._driver = None

        def attach(self, other):
            self._driver = other

        @property
        def is_attached(self):
            return self._driver is not None

        def sample(self, query_queue, *args, **kwargs):
            query_queue.put(KratosALSample("query_0", parameters={"dummy_settings/alpha": 1.5}))
            query_queue.put(KratosALSample("query_1", parameters={"dummy_settings/alpha": -1.0}))

    class _TinyLearner(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.linear = torch.nn.Linear(1, 1)

        def forward(self, x):
            return self.linear(x)


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestActiveLearningDriverIntegration(KratosUnittest.TestCase):
    """End-to-end: a real physicsnemo Driver labeling through a Kratos solve."""

    def setUp(self):
        self.working_directory = Path("test_al_driver_cases")
        self.log_directory = Path("test_al_driver_logs")
        self.parameters_file = Path("test_al_driver_parameters.json")
        with open(_AUX_FILES / "template" / "ProjectParameters.json") as f:
            self.parameters_file.write_text(f.read())

    def tearDown(self):
        KratosUtilities.DeleteDirectoryIfExisting(str(self.working_directory))
        KratosUtilities.DeleteDirectoryIfExisting(str(self.log_directory))
        KratosUtilities.DeleteFileIfExisting(str(self.parameters_file))
        KratosUtilities.DeleteFileIfExisting("test_al_driver_metrology.json")

    def _CreateBackend(self):
        return InProcessBackend(Kratos.Parameters("""{
            "project_parameters_file" : "test_al_driver_parameters.json",
            "analysis_stage_module"   : "dummy_analysis",
            "working_directory"       : "test_al_driver_cases",
            "model_part_name"         : "Main",
            "output_field_specs"      : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ]
        }"""))

    def test_OneActiveLearningStep(self):
        label_strategy = kratos_label_strategy.CreateKratosLabelStrategy(
            self._CreateBackend(), provides_fields={"PRESSURE__node_historical"})

        driver_config = DriverConfig(
            batch_size=1,
            max_active_learning_steps=1,
            skip_training=True,
            skip_metrology=True,
            checkpoint_interval=0,
            checkpoint_on_labeling=False,
            root_log_dir=str(self.log_directory))
        strategies_config = StrategiesConfig(
            query_strategies=[_StubQueryStrategy()],
            queue_cls=queue.Queue,
            label_strategy=label_strategy)
        # Labeled samples are appended to the training pool by the driver; a
        # plain list satisfies the DataPool protocol.
        train_datapool = []
        training_config = TrainingConfig(
            train_datapool=train_datapool,
            max_training_epochs=1)

        driver = Driver(
            config=driver_config,
            learner=_TinyLearner(),
            strategies_config=strategies_config,
            training_config=training_config)
        driver.active_learning_step()

        # The query -> label -> serialize contract executed: both samples
        # were labeled by an actual Kratos solve and the driver appended them
        # to the training pool.
        labeled = train_datapool
        self.assertEqual(len(labeled), 2)
        by_id = {sample.sample_id: sample for sample in labeled}
        self.assertTrue(numpy.allclose(
            by_id["query_0"].fields["PRESSURE__node_historical"], 1.5 * numpy.arange(5.0)))
        self.assertTrue(numpy.allclose(
            by_id["query_1"].fields["PRESSURE__node_historical"], -1.0 * numpy.arange(5.0)))
        self.assertEqual(label_strategy.failed_samples, 0)

    def test_FullStrategyStackStep(self):
        # The application's own query, label AND metrology strategies wired
        # into one real physicsnemo Driver step.
        query_strategy = query_strategies.CreateSolverResidualStrategy(
            Kratos.Parameters("""{ "max_samples": 2, "candidate_pool_size": 3 }"""),
            lambda n: [{"dummy_settings/alpha": float(v)} for v in (0.5, 3.0, -2.0)],
            lambda candidate: abs(candidate["dummy_settings/alpha"]))
        label_strategy = kratos_label_strategy.CreateKratosLabelStrategy(
            self._CreateBackend(), provides_fields={"PRESSURE__node_historical"})
        metrology_strategy = metrology.CreateValidationMetricsMetrology(
            Kratos.Parameters("""{
                "metrics"     : ["mse", "max_abs_error"],
                "output_file" : "test_al_driver_metrology.json"
            }"""),
            lambda: {"PRESSURE": (numpy.ones(5), numpy.zeros(5))})

        driver_config = DriverConfig(
            batch_size=1,
            max_active_learning_steps=1,
            skip_training=True,
            skip_metrology=False,
            checkpoint_interval=0,
            checkpoint_on_labeling=False,
            root_log_dir=str(self.log_directory))
        strategies_config = StrategiesConfig(
            query_strategies=[query_strategy],
            queue_cls=queue.Queue,
            label_strategy=label_strategy,
            metrology_strategies=[metrology_strategy])
        train_datapool = []
        training_config = TrainingConfig(
            train_datapool=train_datapool,
            max_training_epochs=1)

        driver = Driver(
            config=driver_config,
            learner=_TinyLearner(),
            strategies_config=strategies_config,
            training_config=training_config)
        driver.active_learning_step()

        # The residual strategy picked the two highest |alpha| candidates and
        # Kratos labeled them.
        self.assertEqual(len(train_datapool), 2)
        alphas = sorted(sample.parameters["dummy_settings/alpha"] for sample in train_datapool)
        self.assertEqual(alphas, [-2.0, 3.0])
        for sample in train_datapool:
            self.assertTrue(numpy.allclose(
                sample.fields["PRESSURE__node_historical"],
                sample.parameters["dummy_settings/alpha"] * numpy.arange(5.0)))
        self.assertEqual(label_strategy.failed_samples, 0)

        # The metrology phase computed and serialized one record.
        self.assertEqual(len(metrology_strategy.records), 1)
        self.assertAlmostEqual(metrology_strategy.records[0]["PRESSURE"]["mse"], 1.0, places=12)
        self.assertTrue(Path("test_al_driver_metrology.json").is_file())


if __name__ == '__main__':
    KratosUnittest.main()
