"""Process validating ML-predicted fields against reference Kratos fields.

Per configured interval, computes physicsnemo.metrics.general metrics between
a predicted and a reference field, logs them, and writes the accumulated
history as a JSON report in ExecuteFinalize. This provides model-vs-solver
benchmarking with core physicsnemo only.

physicsnemo/torch are imported lazily at first evaluation.
"""

import json
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor

SUPPORTED_METRICS = ("mse", "rmse", "max_abs_error", "wasserstein",
                     "relative_l2", "weighted_mse", "weighted_rmse",
                     "relative_mse", "histogram_l1", "entropy_difference")
DISTRIBUTION_METRICS = ("histogram_l1", "entropy_difference")
SPECTRAL_METRICS = ("power_spectrum_relative_l2", "high_wavenumber_energy_ratio")
WEIGHTED_METRICS = ("weighted_mse", "weighted_rmse")
ENSEMBLE_METRICS = ("crps", "kcrps")
CALIBRATION_METRICS = ("coverage", "nll", "sharpness", "calibration_error")


def _TryImportPhysicsNemoMetrics():
    try:
        from physicsnemo.metrics.general import mse as mse_module
        from physicsnemo.metrics.general import wasserstein as wasserstein_module
        return mse_module, wasserstein_module
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.validation_metrics_process requires physicsnemo, which "
            "could not be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def _TryImportPhysicsNemoCrps():
    try:
        from physicsnemo.metrics.general import crps as crps_module
        return crps_module
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.validation_metrics_process requires physicsnemo, which "
            "could not be imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def _TryImportDistributionMetrics():
    try:
        from physicsnemo.metrics.general.histogram import histogram
        from physicsnemo.metrics.general.entropy import entropy_from_counts
        return histogram, entropy_from_counts
    except ImportError as e:
        raise ImportError(
            "The histogram and entropy metrics require physicsnemo, which could not be "
            "imported. Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def _TryImportPowerSpectrum():
    try:
        from physicsnemo.metrics.general.power_spectrum import power_spectrum
        return power_spectrum
    except ImportError as e:
        raise ImportError(
            "The spectral metrics require physicsnemo, which could not be imported. "
            "Install it with e.g. 'pip install nvidia-physicsnemo'.") from e


def _SharedHistograms(predicted, reference, bins: int):
    """Both fields' histograms over ONE set of bin edges.

    Comparing two distributions bin by bin is only meaningful on common
    edges; histogram(a, b) computes edges covering both, and each field is
    then counted into those edges separately.
    """
    histogram, _ = _TryImportDistributionMetrics()
    torch = torch_bridge._TryImportTorch()
    flat_predicted = predicted.reshape(-1).to(torch.float64)
    flat_reference = reference.reshape(-1).to(torch.float64)
    edges, _ = histogram(flat_predicted, flat_reference, bins=bins)
    _, predicted_counts = histogram(flat_predicted, bins=edges)
    _, reference_counts = histogram(flat_reference, bins=edges)
    return edges, predicted_counts.to(torch.float64), reference_counts.to(torch.float64)


def ComputeSpectralMetricValues(predicted_grid, reference_grid, metric_names,
                                high_wavenumber_fraction: float = 0.5):
    """Compares the azimuthally averaged power spectra of two gridded fields.

    The point of these metrics is what pointwise errors miss. A
    superresolution model that returns a smooth, plausible field can score a
    small RMSE while reproducing none of the fine-scale energy it was meant
    to recover; its spectrum shows that at once.

    Args:
        predicted_grid, reference_grid: torch tensors shaped (..., H, W) -
            the last two axes are the plane the spectrum is taken over, and
            any leading axes (channels, a depth axis) are averaged.
        metric_names: from SPECTRAL_METRICS.
        high_wavenumber_fraction: where "high wavenumbers" start, as a
            fraction of the wavenumber bins.

    Returns:
        (values, spectra): {metric_name: float}, and the wavenumbers with
        both averaged spectra as lists, for the report.
    """
    power_spectrum = _TryImportPowerSpectrum()
    torch = torch_bridge._TryImportTorch()

    if predicted_grid.shape != reference_grid.shape:
        raise ValueError(
            f"Shape mismatch between predicted grid {list(predicted_grid.shape)} and "
            f"reference grid {list(reference_grid.shape)}.")
    if predicted_grid.ndim < 2 or min(predicted_grid.shape[-2:]) < 4:
        raise ValueError(
            "A power spectrum needs a plane of at least 4x4 points in the last two "
            f"axes; got shape {list(predicted_grid.shape)}.")

    wavenumbers, predicted_power = power_spectrum(predicted_grid.to(torch.float64))
    _, reference_power = power_spectrum(reference_grid.to(torch.float64))
    predicted_power = predicted_power.reshape(-1, predicted_power.shape[-1]).mean(dim=0)
    reference_power = reference_power.reshape(-1, reference_power.shape[-1]).mean(dim=0)

    split = max(1, int(round(high_wavenumber_fraction * predicted_power.shape[0])))
    values = {}
    for metric in metric_names:
        if metric == "power_spectrum_relative_l2":
            values[metric] = float(
                torch.linalg.vector_norm(predicted_power - reference_power)
                / torch.linalg.vector_norm(reference_power))
        elif metric == "high_wavenumber_energy_ratio":
            predicted_share = predicted_power[split:].sum() / predicted_power.sum()
            reference_share = reference_power[split:].sum() / reference_power.sum()
            values[metric] = float(predicted_share / reference_share)
        else:
            raise ValueError(
                f"Unsupported spectral metric \"{metric}\". "
                f"Supported: {', '.join(SPECTRAL_METRICS)}.")
    spectra = {
        "wavenumbers": [float(k) for k in wavenumbers.reshape(-1)],
        "predicted": [float(v) for v in predicted_power],
        "reference": [float(v) for v in reference_power],
    }
    return values, spectra


def ComputeMetricValues(predicted, reference, metric_names, weights=None, bins: int = 32) -> dict:
    """Computes SUPPORTED_METRICS values between two same-shape torch tensors.

    The single metric implementation shared by ValidationMetricsProcess and
    the active-learning metrology strategy. Returns {metric_name: float}.

    "relative_l2" is ||predicted - reference|| / ||reference|| over the
    flattened fields (implemented locally - physicsnemo.metrics has no
    relative-L2), and "relative_mse" its square, also local: physicsnemo
    2.2 has no relative_mse either, whatever the changelog suggests.

    The two DISTRIBUTION metrics compare what values occur, not where:
    "histogram_l1" is the L1 distance between the two normalized histograms
    on shared bin edges (0 identical, 2 disjoint), and "entropy_difference"
    the absolute difference of their entropies. A surrogate can match a
    field pointwise on average while collapsing its extremes, and these say
    so. Entropies are taken UNNORMALIZED: physicsnemo's normalized
    entropy_from_counts is documented to map onto [0, 1] and returns about
    -0.9 for a delta distribution. The "weighted_*" metrics need a weights tensor whose first
    dimension matches the fields' (e.g. nodal areas): the squared error is
    averaged over each entity's components and weighted-averaged over
    entities.
    """
    mse_module, wasserstein_module = _TryImportPhysicsNemoMetrics()
    torch = torch_bridge._TryImportTorch()

    if predicted.shape != reference.shape:
        raise ValueError(
            f"Shape mismatch between predicted {list(predicted.shape)} and reference "
            f"{list(reference.shape)}.")

    weighted_error2 = None
    if any(metric in WEIGHTED_METRICS for metric in metric_names):
        if weights is None:
            raise ValueError(
                f"Metrics {WEIGHTED_METRICS} need a weights tensor (e.g. nodal areas); "
                "none was given.")
        weights = weights.reshape(-1)
        if weights.shape[0] != predicted.shape[0]:
            raise ValueError(
                f"weights has {weights.shape[0]} entries but the fields have "
                f"{predicted.shape[0]} entities.")
        entity_error2 = (predicted - reference).reshape(predicted.shape[0], -1).square().mean(dim=1)
        weighted_error2 = float((weights * entity_error2).sum() / weights.sum())

    distributions = None
    if any(metric in DISTRIBUTION_METRICS for metric in metric_names):
        distributions = _SharedHistograms(predicted, reference, bins)

    values = {}
    for metric in metric_names:
        if metric == "mse":
            values[metric] = float(mse_module.mse(predicted, reference))
        elif metric == "rmse":
            values[metric] = float(mse_module.rmse(predicted, reference))
        elif metric == "max_abs_error":
            values[metric] = float(torch.max(torch.abs(predicted - reference)))
        elif metric == "wasserstein":
            values[metric] = float(wasserstein_module.wasserstein_from_samples(
                predicted.reshape(-1), reference.reshape(-1)))
        elif metric == "relative_l2":
            values[metric] = float(
                torch.linalg.vector_norm((predicted - reference).reshape(-1))
                / torch.linalg.vector_norm(reference.reshape(-1)))
        elif metric == "weighted_mse":
            values[metric] = weighted_error2
        elif metric == "weighted_rmse":
            values[metric] = weighted_error2 ** 0.5
        elif metric == "relative_mse":
            values[metric] = float(
                (predicted - reference).square().sum() / reference.square().sum())
        elif metric == "histogram_l1":
            _, predicted_counts, reference_counts = distributions
            values[metric] = float(
                (predicted_counts / predicted_counts.sum()
                 - reference_counts / reference_counts.sum()).abs().sum())
        elif metric == "entropy_difference":
            _, entropy_from_counts = _TryImportDistributionMetrics()
            edges, predicted_counts, reference_counts = distributions
            values[metric] = float(abs(
                entropy_from_counts(predicted_counts, edges, normalized=False)
                - entropy_from_counts(reference_counts, edges, normalized=False)))
        else:
            raise ValueError(
                f"Unsupported metric \"{metric}\". Supported: {', '.join(SUPPORTED_METRICS)}.")
    return values


def ComputeCalibrationMetricValues(mean, std, reference, metric_names,
                                   confidence_z: float = 1.96) -> dict:
    """Computes CALIBRATION_METRICS for a predicted (mean, std) against truth.

    These answer the question an error metric cannot: is the model's stated
    uncertainty honest? A surrogate can have an excellent RMSE and still be
    badly calibrated, which is what makes its error bars unusable for
    decisions.

    - "coverage": fraction of references inside mean +/- confidence_z * std.
      For a well-calibrated Gaussian at z = 1.96 this is ~0.95.
    - "calibration_error": |coverage - nominal|, the signed miss turned
      absolute (0 is perfect; large means over- or under-confident).
    - "nll": mean Gaussian negative log-likelihood, which unlike coverage
      punishes both over- and under-confidence continuously.
    - "sharpness": mean predicted std. Meaningful only next to coverage -
      a model can be arbitrarily sharp by being wrong, or trivially
      well-covered by being vague.

    Returns {metric_name: float}.
    """
    torch = torch_bridge._TryImportTorch()

    if tuple(mean.shape) != tuple(reference.shape) or tuple(std.shape) != tuple(reference.shape):
        raise ValueError(
            f"Shape mismatch between mean {list(mean.shape)}, std {list(std.shape)} and "
            f"reference {list(reference.shape)}.")
    if float(std.min()) < 0.0:
        raise ValueError("std carries negative entries; it must be a standard deviation.")

    safe_std = std.clamp_min(1e-12)
    nominal = float(torch.special.erf(torch.tensor(
        confidence_z / (2.0 ** 0.5), dtype=torch.float64)))

    values = {}
    for metric in metric_names:
        if metric == "coverage":
            inside = ((reference >= mean - confidence_z * std)
                      & (reference <= mean + confidence_z * std))
            values[metric] = float(inside.to(torch.float64).mean())
        elif metric == "calibration_error":
            inside = ((reference >= mean - confidence_z * std)
                      & (reference <= mean + confidence_z * std))
            values[metric] = float(abs(float(inside.to(torch.float64).mean()) - nominal))
        elif metric == "nll":
            residual = (reference - mean) / safe_std
            values[metric] = float((0.5 * residual.square()
                                    + safe_std.log()
                                    + 0.5 * float(torch.log(torch.tensor(
                                        2.0 * torch.pi, dtype=torch.float64)))).mean())
        elif metric == "sharpness":
            values[metric] = float(std.mean())
        else:
            raise ValueError(
                f"Unsupported calibration metric \"{metric}\". "
                f"Supported: {', '.join(CALIBRATION_METRICS)}.")
    return values


def ComputeEnsembleMetricValues(ensemble, reference, metric_names) -> dict:
    """Computes ENSEMBLE_METRICS between an (M, ...) prediction ensemble and
    a (...)-shaped reference, via physicsnemo.metrics.general.crps.

    "crps" is the kernel CRPS estimator, "kcrps" the (biased) kernel form;
    both are averaged over all field entries. Returns {metric_name: float}.
    """
    crps_module = _TryImportPhysicsNemoCrps()

    if ensemble.ndim < 2 or ensemble.shape[0] < 2:
        raise ValueError(
            f"ensemble must be (M >= 2, ...); got shape {list(ensemble.shape)}.")
    if tuple(ensemble.shape[1:]) != tuple(reference.shape):
        raise ValueError(
            f"Shape mismatch between ensemble members {list(ensemble.shape[1:])} and "
            f"reference {list(reference.shape)}.")

    values = {}
    for metric in metric_names:
        if metric == "crps":
            values[metric] = float(crps_module.crps(ensemble, reference, dim=0).mean())
        elif metric == "kcrps":
            values[metric] = float(crps_module.kcrps(ensemble, reference, dim=0).mean())
        else:
            raise ValueError(
                f"Unsupported ensemble metric \"{metric}\". "
                f"Supported: {', '.join(ENSEMBLE_METRICS)}.")
    return values


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "ValidationMetricsProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return ValidationMetricsProcess(model, settings["Parameters"])


class ValidationMetricsProcess(Kratos.Process):
    """Evaluates configured metrics every output_interval steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_name"     : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "list_of_comparisons" : [
                {
                    "predicted_variable" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "predicted_location" : "node_historical",
                    "reference_variable" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "reference_location" : "node_historical",
                    "weight_variable"    : "",
                    "weight_location"    : "node_non_historical",
                    "bins"               : 32,
                    "metrics"            : ["mse", "rmse", "max_abs_error"]
                }
            ],
            "cfd_metrics"         : [],
            "uncertainty_comparisons" : [],
            "ensemble_comparisons"    : [],
            "spectral_comparisons"    : [],
            "output_interval"     : 1,
            "output_file"         : "validation_metrics.json"
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for i in range(settings["list_of_comparisons"].size()):
            settings["list_of_comparisons"][i].ValidateAndAssignDefaults(default_settings["list_of_comparisons"][0])

        self.model_part = model[settings["model_part_name"].GetString()]
        self.comparisons = []
        for i in range(settings["list_of_comparisons"].size()):
            comparison = settings["list_of_comparisons"][i]
            metrics = comparison["metrics"].GetStringArray()
            for metric in metrics:
                if metric not in SUPPORTED_METRICS:
                    raise ValueError(
                        f"Unsupported metric \"{metric}\". Supported: {', '.join(SUPPORTED_METRICS)}.")
            weight_variable = comparison["weight_variable"].GetString()
            if any(metric in WEIGHTED_METRICS for metric in metrics) and not weight_variable:
                raise ValueError(
                    f"Metrics {WEIGHTED_METRICS} need a \"weight_variable\" "
                    "(e.g. NODAL_AREA) in the comparison.")
            self.comparisons.append({
                "predicted": (comparison["predicted_variable"].GetString(),
                              comparison["predicted_location"].GetString()),
                "reference": (comparison["reference_variable"].GetString(),
                              comparison["reference_location"].GetString()),
                "weight": ((weight_variable, comparison["weight_location"].GetString())
                           if weight_variable else None),
                "metrics": metrics,
                "bins": comparison["bins"].GetInt(),
            })
        # Optional physicsnemo-cfd registry metrics: each entry names a
        # registered metric plus the Kratos fields feeding its SEMANTIC keys
        # ("pressure", "velocity", ...). The "fields" sub-block is free-form
        # (keyed by those semantic names), so it is validated manually.
        self.cfd_metric_specs = []
        for i in range(settings["cfd_metrics"].size()):
            entry = settings["cfd_metrics"][i]
            field_pairs = {}
            if entry.Has("fields"):
                pair_defaults = Kratos.Parameters("""{
                    "predicted_variable" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "predicted_location" : "node_historical",
                    "reference_variable" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "reference_location" : "node_historical"
                }""")
                for semantic_name in entry["fields"].keys():
                    pair = entry["fields"][semantic_name]
                    pair.ValidateAndAssignDefaults(pair_defaults)
                    field_pairs[semantic_name] = (
                        (pair["predicted_variable"].GetString(), pair["predicted_location"].GetString()),
                        (pair["reference_variable"].GetString(), pair["reference_location"].GetString()))
                entry.RemoveValue("fields")
            entry.ValidateAndAssignDefaults(Kratos.Parameters("""{
                "name"   : "PLEASE_SPECIFY_METRIC_NAME",
                "domain" : "surface"
            }"""))
            if not field_pairs:
                raise ValueError(
                    f"cfd_metrics entry \"{entry['name'].GetString()}\" needs a \"fields\" block "
                    "mapping the metric's semantic keys to Kratos variables.")
            self.cfd_metric_specs.append(
                (entry["name"].GetString(), entry["domain"].GetString(), field_pairs))

        uncertainty_comparison_defaults = Kratos.Parameters("""{
            "mean_variable"      : "PLEASE_SPECIFY_VARIABLE_NAME",
            "mean_location"      : "node_historical",
            "std_variable"       : "PLEASE_SPECIFY_VARIABLE_NAME",
            "std_location"       : "node_non_historical",
            "reference_variable" : "PLEASE_SPECIFY_VARIABLE_NAME",
            "reference_location" : "node_historical",
            "confidence_z"       : 1.96,
            "metrics"            : ["coverage", "nll"]
        }""")
        self.uncertainty_comparisons = []
        for i in range(settings["uncertainty_comparisons"].size()):
            entry = settings["uncertainty_comparisons"][i]
            entry.ValidateAndAssignDefaults(uncertainty_comparison_defaults)
            metrics = entry["metrics"].GetStringArray()
            for metric in metrics:
                if metric not in CALIBRATION_METRICS:
                    raise ValueError(
                        f"Unsupported calibration metric \"{metric}\". "
                        f"Supported: {', '.join(CALIBRATION_METRICS)}.")
            self.uncertainty_comparisons.append({
                "mean": (entry["mean_variable"].GetString(),
                         entry["mean_location"].GetString()),
                "std": (entry["std_variable"].GetString(),
                        entry["std_location"].GetString()),
                "reference": (entry["reference_variable"].GetString(),
                              entry["reference_location"].GetString()),
                "confidence_z": entry["confidence_z"].GetDouble(),
                "metrics": metrics,
            })

        ensemble_comparison_defaults = Kratos.Parameters("""{
            "member_variables"   : [],
            "member_location"    : "node_non_historical",
            "reference_variable" : "PLEASE_SPECIFY_VARIABLE_NAME",
            "reference_location" : "node_historical",
            "metrics"            : ["crps"]
        }""")
        self.ensemble_comparisons = []
        for i in range(settings["ensemble_comparisons"].size()):
            entry = settings["ensemble_comparisons"][i]
            entry.ValidateAndAssignDefaults(ensemble_comparison_defaults)
            members = entry["member_variables"].GetStringArray()
            # crps/kcrps are scoring rules over the FULL ensemble and cannot be
            # recovered from a (mean, std) pair, so the members are named
            # explicitly rather than reduced first
            if len(members) < 2:
                raise ValueError(
                    f"\"ensemble_comparisons\" needs at least 2 \"member_variables\", got "
                    f"{len(members)}: an ensemble metric is undefined for a single member.")
            metrics = entry["metrics"].GetStringArray()
            for metric in metrics:
                if metric not in ENSEMBLE_METRICS:
                    raise ValueError(
                        f"Unsupported ensemble metric \"{metric}\". "
                        f"Supported: {', '.join(ENSEMBLE_METRICS)}.")
            self.ensemble_comparisons.append({
                "members": [(name, entry["member_location"].GetString()) for name in members],
                "reference": (entry["reference_variable"].GetString(),
                              entry["reference_location"].GetString()),
                "metrics": metrics,
            })

        # spectral comparisons: both fields sampled onto a grid, their
        # azimuthally averaged power spectra compared. The template stays
        # separate from the (empty) default so an unconfigured process never
        # inherits a phantom entry.
        spectral_comparison_defaults = Kratos.Parameters("""{
            "predicted_variable"       : "PLEASE_SPECIFY_VARIABLE_NAME",
            "predicted_location"       : "node_historical",
            "reference_variable"       : "PLEASE_SPECIFY_VARIABLE_NAME",
            "reference_location"       : "node_historical",
            "grid_shape"               : [32, 32, 2],
            "bounding_box"             : [],
            "squeeze_axis"             : 2,
            "high_wavenumber_fraction" : 0.5,
            "metrics"                  : ["power_spectrum_relative_l2", "high_wavenumber_energy_ratio"]
        }""")
        self.spectral_comparisons = []
        for i in range(settings["spectral_comparisons"].size()):
            entry = settings["spectral_comparisons"][i]
            entry.ValidateAndAssignDefaults(spectral_comparison_defaults)
            metrics = entry["metrics"].GetStringArray()
            for metric in metrics:
                if metric not in SPECTRAL_METRICS:
                    raise ValueError(
                        f"Unsupported spectral metric \"{metric}\". "
                        f"Supported: {', '.join(SPECTRAL_METRICS)}.")
            grid_shape = tuple(int(n) for n in entry["grid_shape"].GetVector())
            if len(grid_shape) != 3:
                raise ValueError(f"\"grid_shape\" must have three entries, got {grid_shape}.")
            box = entry["bounding_box"].GetVector()
            if len(box) not in (0, 6):
                raise ValueError("\"bounding_box\" must be empty or [x0,y0,z0,x1,y1,z1].")
            squeeze_axis = entry["squeeze_axis"].GetInt()
            if squeeze_axis not in (-1, 0, 1, 2):
                raise ValueError(
                    f"\"squeeze_axis\" must be -1 (off), 0, 1 or 2, got {squeeze_axis}.")
            self.spectral_comparisons.append({
                "predicted": (entry["predicted_variable"].GetString(),
                              entry["predicted_location"].GetString()),
                "reference": (entry["reference_variable"].GetString(),
                              entry["reference_location"].GetString()),
                "grid_shape": grid_shape,
                "bounding_box": ((list(box[:3]), list(box[3:])) if len(box) == 6 else None),
                "squeeze_axis": None if squeeze_axis == -1 else squeeze_axis,
                "high_wavenumber_fraction": entry["high_wavenumber_fraction"].GetDouble(),
                "metrics": metrics,
            })

        self.output_interval = settings["output_interval"].GetInt()
        self.output_file = Path(settings["output_file"].GetString())
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

        self.history = []

    def ExecuteFinalizeSolutionStep(self) -> None:
        step = self.model_part.ProcessInfo[Kratos.STEP]
        if step % self.output_interval != 0:
            return
        self.Evaluate()

    def Evaluate(self) -> None:
        record = {
            "STEP": self.model_part.ProcessInfo[Kratos.STEP],
            "TIME": self.model_part.ProcessInfo[Kratos.TIME],
        }
        for comparison in self.comparisons:
            predicted = self._Gather(*comparison["predicted"])
            reference = self._Gather(*comparison["reference"])
            weights = (self._Gather(*comparison["weight"])
                       if comparison["weight"] is not None else None)
            values = ComputeMetricValues(predicted, reference, comparison["metrics"], weights,
                                         bins=comparison["bins"])

            name = f"{comparison['predicted'][0]}_vs_{comparison['reference'][0]}"
            record[name] = values
            Kratos.Logger.PrintInfo(
                "ValidationMetricsProcess",
                f"step {record['STEP']}: {name}: " +
                ", ".join(f"{k}={v:.6e}" for k, v in values.items()))
        for comparison in self.uncertainty_comparisons:
            mean = self._Gather(*comparison["mean"])
            std = self._Gather(*comparison["std"])
            reference = self._Gather(*comparison["reference"])
            values = ComputeCalibrationMetricValues(
                mean, std, reference, comparison["metrics"], comparison["confidence_z"])

            name = f"calibration_{comparison['mean'][0]}_vs_{comparison['reference'][0]}"
            record[name] = values
            Kratos.Logger.PrintInfo(
                "ValidationMetricsProcess",
                f"step {record['STEP']}: {name}: " +
                ", ".join(f"{k}={v:.6e}" for k, v in values.items()))

        for comparison in self.ensemble_comparisons:
            torch = torch_bridge._TryImportTorch()
            ensemble = torch.stack([self._Gather(*member) for member in comparison["members"]])
            reference = self._Gather(*comparison["reference"])
            values = ComputeEnsembleMetricValues(ensemble, reference, comparison["metrics"])

            name = f"ensemble_{comparison['reference'][0]}"
            record[name] = values
            Kratos.Logger.PrintInfo(
                "ValidationMetricsProcess",
                f"step {record['STEP']}: {name}: " +
                ", ".join(f"{k}={v:.6e}" for k, v in values.items()))

        for comparison in self.spectral_comparisons:
            from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
            torch = torch_bridge._TryImportTorch()

            bounding_box = comparison["bounding_box"]
            if bounding_box is None:
                bounding_box = grid_bridge.ComputeBoundingBox(self.model_part)
            grids = []
            for field in (comparison["predicted"], comparison["reference"]):
                grid, _ = grid_bridge.SampleFieldsOnGrid(
                    self.model_part, [field], comparison["grid_shape"], bounding_box)
                if comparison["squeeze_axis"] is not None:
                    grid = grid.mean(axis=1 + comparison["squeeze_axis"])
                grids.append(torch.from_numpy(grid))
            values, spectra = ComputeSpectralMetricValues(
                grids[0], grids[1], comparison["metrics"],
                comparison["high_wavenumber_fraction"])

            name = f"spectrum_{comparison['predicted'][0]}_vs_{comparison['reference'][0]}"
            record[name] = dict(values, spectra=spectra)
            Kratos.Logger.PrintInfo(
                "ValidationMetricsProcess",
                f"step {record['STEP']}: {name}: " +
                ", ".join(f"{k}={v:.6e}" for k, v in values.items()))

        if self.cfd_metric_specs:
            from KratosMultiphysics.PhysicsNeMoApplication.bridges import cfd_bridge
            for name, domain, field_pairs in self.cfd_metric_specs:
                ground_truth = {}
                predictions = {}
                for semantic_name, (predicted, reference) in field_pairs.items():
                    predictions[semantic_name] = self._Gather(*predicted).numpy()
                    ground_truth[semantic_name] = self._Gather(*reference).numpy()
                value = cfd_bridge.EvaluateCfdMetrics(
                    [(name, domain)], ground_truth, predictions)[name]
                record[f"cfd_{name}"] = value
                Kratos.Logger.PrintInfo(
                    "ValidationMetricsProcess",
                    f"step {record['STEP']}: cfd {name} ({domain}): {value:.6e}")
        self.history.append(record)

    def _Gather(self, variable_name: str, data_location: str):
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        tensor_adaptor = GetTensorAdaptor(self.model_part, data_location, variable)
        return torch_bridge.KratosTensorToTorch(tensor_adaptor).clone()

    def ExecuteFinalize(self) -> None:
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(self.output_file, "w") as f:
            json.dump(self.history, f, indent=4)
        Kratos.Logger.PrintInfo(
            "ValidationMetricsProcess",
            f"Wrote {len(self.history)} record(s) to \"{self.output_file}\".")
