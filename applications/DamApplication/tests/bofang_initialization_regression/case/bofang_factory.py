import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam


def Factory(settings, Model):
    """Instantiate the raw C++ DamBofangConditionTemperatureProcess.

    A dedicated factory is used (instead of the production
    impose_reservoir_temperature_condition_process wrapper) so that the C++
    process lifecycle is exercised directly by the analysis. The production
    wrapper does not forward Process::ExecuteInitialize(), which would mask the
    lifecycle behaviour under investigation.
    """
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    params = settings["Parameters"]
    model_part = Model[params["model_part_name"].GetString()]
    return KratosDam.DamBofangConditionTemperatureProcess(model_part, params)
