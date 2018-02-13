import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import hdf5_output

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PartitionedMultipleMeshTemporalOutputProcess(Model, settings["Parameters"])

class PartitionedMultipleMeshTemporalOutputProcess(KratosMultiphysics.Process):
    """A process for writing partitioned simulation results for multiple meshes to HDF5."""

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "please_specify_model_part_name",
                "file_output_settings" : {},
                "model_part_output_settings" : {},
                "results_settings" : {},
                "output_time_settings": {}
            }
            """)
        settings.ValidateAndAssignDefaults(default_settings)
        model_part = Model[settings["model_part_name"].GetString()]

        hdf5_file_factory = hdf5_output.HDF5ParallelFileFactory(settings["file_output_settings"])
        model_part_output = hdf5_output.PartitionedModelPartOutput(settings["model_part_output_settings"])
        results_output = hdf5_output.PartitionedNodalResultsOutput(settings["results_settings"])

        self._static_output = hdf5_output.StaticOutputProcess(model_part, hdf5_file_factory)
        self._static_output.AddOutput(model_part_output)
        self._static_output.AddOutput(results_output)
        
        self._temporal_output = hdf5_output.TemporalOutputProcess(
            model_part, hdf5_file_factory, settings["output_time_settings"])
        self._temporal_output.AddOutput(model_part_output)
        self._temporal_output.AddOutput(results_output)

    def ExecuteBeforeSolutionLoop(self):
        self._static_output.Execute()
        self._temporal_output.ExecuteBeforeSolutionLoop()

    def ExecuteFinalizeSolutionStep(self):
        self._temporal_output.ExecuteFinalizeSolutionStep()
