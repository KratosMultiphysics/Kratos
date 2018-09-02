import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import hdf5_io
import hdf5_defaults


class FactoryHelper:

    def Execute(self, settings, Model):
        """Return objects needed for constructing a temporal output process."""
        if not isinstance(settings, KratosMultiphysics.Parameters):
            raise Exception("expected input shall be a Parameters object, encapsulating a json string")
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.hdf5_default_settings)
        settings = settings.Clone()
        settings.ValidateAndAssignDefaults(default_settings)
        model_part = Model[settings["model_part_name"].GetString()]
        hdf5_file_factory = self.FileFactory(settings["file_settings"])
        model_part_output = self.ModelPartOutput(settings["model_part_output_settings"])
        nodal_solution_step_output = self.NodalSolutionStepDataOutput(settings["nodal_solution_step_data_settings"])
        element_data_value_output = hdf5_io.ElementDataValueOutput(settings["element_data_value_settings"])
        nodal_data_value_output = hdf5_io.NodalDataValueOutput(settings["nodal_data_value_settings"])
        temporal_output_process = hdf5_io.TemporalOutputProcess(
            model_part, hdf5_file_factory, settings["output_time_settings"], [model_part_output, nodal_solution_step_output, element_data_value_output, nodal_data_value_output])
        return (temporal_output_process, model_part_output, [nodal_solution_step_output, element_data_value_output, nodal_data_value_output])


class SerialFactory:

    def FileFactory(self, file_settings):
        return hdf5_io.HDF5SerialFileFactory(file_settings)

    def ModelPartOutput(self, model_part_output_settings):
        return hdf5_io.ModelPartOutput(model_part_output_settings)

class ParallelFactory:

    def FileFactory(self, file_settings):
        return hdf5_io.HDF5ParallelFileFactory(file_settings)

    def ModelPartOutput(self, model_part_output_settings):
        return hdf5_io.PartitionedModelPartOutput(model_part_output_settings)


class ResultsFactory:

    def NodalSolutionStepDataOutput(self, nodal_solution_step_data_settings):
        return hdf5_io.NodalSolutionStepDataOutput(nodal_solution_step_data_settings)

class PrimalResultsFactory:

    def __init__(self, alpha_bossak):
        self.alpha_bossak = alpha_bossak

    def NodalSolutionStepDataOutput(self, nodal_solution_step_data_settings):
        return hdf5_io.PrimalBossakOutput(nodal_solution_step_data_settings, self.alpha_bossak)

class TemporalOutputFactoryHelper(FactoryHelper, SerialFactory, ResultsFactory):
    pass


class PrimalOutputFactoryHelper(FactoryHelper, SerialFactory, PrimalResultsFactory):

    def __init__(self, alpha_bossak):
        PrimalResultsFactory.__init__(self, alpha_bossak)


class PartitionedTemporalOutputFactoryHelper(FactoryHelper, ParallelFactory, ResultsFactory):
    pass


class PartitionedPrimalOutputFactoryHelper(FactoryHelper, ParallelFactory, PrimalResultsFactory):

    def __init__(self, alpha_bossak):
        PrimalResultsFactory.__init__(self, alpha_bossak)
