'''HDF5 Xdmf operations.

license: HDF5Application/license.txt

Main authors:
    Philipp Bucher
    Michael Andre
'''

# --- Kratos Imports ---
from typing import Any
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5


# in case the h5py-module is not installed (e.g. on clusters) we don't want it to crash the simulation!
# => in such a case the xdmf can be created manually afterwards locally
try:
    import KratosMultiphysics.HDF5Application.xdmf_utils
    from KratosMultiphysics.HDF5Application.core.dataset_generator import DatasetGenerator
    from KratosMultiphysics.HDF5Application.core.dataset_generator import SingleMeshMultiFileSameDatasetGenerator
    from KratosMultiphysics.HDF5Application.core.dataset_generator import GenericDatasetGenerator
    from KratosMultiphysics.HDF5Application.core.dataset_generator import GetDataSetPatterns
    from KratosMultiphysics.HDF5Application.core.dataset_generator import HasTags
    from KratosMultiphysics.HDF5Application.xdmf_utils import WriteDatasetsToXdmf
except ImportError:
    # if we failed to import, then assign a dummy function that does nothing.
    WriteMultifileTemporalAnalysisToXdmf = lambda *args: None
    warn_msg = "XDMF-Writing is not available,\nOnly HDF5 files are written"
    KratosMultiphysics.Logger.PrintWarning(__name__, warn_msg)


class XdmfOutput(KratosMultiphysics.Operation):
    '''Output that creates the xdmf-file for the given h5-files.'''
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "dataset_pattern"         : "",
            "temporal_tag_value_index": 0,
            "single_mesh_datasets"    : true
        }""")

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 parameters: KratosMultiphysics.Parameters,
                 file: KratosHDF5.HDF5File):
        super().__init__()
        self.__model_part = model_part
        self.__file = file
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.__dataset_pattern = parameters["dataset_pattern"].GetString()
        self.__is_single_mesh_datasets = parameters["single_mesh_datasets"].GetBool()
        self.__temporal_tag_value_index = parameters["temporal_tag_value_index"].GetInt()
        file_name_pattern, dataset_pattern = GetDataSetPatterns(self.__dataset_pattern)

        if self.__is_single_mesh_datasets:
            if not HasTags(file_name_pattern, {"<time>": float, "<step>": int}):
                raise RuntimeError(f"Single mesh datasets require at least one tag to be present in the filename.")
            if HasTags(dataset_pattern, {"<time>": float, "<step>": int}):
                raise RuntimeError(f"Single mesh datasets require cannot have any tags in the dataset prefix.")

    def Execute(self) -> None:
        self.__model_part.GetCommunicator().GetDataCommunicator().Barrier()
        # write xdmf only on one rank!
        if self.__model_part.GetCommunicator().MyPID() == 0:
            dataset_generator: DatasetGenerator
            if self.__is_single_mesh_datasets:
                dataset_generator = SingleMeshMultiFileSameDatasetGenerator(self.__dataset_pattern, self.__temporal_tag_value_index)
            else:
                dataset_generator = GenericDatasetGenerator(self.__dataset_pattern, self.__temporal_tag_value_index)
            WriteDatasetsToXdmf(dataset_generator, self.__model_part.FullName() + ".xdmf")

def Create(settings):
    '''Return an operation specified by the setting's 'operation_type'.

    This method is normally not used directly, but rather it is imported
    in core.operations.model_part.Create using the 'module_name' setting.
    '''
    operation_type = settings['operation_type'].GetString()
    if operation_type == 'xdmf_output':
        return XdmfOutput
    else:
        raise ValueError(f'"operation_type" has invalid value "{operation_type}"')
