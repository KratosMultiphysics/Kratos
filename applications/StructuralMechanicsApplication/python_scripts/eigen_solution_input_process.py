from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EigenSolutionInputProcess(Model, settings["Parameters"])

class EigenSolutionInputProcess(KratosMultiphysics.Process):
    """A process for reading eigenvalue and eigenvector results."""

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
            {
                "help"            : "A process for reading eigenvalue and eigenvector results.",
                "model_part_name" : "PLEASE_SPECIFY_MODEL_PART",
                "file_settings" : {
                },
                "prefix" : "/EigenResults"
            }
            """)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)
        self._model_part = Model[self.settings["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        hdf5_file = self._GetFile()
        prefix = self.settings["prefix"].GetString()
        KratosHDF5.ReadDataValueContainer(hdf5_file, prefix, self._model_part.ProcessInfo)
        nodal_io_settings = KratosMultiphysics.Parameters("""
            {
                "list_of_variables": ["EIGENVECTOR_MATRIX"],
                "prefix" : ""
            }
            """)
        nodal_io_settings["prefix"].SetString(prefix)
        nodal_data_value_io = KratosHDF5.HDF5NodalDataValueIO(nodal_io_settings, hdf5_file)
        nodal_data_value_io.ReadNodalResults(self._model_part.Nodes, self._model_part.GetCommunicator())

    def _GetFile(self):
        return KratosHDF5.HDF5FileSerial(self.settings["file_settings"])
