# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

# Import applications
import KratosMultiphysics.HDF5Application as KratosHDF5

def Factory(settings: Kratos.Parameters, model: Kratos.Model):
    if(type(settings) != Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EigenvectorSolutionOutputProcess(model, settings["Parameters"])

class EigenvectorSolutionOutputProcess(Kratos.OutputProcess):
    """A process for writing (measurement) eigenvector results."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
        Kratos.OutputProcess.__init__(self)
        
        default_settings = Kratos.Parameters("""
            {
                "help"            : "This process generates a postprocess file in a HDF5 file for eigenvectors",
                "model_part_name" : "PLEASE_SPECIFY_MODEL_PART",
                "file_settings" : {
                },
                "prefix" : "/EigenResults"
            }
            """)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)
        self._model_part = model[self.settings["model_part_name"].GetString()]

    def PrintOutput(self) -> None:
        pass
    
    def ExecuteFinalize(self):

        self._model_part.ProcessInfo[KratosSI.MEASURED_EIGENVALUE_VECTOR] = self._model_part.ProcessInfo[KratosSI.EIGENVALUE_VECTOR]
        #print(self._model_part.GetValue(KratosSI.MEASURED_EIGENVALUE_VECTOR))

        for node in self._model_part.Nodes:
            node: Kratos.Node
            node.SetValue(KratosSI.MEASURED_EIGENVECTOR_MATRIX, node.GetValue(KratosSI.EIGENVECTOR_MATRIX))

        hdf5_file = self._GetFile()
        prefix = self.settings["prefix"].GetString()
        KratosHDF5.WriteDataValueContainer(hdf5_file, prefix, self._model_part.ProcessInfo)

        nodal_io_settings = Kratos.Parameters("""
            {
                "list_of_variables": ["MEASURED_EIGENVECTOR_MATRIX"],
                "prefix" : ""
            }
            """)
        nodal_io_settings["prefix"].SetString(prefix)
        nodal_data_value_io = KratosHDF5.HDF5NodalDataValueIO(nodal_io_settings, hdf5_file)
        nodal_data_value_io.WriteNodalResults(self._model_part.Nodes)

    def _GetFile(self):
        return KratosHDF5.HDF5File(self.settings["file_settings"])
