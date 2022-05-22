import KratosMultiphysics as Kratos
import KratosMultiphysics.LinearSolversApplication as KratosLSA
import KratosMultiphysics.HDF5Application as KratosHDF5

from pathlib import Path

def Factory(settings, Model):
    if(type(settings) != Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeAdjointStabilizationCoefficients(Model, settings["Parameters"])

class ComputeAdjointStabilizationCoefficients(Kratos.Process):
    def __init__(self, model, settings):
        Kratos.Process.__init__(self)

        default_settings = Kratos.Parameters("""
        {
            "model_part_name"     : "PLEASE_SPECIFIY_MODEL_PART_NAME",
            "output_folder_name"  : "PLEASE_SPECIFY_OUTPUT_FOLDER_NAME",
            "output_variable_name": "PLEASE_SPECIFY_SCALAR_VARIABLE_NAME",
            "echo_level"          : 0
        }
        """)

        settings.ValidateAndAssignDefaults(default_settings)

        self.model = model

        max_svd_process_parameters = Kratos.Parameters("""{
            "model_part_name"     : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "input_variable_name" : "PLEASE_SPECIFY_MATRIX_VARIABLE",
            "output_variable_name": "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "echo_level"          : 0,
            "container_type"      : "elements"
        }""")
        max_svd_process_parameters["model_part_name"].SetString(settings["model_part_name"].GetString())
        max_svd_process_parameters["input_variable_name"].SetString("PRIMAL_STEADY_RESIDUAL_FIRST_DERIVATIVES")
        max_svd_process_parameters["output_variable_name"].SetString(settings["output_variable_name"].GetString())
        max_svd_process_parameters["echo_level"].SetInt(settings["echo_level"].GetInt())

        self.max_svd_process = KratosLSA.MaxSingularValueDecompositionProcess(model, max_svd_process_parameters)
        self.output_folder_name = settings["output_folder_name"].GetString()
        self.model_part_name = settings["model_part_name"].GetString()
        self.output_variable_name = settings["output_variable_name"].GetString()
        self.echo_level = settings["echo_level"].GetInt()

    def ExecuteInitializeSolutionStep(self):
        model_part = self.model[self.model_part_name]
        step = model_part.ProcessInfo[Kratos.STEP]

        file_name_path = Path(self.output_folder_name) / "{:s}-{:d}.h5".format(self.model_part_name, step)

        hdf5_file_settings = Kratos.Parameters("""{
            "file_name"       : "",
            "file_access_mode": ""
        }""")
        hdf5_file_settings["file_name"].SetString(str(file_name_path))

        hdf5_input_parameters = Kratos.Parameters("""
        {
            "prefix": "/ResultsData",
            "list_of_variables":[]
        }""")
        hdf5_input_parameters["list_of_variables"].SetStringArray([self.output_variable_name])

        # if existing file is found, then read and populate data from the file
        recalculate_coefficients = True
        if file_name_path.is_file():
            try:
                hdf5_file_settings["file_access_mode"].SetString("read_only")
                h5_file = KratosHDF5.HDF5FileSerial(hdf5_file_settings)
                element_io = KratosHDF5.HDF5ElementDataValueIO(hdf5_input_parameters, h5_file)
                element_io.ReadElementResults(model_part.Elements, model_part.GetCommunicator())

                if self.echo_level > 0:
                    Kratos.Logger.PrintInfo(self.__class__.__name__, "Using existing stabilization coefficient calculation at {:s}.".format(str(file_name_path)))
                recalculate_coefficients = False
            except:
                if self.echo_level > 0:
                    Kratos.Logger.PrintInfo(self.__class__.__name__, "Corrupted existing stabilization coefficient calculation at {:s}.".format(str(file_name_path)))

        if recalculate_coefficients:
            file_name_path.parent.mkdir(parents=True, exist_ok=True)
            # no existing data is found (or existing data is corrupt). then calculate and output it to file

            # calculate maximum svd
            self.max_svd_process.Execute()

            # store the values
            hdf5_file_settings["file_access_mode"].SetString("exclusive")
            h5_file = KratosHDF5.HDF5FileSerial(hdf5_file_settings)
            element_io = KratosHDF5.HDF5ElementDataValueIO(hdf5_input_parameters, h5_file)
            element_io.WriteElementResults(model_part.Elements)

            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Stored stabilization coefficient calculation at {:s}.".format(str(file_name_path)))


