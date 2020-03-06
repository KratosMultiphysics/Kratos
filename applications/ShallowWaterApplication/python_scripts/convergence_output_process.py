# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

# Other imports
import h5py
import time
import numpy as np

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ConvergenceOutputProcess(model, settings["Parameters"])

class ConvergenceOutputProcess(KM.Process):

    def __init__(self, model, settings):
        super(ConvergenceOutputProcess, self).__init__()

        default_settings = KM.Parameters("""
            {
                "model_part_name"            : "model_part",
                "file_name"                  : "output_file",
                "analysis_label"             : "label",
                "analysis_attributes"        : {},
                "convergence_variables_list" : [],
                "weighting_variable"         : "NODAL_AREA"
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]

        self.f = h5py.File(self.settings["file_name"].GetString() + ".hdf5", 'a') # 'a' means append mode
        self.f.flush() # Since this process will reuse the existing files, we need to ensure the file construction is finished. Otherwise problems may arise if the first analysis didn't finish.

        self.variables = [KM.KratosGlobals.GetVariable(v) for v in self.settings["convergence_variables_list"].GetStringArray()]
        self.weight_variable = KM.KratosGlobals.GetVariable(settings["weighting_variable"].GetString())

    def ExecuteBeforeSolutionLoop(self):
        if self.weight_variable == KM.NODAL_AREA:
            KM.CalculateNonHistoricalNodalAreaProcess(self.model_part).Execute()
        self.dset = self._GetDataset()
        self.start_time = time.time()

    @staticmethod
    def IsOutputStep():
        return False

    @staticmethod
    def PrintOutput():
        pass

    def ExecuteFinalize(self):
        self._WriteAverageError()

    def Check(self):
        for variable in self.variables:
            if not isinstance(variable, KM.DoubleVariable) and not isinstance(variable, KM.Array1DComponentVariable):
                raise Exception("This process is expecting only double or component variables")

        if not isinstance(self.weight_variable, KM.DoubleVariable):
            raise Exception("The weighting variable is expected to be a double variable")

    def _WriteAttributes(self, dset):
        for key, param in self.settings["analysis_attributes"].items():
            if param.IsBool():
                value = param.GetBool()
            elif param.IsInt():
                value = param.GetInt()
            elif param.IsDouble():
                value = param.GetDouble()
            elif param.IsString():
                value = param.GetString()
            else:
                msg = "Unknown type for " + str(param)
                msg = msg.rstrip() + " with key : \"" + str(key) + "\""
                raise Exception(msg)
            dset.attrs[key] = value

    def _CheckAttributes(self, attributes):
        match = None
        for key, param in self.settings["analysis_attributes"].items():
            match = False
            if param.IsBool():
                value = param.GetBool()
            elif param.IsInt():
                value = param.GetInt()
            elif param.IsDouble():
                value = param.GetDouble()
            elif param.IsString():
                value = param.GetString()
            else:
                msg = "Unknown type for " + str(param)
                msg = msg.rstrip() + " with key : \"" + str(key) + "\""
                raise Exception(msg)

            if attributes.get(key) is not None:
                if attributes[key] == value:
                    match = True
            if not match:
                break

        return match

    def _GetDataset(self):
        for name, data in self.f.items():
            if self._CheckAttributes(data.attrs):
                return data

        dset_name = "analysis_{:03d}".format(self.f.items().__len__())
        header = self._GetHeaderDtype()
        dset = self.f.create_dataset(dset_name, (0,), maxshape=(None,), chunks=True, dtype=header)
        self._WriteAttributes(dset)
        return dset

    def _GetHeaderDtype(self):
        header = [
            ("label", h5py.special_dtype(vlen=str)),
            ("num_nodes", np.uint32),
            ("num_elems", np.uint32),
            ("time_step", np.float),
            ("computational_time", np.float)]

        for variable in self.variables:
            header.append((variable.Name(), np.float))

        return np.dtype(header)

    def _WriteAverageError(self):
        elapsed_time = time.time() - self.start_time
        case_data = [
            self.settings["analysis_label"].GetString(),
            self.model_part.NumberOfNodes(),
            self.model_part.NumberOfElements(),
            self.model_part.ProcessInfo[KM.DELTA_TIME],
            elapsed_time]

        for variable in self.variables:
            value = SW.ShallowWaterUtilities().RootMeanSquareNonHistorical(variable, self.weight_variable, self.model_part.Nodes)
            case_data.append(value)

        case_idx = self.dset.len()
        self.dset.resize((case_idx+1,))
        self.dset[case_idx] = tuple(case_data)
