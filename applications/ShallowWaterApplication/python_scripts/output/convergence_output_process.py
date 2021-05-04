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
        super().__init__()

        default_settings = KM.Parameters("""
            {
                "model_part_name"            : "model_part",
                "file_name"                  : "output_file",
                "printing_times"             : [],
                "analysis_label"             : "label",
                "analysis_attributes"        : {},
                "convergence_variables_list" : [],
                "low_corner"                 : [],
                "high_corner"                : []
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]

        self.f = h5py.File(self.settings["file_name"].GetString() + ".hdf5", 'a') # 'a' means append mode
        self.f.flush() # Since this process will reuse the existing files, we need to ensure the file construction is finished. Otherwise problems may arise if the first analysis didn't finish.

        self.variables = [KM.KratosGlobals.GetVariable(v) for v in self.settings["convergence_variables_list"].GetStringArray()]

        # Initialize output control variables
        self.printing_times = self.settings["printing_times"].GetVector()
        self.is_printed = [False] * len(self.printing_times)

        if self.settings["low_corner"].GetVector().Size() == 0:
            self.integrate_over_all_the_domain = True
        else:
            self.integrate_over_all_the_domain = False

    def ExecuteBeforeSolutionLoop(self):
        self.dset = self._GetDataset()
        self.start_time = time.time()

    def IsOutputStep(self):
        """Check if the current time step is near enough to the specified printing times."""
        time = self.model_part.ProcessInfo.GetValue(KM.TIME)
        for i in range(len(self.printing_times)):
            if time >= self.printing_times[i] and not self.is_printed[i]:
                self.is_printed[i] = True
                return True
        return False

    def PrintOutput(self):
        self._WriteAverageError()

    def Check(self):
        for variable in self.variables:
            if not isinstance(variable, KM.DoubleVariable):
                raise Exception("This process is expecting only double or component variables")
        
        low_corner = self.settings["low_corner"].GetVector()
        high_corner = self.settings["high_corner"].GetVector()
        if not low_corner.Size() == high_corner.Size():
            raise Exception("The low and high corners does not have the same dimension")

        if low_corner.Size() == 0:
            pass
        elif low_corner.Size() == 2:
            self.settings["low_corner"].Append(0.0)
            self.settings["high_corner"].Append(0.0)
        elif low_corner.Size() == 3:
            pass
        else:
            raise Exception("The corners must be specified with 2 or 3 coordinates")

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
            ("time", np.float),
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
            self.model_part.ProcessInfo[KM.TIME],
            elapsed_time]

        if not self.integrate_over_all_the_domain:
            low_corner = KM.Point(self.settings["low_corner"].GetVector())
            high_corner = KM.Point(self.settings["high_corner"].GetVector())

        for variable in self.variables:
            if self.integrate_over_all_the_domain:
                value = SW.ShallowWaterUtilities().ComputeL2NormNonHistorical(self.model_part, variable)
            else:
                value = SW.ShallowWaterUtilities().ComputeL2NormNonHistorical(self.model_part, variable, low_corner, high_corner)
            case_data.append(value)

        case_idx = self.dset.len()
        self.dset.resize((case_idx+1,))
        self.dset[case_idx] = tuple(case_data)
