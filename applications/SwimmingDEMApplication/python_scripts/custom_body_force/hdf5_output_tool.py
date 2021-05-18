# Importing the Kratos Library
import KratosMultiphysics

import h5py
import numpy as np

class Hdf5OutputTool:
    def __init__(self, model, settings, ):

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"   : "main_model_part",
                "file_name"         : "output_file",
                "vortex_dataset_id" : 0,
                "framework"         : "please_specify_it"    
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)
        self.settings = settings

        self.model_part = model[self.settings["model_part_name"].GetString()]

        case_dtype = np.dtype([("framework", h5py.special_dtype(vlen=str)),
                               ("num_nodes", np.uint32),
                               ("num_elems", np.uint32),
                               ("time_step", np.float),
                               ("rel_error", np.float)])

        self.f = h5py.File(self.settings["file_name"].GetString() + ".hdf5", 'a') # 'a' means append mode

        vortex_name = "vortex_{:03d}".format(self.settings["vortex_dataset_id"].GetInt())
        if vortex_name in self.f:
            self.vortex_dset = self.f[vortex_name]
            self.dataset_is_initialized = True
        else:
            self.vortex_dset = self.f.create_dataset(vortex_name, (0,), maxshape=(None,), chunks=True, dtype=case_dtype)
            self.dataset_is_initialized = False

    def WriteBodyForceAttributes(self, body_force_settings):
        for attr, value in body_force_settings.items():
            if self.dataset_is_initialized:
                if self.vortex_dset.attrs[attr] != value.GetDouble():
                    self.vortex_dset.attrs["Warning"] = "There are several vortex definitions in this dataset"
            self.vortex_dset.attrs[attr] = value.GetDouble()

    def WriteAverageRelativeError(self, rel_err):
        case_data = (self.settings["framework"].GetString(),
                     self.model_part.Nodes.__len__(),
                     self.model_part.Elements.__len__(),
                     self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME],
                     rel_err)

        case_idx = self.vortex_dset.len()
        self.vortex_dset.resize((case_idx+1,))
        self.vortex_dset[case_idx] = case_data
