# Importing the Kratos Library
import KratosMultiphysics

import h5py
import numpy as np
import datetime as dt

class Hdf5OutputTool:
    def __init__(self, model, settings):

        default_output_settings = KratosMultiphysics.Parameters("""
            {
                "file_name" : "output_file",
                "case_id"   : 0,
                "framework" : "please_specify_it",
                "h"         : 0.0,
                "courant"   : 0.0       
            }
            """
            )

        settings["output_parameters"].ValidateAndAssignDefaults(default_output_settings)
        self.settings = settings

        self.model_part = model[settings["model_part_name"].GetString()]

        self.f = h5py.File(self.settings["output_parameters"]["file_name"].GetString() + ".hdf5", 'a') # 'a' means append mode

        case_name = "case_{:03d}".format(self.settings["output_parameters"]["case_id"].GetInt())
        self.case = self.f.create_group(case_name)

    def WriteDiscretizationAttributes(self):
        d = self.case.create_group("discretization")
        d.attrs["framework"] = self.settings["output_parameters"]["framework"].GetString()
        d.attrs["h"] = self.settings["output_parameters"]["h"].GetDouble()
        d.attrs["nodes"] = self.model_part.Nodes.__len__()
        d.attrs["courant"] = self.settings["output_parameters"]["courant"].GetDouble()
        d.attrs["dt"] = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

    def WriteBodyForceAttributes(self):
        b = self.case.create_group("body_force")
        for attr, value in self.settings["benchmark_parameters"].items():
            b.attrs[attr] = value.GetDouble()
    
    def WriteAverageRelativeError(self, err):
        self.case.attrs["relative_error"] = err
