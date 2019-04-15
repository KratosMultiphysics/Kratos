import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

import os

## This process imports the topographic model part and initializes the computing one
class ModelImportUtilities:
    def __init__(self, model, settings):
        default_settings = KM.Parameters("""
        {
            "input_type"                  : "mdpa",
            "input_filename"              : "unknown_name",
            "topographic_model_part_name" : "topographic_model_part",
            "computing_model_part_name"   : "computing_model_part",
            "use_different_model_parts"   : true,
            "create_regular_grid"         : true
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)

        self.settings = settings
        self.model = model

        self.topographic_model_part = self._GetOrCreateModelPart(settings["topographic_model_part_name"].GetString())
        self.computing_model_part = self._GetOrCreateModelPart(settings["computing_model_part_name"].GetString())


    def ImportModelPart(self):
        self._ImportTopographicModelPart()
        self._InitializeComputingModelPart()


    def _ImportTopographicModelPart(self):
        KM.Logger.PrintInfo("::[ModelImportUtilities]::", "Reading model part.")
        input_type = self.settings["input_type"].GetString()

        if (input_type == "mdpa"):
            problem_path = os.getcwd()
            input_filename = self.settings["input_filename"].GetString()

            # Setting some mdpa-import-related flags
            import_flags = KM.ModelPartIO.READ
            skip_timer = True
            import_flags = KM.ModelPartIO.SKIP_TIMER|import_flags

            # Import model part from mdpa file.
            KM.Logger.PrintInfo("::[ModelImportUtilities]::", "Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa")
            KM.ModelPartIO(input_filename, import_flags).ReadModelPart(self.topographic_model_part)
            KM.Logger.PrintInfo("::[ModelImportUtilities]::", "Finished reading model part from mdpa file.")

        else:
            raise Exception("Other model part input options are not implemented.")

        KM.Logger.PrintInfo("ModelPart", self.topographic_model_part)
        KM.Logger.PrintInfo("::[PythonSolver]:: ", "Finished reading model part.")

    def _InitializeComputingModelPart(self):
        if self.settings["use_different_model_parts"].GetBool():
            if self.settings["create_regular_grid"].GetBool():
                raise Exception("Regular grid not yet implemented.")
            else:
                SW.ReplicateModelPartUtility(
                    self.topographic_model_part,
                    self.computing_model_part
                ).Replicate()
        else:
            self.computing_model_part = self.topographic_model_part

    def _GetOrCreateModelPart(self, name):
        if self.model.HasModelPart(name):
            model_part = self.model[name]
        else:
            model_part = self.model.CreateModelPart(name)
        return model_part
