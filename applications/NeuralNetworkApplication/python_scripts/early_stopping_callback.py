import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.callback_class import CallbackClass
import tensorflow.keras.callbacks as callbacks

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EarlyStoppingCallback(settings)


class EarlyStoppingCallback(CallbackClass):

    def __init__(self, settings):

        default_settings = KM.Parameters("""{
            "monitor"               : "val_loss",
            "min_delta"             : 0,
            "patience"              : 0,
            "verbose"               : 0,
            "mode"                  : "auto",
            "baseline"              : "",
            "restore_best_weights" : false
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.monitor = settings["monitor"].GetString()
        self.min_delta = settings["min_delta"].GetDouble()
        self.patience = settings["patience"].GetInt()
        self.verbose = settings["verbose"].GetInt()
        self.mode = settings["mode"].GetString()
        self.restore_best_weights = settings["restore_best_weights"].GetBool()
        if not settings["baseline"].GetString() == "":
            self.baseline = settings["baseline"].GetDouble()
        else:
            self.baseline = None

    def Build(self, hp = None):
        return callbacks.EarlyStopping(monitor = self.monitor, min_delta = self.min_delta, patience = self.patience,
        verbose = self.verbose, mode = self.mode, baseline = self.baseline, restore_best_weights = self.restore_best_weights)