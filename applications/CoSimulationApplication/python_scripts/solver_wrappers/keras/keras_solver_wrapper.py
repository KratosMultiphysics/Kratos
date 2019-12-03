from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import numpy as np
# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from tensorflow.keras.models import load_model


def Create(settings, solver_name):
    return KerasSolverWrapper(settings, solver_name)

class KerasSolverWrapper(CoSimulationSolverWrapper):
    """ This class implements a wrapper for an SDof solver to be used in CoSimulation
    """
    def __init__(self, settings, solver_name):
        super(KerasSolverWrapper, self).__init__(settings, solver_name)

        input_file_name = self.settings["solver_wrapper_settings"]["input_file"].GetString()
        
        self.keras = load_model(input_file_name)
        self.mp = self.model.CreateModelPart("Keras")
        self.mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.mp.AddNodalSolutionStepVariable(KM.FORCE)
        mdpa_file_name = self.settings["solver_wrapper_settings"]["mdpa_file_name"].GetString()
        KM.ModelPartIO(mdpa_file_name).ReadModelPart(self.mp)

        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 2

        print("Used keras model is",self.keras.summary())

    def AdvanceInTime(self, current_time):
        return 0.0

    def SolveSolutionStep(self):

        print("Receiving data into the keras model")

        loads = self.GetInterfaceData("load").GetData()
        loads = loads.reshape(loads.shape[0],1)

        disps = self.keras(loads.T)

        disp = disps.numpy()

        self.GetInterfaceData("disp").SetData(disp[0])

        print("Predicted data using Keras model")



