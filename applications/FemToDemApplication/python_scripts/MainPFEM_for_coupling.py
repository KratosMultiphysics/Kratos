from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import time as timer
import KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis as PfemFluidDynamicsAnalysis
import os

def Wait():
    input("Alejandro -> Press Something")

class MainPFEM_for_coupling_solution(PfemFluidDynamicsAnalysis.PfemFluidDynamicsAnalysis):
    """
    The derived class for the PfemFluidDynamicsAnalysis
    """

    def __init__(self, model, parameters):
        """
        The constructor of the MainPFEM_for_coupling_solution-Object.
        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        parameters -- The ProjectParameters used
        """

        # We change the name of the inputs in order to be separated from FEMDEM
        problem_name = parameters["problem_data"]["problem_name"].GetString()
        parameters["problem_data"]["problem_name"].SetString("PFEM" + problem_name)
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("PFEM" + problem_name)

        super(MainPFEM_for_coupling_solution, self).__init__(model, parameters)

    def Initialize(self):
        super(MainPFEM_for_coupling_solution, self).Initialize()
        for node in self.main_model_part.Nodes:
            print(node.Id)
        Wait()