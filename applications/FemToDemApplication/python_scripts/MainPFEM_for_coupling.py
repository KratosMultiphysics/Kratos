from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import time as timer
import KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis as PfemFluidDynamicsAnalysis
import os
from importlib import import_module

def Wait():
    input("Alejandro -> Press Something")

#============================================================================================================================
class MainPFEM_for_coupling_solution(PfemFluidDynamicsAnalysis.PfemFluidDynamicsAnalysis):
    """
    The derived class for the PfemFluidDynamicsAnalysis
    """
#============================================================================================================================
#============================================================================================================================
    def __init__(self, model, FEM_model_part, parameters):
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

        self.FEM_model_part = FEM_model_part
        super(MainPFEM_for_coupling_solution, self).__init__(model, parameters)

#============================================================================================================================
    def SetCustomGraphicalOutput(self, custom_parameters):
        """This function sets the settings for the graphical
        output
        """
        if custom_parameters.Has("output_configuration"):
            from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_gid_output_process import GiDOutputProcess
            self.output_settings = custom_parameters["output_configuration"]
            self.output_settings["result_file_configuration"].RemoveValue("nodal_results")
            self.output_settings["result_file_configuration"].AddValue("nodal_results", self.project_parameters["output_configuration"]["result_file_configuration"]["nodal_results"])
            self.output_settings["result_file_configuration"].RemoveValue("gauss_point_results")
            self.output_settings["result_file_configuration"].AddValue("gauss_point_results", self.project_parameters["output_configuration"]["result_file_configuration"]["gauss_point_results"])
            return GiDOutputProcess(self.post_process_model_part,
                                    self.problem_name,
                                    self.output_settings)
        else:
            return (KratosMultiphysics.Process())

#============================================================================================================================
    def _CreateSolver(self):
        """Create the solver
        """
        python_module_name = "KratosMultiphysics.FemToDemApplication"
        full_module_name = python_module_name + "." + "pfem_fluid_solver_for_coupling"
        solver_module = import_module(full_module_name)
        solver = solver_module.CreateSolver(self.model, self.FEM_model_part, self.project_parameters["solver_settings"])
        return solver
        
#============================================================================================================================