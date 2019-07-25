from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
import KratosMultiphysics.ChimeraApplication as KratosChimera
import KratosMultiphysics.ChimeraApplication.python_solvers_wrapper_fluid_chimera

class FluidChimeraAnalysis(FluidDynamicsAnalysis):
    '''Main script for fluid chimera simulations using the navier stokes family of python solvers.'''

    def __init__(self,model,parameters):
        # Deprecation warnings
        self.parameters = parameters
        if self.parameters["solver_settings"].Has("chimera_parts"):
            self.chimera_params = self.parameters["solver_settings"]["chimera_parts"].Clone()
            self.parameters["solver_settings"].RemoveValue("chimera_parts")
        else:
            raise Exception("The \"solver_settings\" should have the entry \"chimera_parts\" ")

        self.chimera_echo_lvl = 0
        if self.parameters["solver_settings"].Has("chimera_echo_level"):
            self.chimera_echo_lvl = self.parameters["solver_settings"]["chimera_echo_level"].GetInt()
            self.parameters["solver_settings"].RemoveValue("chimera_echo_level")
        else:
            self.chimera_echo_lvl = self.parameters["solver_settings"]["echo_level"].GetInt()


        if self.parameters["solver_settings"].Has("internal_parts_for_chimera"):
            self.chimera_internal_parts = self.parameters["solver_settings"]["internal_parts_for_chimera"].Clone()
            self.parameters["solver_settings"].RemoveValue("internal_parts_for_chimera")

        self.reformulate_every_step = False
        if self.parameters["solver_settings"].Has("reformulate_chimera_every_step"):
            self.reformulate_every_step = self.parameters["solver_settings"]["reformulate_chimera_every_step"].GetBool()
            self.parameters["solver_settings"].RemoveValue("reformulate_chimera_every_step")
            # Setting reform dofs every step to true "reform_dofs_at_each_step": false,
            if not self.parameters["solver_settings"].Has("reform_dofs_at_each_step"):
                self.parameters["solver_settings"].AddEmptyValue("reform_dofs_at_each_step")
            self.parameters["solver_settings"]["reform_dofs_at_each_step"].SetBool(True)

        # Import parallel modules if needed
        # has to be done before the base-class constuctor is called (in which the solver is constructed)
        if (parameters["problem_data"]["parallel_type"].GetString() == "MPI"):
            raise Exception("MPI-Chimera is not implemented yet")

        super(FluidChimeraAnalysis,self).__init__(model,parameters)

    def Initialize(self):
        super(FluidChimeraAnalysis,self).Initialize()
        self.__SetChimeraInternalPartsFlag()

    def _CreateProcesses(self, parameter_name, initialization_order):

        list_of_processes = super(FluidChimeraAnalysis,self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            if self.parameters["solver_settings"].Has("fluid_solver_settings"):
                self.solver_settings =  self.parameters["solver_settings"]["fluid_solver_settings"]
            else:
                self.solver_settings = self.parameters["solver_settings"]

            main_model_part = self.model[self.solver_settings["model_part_name"].GetString()]
            domain_size = main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            solver_type = self.parameters["solver_settings"]["solver_type"].GetString()

            if domain_size == 2:
                if(solver_type == "Monolithic" or solver_type == "monolithic"):
                    self.chimera_process = KratosChimera.ApplyChimeraProcessMonolithic2d(main_model_part,self.chimera_params)
                elif (solver_type == "fractional_step" or solver_type == "FractionalStep"):
                    self.chimera_process = KratosChimera.ApplyChimeraProcessFractionalStep2d(main_model_part,self.chimera_params)
            else:
                if(solver_type == "Monolithic" or solver_type == "monolithic"):
                    self.chimera_process = KratosChimera.ApplyChimeraProcessMonolithic3d(main_model_part,self.chimera_params)
                elif (solver_type == "fractional_step" or solver_type == "FractionalStep"):
                    self.chimera_process = KratosChimera.ApplyChimeraProcessFractionalStep3d(main_model_part,self.chimera_params)

            self.chimera_process.SetEchoLevel(self.chimera_echo_lvl)
            self.chimera_process.SetReformulateEveryStep(self.reformulate_every_step)

        return list_of_processes

    def __SetChimeraInternalPartsFlag(self):
        for mp_name in self.chimera_internal_parts:
            KratosMultiphysics.VariableUtils().SetFlag(KratosChimera.CHIMERA_INTERNAL_BOUNDARY, True,  self.model[mp_name.GetString()].Nodes)

    def InitializeSolutionStep(self):
        super(FluidChimeraAnalysis,self).InitializeSolutionStep()
        self.chimera_process.ExecuteInitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.chimera_process.ExecuteFinalizeSolutionStep()
        super(FluidChimeraAnalysis,self).FinalizeSolutionStep()


    def _CreateSolver(self):
        return KratosChimera.python_solvers_wrapper_fluid_chimera.CreateSolver(self.model, self.project_parameters)


    def _GetSimulationName(self):
        return "Fluid Chimera Analysis"

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_chimera_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_chimera_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = FluidChimeraAnalysis(model,parameters)
    simulation.Run()
