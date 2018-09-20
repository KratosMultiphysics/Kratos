from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.ChimeraApplication as KratosChimera

from fluid_dynamics_analysis import FluidDynamicsAnalysis

class FluidChimeraAnalysis(FluidDynamicsAnalysis):
    '''Main script for fluid chimera simulations using the navier_stokes family of python solvers.'''

    def __init__(self,model,parameters):
        # Deprecation warnings
        self.parameters = parameters
        solver_settings = parameters["solver_settings"]
        if not solver_settings.Has("domain_size"):
            Kratos.Logger.PrintInfo("FluidChimeraAnalysis", "Using the old way to pass the domain_size, this will be removed!")
            solver_settings.AddEmptyValue("domain_size")
            solver_settings["domain_size"].SetInt(parameters["problem_data"]["domain_size"].GetInt())

        if not solver_settings.Has("model_part_name"):
            Kratos.Logger.PrintInfo("FluidChimeraAnalysis", "Using the old way to pass the model_part_name, this will be removed!")
            solver_settings.AddEmptyValue("model_part_name")
            solver_settings["model_part_name"].SetString(parameters["problem_data"]["model_part_name"].GetString())

        # Import parallel modules if needed
        # has to be done before the base-class constuctor is called (in which the solver is constructed)
        if (parameters["problem_data"]["parallel_type"].GetString() == "MPI"):
            raise Exception("MPI-Chimera is not implemented yet")

        super(FluidChimeraAnalysis,self).__init__(model,parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):

        list_of_processes = super(FluidChimeraAnalysis,self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            # Adding the Chimera Process to the list of processes
            chimera_params = self.project_parameters["chimera"][0]
            main_model_part = self.model["MainModelPart"]
            domain_size = main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]

            print("domain size is",domain_size)

            solver_type = self.parameters["solver_settings"]["solver_type"].GetString()

            if domain_size == 2:
                if(solver_type == "Monolithic" or solver_type == "monolithic"):
                    self.ChimeraProcess = KratosChimera.ApplyChimeraProcessMonolithic2d(main_model_part,chimera_params)
                elif (solver_type == "fractional_step" or solver_type == "FractionalStep"):
                    self.ChimeraProcess = KratosChimera.ApplyChimeraProcessFractionalStep2d(main_model_part,chimera_params)
            else:
                if(solver_type == "Monolithic" or solver_type == "monolithic"):
                    self.ChimeraProcess = KratosChimera.ApplyChimeraProcessMonolithic3d(main_model_part,chimera_params)
                elif (solver_type == "fractional_step" or solver_type == "FractionalStep"):
                    self.ChimeraProcess = KratosChimera.ApplyChimeraProcessFractionalStep3d(main_model_part,chimera_params)

            #list_of_processes.insert(0, ChimeraProcess)
            #print(list_of_processes)
            #err
        return list_of_processes

    ''' def Initialize(self):
        super(FluidChimeraAnalysis,self).Initialize()
        chimera_params = self.project_parameters["chimera"][0]
        main_model_part = self.model["MainModelPart"]
        self.ChimeraProcess = KratosChimera.ApplyChimeraProcessMonolithic2d(main_model_part,chimera_params)
    '''

    def InitializeSolutionStep(self):
        print(" chimera stage")
        self.ChimeraProcess.ExecuteInitializeSolutionStep()
        super(FluidChimeraAnalysis,self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(FluidChimeraAnalysis,self).FinalizeSolutionStep()
        self.ChimeraProcess.ExecuteFinalizeSolutionStep()


    def _CreateSolver(self):
        import python_solvers_wrapper_fluid_chimera
        return python_solvers_wrapper_fluid_chimera.CreateSolver(self.model, self.project_parameters)


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
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = FluidChimeraAnalysis(model,parameters)
    simulation.Run()
