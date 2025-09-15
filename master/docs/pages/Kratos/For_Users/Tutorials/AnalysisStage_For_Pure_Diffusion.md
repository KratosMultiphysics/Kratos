---
title: Analysis stage for pure diffusion problem
keywords: 
tags: [Tutorial-Analysis-stage-for-pure-diffusion-problem.md]
sidebar: kratos_for_users
summary: 
---

# Overview

In order to create a analysis stage we should follow the structure designed by the *Kratos* team. The base python interface can be found in *Kratos/kratos/python_scripts/analysis_stage.py*.

# Approach

Following that very same structure we just need to define the following methods:
1. The constructor (always needed). As a derived one will call first the base class constructor.
2. `_CreateSolver` In order to create the solver
3. `_CreateProcesses` this one will create the list of processes to be executed
4. `_SetUpGiDOutput` To postprocess the problem
5. `_GetSimulationName` Minor method to have a nice and consistent name on the information printed

# Code 

```python


# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.MyLaplacianApplication as Poisson

# Importing the solvers (if available)
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("ExternalSolversApplication", "not imported")

# Importing the base class
from analysis_stage import AnalysisStage

# Other imports
import sys

class PoissonAnalysis(AnalysisStage):
    """
    This class is the main-script of the Poisson put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Calling the constructor of the base class
        super(PoissonAnalysis, self).__init__(model, project_parameters)

        ## Import parallel modules if needed
        if (self.parallel_type == "MPI"):
            import KratosMultiphysics.MetisApplication as MetisApplication
            import KratosMultiphysics.TrilinosApplication as TrilinosApplication

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return pure_diffusion_solver.PureDiffusionSolver(self.model, self.project_parameters)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        """
        list_of_processes = super(PoissonAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["constraints_process_list", "fluxes_process_list", "list_other_processes"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                KratosMultiphysics.Logger.PrintInfo("PoissonAnalysis", "Using the old way to create the processes, this will be removed!")
                from process_factory import KratosProcessFactory
                factory = KratosProcessFactory(self.model)
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        list_of_processes += factory.ConstructListOfProcesses(self.project_parameters[process_name])
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not alowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                #KratosMultiphysics.Logger.PrintInfo("PoissonAnalysis", "Using the old way to create the gid-output, this will be removed!")
                gid_output= self._SetUpGiDOutput()
                list_of_processes += [gid_output,]
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _SetUpGiDOutput(self):
        '''Initialize a GiD output instance'''
        self.__CheckForDeprecatedGiDSettings()
        if self.parallel_type == "OpenMP":
            from gid_output_process import GiDOutputProcess as OutputProcess
        elif self.parallel_type == "MPI":
            from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

        gid_output = OutputProcess(self._GetSolver().GetComputingModelPart(),
                                   self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                   self.project_parameters["output_configuration"])

        return gid_output

    def _GetSimulationName(self):
        return "::[Poisson Simulation]:: "

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 poisson_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 poisson_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = PoissonAnalysis(model, parameters)
    simulation.Run()
```

We save the file to the folder python_scripts with the name `poisson_analysis`.
