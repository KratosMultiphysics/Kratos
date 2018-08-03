from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as Kratos
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.ShallowWaterApplication as Shallow

from analysis_stage import AnalysisStage

class ShallowWaterAnalysis(AnalysisStage):
    ''' Main script for shallow water simulations '''

    def _CreateSolver(self):
        solver_module = __import__(self.project_parameters["solver_settings"]["solver_type"].GetString())
        solver = solver_module.CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver

    def _GetOrderOfProcessesInitialization(self):
        return ["bathymetry_process_list",
                "initial_conditions_process_list",
                "boundary_conditions_process_list"]

    def _GetSimulationName(self):
        return "Shallow Water Analysis"

    # DANGER: look this code
    def SetUpModel(self):
        '''Initialize the model part for the problem and other general model data.'''
        
        ## Defining variables ----------------------------------------------------------------------------------------
        gravity             = self.ProjectParameters["problem_data"]["gravity"].GetDouble()
        time_scale          = self.ProjectParameters["problem_data"]["time_scale"].GetString()
        water_height_scale  = self.ProjectParameters["problem_data"]["water_height_scale"].GetString()

        # Time unit converter
        if   time_scale == "seconds":
            time_unit_converter =     1
        elif time_scale == "minutes":
            time_unit_converter =    60
        elif time_scale == "hours":
            time_unit_converter =  3600
        elif time_scale == "days":
            time_unit_converter = 86400
        else:
            raise Exception("unknown time scale")

        # Water height unit converter
        if   water_height_scale == "meters":
            water_height_unit_converter = 1.0
        elif water_height_scale == "millimeters":
            water_height_unit_converter = 0.001
        else:
            raise Exception("unknown water height scale")

        ## Model part ------------------------------------------------------------------------------------------------

        # Defining the model part
        # self.main_model_part = Kratos.ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(Kratos.GRAVITY_Z, gravity * time_unit_converter**2)
        self.main_model_part.ProcessInfo.SetValue(Shallow.TIME_UNIT_CONVERTER, time_unit_converter)
        self.main_model_part.ProcessInfo.SetValue(Shallow.WATER_HEIGHT_UNIT_CONVERTER, water_height_unit_converter)

        # Solver construction (main settings methods are located in the solver_module)
        solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
        self.solver = solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])

        self.solver.AddVariables()
        self.solver.ImportModelPart()
        self.solver.AddDofs()

        # Fill a Model instance using input
        self.model = Kratos.Model()
        self.model.AddModelPart(self.main_model_part)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(ShallowWaterAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                gid_output = self._SetUpGiDOutput()
                list_of_processes += [gid_output,]
        
        return list_of_processes

    def _SetUpGiDOutput(self):
        '''Initialize self.output as a GiD output instance.'''
        if self.parallel_type == "OpenMP":
            from gid_output_process import GiDOutputProcess as OutputProcess
        elif self.parallel_type == "MPI":
            from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

        output = OutputProcess(self._GetSolver().GetComputingModelPart(),
                                self.project_parameters["problem_data"]["problem_name"].GetString(),
                                self.project_parameters["output_configuration"])

        return output


if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    model = Kratos.Model()
    ShallowWaterAnalysis(model, project_parameters_file_name).Run()