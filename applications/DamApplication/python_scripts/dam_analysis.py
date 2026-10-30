# Time monitoring
import time as timer
print(timer.ctime())
initial_time = timer.perf_counter()

import os

# Importing Kratos Core, Applications and Dependencies
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.DamApplication as KratosDam

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

from importlib import import_module

class DamAnalysis(AnalysisStage):
    """Main-script of the DamApplication put in a class."""
    def __init__(self, model, project_parameters):
        """DamApplication analysis."""
        self.model = model
        self.project_parameters = project_parameters
        self.DefineParallelType()
        self.DefineVariables()

        self.CreateModelPart()

        if self.add_previous_results:
            if self.type_of_results == "Mechanical":
                self.model_mechanical = KratosMultiphysics.Model()
                self.CreatePostModelPartMechanical()
            elif self.type_of_results == "Thermal":
                self.model_thermal = KratosMultiphysics.Model()
                self.CreatePostModelPartThermal()
            else:
                self.model_mechanical = KratosMultiphysics.Model()
                self.model_thermal = KratosMultiphysics.Model()
                self.CreatePostModelPartMechanical()
                self.CreatePostModelPartThermal()

        self.SetSolver()

    def DefineParallelType(self):
        KratosMultiphysics.ParallelUtilities.SetNumThreads(self.project_parameters["problem_data"]["number_of_threads"].GetInt())
        print("OpenMP parallel configuration. OMP_NUM_THREADS =",KratosMultiphysics.ParallelUtilities.GetNumThreads())

    def DefineVariables(self):
        self.domain_size = self.project_parameters["problem_data"]["domain_size"].GetInt()
        self.problem_name = self.project_parameters["problem_data"]["problem_name"].GetString()
        self.problem_path = os.getcwd()
        self.echo_level = self.project_parameters["solver_settings"]["echo_level"].GetInt()
        self.buffer_size = self.project_parameters["solver_settings"]["buffer_size"].GetInt()
        if self.project_parameters["problem_data"].Has("consider_selfweight"):
            if self.project_parameters["problem_data"]["consider_selfweight"].GetBool():
                raise Exception(
                    "The legacy SelfWeight workflow has been removed. "
                    "'problem_data.consider_selfweight' must not be true; "
                    "use standard gravity/body-force loading instead."
                )
        self.consider_construction = self.project_parameters["problem_data"]["consider_construction"].GetBool()
        self.use_streamline_utility = self.project_parameters["problem_data"]["streamlines_utility"].GetBool()
        self.delta_time = self.project_parameters["problem_data"]["time_step"].GetDouble()
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()
        self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
        self.tol = self.delta_time*1.0e-10
        self.time_scale = self.project_parameters["problem_data"]["time_scale"].GetString()
        self.save_intermediate_mechanical_variables = self.project_parameters["transfer_results_process"]["save_intermediate_mechanical_variables"].GetBool()
        self.save_intermediate_thermal_variables = self.project_parameters["transfer_results_process"]["save_intermediate_thermal_variables"].GetBool()
        self.save_intermediate_variables_step = self.project_parameters["transfer_results_process"]["save_intermediate_variables_step"].GetInt()
        self.save_final_mechanical_variables = self.project_parameters["transfer_results_process"]["save_final_mechanical_variables"].GetBool()
        self.save_final_thermal_variables = self.project_parameters["transfer_results_process"]["save_final_thermal_variables"].GetBool()
        self.add_previous_results = self.project_parameters["transfer_results_process"]["add_previous_results"].GetBool()
        self.type_of_results = self.project_parameters["transfer_results_process"]["type_of_results"].GetString()
        self.add_displacement = self.project_parameters["transfer_results_process"]["add_displacement"].GetBool()
        self.add_stress = self.project_parameters["transfer_results_process"]["add_stress"].GetBool()
        self.add_temperature = self.project_parameters["transfer_results_process"]["add_temperature"].GetBool()
        self.add_reference_temperature = self.project_parameters["transfer_results_process"]["add_reference_temperature"].GetBool()

        # Time Units Converter
        if(self.time_scale=="Weeks"):               # Factor to pass from weeks to seconds
            self.time_unit_converter = 604800.0
        elif(self.time_scale=="Days"):               # Factor to pass from days to seconds
            self.time_unit_converter = 86400.0
        elif(self.time_scale=="Hours"):              # Factor to pass from hours to seconds
            self.time_unit_converter = 3600.0
        else:                                       # No changes
            self.time_unit_converter = 1.0

        # Update time variables
        self.start_time = self.time
        self.delta_time = self.delta_time * self.time_unit_converter
        self.end_time = self.end_time * self.time_unit_converter
        self.time = self.time * self.time_unit_converter
        self.tol = self.tol * self.time_unit_converter

    def CreateModelPart(self):
        self.main_model_part = self.model.CreateModelPart(self.project_parameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.domain_size)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.time)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.delta_time)
        self.main_model_part.ProcessInfo.SetValue(KratosDam.TIME_UNIT_CONVERTER, self.time_unit_converter)

    def CreatePostModelPartMechanical(self):
        self.file_name_mechanical = self.project_parameters["transfer_results_process"]["file_name_mechanical"].GetString()
        self.post_model_part_mechanical = self.model_mechanical.CreateModelPart(self.file_name_mechanical)
        self.post_model_part_mechanical.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.domain_size)
        self.post_model_part_mechanical.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.time)
        self.post_model_part_mechanical.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.delta_time)
        self.post_model_part_mechanical.ProcessInfo.SetValue(KratosDam.TIME_UNIT_CONVERTER, self.time_unit_converter)

    def CreatePostModelPartThermal(self):
        self.file_name_thermal = self.project_parameters["transfer_results_process"]["file_name_thermal"].GetString()
        self.post_model_part_thermal = self.model_thermal.CreateModelPart(self.file_name_thermal)
        self.post_model_part_thermal.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.domain_size)
        self.post_model_part_thermal.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.time)
        self.post_model_part_thermal.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.delta_time)
        self.post_model_part_thermal.ProcessInfo.SetValue(KratosDam.TIME_UNIT_CONVERTER, self.time_unit_converter)

    def SetSolver(self):
        python_module_name = "KratosMultiphysics.DamApplication"
        full_module_name = python_module_name + "." + self.project_parameters["solver_settings"]["solver_type"].GetString()
        solver_module = import_module(full_module_name)
        self.solver = solver_module.CreateSolver(self.main_model_part, self.project_parameters["solver_settings"])

    def Run(self):
        self.Initialize()

        self.RunMainTemporalLoop()

        self.Finalize()

    def Initialize(self):

        self.solver.AddVariables() # Add problem variables
        self.solver.ImportModelPart() # Read model_part (note: the buffer_size is set here)
        self.solver.AddDofs() # Add degrees of freedom

        #DamModel = KratosMultiphysics.Model() # Creation of Kratos model
        #DamModel.AddModelPart(self.main_model_part)

        # Print model_part and properties
        if(self.echo_level > 1):
            print(self.main_model_part)
            for self.properties in self.main_model_part.Properties:
                print(self.properties)

        # Initialize GiD I/O
        computing_model_part = self.solver.GetComputingModelPart()

        ## Initialize Construction Utility
        if self.consider_construction:
            thermal_computing_model_part = self.solver.GetComputingThermalModelPart()
            from KratosMultiphysics.DamApplication import dam_construction_utility
            self.construction_utilities = dam_construction_utility.DamConstructionUtility(computing_model_part, thermal_computing_model_part, self.project_parameters["construction_process"])
            self.construction_utilities.Initialize()

        # Construct processes to be applied
        import KratosMultiphysics.process_factory
        self.list_of_processes = KratosMultiphysics.process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.project_parameters["constraints_process_list"] )
        self.list_of_processes += KratosMultiphysics.process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.project_parameters["loads_process_list"] )
        self.list_of_processes += KratosMultiphysics.process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.project_parameters["temperature_by_device_list"] )

        if self.project_parameters.Has("output_device_list"):
            self.list_of_output_processes = KratosMultiphysics.process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.project_parameters["output_device_list"] )
        else:
            self.list_of_output_processes = []

        # Construct standard output processes from the optional output_processes block
        self.list_of_standard_output_processes = []
        if self.project_parameters.Has("output_processes"):
            output_processes = self.project_parameters["output_processes"]
            for _, process_list in output_processes.items():
                self.list_of_standard_output_processes += KratosMultiphysics.process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(process_list)

        # A standard GiD process takes precedence over the legacy output_configuration
        from KratosMultiphysics.gid_output_process import GiDOutputProcess
        has_standard_gid_output = any(
            isinstance(output_process, GiDOutputProcess)
            for output_process in self.list_of_standard_output_processes
        )
        use_legacy_gid_output = (
            self.project_parameters.Has("output_configuration") and
            not has_standard_gid_output
        )
        if has_standard_gid_output and self.project_parameters.Has("output_configuration"):
            KratosMultiphysics.Logger.PrintWarning(
                "DamAnalysis",
                "Both 'output_configuration' and a standard GiD process in 'output_processes' "
                "are present. The standard GiD process is used and the legacy "
                "'output_configuration' is ignored."
            )

        # Standard GiD must clean previous post files before its ExecuteInitialize below.
        # Legacy GiD keeps its original cleanup point, immediately before construction.
        # In particular, VTK-only runs must not delete existing GiD results.
        if has_standard_gid_output:
            from KratosMultiphysics.DamApplication import dam_cleaning_utility
            dam_cleaning_utility.CleanPreviousFiles(self.problem_path)

        self.gid_output = None

        # Print list of constructed processes
        if(self.echo_level>1):
            for self.process in self.list_of_processes:
                print(self.process)

            for self.output_process in self.list_of_output_processes:
                print(self.output_process)

            for self.output_process in self.list_of_standard_output_processes:
                print(self.output_process)

        # Initialize processes
        for self.process in self.list_of_processes:
            self.process.ExecuteInitialize()

        for self.output_process in self.list_of_output_processes:
            self.output_process.ExecuteInitialize()

        for self.output_process in self.list_of_standard_output_processes:
            self.output_process.ExecuteInitialize()

        # Set TIME and DELTA_TIME and fill the previous steps of the buffer with the initial conditions
        self.time = self.time - (self.buffer_size-1)*self.delta_time
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.time)
        for step in range(self.buffer_size-1):
            self.time = self.time + self.delta_time
            self.main_model_part.CloneTimeStep(self.time)


        ## Initialize Mapping Variables Utility
        if self.add_previous_results:

            from KratosMultiphysics.DamApplication import dam_mapping_variables_utility
            dam_mapping_variables_utility = dam_mapping_variables_utility.MappingVariablesUtility(self.domain_size)

            if self.type_of_results == "Mechanical":
                self.post_model_part_mechanical.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
                self.post_model_part_mechanical.AddNodalSolutionStepVariable(KratosPoro.NODAL_CAUCHY_STRESS_TENSOR)
                self.aux_file_name_mechanical = self.file_name_mechanical.replace('.mdpa','')
                KratosMultiphysics.ModelPartIO(self.aux_file_name_mechanical).ReadModelPart(self.post_model_part_mechanical)

                dam_mapping_variables_utility.AddPreviousModelPartMechanical(self.main_model_part,self.post_model_part_mechanical,self.add_displacement,self.add_stress)

            if self.type_of_results == "Thermal":
                self.post_model_part_thermal.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
                self.post_model_part_thermal.AddNodalSolutionStepVariable(KratosDam.NODAL_REFERENCE_TEMPERATURE)
                self.aux_file_name_thermal = self.file_name_thermal.replace('.mdpa','')
                KratosMultiphysics.ModelPartIO(self.aux_file_name_thermal).ReadModelPart(self.post_model_part_thermal)

                dam_mapping_variables_utility.AddPreviousModelPartThermal(self.main_model_part,self.post_model_part_thermal,self.add_temperature,self.add_reference_temperature)

            if self.type_of_results == "Thermo-Mechanical":
                self.post_model_part_mechanical.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
                self.post_model_part_mechanical.AddNodalSolutionStepVariable(KratosPoro.NODAL_CAUCHY_STRESS_TENSOR)
                self.aux_file_name_mechanical = self.file_name_mechanical.replace('.mdpa','')
                KratosMultiphysics.ModelPartIO(self.aux_file_name_mechanical).ReadModelPart(self.post_model_part_mechanical)

                self.post_model_part_thermal.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
                self.post_model_part_thermal.AddNodalSolutionStepVariable(KratosDam.NODAL_REFERENCE_TEMPERATURE)
                self.aux_file_name_thermal = self.file_name_thermal.replace('.mdpa','')
                KratosMultiphysics.ModelPartIO(self.aux_file_name_thermal).ReadModelPart(self.post_model_part_thermal)

                dam_mapping_variables_utility.AddPreviousModelPartThermoMechanical(self.main_model_part,self.post_model_part_mechanical,self.post_model_part_thermal,self.add_displacement,self.add_stress,self.add_temperature,self.add_reference_temperature)

        if use_legacy_gid_output:
            from KratosMultiphysics.DamApplication import dam_cleaning_utility
            dam_cleaning_utility.CleanPreviousFiles(self.problem_path) # Clean previous post files
            output_settings = self.project_parameters["output_configuration"]
            from KratosMultiphysics.DamApplication.gid_dam_output_process import GiDDamOutputProcess
            self.gid_output = GiDDamOutputProcess(computing_model_part,
                                                  self.problem_name,
                                                  self.start_time,
                                                  output_settings)
            self.gid_output.ExecuteInitialize()

        self.solver.Initialize() # Initialize the solver

        if self.consider_construction:
            # Execute initialize solution
            self.construction_utilities.BeforeSolutionLoop()

        # ExecuteBeforeSolutionLoop
        for self.process in self.list_of_processes:
            self.process.ExecuteBeforeSolutionLoop()

        for self.output_process in self.list_of_output_processes:
            self.output_process.ExecuteBeforeSolutionLoop()

        for self.output_process in self.list_of_standard_output_processes:
            self.output_process.ExecuteBeforeSolutionLoop()

        if self.gid_output is not None:
            self.gid_output.ExecuteBeforeSolutionLoop() # Set results when they are written in a single file

        # Initialize streamlines_output_utility
        self.UseStreamlineUtility = False
        if self.use_streamline_utility and self.domain_size==3:
            self.UseStreamlineUtility = True
            from KratosMultiphysics.DamApplication import streamlines_output_utility
            self.streamline_utility = streamlines_output_utility.StreamlinesOutputUtility(self.domain_size)

        if (self.echo_level > 1):
            f = open("ProjectParametersOutput.json", 'w')
            f.write(self.project_parameters.PrettyPrintJsonString())
            f.close()

    def RunMainTemporalLoop(self):

        while( (self.time+self.tol) <= self.end_time ):

            # Update temporal variables
            self.time = self.time + self.delta_time
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.delta_time)
            self.main_model_part.CloneTimeStep(self.time)

            if self.add_previous_results:
                if self.type_of_results == "Mechanical":
                    self.post_model_part_mechanical.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.delta_time)
                if self.type_of_results == "Thermal":
                    self.post_model_part_thermal.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.delta_time)
                if self.type_of_results == "Thermo-Mechanical":
                    self.post_model_part_mechanical.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.delta_time)
                    self.post_model_part_thermal.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.delta_time)

            if self.consider_construction:
                # Execute initialize solution
                self.construction_utilities.InitializeSolutionStep()

            # Update imposed conditions
            for self.process in self.list_of_processes:
                self.process.ExecuteInitializeSolutionStep()

            for self.output_process in self.list_of_output_processes:
                self.output_process.ExecuteInitializeSolutionStep()

            # ExecuteInitializeSolutionStep
            for self.output_process in self.list_of_standard_output_processes:
                self.output_process.ExecuteInitializeSolutionStep()

            if self.gid_output is not None:
                self.gid_output.ExecuteInitializeSolutionStep()

            self.solver.Solve() # Solve step

            # streamlines_output_utility
            if self.UseStreamlineUtility:
                self.streamline_utility.ComputeOutputStep(self.main_model_part, self.domain_size)

            if self.gid_output is not None:
                self.gid_output.ExecuteFinalizeSolutionStep()

            for self.output_process in self.list_of_standard_output_processes:
                self.output_process.ExecuteFinalizeSolutionStep()

            for self.process in self.list_of_processes:
                self.process.ExecuteFinalizeSolutionStep()

            for self.output_process in self.list_of_output_processes:
                self.output_process.ExecuteFinalizeSolutionStep()

            standard_output_flags = [
                (output_process, output_process.IsOutputStep())
                for output_process in self.list_of_standard_output_processes
            ]
            output_device_flags = [
                (output_process, output_process.IsOutputStep())
                for output_process in self.list_of_output_processes
            ]
            legacy_gid_output_step = (
                self.gid_output is not None and self.gid_output.IsOutputStep()
            )
            has_output_step = (
                legacy_gid_output_step or
                any(flag for _, flag in standard_output_flags) or
                any(flag for _, flag in output_device_flags)
            )

            if has_output_step:
                for self.process in self.list_of_processes:
                    self.process.ExecuteBeforeOutputStep()

                for self.output_process in self.list_of_output_processes:
                    self.output_process.ExecuteBeforeOutputStep()

            # add previous results utility
            # The smoothing schemes recompute the total nodal stress from the
            # current Gauss-point state. If the displacement was also mapped,
            # adding the imported initial stress again would double count it.
            # UP/P schemes do not recompute nodal stress and keep the imported
            # contribution for output.
            recomputes_nodal_cauchy_stress = getattr(
                self.solver, "_recomputes_nodal_cauchy_stress", False)
            if (self.add_previous_results and self.add_stress and
                    not (self.add_displacement and recomputes_nodal_cauchy_stress)):
                if self.type_of_results == "Mechanical" or self.type_of_results == "Thermo-Mechanical":
                    from KratosMultiphysics.DamApplication import transfer_initial_stress_utility
                    self.transfer_utility = transfer_initial_stress_utility.TransferInitialStressToMainModelPartUtility()
                    self.transfer_utility.TransferInitialStress(self.main_model_part, self.domain_size)

            if has_output_step:
                for output_process, _ in standard_output_flags:
                    output_process.ExecuteBeforeOutputStep()

                if self.gid_output is not None:
                    self.gid_output.ExecuteBeforeOutputStep()

                # Write legacy GiD results
                if legacy_gid_output_step:
                    self.PrintOutput()

                for output_process, flag in standard_output_flags:
                    if flag:
                        output_process.PrintOutput()

                for output_process, flag in output_device_flags:
                    if flag:
                        output_process.PrintOutput()

                for self.process in self.list_of_processes:
                    self.process.ExecuteAfterOutputStep()

                for self.output_process in self.list_of_output_processes:
                    self.output_process.ExecuteAfterOutputStep()

                for output_process, _ in standard_output_flags:
                    output_process.ExecuteAfterOutputStep()

                if self.gid_output is not None:
                    self.gid_output.ExecuteAfterOutputStep()

            if self.consider_construction:
                #  After initialize solution
                self.construction_utilities.AfterOutputStep()

            # Save results at any time in an auxiliary .mdpa
            if self.time == (self.save_intermediate_variables_step*self.time_unit_converter):
                from KratosMultiphysics.DamApplication import save_variables_utility
                self.save_utility = save_variables_utility.SaveVariablesUtility
                if self.save_intermediate_mechanical_variables:
                    self.save_utility.SaveMechanicalVariables(self.problem_name, self.project_parameters, self.main_model_part, self.save_intermediate_variables_step)
                if self.save_intermediate_thermal_variables:
                    self.save_utility.SaveThermalVariables(self.problem_name, self.project_parameters, self.main_model_part, self.save_intermediate_variables_step)

    def PrintOutput(self):
        if self.gid_output is not None:
            self.gid_output.PrintOutput()

    def Finalize(self):
        # Save final results in an auxiliary .mdpa
        from KratosMultiphysics.DamApplication import save_variables_utility
        self.save_utility = save_variables_utility.SaveVariablesUtility
        if self.save_final_mechanical_variables:
            self.save_utility.SaveFinalMechanicalVariables(self.problem_name, self.project_parameters, self.main_model_part)

        if self.save_final_thermal_variables:
            self.save_utility.SaveFinalThermalVariables(self.problem_name, self.project_parameters, self.main_model_part)

        if self.gid_output is not None:
            self.gid_output.ExecuteFinalize() # Finalizing output files

        for self.output_process in self.list_of_standard_output_processes:
            self.output_process.ExecuteFinalize()

        for self.process in self.list_of_processes:
            self.process.ExecuteFinalize()

        for self.output_process in self.list_of_output_processes:
            self.output_process.ExecuteFinalize()

        # Finalizing strategy
        #if self.parallel_type == "OpenMP":
        self.solver.Clear()

        # Time control
        print("Analysis Completed. Elapsed Time = %.3f" % (timer.perf_counter() - initial_time)," seconds.")
        print(timer.ctime())

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 dam_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 dam_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = DamAnalysis(model, parameters)
    simulation.Run()
