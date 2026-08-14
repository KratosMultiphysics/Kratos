import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.DamApplication as KratosDam

from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis


class InstrumentedDamAnalysis(DamAnalysis):
    """DamAnalysis replica that records the Bofang lifecycle.

    It reproduces byte-for-byte the DamAnalysis::Initialize / temporal-loop
    flow but inserts recorder calls at the lifecycle checkpoints that matter
    for the Bofang temperature assignment.

    Parameters
    ----------
    legacy_process_initialize : bool
        If True the pre-solver-initialize process call is
        ``process.ExecuteInitialize()`` (behaviour immediately before PR
        #13472). If False it is ``process.ExecuteBeforeSolutionLoop()``
        (current behaviour).
    """

    def __init__(self, model, project_parameters, recorder, legacy_process_initialize=False):
        super().__init__(model, project_parameters)
        self.recorder = recorder
        self.legacy_process_initialize = legacy_process_initialize
        self.first_step_handled = False

    def Initialize(self):
        self.solver.AddVariables()  # Add problem variables
        self.solver.ImportModelPart()  # Read model_part (note: the buffer_size is set here)
        self.solver.AddDofs()  # Add degrees of freedom

        if (self.echo_level > 1):
            print(self.main_model_part)
            for self.properties in self.main_model_part.Properties:
                print(self.properties)

        computing_model_part = self.solver.GetComputingModelPart()

        if self.consider_construction:
            thermal_computing_model_part = self.solver.GetComputingThermalModelPart()
            from KratosMultiphysics.DamApplication import dam_construction_utility
            self.construction_utilities = dam_construction_utility.DamConstructionUtility(
                computing_model_part, thermal_computing_model_part,
                self.project_parameters["construction_process"])
            self.construction_utilities.Initialize()

        import KratosMultiphysics.process_factory
        self.list_of_processes = KratosMultiphysics.process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(
            self.project_parameters["constraints_process_list"])
        self.list_of_processes += KratosMultiphysics.process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(
            self.project_parameters["loads_process_list"])
        self.list_of_processes += KratosMultiphysics.process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(
            self.project_parameters["temperature_by_device_list"])

        self.list_of_output_processes = KratosMultiphysics.process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses(
            self.project_parameters["output_device_list"])

        if (self.echo_level > 1):
            for self.process in self.list_of_processes:
                print(self.process)
            for self.output_process in self.list_of_output_processes:
                print(self.output_process)

        # [lifecycle] model imported and processes constructed
        self.recorder.record_lifecycle_stage("after_process_construction")

        # Initialize processes
        if self.legacy_process_initialize:
            for self.process in self.list_of_processes:
                self.process.ExecuteInitialize()
        else:
            for self.process in self.list_of_processes:
                self.process.ExecuteBeforeSolutionLoop()

        # [lifecycle] immediately after the pre-solver-initialize process call
        self.recorder.record_lifecycle_stage("after_process_initialize")

        for self.output_process in self.list_of_output_processes:
            self.output_process.ExecuteInitialize()

        # Set TIME and DELTA_TIME and fill the previous steps of the buffer
        self.time = self.time - (self.buffer_size - 1) * self.delta_time
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.time)
        for step in range(self.buffer_size - 1):
            self.time = self.time + self.delta_time
            self.main_model_part.CloneTimeStep(self.time)

        if self.add_previous_results:
            from KratosMultiphysics.DamApplication import dam_mapping_variables_utility
            dam_mapping_variables_utility = dam_mapping_variables_utility.MappingVariablesUtility(self.domain_size)
            if self.type_of_results == "Mechanical":
                self.post_model_part_mechanical.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
                self.post_model_part_mechanical.AddNodalSolutionStepVariable(KratosPoro.NODAL_CAUCHY_STRESS_TENSOR)
                self.aux_file_name_mechanical = self.file_name_mechanical.replace('.mdpa', '')
                KratosMultiphysics.ModelPartIO(self.aux_file_name_mechanical).ReadModelPart(self.post_model_part_mechanical)
                dam_mapping_variables_utility.AddPreviousModelPartMechanical(
                    self.main_model_part, self.post_model_part_mechanical, self.add_displacement, self.add_stress)
            if self.type_of_results == "Thermal":
                self.post_model_part_thermal.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
                self.post_model_part_thermal.AddNodalSolutionStepVariable(KratosDam.NODAL_REFERENCE_TEMPERATURE)
                self.aux_file_name_thermal = self.file_name_thermal.replace('.mdpa', '')
                KratosMultiphysics.ModelPartIO(self.aux_file_name_thermal).ReadModelPart(self.post_model_part_thermal)
                dam_mapping_variables_utility.AddPreviousModelPartThermal(
                    self.main_model_part, self.post_model_part_thermal, self.add_temperature, self.add_reference_temperature)
            if self.type_of_results == "Thermo-Mechanical":
                self.post_model_part_mechanical.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
                self.post_model_part_mechanical.AddNodalSolutionStepVariable(KratosPoro.NODAL_CAUCHY_STRESS_TENSOR)
                self.aux_file_name_mechanical = self.file_name_mechanical.replace('.mdpa', '')
                KratosMultiphysics.ModelPartIO(self.aux_file_name_mechanical).ReadModelPart(self.post_model_part_mechanical)
                self.post_model_part_thermal.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
                self.post_model_part_thermal.AddNodalSolutionStepVariable(KratosDam.NODAL_REFERENCE_TEMPERATURE)
                self.aux_file_name_thermal = self.file_name_thermal.replace('.mdpa', '')
                KratosMultiphysics.ModelPartIO(self.aux_file_name_thermal).ReadModelPart(self.post_model_part_thermal)
                dam_mapping_variables_utility.AddPreviousModelPartThermoMechanical(
                    self.main_model_part, self.post_model_part_mechanical, self.post_model_part_thermal,
                    self.add_displacement, self.add_stress, self.add_temperature, self.add_reference_temperature)

        output_settings = self.project_parameters["output_configuration"]
        from KratosMultiphysics.DamApplication import dam_cleaning_utility
        dam_cleaning_utility.CleanPreviousFiles(self.problem_path)
        from KratosMultiphysics.DamApplication.gid_dam_output_process import GiDDamOutputProcess
        self.gid_output = GiDDamOutputProcess(computing_model_part, self.problem_name, self.start_time, output_settings)
        self.gid_output.ExecuteInitialize()

        # [lifecycle] immediately before solver.Initialize()
        self.recorder.record_lifecycle_stage("before_solver_initialize")

        self.solver.Initialize()  # Initialize the solver

        # [lifecycle] immediately after solver.Initialize()
        self.recorder.record_lifecycle_stage("after_solver_initialize")

        if self.consider_construction:
            self.construction_utilities.BeforeSolutionLoop()

        # [lifecycle] immediately before ExecuteBeforeSolutionLoop()
        self.recorder.record_lifecycle_stage("before_execute_before_solution_loop")

        # ExecuteBeforeSolutionLoop
        for self.process in self.list_of_processes:
            self.process.ExecuteBeforeSolutionLoop()

        # [lifecycle] immediately after ExecuteBeforeSolutionLoop()
        self.recorder.record_lifecycle_stage("after_execute_before_solution_loop")

        for self.output_process in self.list_of_output_processes:
            self.output_process.ExecuteBeforeSolutionLoop()

        self.gid_output.ExecuteBeforeSolutionLoop()

        self.UseStreamlineUtility = False
        if self.use_streamline_utility and self.domain_size == 3:
            self.UseStreamlineUtility = True
            from KratosMultiphysics.DamApplication import streamlines_output_utility
            self.streamline_utility = streamlines_output_utility.StreamlinesOutputUtility(self.domain_size)

        if (self.echo_level > 1):
            f = open("ProjectParametersOutput.json", 'w')
            f.write(self.project_parameters.PrettyPrintJsonString())
            f.close()

        # [results] initial state at t0 after full initialization
        self.recorder.record_step()

    def RunMainTemporalLoop(self):
        while (self.time + self.tol) <= self.end_time:
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
                self.construction_utilities.InitializeSolutionStep()

            if not self.first_step_handled:
                # [lifecycle] at the start of the first InitializeSolutionStep()
                self.recorder.record_lifecycle_stage("start_first_initialize_solution_step")

            # Update imposed conditions
            for self.process in self.list_of_processes:
                self.process.ExecuteInitializeSolutionStep()

            for self.output_process in self.list_of_output_processes:
                self.output_process.ExecuteInitializeSolutionStep()

            if not self.first_step_handled:
                # [lifecycle] immediately before solving the first step
                self.recorder.record_lifecycle_stage("before_first_solve")

            self.gid_output.ExecuteInitializeSolutionStep()
            self.solver.Solve()  # Solve step

            if not self.first_step_handled:
                # [lifecycle] immediately after solving the first step
                self.recorder.record_lifecycle_stage("after_first_solve")
                self.first_step_handled = True

            self.recorder.record_step()

            if self.UseStreamlineUtility:
                self.streamline_utility.ComputeOutputStep(self.main_model_part, self.domain_size)
            self.gid_output.ExecuteFinalizeSolutionStep()

            for self.process in self.list_of_processes:
                self.process.ExecuteFinalizeSolutionStep()

            for self.output_process in self.list_of_output_processes:
                self.output_process.ExecuteFinalizeSolutionStep()

            for self.process in self.list_of_processes:
                self.process.ExecuteBeforeOutputStep()

            for self.output_process in self.list_of_output_processes:
                self.output_process.ExecuteBeforeOutputStep()

            if self.consider_selfweight:
                self.transfer_utility.Transfer(self.self_weight_model_part, self.main_model_part, self.domain_size)

            if self.add_previous_results and self.add_stress:
                if self.type_of_results == "Mechanical" or self.type_of_results == "Thermo-Mechanical":
                    from KratosMultiphysics.DamApplication import transfer_selfweight_stress_utility
                    self.transfer_utility = transfer_selfweight_stress_utility.TransferSelfweightStressToMainModelPartUtility()
                    self.transfer_utility.TransferInitialStress(self.main_model_part, self.domain_size)

            if self.gid_output.IsOutputStep():
                self.PrintOutput()

            for self.process in self.list_of_processes:
                self.process.ExecuteAfterOutputStep()

            for self.output_process in self.list_of_output_processes:
                self.output_process.ExecuteAfterOutputStep()
                if self.output_process.IsOutputStep():
                    self.output_process.PrintOutput()

            if self.consider_construction:
                self.construction_utilities.AfterOutputStep()

            if self.time == (self.save_intermediate_variables_step * self.time_unit_converter):
                from KratosMultiphysics.DamApplication import save_variables_utility
                self.save_utility = save_variables_utility.SaveVariablesUtility
                if self.save_intermediate_mechanical_variables:
                    self.save_utility.SaveMechanicalVariables(
                        self.problem_name, self.project_parameters, self.main_model_part,
                        self.save_intermediate_variables_step)
                if self.save_intermediate_thermal_variables:
                    self.save_utility.SaveThermalVariables(
                        self.problem_name, self.project_parameters, self.main_model_part,
                        self.save_intermediate_variables_step)

    def Finalize(self):
        from KratosMultiphysics.DamApplication import save_variables_utility
        self.save_utility = save_variables_utility.SaveVariablesUtility
        if self.save_final_mechanical_variables:
            self.save_utility.SaveFinalMechanicalVariables(self.problem_name, self.project_parameters, self.main_model_part)
        if self.save_final_thermal_variables:
            self.save_utility.SaveFinalThermalVariables(self.problem_name, self.project_parameters, self.main_model_part)

        self.gid_output.ExecuteFinalize()

        for self.process in self.list_of_processes:
            self.process.ExecuteFinalize()

        for self.output_process in self.list_of_output_processes:
            self.output_process.ExecuteFinalize()

        self.solver.Clear()


def RunInstrumentedAnalysis(project_parameters_file, output_dir, variant,
                            monitored_node_ids, legacy_process_initialize=False):
    project_parameters_file = os.path.abspath(project_parameters_file)
    case_dir = os.path.dirname(project_parameters_file)
    sys.path.insert(0, case_dir)
    os.chdir(case_dir)

    with open(project_parameters_file, 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()

    model_part_name = parameters["problem_data"]["model_part_name"].GetString()
    from instrumentation.bofang_lifecycle_recorder import BofangLifecycleRecorder
    recorder = BofangLifecycleRecorder(output_dir, variant, model, model_part_name, monitored_node_ids)

    simulation = InstrumentedDamAnalysis(model, parameters, recorder,
                                         legacy_process_initialize=legacy_process_initialize)
    simulation.Run()
    print("Instrumented analysis completed")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 5:
        print("Usage: instrumented_dam_analysis.py <ProjectParameters.json> <output_dir> <variant> <node_id,...> [legacy]")
        sys.exit(1)
    proj = sys.argv[1]
    out_dir = sys.argv[2]
    variant = sys.argv[3]
    node_ids = [int(x) for x in sys.argv[4].split(",")]
    legacy = len(sys.argv) > 5 and sys.argv[5] == "legacy"
    RunInstrumentedAnalysis(proj, out_dir, variant, node_ids, legacy_process_initialize=legacy)
