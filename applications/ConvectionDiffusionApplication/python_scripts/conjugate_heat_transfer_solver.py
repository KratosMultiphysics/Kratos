from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    return ConjugateHeatTransferSolver(main_model_part, custom_settings)

class ConjugateHeatTransferSolver(PythonSolver):

    def _validate_settings(self, custom_settings):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "conjugate_heat_transfer",
            "domain_size": -1,
            "echo_level": 0,
            "fluid_domain_solver_settings": {
                "solver_type": "ThermallyCoupled",
                "domain_size": -1,
                "echo_level": 1,
                "fluid_solver_settings": {
                    "solver_type": "Monolithic",
                    "model_import_settings": {
                        "input_type": "mdpa",
                        "input_filename": "unknown_name"
                    }
                },
                "thermal_solver_settings":{
                    "model_part_name": "FluidThermalModelPart",
                    "solver_type": "Transient",
                    "analysis_type": "linear",
                    "model_import_settings": {
                        "input_type": "use_input_model_part"
                    },
                    "material_import_settings": {
                        "materials_filename": "ThermicMaterialsFluid.json"
                    }
                }
            },
            "solid_domain_solver_settings":{
                "solid_solver_settings": {
                },
                "thermal_solver_settings": {
                    "model_part_name": "SolidThermalModelPart",
                    "solver_type": "Transient",
                    "analysis_type": "linear",
                    "model_import_settings": {
                        "input_type": "mdpa",
                        "input_filename": "unknown_name"
                    },
                    "material_import_settings": {
                        "materials_filename": "ThermicMaterialsSolid.json"
                    }
                }
            },
            "coupling_settings":{
                "max_iteration": 10,
                "temperature_relative_tolerance": 1e-5,
                "dirichlet_coupling_interface": "fluid",
                "variable_redistribution_settings": {
                    "absolute_tolerance": 1.0e-9,
                    "max_iterations": 200
                },
                "mappers_settings": {
                    "echo_level": 0,
                    "distance_threshold": 1.0e+24,
                    "absolute_convergence_tolerance": 1.0e-9,
                    "relative_convergence_tolerance": 1.0e-7,
                    "max_number_iterations": 10,
                    "integration_order": 2,
                    "search_parameters": {
                        "allocation_size": 1000,
                        "bucket_size": 4,
                        "search_factor": 1.0
                    }
                },
                "convergence_accelerator_settings": {
                    "solver_type": "Relaxation",
                    "acceleration_type": "Aitken",
                    "w_0": 0.5
                },
                "fluid_interfaces_list": [],
                "solid_interfaces_list": []
            }
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Validate fluid and solid domains settings
        self.settings["fluid_domain_solver_settings"].ValidateAndAssignDefaults(default_settings["fluid_domain_solver_settings"])
        self.settings["solid_domain_solver_settings"].ValidateAndAssignDefaults(default_settings["solid_domain_solver_settings"])

        ## Validate coupling settings
        self.settings["coupling_settings"].ValidateAndAssignDefaults(default_settings["coupling_settings"])

    def __init__(self, model, custom_settings):
        super(ConjugateHeatTransferSolver, self).__init__(model, custom_settings)
        ## Validate custom settings against defaults
        self._validate_settings(custom_settings)

        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        ## Set the fluid dynamics solver
        from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_domain_solver_settings"]["fluid_solver_settings"], "OpenMP")

        # Set the fluid and solid heat solvers
        from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion
        self.fluid_thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model, self.settings["fluid_domain_solver_settings"]["thermal_solver_settings"], "OpenMP")
        self.solid_thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model, self.settings["solid_domain_solver_settings"]["thermal_solver_settings"], "OpenMP")

    def AddVariables(self):
        self.fluid_solver.AddVariables()
        self.fluid_thermal_solver.AddVariables()

        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.AUX_FLUX)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.AUX_TEMPERATURE)
        self.fluid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.AUX_FLUX)
        self.fluid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.fluid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.AUX_TEMPERATURE)
        self.fluid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SCALAR_INTERFACE_RESIDUAL)

        # Temporary container for un-relaxed temperature
        KratosMultiphysics.MergeVariableListsUtility().Merge(
            self.fluid_solver.main_model_part,
            self.fluid_thermal_solver.main_model_part)

        self.solid_thermal_solver.AddVariables()
        self.solid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.AUX_FLUX)
        self.solid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.solid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.AUX_TEMPERATURE)
        self.solid_thermal_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SCALAR_INTERFACE_RESIDUAL)

    def ImportModelPart(self):
        # Check that both thermal solvers have a different model part name. If
        # both model part names coincide the solver will fail to acces them. This
        # is the case if the default one in the convection diffusion is taken.
        fluid_thermal_model_part_name = self.settings["fluid_domain_solver_settings"]["thermal_solver_settings"]["model_part_name"].GetString()
        solid_thermal_model_part_name = self.settings["solid_domain_solver_settings"]["thermal_solver_settings"]["model_part_name"].GetString()
        if fluid_thermal_model_part_name == solid_thermal_model_part_name:
            err_msg = "\nFluid thermal solver settings model_part_name and solid thermal solver settings model_part_name can not coincide.\n"
            err_msg += "- fluid model_part_name: " + fluid_thermal_model_part_name + "\n"
            err_msg += "- solid model_part_name: " + solid_thermal_model_part_name + "\n"
            err_msg += "Provide different model_part_names in the JSON settings file."
            raise Exception(err_msg)

        # Import the fluid domain in the fluid dynamics solver
        self.fluid_solver.ImportModelPart()

        # In order to consider the buoyancy effects, the nodes in the fluid model part must
        # be shared with the nodes in the fluid thermal model part. To do that, we use the modeler
        # Save the convection diffusion settings
        convection_diffusion_settings = self.fluid_thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)

        # Here the fluid model part is cloned to be thermal model part so that the nodes are shared
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        if(self.domain_size == 2):
            modeler.GenerateModelPart(self.fluid_solver.main_model_part,
                                      self.fluid_thermal_solver.main_model_part,
                                      "Element2D3N",
                                      "LineCondition2D2N")
        else:
            modeler.GenerateModelPart(self.fluid_solver.main_model_part,
                                      self.fluid_thermal_solver.main_model_part,
                                      "Element3D4N",
                                      "SurfaceCondition3D3N")

        # Set the saved convection diffusion settings to the new thermal model part
        self.fluid_thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)

        # Confirm that the buffer size in the shared nodes is the maximum required one
        fluid_solver_buffer = self.fluid_solver.GetMinimumBufferSize()
        fluid_thermal_solver_buffer = self.fluid_thermal_solver.GetMinimumBufferSize()
        if (fluid_solver_buffer != fluid_thermal_solver_buffer):
            max_buffer_size = max(fluid_solver_buffer, fluid_thermal_solver_buffer)
            self.fluid_solver.min_buffer_size = max_buffer_size
            self.fluid_thermal_solver.min_buffer_size = max_buffer_size
            if self.settings["echo_level"].GetInt() >= 1:
                warning_msg = "Fluid solver and fluid thermal solver have different buffer size:\n"
                warning_msg += " - Fluid solver buffer size: " + str(fluid_solver_buffer) + "\n"
                warning_msg += " - Fluid thermal solver buffer size: " + str(fluid_thermal_solver_buffer) + "\n"
                warning_msg += "Setting buffer size equal to " + str(max_buffer_size) + " in both solvers."
                KratosMultiphysics.Logger.PrintWarning("::[ConjugateHeatTransferSolver]::", warning_msg)

        # Import the solid domain
        self.solid_thermal_solver.ImportModelPart()

    def AddDofs(self):
        (self.fluid_solver).AddDofs()
        (self.fluid_thermal_solver).AddDofs()
        (self.solid_thermal_solver).AddDofs()

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()

    def ComputeDeltaTime(self):
        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_thermal_fluid = self.fluid_thermal_solver.GetMinimumBufferSize()
        buffer_size_thermal_solid = self.solid_thermal_solver.GetMinimumBufferSize()

        return max(buffer_size_fluid, buffer_size_thermal_fluid, buffer_size_thermal_solid)

    def Initialize(self):
        # Check that the reactions are computed in the Dirichlet side
        if self.settings["coupling_settings"]["dirichlet_coupling_interface"].GetString() == "fluid":
            if not (self.fluid_thermal_solver).settings["compute_reactions"].GetBool():
                (self.fluid_thermal_solver).settings["compute_reactions"].SetBool(True)
        else:
            if not (self.solid_thermal_solver).settings["compute_reactions"].GetBool():
                (self.solid_thermal_solver).settings["compute_reactions"].SetBool(True)

        # Initialize the fluid and solid solvers
        (self.fluid_solver).Initialize()
        (self.fluid_thermal_solver).Initialize()
        (self.solid_thermal_solver).Initialize()

        # Create the fluid and solid interface mapper
        self._set_up_mappers()

        # Coupling utility initialization
        # The _get_convergence_accelerator is supposed to construct the convergence accelerator in here
        self._get_convergence_accelerator().Initialize()

    def Clear(self):
        (self.fluid_solver).Clear()
        (self.fluid_thermal_solver).Clear()
        (self.solid_thermal_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.fluid_thermal_solver).Check()
        (self.solid_thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.fluid_thermal_solver).SetEchoLevel(level)
        (self.solid_thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        #The cloning is done ONLY ONCE since the nodes are shared between the fluid and thermal solvers
        new_time = self.fluid_solver.AdvanceInTime(current_time)

        # Do the time advance in the solid thermal solver
        self.solid_thermal_solver.main_model_part.CloneTimeStep(new_time)
        self.solid_thermal_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def PrepareModelPart(self):
        self.fluid_solver.PrepareModelPart()
        # TODO: CHECK THIS (if we switch the order in thermal solvers the solver breaks)
        self.solid_thermal_solver.PrepareModelPart()
        self.fluid_thermal_solver.PrepareModelPart()

        self._set_up_dirichlet_coupling_boundary()

    def InitializeSolutionStep(self):
        if self._time_buffer_is_initialized():
            self.fluid_solver.InitializeSolutionStep()
            self.fluid_thermal_solver.InitializeSolutionStep()
            self.solid_thermal_solver.InitializeSolutionStep()
            self._get_convergence_accelerator().InitializeSolutionStep()

    def Predict(self):
        if self._time_buffer_is_initialized():
            self.fluid_solver.Predict()
            self.fluid_thermal_solver.Predict()
            self.solid_thermal_solver.Predict()

    def SolveSolutionStep(self):
        if self._time_buffer_is_initialized():
            max_iteration = self.settings["coupling_settings"]["max_iteration"].GetInt()
            temp_rel_tol = self.settings["coupling_settings"]["temperature_relative_tolerance"].GetDouble()
            redistribution_tolerance = self.settings["coupling_settings"]["variable_redistribution_settings"]["absolute_tolerance"].GetDouble()
            redistribution_max_iterations = self.settings["coupling_settings"]["variable_redistribution_settings"]["max_iterations"].GetInt()

            # Solve the buoyancy solver
            self.fluid_solver.SolveSolutionStep()

            # Interface temperature prediction
            self._temperature_coupling_prediction()

            # Initialize iteration value vector
            self._initialize_iteration_value_vector()

            # Couple the solid and fluid thermal problems
            iteration = 0
            KratosMultiphysics.Logger.PrintInfo("::[ConjugateHeatTransferSolver]::", "Starting non-linear temperature coupling")
            while iteration < max_iteration:
                # Initialize non-linear iteration
                iteration += 1
                self.solid_thermal_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = iteration
                self.fluid_thermal_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = iteration
                self._get_convergence_accelerator().InitializeNonLinearIteration()

                # Solve Dirichlet side to get reactions from fluid domain
                # self.fluid_thermal_solver.SolveSolutionStep()
                self._get_dirichlet_interface_thermal_solver().SolveSolutionStep()

                # Map reactions to the solid interface. Note that we first call the redistribution utility to convert the point values to distributed ones
                KratosMultiphysics.VariableRedistributionUtility.DistributePointValues(
                    self._get_dirichlet_coupling_interface(),
                    KratosMultiphysics.REACTION_FLUX,
                    KratosConvDiff.AUX_FLUX,
                    redistribution_tolerance,
                    redistribution_max_iterations)

                self.flux_mapper.Execute()

                # Solve Neumann side to get temperature values from fluid solid
                # self.solid_thermal_solver.Solve()
                self._get_neumann_interface_thermal_solver().Solve()

                # Map back the Neumann domain obtained temperature
                self.temp_mapper.Execute()

                # Compute the interface residual
                temp_residual = KratosMultiphysics.Vector(self._get_partitioned_FSI_utilities().GetInterfaceResidualSize(self._get_dirichlet_coupling_interface()))
                self._get_partitioned_FSI_utilities().ComputeInterfaceResidualVector(
                    self._get_dirichlet_coupling_interface(),
                    KratosMultiphysics.TEMPERATURE,
                    KratosConvDiff.AUX_TEMPERATURE,
                    KratosMultiphysics.SCALAR_INTERFACE_RESIDUAL,
                    temp_residual,
                    "nodal",
                    KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM) #TODO: Rename this variable

                # Residual computation
                rel_res_norm = self._get_dirichlet_coupling_interface().ProcessInfo[KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM] / len(self._get_dirichlet_coupling_interface().Nodes)
                KratosMultiphysics.Logger.PrintInfo("::[ConjugateHeatTransferSolver]::", "Iteration: " + str(iteration) + " Relative residual: " + str(rel_res_norm))

                # Perform the convergence accelerator solution update
                self._get_convergence_accelerator().UpdateSolution(temp_residual, self.iteration_value)

                # Update the interface with the corrected values
                self._get_partitioned_FSI_utilities().UpdateInterfaceValues(
                    self._get_dirichlet_coupling_interface(),
                    KratosMultiphysics.TEMPERATURE,
                    self.iteration_value)

                # Finalize the congergence accelerator iteration
                self._get_convergence_accelerator().FinalizeNonLinearIteration()

                # Check convergence
                if rel_res_norm <= temp_rel_tol:
                    KratosMultiphysics.Logger.PrintInfo("::[ConjugateHeatTransferSolver]::", "Converged in " + str(iteration) + " iterations.")
                    break
                elif iteration == max_iteration:
                    KratosMultiphysics.Logger.PrintInfo("::[ConjugateHeatTransferSolver]::", "Did not converge in " + str(iteration) + " iterations.")

    def FinalizeSolutionStep(self):
        if self._time_buffer_is_initialized():
            self.fluid_solver.FinalizeSolutionStep()
            self.fluid_thermal_solver.FinalizeSolutionStep()
            self.solid_thermal_solver.FinalizeSolutionStep()
            self._get_convergence_accelerator().FinalizeSolutionStep()

    def _set_up_dirichlet_coupling_boundary(self):
        # Run the solid interfaces list to fix the temperature DOFs.
        for node in self._get_dirichlet_coupling_interface().Nodes:
            node.Fix(KratosMultiphysics.TEMPERATURE)

    def _set_up_mappers(self):
        # Set mappers settings
        mappers_settings = self.settings["coupling_settings"]["mappers_settings"]

        flux_mapper_parameters = KratosMultiphysics.Parameters(r'''{
            "mapping_coefficient": -1.0,
            "origin_variable": "AUX_FLUX",
            "destination_variable": "FACE_HEAT_FLUX"
        }''')
        flux_mapper_parameters.AddValue("echo_level", mappers_settings["echo_level"])
        flux_mapper_parameters.AddValue("integration_order", mappers_settings["integration_order"])
        flux_mapper_parameters.AddValue("distance_threshold", mappers_settings["distance_threshold"])
        flux_mapper_parameters.AddValue("max_number_iterations", mappers_settings["max_number_iterations"])
        flux_mapper_parameters.AddValue("absolute_convergence_tolerance", mappers_settings["absolute_convergence_tolerance"])
        flux_mapper_parameters.AddValue("relative_convergence_tolerance", mappers_settings["relative_convergence_tolerance"])
        flux_mapper_parameters.AddValue("search_parameters", mappers_settings["search_parameters"])

        temp_mapper_parameters = KratosMultiphysics.Parameters(r'''{
            "origin_variable": "TEMPERATURE",
            "destination_variable": "AUX_TEMPERATURE"
        }''')
        temp_mapper_parameters.AddValue("echo_level", mappers_settings["echo_level"])
        temp_mapper_parameters.AddValue("integration_order", mappers_settings["integration_order"])
        temp_mapper_parameters.AddValue("distance_threshold", mappers_settings["distance_threshold"])
        temp_mapper_parameters.AddValue("max_number_iterations", mappers_settings["max_number_iterations"])
        temp_mapper_parameters.AddValue("absolute_convergence_tolerance", mappers_settings["absolute_convergence_tolerance"])
        temp_mapper_parameters.AddValue("relative_convergence_tolerance", mappers_settings["relative_convergence_tolerance"])
        temp_mapper_parameters.AddValue("search_parameters", mappers_settings["search_parameters"])

        # Create flux mapper
        self.flux_mapper = KratosMultiphysics.SimpleMortarMapperProcess(
            self._get_dirichlet_coupling_interface(),
            self._get_neumann_coupling_interface(),
            flux_mapper_parameters)

        # Create temperature mapper
        self.temp_mapper = KratosMultiphysics.SimpleMortarMapperProcess(
            self._get_neumann_coupling_interface(),
            self._get_dirichlet_coupling_interface(),
            temp_mapper_parameters)

    # This method returns the fluid thermal interface
    def _get_fluid_thermal_interface(self):
        if self.settings["coupling_settings"]["fluid_interfaces_list"].size() > 1:
            raise Exception("More than one fluid interfaces is not supported yet.")

        fluid_int_name = self.settings["coupling_settings"]["fluid_interfaces_list"][0].GetString()
        return self.model.GetModelPart(fluid_int_name)

    # This method returns the solid thermal interface
    def _get_solid_thermal_interface(self):
        if self.settings["coupling_settings"]["solid_interfaces_list"].size() > 1:
            raise Exception("More than one solid interfaces is not supported yet.")

        solid_int_name = self.settings["coupling_settings"]["solid_interfaces_list"][0].GetString()
        return self.model.GetModelPart(solid_int_name)

    # This method returns the Dirichlet coupling interface model part
    def _get_dirichlet_coupling_interface(self):
        if self.settings["coupling_settings"]["dirichlet_coupling_interface"].GetString() == "fluid":
            return self._get_fluid_thermal_interface()
        else:
            return self._get_solid_thermal_interface()

    # This method returns the Neumann coupling interface model part
    def _get_neumann_coupling_interface(self):
        if self.settings["coupling_settings"]["dirichlet_coupling_interface"].GetString() == "fluid":
            return self._get_solid_thermal_interface()
        else:
            return self._get_fluid_thermal_interface()

    # This method returns the domain thermal solver corresponding to the Dirichlet coupling interface
    def _get_dirichlet_interface_thermal_solver(self):
        if self.settings["coupling_settings"]["dirichlet_coupling_interface"].GetString() == "fluid":
            return self.fluid_thermal_solver
        else:
            return self.solid_thermal_solver

    # This method returns the domain thermal solver corresponding to the Neumann coupling interface
    def _get_neumann_interface_thermal_solver(self):
        if self.settings["coupling_settings"]["dirichlet_coupling_interface"].GetString() == "fluid":
            return self.solid_thermal_solver
        else:
            return self.fluid_thermal_solver

    # This method returns the convergence accelerator.
    # If it is not created yet, it calls the _create_convergence_accelerator first
    def _get_convergence_accelerator(self):
        if not hasattr(self, '_convergence_accelerator'):
            self._convergence_accelerator = self._create_convergence_accelerator()
        return self._convergence_accelerator

    # This method constructs the convergence accelerator coupling utility
    def _create_convergence_accelerator(self):
        conv_acc_parameters = self.settings["coupling_settings"]["convergence_accelerator_settings"]
        convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(conv_acc_parameters)
        KratosMultiphysics.Logger.PrintInfo("::[ConjugateHeatTransferSolver]::", "Convergence accelerator construction finished.")
        return convergence_accelerator

    # This method returns the partitioned FSI utilities class
    def _get_partitioned_FSI_utilities(self):
        if (self.domain_size == 2):
            return KratosFSI.PartitionedFSIUtilitiesDouble2D()
        else:
            return KratosFSI.PartitionedFSIUtilitiesDouble3D()

    def _initialize_iteration_value_vector(self):
        # Initialize the iteration value for the residual computation
        self.iteration_value = KratosMultiphysics.Vector(
            self._get_partitioned_FSI_utilities().GetInterfaceResidualSize(self._get_dirichlet_coupling_interface()))
        i = 0
        for node in self._get_dirichlet_coupling_interface().Nodes:
            self.iteration_value[i] = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            i += 1

    def _temperature_coupling_prediction(self):
        # Create temperature prediction mapper
        mappers_settings = self.settings["coupling_settings"]["mappers_settings"]
        temp_pred_mapper_parameters = KratosMultiphysics.Parameters(r'''{
            "origin_variable": "TEMPERATURE",
            "destination_variable": "TEMPERATURE"
        }''')
        temp_pred_mapper_parameters.AddValue("echo_level", mappers_settings["echo_level"])
        temp_pred_mapper_parameters.AddValue("integration_order", mappers_settings["integration_order"])
        temp_pred_mapper_parameters.AddValue("distance_threshold", mappers_settings["distance_threshold"])
        temp_pred_mapper_parameters.AddValue("max_number_iterations", mappers_settings["max_number_iterations"])
        temp_pred_mapper_parameters.AddValue("absolute_convergence_tolerance", mappers_settings["absolute_convergence_tolerance"])
        temp_pred_mapper_parameters.AddValue("relative_convergence_tolerance", mappers_settings["relative_convergence_tolerance"])
        temp_pred_mapper_parameters.AddValue("search_parameters", mappers_settings["search_parameters"])

        # Map TEMPERATURE from the Neumann side to the Dirichlet side
        KratosMultiphysics.SimpleMortarMapperProcess(
            self._get_neumann_coupling_interface(),
            self._get_dirichlet_coupling_interface(),
            temp_pred_mapper_parameters).Execute()

    def _time_buffer_is_initialized(self):
        # Get current step counter. Note that the buoyancy and thermal fluid domain main_model_part
        # share the same ProcessInfo(), so the STEP counter is always the same.
        fluid_step = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        solid_step = self.solid_thermal_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

        if (fluid_step != solid_step):
            err_msg = "Fluid domain and solid domain time steps do not coincide\n"
            err_msg += "-Fluid domain time step: " + str(fluid_step) + "\n"
            err_msg += "-Solid domain time step: " + str(solid_step) + "\n"
            err_msg += "Check AdvanceInTime() method."
            raise Exception(err_msg)

        # We always have one extra old step (step 0, read from input)
        return fluid_step + 1 >= self.GetMinimumBufferSize()