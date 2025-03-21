# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtilities

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.FluidDynamicsHydraulicsApplication as KratosHydraulics

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

import math
from pathlib import Path

def CreateSolver(model, custom_settings):
    return NavierStokesTwoFluidsHydraulicFractionalSolver(model, custom_settings)

class NavierStokesTwoFluidsHydraulicFractionalSolver(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "two_fluid_hydraulic_fractional",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "maximum_iterations": 7,
            "echo_level": 0,
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": false,
            "consider_periodic_conditions": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "assign_neighbour_elements_to_conditions": true,
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-2,
                "maximum_delta_time"  : 1.0,
                "time_step"           : 0.0
            },
            "move_mesh_flag": false,
            "formulation": {
                "dynamic_tau": 1.0,
                "mass_source":true
            },
            "artificial_viscosity": false,
            "artificial_visocosity_settings":{
                "limiter_coefficient": 1000
            },
            "energy_measurement":true,
            "file_name" : "energy.txt",
            "time_scheme": "bdf2",
            "fractional_splitting_settings":{
                "element_type" : "ns_fractional_velocity_convection"
             },
            "levelset_convection_settings": {
                "max_CFL" : 1.0,
                "max_substeps" : 0,
                "eulerian_error_compensation" : false,
                "element_type" : "levelset_convection_supg",
                "element_settings" : {
                    "dynamic_tau" : 1.0,
                    "tau_nodal":true
                }
            },
            "distance_reinitialization_type" :"variational",
            "distance_reinitialization_settings":{
            },
            "distance_modification_settings": {
                "model_part_name": "",
                "distance_threshold": 1e-5,
                "continuous_distance": true,
                "check_at_each_time_step": true,
                "avoid_almost_empty_elements": false,
                "deactivate_full_negative_elements": false
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesTwoFluidsHydraulicFractionalSolver, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        super().__init__(model,custom_settings)

        # TODO: At this moment only thd bdf2 is available but in a future alpha fractional element will be merged.
        self.element_name = "TwoFluidNavierStokesFractional"
        self.min_buffer_size = 3
        self.condition_name = "TwoFluidNavierStokesWallCondition"
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = True

        # Set the levelset characteristic variables and add them to the convection settings
        # These are required to be set as some of the auxiliary processes admit user-defined variables
        self._levelset_variable = KratosMultiphysics.DISTANCE
        self._levelset_gradient_variable = KratosMultiphysics.DISTANCE_GRADIENT
        self._levelset_convection_variable = KratosMultiphysics.VELOCITY
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_variable_name").SetString("DISTANCE")
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_gradient_variable_name").SetString("DISTANCE_GRADIENT")
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_convection_variable_name").SetString("VELOCITY")
        self.settings["levelset_convection_settings"].AddEmptyValue("convection_model_part_name").SetString("LevelSetConvectionModelPart")

        self.settings["fractional_splitting_settings"].AddEmptyValue("model_part_name").SetString(self.main_model_part.Name )


        dynamic_tau = self.settings["formulation"]["dynamic_tau"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, dynamic_tau)

        self.artificial_viscosity = self.settings["artificial_viscosity"].GetBool()
        if self.artificial_viscosity:
            self.artificial_limiter_coefficient = self.settings["artificial_visocosity_settings"]["limiter_coefficient"].GetDouble()
        self._reinitialization_type = self.settings["distance_reinitialization_type"].GetString()
        # Note that this will be computed only once in the first InitializeSolutionStep call
        self.__initial_water_system_volume = None

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesTwoFluidsHydraulicFractionalSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.OUTLET_NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE) # Distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT) # Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.AUXILIAR_VECTOR_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.FRACTIONAL_VELOCITY)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def AddDofs(self):
        super().AddDofs()

        dofs_and_reactions_to_add = []
        dofs_and_reactions_to_add.append("FRACTIONAL_VELOCITY_X")
        dofs_and_reactions_to_add.append("FRACTIONAL_VELOCITY_Y")
        dofs_and_reactions_to_add.append("FRACTIONAL_VELOCITY_Z")

        KratosMultiphysics.VariableUtils.AddDofsList( dofs_and_reactions_to_add, self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "DOF added correctly.")

    def Initialize(self):
        computing_model_part = self.GetComputingModelPart()

        for prop in self.main_model_part.Properties:
            print(f"Property ID: {prop.Id}")
        # Calculate boundary normals
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(
            computing_model_part,
            computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Finding nodal and elemental neighbors
        data_communicator = computing_model_part.GetCommunicator().GetDataCommunicator()
        neighbour_search = KratosMultiphysics.FindGlobalNodalNeighboursProcess(
            data_communicator,
            computing_model_part)
        neighbour_search.Execute()
        elemental_neighbour_search = KratosMultiphysics.GenericFindElementalNeighboursProcess(computing_model_part)
        elemental_neighbour_search.Execute()

        # Set and Initialize the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

        # Set nodal properties after setting distance(level-set).
        self._SetNodalProperties()

        # Initialize the distance correction process
        self._GetDistanceModificationProcess().ExecuteInitialize()
        self._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()


        # Instantiate the level set convection process
        # Note that is is required to do this in here in order to validate the defaults and set the corresponding distance gradient flag
        # Note that the nodal gradient of the distance is required either for the eulerian BFECC limiter or by the algebraic element antidiffusivity
        self._GetLevelSetConvectionProcess()
        self._GetNSFractionalSplittingProcess()

        # Set the flag for the mass loss correction
        if self.settings["formulation"].Has("mass_source"):
            self.mass_source = self.settings["formulation"]["mass_source"].GetBool()

        # Initialize non historical artificial velocity to zero
        if self.artificial_viscosity:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.ARTIFICIAL_DYNAMIC_VISCOSITY, self.main_model_part.Elements)

        self.domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

        self.previous_dt = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self.energy_process_activation = self.settings["energy_measurement"].GetBool(
        )
        if self.energy_process_activation:
            self.post_file_name = self.settings["file_name"].GetString()
            self.my_energy_process = KratosCFD.EnergyCheckProcess(
                self.main_model_part, self.domain_size, self.post_file_name)


    def Check(self):
        super().Check()
        # Check if Inlet and Outlet boundary conditions are defined
        self._HydraulicBoundaryConditionCheck(KratosMultiphysics.INLET,"INLET")
        self._HydraulicBoundaryConditionCheck(KratosMultiphysics.OUTLET,"OUTLET")

    def InitializeSolutionStep(self):

        # Inlet and outlet water discharge is calculated for current time step, first discharge and the considering the time step inlet and outlet volume is calculated
        if self.mass_source:
            self._ComputeStepInitialWaterVolume()

        # Recompute the BDF2 coefficients
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

        # STEP I: NS Fractional part 1
        # Perform the pure convection of the fractional velocity which corresponds to the first part of the NS fractional splitting.
        self.__PerformNSFractionalSplitting()
        self.vectorial_convection_iterations = self.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER]
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Navier Stokes fractional convection part is performed.")

        # STEP II: Convect the free surface according to the fractional velocity
        # Before doing this second step, the fractional velocity data is copied to the velocity data since the level set convection process takes velocity variable as convection variable.
        # And the previous previous velocity is copied in an auxiliar variable
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosCFD.FRACTIONAL_VELOCITY,KratosMultiphysics.VELOCITY, self.main_model_part, self.main_model_part, 0, 0)
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.VELOCITY,KratosCFD.AUXILIAR_VECTOR_VELOCITY, self.main_model_part, self.main_model_part, 0, 0)

        self.__PerformLevelSetConvection()
        self.levelset_iterations = self.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER]
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Level-set convection is performed.")

        # After the convection process, the velocity is copied back to the original state.
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosCFD.AUXILIAR_VECTOR_VELOCITY,KratosMultiphysics.VELOCITY, self.main_model_part, self.main_model_part, 0, 0)

        # Perform distance correction to prevent ill-conditioned cuts
        self._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        # Update the DENSITY and DYNAMIC_VISCOSITY values according to the new level-set
        self._SetNodalProperties()

        # Accumulative water volume error ratio due to level set. Adding source term
        self._ComputeVolumeError()

        # Calculate residual-based artificial viscosity
        if self.artificial_viscosity:
            self.__CalculateArtificialViscosity()

        # Initialize the solver current step
        self._GetSolutionStrategy().InitializeSolutionStep()

        # We set this value at every time step as other processes/solvers also use them
        # Note that this is required as the convection processes may set a different value (this is the one to be used in the Navier-Stokes element)
        dynamic_tau = self.settings["formulation"]["dynamic_tau"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, dynamic_tau)

    def FinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Mass and momentum conservation equations are solved.")

        # Recompute the distance field according to the new level-set position
        if self._reinitialization_type != "none":
            self._GetDistanceReinitializationProcess().Execute()
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Redistancing process is finished.")

        # Prepare distance correction for next step
        self._GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()

        # FinalizeSolutionStep of Navier-Stokes strategy
        self._GetSolutionStrategy().FinalizeSolutionStep()

        if self.energy_process_activation:
            self.my_energy_process.Execute()

    def _ComputeStepInitialWaterVolume(self):

        # This function calculates the theoretical water volume at each time step.
        # Reminder: Despite adding the source term to both air and water, the absolute volume error
        # is referenced to the water volume, since what is lost from water is gained by air and vice versa.

        # Here the initial water volume of the system is calculated without considering inlet and outlet flow rate
        if self.__initial_water_system_volume is None:
            self.__initial_water_system_volume = KratosCFD.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(self.GetComputingModelPart())

        current_dt = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        # Calculate the inlet and outlet volume discharges
        water_outlet_discharge = KratosCFD.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(self.GetComputingModelPart(), KratosMultiphysics.OUTLET)
        water_inlet_discharge = KratosCFD.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(self.GetComputingModelPart(), KratosMultiphysics.INLET)
        KratosMultiphysics.Logger.PrintWarning( self.__class__.__name__, str(water_inlet_discharge) + "m3/s is the inlet discharge")
        inlet_water_volume = -current_dt * water_inlet_discharge
        outlet_water_volume = current_dt * water_outlet_discharge
        system_water_volume = inlet_water_volume + self.__initial_water_system_volume - outlet_water_volume

        # System water volume is calculated for current time step considering inlet and outlet discharge.
        self.__initial_water_system_volume = system_water_volume

    def _ComputeVolumeError(self):
        # In this function, the volume of the cut elements is calculated,
        # corresponding to the portions of water and air volumes, as this is the domain where the source term will be added.
        # Meanwhile, the absolute error is calculated within the water domain

        if self.mass_source:
            water_volume_after_transport = KratosCFD.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(self.GetComputingModelPart())
            water_cut_volume_after_transport = KratosCFD.FluidAuxiliaryUtilities.CalculateFluidCutElementsNegativeVolume(self.GetComputingModelPart())
            absolute_water_volume_error = water_volume_after_transport -  self.__initial_water_system_volume
            water_cut_element_error = absolute_water_volume_error /water_cut_volume_after_transport
            air_cut_volume_after_transport = KratosCFD.FluidAuxiliaryUtilities.CalculateFluidCutElementsPositiveVolume(self.GetComputingModelPart())
            air_cut_element_error = absolute_water_volume_error / air_cut_volume_after_transport
        else:
            water_cut_element_error = 0.0
            air_cut_element_error =0.0

        self.main_model_part.ProcessInfo.SetValue(KratosCFD.WATER_VOLUME_ERROR, water_cut_element_error)
        self.main_model_part.ProcessInfo.SetValue(KratosCFD.AIR_VOLUME_ERROR, air_cut_element_error)

    def __CalculateArtificialViscosity(self):
        properties_1 = self.main_model_part.Properties[1]
        water_dynamic_viscosity_max = self.artificial_limiter_coefficient * properties_1.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        KratosHydraulics.HydraulicFluidAuxiliaryUtilities.CalculateNonIntersectedElementsArtificialViscosity(self.main_model_part, water_dynamic_viscosity_max)

    def __PerformLevelSetConvection(self):
        # Solve the levelset convection problem
        self._GetLevelSetConvectionProcess().Execute()

    def __PerformNSFractionalSplitting(self):
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        velocity_components = [KratosMultiphysics.VELOCITY_X,KratosMultiphysics.VELOCITY_Y,KratosMultiphysics.VELOCITY_Z]
        fractional_velocity_components = [KratosCFD.FRACTIONAL_VELOCITY_X, KratosCFD.FRACTIONAL_VELOCITY_Y, KratosCFD.FRACTIONAL_VELOCITY_Z]
        # Trasfer velocity node fixity to fractional velocity
        for i in range(domain_size):
            self.VelocityBoundaryConditionFractional(fractional_velocity_components[i], velocity_components[i])
        self._GetNSFractionalSplittingProcess().Execute()
        # Trasfer velocity slip condition to fractional velocity
        self.__SlipConditonFractionalFixity()
    # TODO: Remove those methods as soon as a new  hydraulic slip process is done.
    def __SlipConditonFractionalFixity(self):
        for node in self.GetComputingModelPart().Nodes:
            if node.Is(KratosMultiphysics.SLIP):
                n= node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
                n/= math.sqrt(n[0]**2+n[1]**2+n[2]**2)
                v= node.GetSolutionStepValue(KratosCFD.FRACTIONAL_VELOCITY)
                v_prooj = self.DotProduct(v,n)
                v-= v_prooj*n
                node.SetSolutionStepValue(KratosCFD.FRACTIONAL_VELOCITY,v)
    def VelocityBoundaryConditionFractional(self, fractional_velocity_componentes, velocity_components):
        for node in self.GetComputingModelPart().Nodes:
            if node.IsFixed(velocity_components):
                v_fix = node.GetSolutionStepValue(velocity_components)
                node.SetSolutionStepValue(fractional_velocity_componentes,v_fix)
                node.Fix(fractional_velocity_componentes)

    def DotProduct(self,A, B):
        result = 0
        for i, j in zip(A, B):
            result += i*j
        return result

    # TODO: Remove this method as soon as the subproperties are available
    def _SetPhysicalProperties(self):

        warn_msg  = '\nThe materials import mechanism used in the two fluids solver is DEPRECATED!\n'
        warn_msg += 'It will be removed to use the base fluid_solver.py one as soon as the subproperties are available.\n'
        KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)

        # Check if the fluid properties are provided using a .json file
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            data_comm = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator() # only using the global comm as the Communicators are not yet created when running in MPI. Hotfix since this method will disappear completely when using subproperties!

            def GetAuxMaterialsFileName(mat_file_name, prop_id):
                p_mat_file_name = Path(mat_file_name)
                new_stem = "{}_p{}".format(p_mat_file_name.stem, prop_id)
                return str(p_mat_file_name.with_name(new_stem).with_suffix(p_mat_file_name.suffix))

            with open(materials_filename,'r') as materials_file:
                materials = KratosMultiphysics.Parameters(materials_file.read())

            if data_comm.Rank() == 0:
                # Create and read an auxiliary materials file for each one of the fields (only on one rank)
                for i_material in materials["properties"]:
                    aux_materials = KratosMultiphysics.Parameters()
                    aux_materials.AddEmptyArray("properties")
                    aux_materials["properties"].Append(i_material)

                    aux_materials_filename = GetAuxMaterialsFileName(materials_filename, i_material["properties_id"].GetInt())
                    with open(aux_materials_filename,'w') as aux_materials_file:
                        aux_materials_file.write(aux_materials.WriteJsonString())

            data_comm.Barrier()

            # read the files on all ranks
            for i_material in materials["properties"]:
                aux_materials_filename = GetAuxMaterialsFileName(materials_filename, i_material["properties_id"].GetInt())
                aux_material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
                aux_material_settings["Parameters"]["materials_filename"].SetString(aux_materials_filename)
                KratosMultiphysics.ReadMaterialsUtility(aux_material_settings, self.model)

            data_comm.Barrier()

            if data_comm.Rank() == 0:
                # remove aux files after every rank read them
                for i_material in materials["properties"]:
                    aux_materials_filename = GetAuxMaterialsFileName(materials_filename, i_material["properties_id"].GetInt())
                    KratosUtilities.DeleteFileIfExisting(aux_materials_filename)

            materials_imported = True
        else:
            materials_imported = False

        # If the element uses nodal material properties, transfer them to the nodes
        if self.element_has_nodal_properties:
            self._SetNodalProperties()

        return materials_imported

    def _SetNodalProperties(self):
        # Get fluid 1 and 2 properties
        properties_1 = self.main_model_part.Properties[1]
        properties_2 = self.main_model_part.Properties[2]

        rho_1 = properties_1.GetValue(KratosMultiphysics.DENSITY)
        rho_2 = properties_2.GetValue(KratosMultiphysics.DENSITY)
        mu_1 = properties_1.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        mu_2 = properties_2.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)

        # Check fluid 1 and 2 properties
        if rho_1 <= 0.0:
            raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho_1, properties_1.Id))
        if rho_2 <= 0.0:
            raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho_2, properties_2.Id))
        if mu_1 <= 0.0:
            raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(mu_1, properties_1.Id))
        if mu_2 <= 0.0:
            raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(mu_2, properties_2.Id))

        # Transfer density and (dynamic) viscostity to the nodes
        for node in self.main_model_part.Nodes:
            if node.GetSolutionStepValue(self._levelset_variable) <= 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho_1)
                node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, mu_1)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho_2)
                node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, mu_2)

    def _GetRedistancingLinearSolver(self):
        # A linear solver configured specifically for distance re-initialization process
        if not hasattr(self, '_redistancing_linear_solver'):
            self._redistancing_linear_solver = self._CreateLinearSolver() # TODO: add customized configuration
        return self._redistancing_linear_solver

    def _GetLevelsetLinearSolver(self):
        # A linear solver configured specifically for the level-set convection process
        if not hasattr(self, '_levelset_linear_solver'):
            self._levelset_linear_solver = self._CreateLinearSolver() # TODO: add customized configuration
        return self._levelset_linear_solver

    def _GetLevelSetConvectionProcess(self):
        if not hasattr(self, '_level_set_convection_process'):
            self._level_set_convection_process = self._CreateLevelSetConvectionProcess()
        return self._level_set_convection_process

    def _GetNSFractionalSplittingProcess(self):
        if not hasattr(self, '_ns_fractional_splitting_process'):
            self._ns_fractional_splitting_process = self._CreateFractionalNSplittingProcess()
        return self._ns_fractional_splitting_process

    def _GetDistanceReinitializationProcess(self):
        if not hasattr(self, '_distance_reinitialization_process'):
            self._distance_reinitialization_process = self._CreateDistanceReinitializationProcess()
        return self._distance_reinitialization_process

    def _GetDistanceModificationProcess(self):
        if not hasattr(self, '_distance_modification_process'):
            self._distance_modification_process = self.__CreateDistanceModificationProcess()
        return self._distance_modification_process

    def _CreateLevelSetConvectionProcess(self):
        # Construct the level set convection process
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        computing_model_part = self.GetComputingModelPart()
        linear_solver = self._GetLevelsetLinearSolver()
        levelset_convection_settings = self.settings["levelset_convection_settings"]
        if domain_size == 2:
            level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess2D(
                computing_model_part,
                linear_solver,
                levelset_convection_settings)
        else:
            level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess3D(
                computing_model_part,
                linear_solver,
                levelset_convection_settings)

        return level_set_convection_process

    def _CreateFractionalNSplittingProcess(self):

        # Construct the level set convection process
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        linear_solver = self._GetLevelsetLinearSolver()
        fractional_splitting_settings = self.settings["fractional_splitting_settings"]
        if domain_size == 2:
            fractional_splitting_process = KratosCFD.TwoFluidNavierStokesFractionalConvectionProcess2D(
            self.model, linear_solver, fractional_splitting_settings)
        else:
            fractional_splitting_process = KratosCFD.TwoFluidNavierStokesFractionalConvectionProcess3D(
             self.model, linear_solver, fractional_splitting_settings)
        return fractional_splitting_process


    def _CreateDistanceReinitializationProcess(self):
        # Construct the variational distance calculation process
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        distance_reinitialization_settings = self.settings["distance_reinitialization_settings"]
        distance_reinitialization_settings.AddEmptyValue("model_part_name").SetString(self.main_model_part.Name)
        if (self._reinitialization_type == "variational"):
            linear_solver = self._GetRedistancingLinearSolver()
            if domain_size == 2:
                distance_reinitialization_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                    self.model,
                    linear_solver,
                    distance_reinitialization_settings)
            else:
                distance_reinitialization_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                    self.model,
                    linear_solver,
                    distance_reinitialization_settings)

        elif (self._reinitialization_type == "parallel"):
            if domain_size == 2:
                distance_reinitialization_process = KratosMultiphysics.ParallelDistanceCalculationProcess2D(
                    self.main_model_part,
                    distance_reinitialization_settings)
            else:
                distance_reinitialization_process = KratosMultiphysics.ParallelDistanceCalculationProcess3D(
                    self.main_model_part,
                    distance_reinitialization_settings)
        elif (self._reinitialization_type == "none"):
                KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Redistancing is turned off.")
        else:
            raise Exception("Please use a valid distance reinitialization type or set it as \'none\'. Valid types are: \'variational\' and \'parallel\'.")

        return distance_reinitialization_process

    def __CreateDistanceModificationProcess(self):
        # Set suitable distance correction settings for free-surface problems
        # Note that the distance modification process is applied to the computing model part
        distance_modification_settings = self.settings["distance_modification_settings"]
        distance_modification_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["distance_modification_settings"])
        distance_modification_settings["model_part_name"].SetString(self.GetComputingModelPart().FullName())

        # Check user provided settings
        if not distance_modification_settings["continuous_distance"].GetBool():
            distance_modification_settings["continuous_distance"].SetBool(True)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'continuous_distance\' is \'False\'. Setting to \'True\'.")
        if not distance_modification_settings["check_at_each_time_step"].GetBool():
            distance_modification_settings["check_at_each_time_step"].SetBool(True)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'check_at_each_time_step\' is \'False\'. Setting to \'True\'.")
        if distance_modification_settings["avoid_almost_empty_elements"].GetBool():
            distance_modification_settings["avoid_almost_empty_elements"].SetBool(False)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'avoid_almost_empty_elements\' is \'True\'. Setting to \'False\' to avoid modifying the distance sign.")
        if distance_modification_settings["deactivate_full_negative_elements"].GetBool():
            distance_modification_settings["deactivate_full_negative_elements"].SetBool(False)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'deactivate_full_negative_elements\' is \'True\'. Setting to \'False\' to avoid deactivating the negative volume (e.g. water).")

        # Create and return the distance correction process
        return KratosCFD.DistanceModificationProcess(
            self.model,
            distance_modification_settings)

    def _HydraulicBoundaryConditionCheck(self,boundary,name):
        # Check if there are inlet and outlet
        computing_model_part = self.GetComputingModelPart()
        not_boundary_nodes=any([node.Is(boundary) for node in computing_model_part.Nodes])
        if not not_boundary_nodes:
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, name +" condition is not defined in the model part.")
    def SolveSolutionStep(self):
        is_converged = super().SolveSolutionStep()
        self.ns_iterations = self.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER]

        self.old_iterations = max(self.levelset_iterations,self.vectorial_convection_iterations,self.ns_iterations)

        return is_converged

    def _ComputeDeltaTime(self):
        dt = super()._ComputeDeltaTime()
        self.dt = max(dt,self.previous_dt)
        if dt > self.previous_dt:
            self.dt =self.previous_dt
        step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        if step>2.0:
            maximum_iterations = self.settings["maximum_iterations"].GetDouble()
            iterations_ratio = self.old_iterations/maximum_iterations
            alpha = 0.8
            if iterations_ratio < 0.5:
                self.previous_dt /= alpha
                self.dt = self.previous_dt
            elif iterations_ratio > 0.9:
                self.previous_dt *=alpha
                self.dt = self.previous_dt
        return self.dt
