# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtilities
import math
# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver


from pathlib import Path

def CreateSolver(model, custom_settings):
    return NavierStokesTwoFluidsHydraulicSolver(model, custom_settings)

class NavierStokesTwoFluidsHydraulicSolver(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "two_fluids_hydraulic",
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
            "time_scheme": "NEEDS TO BE FIXED --> PROBLEMS CONVERGENCE CRITERIO",
            "eulerian_fm_ale": true,
            "eulerian_fm_ale_settings":{
                "max_CFL" : 1.0,
                "max_substeps" : 0,
                "eulerian_error_compensation" : false,
                "element_type" : "levelset_convection_supg",
                "element_settings" : {
                    "dynamic_tau" : 1.0,
                    "tau_nodal":true
                }
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
            "distance_reinitialization": "variational",
            "parallel_redistance_max_layers" : 25,
            "distance_smoothing": false,
            "distance_smoothing_coefficient": 1.0,
            "distance_modification_settings": {
                "model_part_name": "",
                "distance_threshold": 1e-5,
                "continuous_distance": true,
                "check_at_each_time_step": true,
                "avoid_almost_empty_elements": false,
                "deactivate_full_negative_elements": false
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesTwoFluidsHydraulicSolver, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        super().__init__(model,custom_settings)

        self.element_name = "TwoFluidNavierStokesAlphaMethod"
        self.rho_inf = 0.0
        self.main_model_part.ProcessInfo.SetValue(KratosCFD.SPECTRAL_RADIUS_LIMIT, self.rho_inf)
        self.condition_name = "TwoFluidNavierStokesWallCondition"
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = True
        self.min_buffer_size = 2

        # Set the levelset characteristic variables and add them to the convection settings
        # These are required to be set as some of the auxiliary processes admit user-defined variables
        self._levelset_variable = KratosMultiphysics.DISTANCE
        self._levelset_gradient_variable = KratosMultiphysics.DISTANCE_GRADIENT
        self._levelset_convection_variable = KratosMultiphysics.VELOCITY
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_variable_name").SetString("DISTANCE")
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_gradient_variable_name").SetString("DISTANCE_GRADIENT")
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_convection_variable_name").SetString("VELOCITY")
        self.settings["levelset_convection_settings"].AddEmptyValue("convection_model_part_name").SetString("LevelSetConvectionModelPart")



        self.eulerian_fm_ale = self.settings["eulerian_fm_ale"].GetBool()
        if self.eulerian_fm_ale:
            self.fm_ale_variable = KratosCFD.CONVECTION_SCALAR
            self.eulerian_gradient = KratosCFD.CONVECTION_SCALAR_GRADIENT
            self.eulerian_convection_var = KratosCFD.CONVECTION_VELOCITY
            self.settings["eulerian_fm_ale_settings"].AddEmptyValue("levelset_variable_name").SetString("CONVECTION_SCALAR")
            self.settings["eulerian_fm_ale_settings"].AddEmptyValue("levelset_gradient_variable_name").SetString("CONVECTION_SCALAR_GRADIENT")
            self.settings["eulerian_fm_ale_settings"].AddEmptyValue("levelset_convection_variable_name").SetString("CONVECTION_VELOCITY")
            self.settings["eulerian_fm_ale_settings"].AddEmptyValue("convection_model_part_name").SetString("EulerianFMALEModelPart")



        dynamic_tau = self.settings["formulation"]["dynamic_tau"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, dynamic_tau)

        self.artificial_viscosity = self.settings["artificial_viscosity"].GetBool()
        if self.artificial_viscosity:
            self.artificial_limiter_coefficient = self.settings["artificial_visocosity_settings"]["limiter_coefficient"].GetDouble(
            )

        self._reinitialization_type = self.settings["distance_reinitialization"].GetString()

        # Initialize the system volue to None
        # Note that this will be computed only once in the first InitializeSolutionStep call
        self.__initial_system_volume = None

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesTwoFluidsHydraulicSolver finished.")

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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE) # Distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT) # Distance gradient nodal values

        if self.eulerian_fm_ale:
            # Auxiliary variable to store the historical scalar to be convected
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CONVECTION_SCALAR)
            # Auxiliary variable to store the velocity to be used in the historical data convection
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CONVECTION_VELOCITY)
            # Auxiliary variable to store the gradient of the historical scalar to be convected
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CONVECTION_SCALAR_GRADIENT)
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.AUXILIAR_VECTOR_VELOCITY)



        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def Initialize(self):
        computing_model_part = self.GetComputingModelPart()
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

        elemental_neighbour_search = KratosMultiphysics.GenericFindElementalNeighboursProcess(
            computing_model_part)
        elemental_neighbour_search.Execute()

        # Set and initialize the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        # Set nodal properties after setting distance(level-set).
        self._SetNodalProperties()

        # Initialize the distance correction process
        self._GetDistanceModificationProcess().ExecuteInitialize()
        self._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        # Instantiate the level set convection process
        # Note that is is required to do this in here in order to validate the defaults and set the corresponding distance gradient flag
        # Note that the nodal gradient of the distance is required either for the eulerian BFECC limiter or by the algebraic element antidiffusivity
        self._GetLevelSetConvectionProcess()

        # Instantiate the eulerian historical data convection
        if self.eulerian_fm_ale:
            self._GetEulerianFmAleProcess()

        if self.settings["formulation"].Has("mass_source"):
            self.mass_source = self.settings["formulation"]["mass_source"].GetBool()

        # Non historical variable are initilized in order to avoid memory problems
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.ACCELERATION, self.main_model_part.Nodes)
        if self.artificial_viscosity:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosCFD.ARTIFICIAL_DYNAMIC_VISCOSITY, self.main_model_part.Elements)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")
    def Check(self):
        super().Check()
        # Check if Inlet and Outlet boundary conditions are defined
        self._HydraulicBoundaryConditionCheck(KratosMultiphysics.INLET,"INLET")
        self._HydraulicBoundaryConditionCheck(KratosMultiphysics.OUTLET,"OUTLET")

    def InitializeSolutionStep(self):

        # Inlet and outlet water discharge is calculated for current time step, first discharge and the considering the time step inlet and outlet volume is calculated
        if self.mass_source:
            self._ComputeStepInitialWaterVolume()

        # Perform the convection of the historical database (Eulerian FM-ALE)
        if self.eulerian_fm_ale:
            self.__PerformEulerianFmAleVelocity()
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "FM-Lagrangian method is performed.")

        # Perform the level-set convection according to the previous step velocity
        self.__PerformLevelSetConvection()
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Level-set convection is performed.")

        # Perform distance correction to prevent ill-conditioned cuts
        self._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        # Update the DENSITY and DYNAMIC_VISCOSITY values according to the new level-set
        self._SetNodalProperties()

        # Accumulative water volume error ratio due to level set. Adding source term
        self._ComputeVolumeError()

        if self.artificial_viscosity:
            self.__CalculateArtificialViscosity()

        # Initialize the solver current step
        self._GetSolutionStrategy().InitializeSolutionStep()

        # We set this value at every time step as other processes/solvers also use them
        # Note that this is required as the convection processes may set a different value (this is the one to be used in the Navier-Stokes element)
        dynamic_tau = self.settings["formulation"]["dynamic_tau"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, dynamic_tau)

    # TODO:This is a test for do nothing outlet condition.

    # def SolveSolutionStep(self):
    #     n_it =5

    #     for k in range(n_it):
    #         is_converged = self._GetSolutionStrategy().SolveSolutionStep()
    #         for node in self.GetComputingModelPart().Nodes:
    #             if node.Is(KratosMultiphysics.OUTLET):
    #                 p=node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
    #                 node.SetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,p)
    #         if not is_converged:
    #             msg = "Fluid solver did not converge for step " + \
    #                 str(
    #                     self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
    #             msg += "corresponding to time " + \
    #                 str(
    #                     self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
    #             KratosMultiphysics.Logger.PrintWarning(
    #                 self.__class__.__name__, msg)
    #         return is_converged

    def FinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Mass and momentum conservation equations are solved.")

        # Recompute the distance field according to the new level-set position
        if self._reinitialization_type != "none":
            self._GetDistanceReinitializationProcess().Execute()
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Redistancing process is finished.")

        # Prepare distance correction for next step
        self._GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()

        # Finalize the solver current step
        self._GetSolutionStrategy().FinalizeSolutionStep()

        # Acceleration for generalised alpha time integration method.
        self.__CalculateTimeIntegrationAcceleration()

    def _ComputeStepInitialWaterVolume(self):

        # Here the initial water volume of the system is calculated without considering inlet and outlet flow rate
        if self.__initial_system_volume is None:
            self.__initial_system_volume = KratosCFD.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(
                self.GetComputingModelPart())

        # Calculate the inlet and outlet volume discharges
        outlet_discharge = KratosCFD.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(self.GetComputingModelPart(), KratosMultiphysics.OUTLET)
        inlet_discharge = KratosCFD.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(self.GetComputingModelPart(), KratosMultiphysics.INLET)
        current_dt = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        inlet_volume = -current_dt * inlet_discharge
        outlet_volume = current_dt * outlet_discharge

        # System water volume is calculated for current time step considering inlet and outlet discharge.
        system_volume = inlet_volume + self.__initial_system_volume - outlet_volume
        self.__initial_system_volume = system_volume



    def _ComputeVolumeError(self):
        if self.mass_source:
            water_volume_after_transport = KratosCFD.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(self.GetComputingModelPart())

            volume_error = (water_volume_after_transport -  self.__initial_system_volume) /  self.__initial_system_volume
        else:
            volume_error=0

        self.main_model_part.ProcessInfo.SetValue(KratosCFD.VOLUME_ERROR, volume_error)



    def __CalculateTimeIntegrationAcceleration(self):
        # Acceleration for generalised alpha time integration method.
        alpha_m=0.5*((3-self.rho_inf)/(1+self.rho_inf))
        alpha_f=1/(1+self.rho_inf)
        gamma= 0.5+alpha_m-alpha_f
        dt=self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        for node in self.GetComputingModelPart().Nodes:
            vn=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,1)
            v = node.GetSolutionStepValue(
                KratosMultiphysics.VELOCITY, 0)
            acceleration_n=node.GetValue(KratosMultiphysics.ACCELERATION)
            acceleration_alpha_method=(v-vn)/(gamma*dt)+((gamma-1)/gamma)*acceleration_n
            node.SetValue(KratosMultiphysics.ACCELERATION,acceleration_alpha_method)

    def __CalculateArtificialViscosity(self):

        properties_1 = self.main_model_part.Properties[1]
        water_dynamic_viscosity_max = self.artificial_limiter_coefficient * properties_1.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        for element in self.GetComputingModelPart().Elements:
            elem_artificial_viscosity = element.Calculate(KratosCFD.ARTIFICIAL_DYNAMIC_VISCOSITY, self.main_model_part.ProcessInfo)

            if elem_artificial_viscosity > water_dynamic_viscosity_max:
                elem_artificial_viscosity = water_dynamic_viscosity_max
            nodes = element.GetNodes()
            neg_nodes=0
            pos_nodes=0
            for node in nodes:
                distance=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                if distance>0:
                    pos_nodes += 1
                else:
                    neg_nodes+=1
            if neg_nodes>0 and pos_nodes>0:
                elem_artificial_viscosity= 0.0
            # elif pos_nodes==3:
            #     elem_artificial_viscosity = 0.0

            element.SetValue(KratosCFD.ARTIFICIAL_DYNAMIC_VISCOSITY, elem_artificial_viscosity)

    def __PerformLevelSetConvection(self):

        # Solve the levelset convection problem
        self._GetLevelSetConvectionProcess().Execute()

    def __PerformEulerianFmAleVelocity(self):

        # Solve the historical data convection problem
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        velocity_components = [KratosMultiphysics.VELOCITY_X,KratosMultiphysics.VELOCITY_Y,KratosMultiphysics.VELOCITY_Z]
        auxiliar_velocity_componentes = [KratosCFD.AUXILIAR_VECTOR_VELOCITY_X, KratosCFD.AUXILIAR_VECTOR_VELOCITY_Y, KratosCFD.AUXILIAR_VECTOR_VELOCITY_Z]
        mesh_var = [KratosMultiphysics.MESH_VELOCITY_X,KratosMultiphysics.MESH_VELOCITY_Y,KratosMultiphysics.MESH_VELOCITY_Z]

        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.VELOCITY,self.eulerian_convection_var, self.main_model_part, self.main_model_part, 0, 0)
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(KratosMultiphysics.VELOCITY,self.eulerian_convection_var, self.main_model_part, self.main_model_part, 1, 1)

        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(self.eulerian_gradient, self.main_model_part.Nodes)

        for i in range(domain_size):

            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(velocity_components[i], self.fm_ale_variable, self.main_model_part, self.main_model_part, 0, 0)
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(velocity_components[i], self.fm_ale_variable, self.main_model_part, self.main_model_part, 1,1)
            self._GetEulerianFmAleProcess().Execute()
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(self.fm_ale_variable, auxiliar_velocity_componentes[i], self.main_model_part, self.main_model_part, 0, 0)
            self.__CorrectVelocityHistory(velocity_components[i], mesh_var[i])
        self.__SlipConditionFixity()

    def __CorrectVelocityHistory(self,velocity_components, mesh_variable):
        for node in self.GetComputingModelPart().Nodes:
                velocity_history_corrected = node.GetSolutionStepValue(self.fm_ale_variable)
                if node.Is(KratosMultiphysics.SLIP):
                    pass

                elif not node.IsFixed(velocity_components):
                    node.SetSolutionStepValue(velocity_components, velocity_history_corrected)
                    node.SetSolutionStepValue(velocity_components, 1, velocity_history_corrected)
                    node.SetSolutionStepValue(mesh_variable, velocity_history_corrected)
                    node.SetSolutionStepValue(mesh_variable,1, velocity_history_corrected)

    def DotProduct(self,A, B):
        result = 0
        for i, j in zip(A, B):
            result += i*j
        return result

    def __SlipConditionFixity(self):
        for node in self.GetComputingModelPart().Nodes:
            if node.Is(KratosMultiphysics.SLIP):
                    n = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
                    n/=math.sqrt(n[0]**2+n[1]**2+n[2]**2)
                    # n /= math.sqrt(n[0]**2+n[1]**2)
                    v = node.GetSolutionStepValue(KratosCFD.AUXILIAR_VECTOR_VELOCITY)
                    v_prooj = self.DotProduct(v,n)
                    v -= v_prooj * n
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, v)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 1, v)
                    node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY, v)
                    node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY, 1,v)

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

    def _GetEulerianFmAleProcess(self):
        if not hasattr(self, '_eulerian_fm_convection_process'):
            self._eulerian_fm_convection_process = self._CreateFmAleConvectionProcess()
        return self._eulerian_fm_convection_process

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

    def _CreateFmAleConvectionProcess(self):
        # Construct the level set convection process
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        computing_model_part = self.GetComputingModelPart()
        linear_solver = self._GetLevelsetLinearSolver()
        eulerian_fm_ale_settings = self.settings["eulerian_fm_ale_settings"]
        if domain_size == 2:
            eulerian_fm_ale_process = KratosMultiphysics.LevelSetConvectionProcess2D(
                computing_model_part,
                linear_solver,
                eulerian_fm_ale_settings)
        else:
            eulerian_fm_ale_process = KratosMultiphysics.LevelSetConvectionProcess3D(
                computing_model_part,
                linear_solver,
                eulerian_fm_ale_settings)

        return eulerian_fm_ale_process

    def _CreateDistanceReinitializationProcess(self):
        # Construct the variational distance calculation process
        if (self._reinitialization_type == "variational"):
            maximum_iterations = 2 #TODO: Make this user-definable
            linear_solver = self._GetRedistancingLinearSolver()
            computing_model_part = self.GetComputingModelPart()
            if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                distance_reinitialization_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                    computing_model_part,
                    linear_solver,
                    maximum_iterations,
                    KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE)
            else:
                distance_reinitialization_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                    computing_model_part,
                    linear_solver,
                    maximum_iterations,
                    KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE)

        elif (self._reinitialization_type == "parallel"):
            #TODO: move all this to solver settings
            layers = self.settings["parallel_redistance_max_layers"].GetInt()
            parallel_distance_settings = KratosMultiphysics.Parameters("""{
                "max_levels" : 25,
                "max_distance" : 1.0,
                "calculate_exact_distances_to_plane" : true
            }""")
            parallel_distance_settings["max_levels"].SetInt(layers)
            if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                distance_reinitialization_process = KratosMultiphysics.ParallelDistanceCalculationProcess2D(
                    self.main_model_part,
                    parallel_distance_settings)
            else:
                distance_reinitialization_process = KratosMultiphysics.ParallelDistanceCalculationProcess3D(
                    self.main_model_part,
                    parallel_distance_settings)
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

    #TODO: Temporary solution for fixed alpha scheme
    def _CreateScheme(self):
        # "Fake" scheme for those cases in where the element manages the time integration
        # It is required to perform the nodal update once the current time step is solved
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(
            domain_size,
            domain_size + 1)

        return scheme

    def _HydraulicBoundaryConditionCheck(self,boundary,name):
        # Check if there are inlet and outlet
        computing_model_part = self.GetComputingModelPart()
        not_boundary_nodes=any([node.Is(boundary) for node in computing_model_part.Nodes])
        if not not_boundary_nodes:
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, name +" condition is not defined in the model part.")



