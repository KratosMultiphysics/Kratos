from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import main_script as DEM_algorithm
import swimming_DEM_procedures as SDP

sys.path.insert(0,'')
import DEM_explicit_solver_var as DEM_parameters
BaseAlgorithm = DEM_algorithm.Solution

class Algorithm(BaseAlgorithm):
    def __init__(self, pp):
        self.StartTimer()
        self.pp = pp
        self.SetBetaParamters()
        self.SetDoSolveDEMVariable()

    def CreateParts(self):
        # Order must be respected here
        BaseAlgorithm.__init__(self)
        # defining a fluid model
        self.all_model_parts.Add(ModelPart("FluidPart"))
        # defining a model part for the mixed part
        self.all_model_parts.Add(ModelPart("MixedPart"))

    def StartTimer(self):
        import time as timer
        self.timer = timer
        self.simulation_start_time = timer.time()

    def SetBetaParamters(self): # These are input parameters that have not yet been transferred to the interface
        # import the configuration data as read from the GiD
        self.main_path = os.getcwd()
        self.pp.main_path = os.getcwd()



        ##############################################################################
        #                                                                            #
        #    INITIALIZE                                                              #
        #                                                                            #
        ##############################################################################

        #G
        self.pp.CFD_DEM = DEM_parameters
        self.pp.CFD_DEM.fluid_already_calculated = 0
        self.pp.CFD_DEM.recovery_echo_level = 1
        self.pp.CFD_DEM.gradient_calculation_type = 1
        self.pp.CFD_DEM.pressure_grad_recovery_type = 1
        self.pp.CFD_DEM.store_full_gradient = 0
        self.pp.CFD_DEM.laplacian_calculation_type = 0
        self.pp.CFD_DEM.do_search_neighbours = False
        self.pp.CFD_DEM.faxen_terms_type = 0
        self.pp.CFD_DEM.material_acceleration_calculation_type = 1
        self.pp.CFD_DEM.faxen_force_type = 0
        self.pp.CFD_DEM.vorticity_calculation_type = 5
        self.pp.CFD_DEM.print_FLUID_VEL_PROJECTED_RATE_option = 0
        self.pp.CFD_DEM.print_MATERIAL_FLUID_ACCEL_PROJECTED_option = True
        self.pp.CFD_DEM.basset_force_type = 0
        self.pp.CFD_DEM.print_BASSET_FORCE_option = 1
        self.pp.CFD_DEM.basset_force_integration_type = 2
        self.pp.CFD_DEM.n_init_basset_steps = 0
        self.pp.CFD_DEM.time_steps_per_quadrature_step = 1
        self.pp.CFD_DEM.delta_time_quadrature = self.pp.CFD_DEM.time_steps_per_quadrature_step * self.pp.CFD_DEM.MaxTimeStep
        self.pp.CFD_DEM.quadrature_order = 2
        self.pp.CFD_DEM.time_window = 0.1
        self.pp.CFD_DEM.number_of_exponentials = 10
        self.pp.CFD_DEM.number_of_quadrature_steps_in_window = int(self.pp.CFD_DEM.time_window / self.pp.CFD_DEM.delta_time_quadrature)
        self.pp.CFD_DEM.print_steps_per_plot_step = 1
        self.pp.CFD_DEM.PostCationConcentration = False
        self.pp.CFD_DEM.do_impose_flow_from_field = True
        self.pp.CFD_DEM.print_MATERIAL_ACCELERATION_option = True
        self.pp.CFD_DEM.print_FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED_option = False
        self.pp.CFD_DEM.print_VELOCITY_GRADIENT_option = 1
        self.pp.CFD_DEM.print_VORTICITY_option = 1
        self.pp.CFD_DEM.print_MATERIAL_ACCELERATION_option = True
        # Making the fluid step an exact multiple of the DEM step
        self.pp.Dt = int(self.pp.Dt / self.pp.CFD_DEM.MaxTimeStep) * self.pp.CFD_DEM.MaxTimeStep
        self.pp.viscosity_modification_type = 0.0
        self.domain_size = 3

        # defining and adding imposed porosity fields
        import swimming_DEM_procedures as SDP
        self.pp.fluid_fraction_fields = []
        field1 = SDP.FluidFractionFieldUtility.LinearField(0.0,
                                                          [0.0, 0.0, 0.0],
                                                          [-1.0, -1.0, 0.15],
                                                          [1.0, 1.0, 0.3])
        from math import pi
        self.pp.CFD_DEM.fluid_domain_volume = 0.5 ** 2 * 2 * pi # write down the volume you know it has

        self.pp.fluid_fraction_fields.append(field1)

    def SetDoSolveDEMVariable(self):
        self.pp.do_solve_dem = not self.pp.CFD_DEM.flow_in_porous_DEM_medium_option

    def SetCustomBetaParamters(self, dictionary):
        if len(dictionary) == 0:
            return
        else: # assign the specified values to the specified variables
            var_names = [k for k in dictionary.keys()]
            var_values = [k for k in dictionary.values()]
            for name, value in zip(var_names, var_values):
                globals()['self.pp.CFD_DEM.' + name] = value

    def SetUpResultsDatabase(self):
        pass

    def ReadModelParts(self, starting_node_Id = 0, starting_elem_Id = 0, starting_cond_Id = 0):
        self.ReadFluidModelPart()
        fluid_mp = self.all_model_parts.Get('FluidPart')
        max_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(fluid_mp)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(fluid_mp)
        max_cond_Id = self.creator_destructor.FindMaxConditionIdInModelPart(fluid_mp)
        self.ReadDEMModelParts(max_node_Id + 1, max_elem_Id + 1, max_cond_Id + 1)

    def ReadFluidModelPart(self):
        os.chdir(self.pp.main_path)
        model_part_io_fluid = ModelPartIO(self.pp.problem_name)
        model_part_io_fluid.ReadModelPart(self.all_model_parts.Get('FluidPart'))

    def ReadDEMModelParts(self, starting_node_Id = 0, starting_elem_Id = 0, starting_cond_Id = 0):
        BaseAlgorithm.ReadModelParts(self, starting_node_Id, starting_elem_Id, starting_cond_Id)

    def Initialize(self):
        self.fluid_solver.Initialize()
        BaseAlgorithm.Initialize(self)

    def AddExtraVariables(self):
        import variables_management as vars_man
        fluid_model_part = self.all_model_parts.Get('FluidPart')

                # building lists of variables for which memory is to be allocated
        # TEMPORARY, HORRIBLE !!!
        vars_man.ConstructListsOfVariables(self.pp)
        #_____________________________________________________________________________________________________________________________________
        #
        #                               F L U I D    B L O C K    B E G I N S
        #_____________________________________________________________________________________________________________________________________

        # defining variables to be used
        # GID IO IS NOT USING THIS NOW. TO BE REMOVED ONCE THE "PRINT IN POINTS"
        # CODE IS NOT USING IT

        variables_dictionary = {"PRESSURE"   : PRESSURE,
                                "VELOCITY"   : VELOCITY,
                                "MU"         : MU,         #    MOD.
                                "BUOYANCY"   : BUOYANCY,   #    MOD.
                                "DRAG_FORCE" : DRAG_FORCE,  #    MOD.
                                "LIFT_FORCE" : LIFT_FORCE} #    MOD.

        fluid_model_part = self.all_model_parts.Get('FluidPart')

        if "REACTION" in self.pp.nodal_results:
            fluid_model_part.AddNodalSolutionStepVariable(REACTION)
        if "DISTANCE" in self.pp.nodal_results:
            fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

        # importing the solvers needed
        SolverSettings = self.pp.FluidSolverConfiguration
        self.solver_module = import_solver(SolverSettings)

        # caution with breaking up this block (memory allocation)! {
        self.solver_module.AddVariables(fluid_model_part, SolverSettings)
        vars_man.AddNodalVariables(fluid_model_part, self.pp.fluid_vars)  #     MOD.
        # }



        # Creating necessary directories
        [post_path, data_and_results, graphs_path, MPI_results] = self.procedures.CreateDirectories(str(self.main_path), str(self.pp.CFD_DEM.problem_name))

        #_____________________________________________________________________________________________________________________________________
        #
        #                               F L U I D    B L O C K    E N D S
        #_____________________________________________________________________________________________________________________________________

        # Add variables

        vars_man.AddNodalVariables(self.spheres_model_part, self.pp.dem_vars)
        vars_man.AddNodalVariables(self.rigid_face_model_part, self.pp.rigid_faces_vars)
        vars_man.AddNodalVariables(self.DEM_inlet_model_part, self.pp.inlet_vars)
        # adding extra process info variables
        vars_man.AddingExtraProcessInfoVariables(self.pp, fluid_model_part, self.spheres_model_part)

        fluid_model_part.SetBufferSize(3)
        self.solver_module.AddDofs(fluid_model_part, SolverSettings)
        SDP.AddExtraDofs(self.pp, fluid_model_part, self.spheres_model_part, self.cluster_model_part, self.DEM_inlet_model_part)

        os.chdir(self.main_path)

        self.KRATOSprint("\nInitializing Problem...")

    def DEMSolve(self, time = 'None'): # time is passed in case it is needed
        self.solver.Solve()

    def PerformZeroStepInitializations(self):
        pass

    def PerformInitialDEMStepOperations(self, time = None):
        pass

    def SetSolverStrategy(self):
        import swimming_sphere_strategy as SolverStrategy
        return SolverStrategy

    def SelectScheme(self):
        scheme = BaseAlgorithm.SelectScheme(self)
        if scheme == None:
            if self.pp.CFD_DEM.IntegrationScheme == 'Hybrid_Bashforth':
                return HybridBashforthScheme()
            elif self.pp.CFD_DEM.IntegrationScheme == "TerminalVelocityScheme":
                return TerminalVelocityScheme()
            else:
                return None
        else:
            return scheme

    def SetSolver(self):
        return self.solver_strategy.SwimmingStrategy(self.all_model_parts, self.creator_destructor, self.dem_fem_search, self.scheme, self.pp.CFD_DEM, self.procedures)


    def TellTime(self, time):
        print("\n", "TIME = ", time)
        print('ELAPSED TIME = ', self.timer.time() - self.simulation_start_time)
        print()
        sys.stdout.flush()

    def TellFinalSummary(self, time, step, DEM_step):
        print()
        print()
        print('*************************************************************')
        print("CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.")
        simulation_elapsed_time = self.timer.time() - self.simulation_start_time
        print("Elapsed time: " + "%.5f"%(simulation_elapsed_time) + " s ")
        print("per fluid time step: " + "%.5f"%(simulation_elapsed_time / step) + " s ")
        print("per DEM time step: " + "%.5f"%(simulation_elapsed_time / DEM_step) + " s")
        print('*************************************************************')
        print()
        sys.stdout.flush()

    def GetFluidSolveCounter(self):
        return SDP.Counter(is_dead = (self.pp.CFD_DEM.drag_force_type == 9))

    def GetEmbeddedCounter(self):
        return SDP.Counter(1, 3, self.pp.CFD_DEM.embedded_option)  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 1, self.pp.CFD_DEM.coupling_level_type > 1)

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 1, self.pp.CFD_DEM.coupling_level_type or self.pp.CFD_DEM.print_PRESSURE_GRADIENT_option)

    def GetStationarityCounter(self):
        return SDP.Counter(self.pp.CFD_DEM.time_steps_per_stationarity_step, 1, self.pp.CFD_DEM.stationary_problem_option)

    def GetPrintCounter(self):
        return SDP.Counter(1, 1, 10) # still unused

    def GetDebugInfo(self):
        return SDP.Counter(self.pp.CFD_DEM.debug_tool_cycle, 1, self.pp.CFD_DEM.print_debug_info_option)

    def GetParticlesResultsCounter(self):
        return SDP.Counter(self.pp.CFD_DEM.print_particles_results_cycle, 1, self.pp.CFD_DEM.print_particles_results_option)

    def HistoryForceQuadratureCounter(self):
        return SDP.Counter(self.pp.CFD_DEM.time_steps_per_quadrature_step, 1, self.pp.CFD_DEM.basset_force_type)

    def GetRunCode(self):
        return SDP.CreateRunCode(self.pp)

    def ActivateTurbulenceModel(self):
        fluid_model_part = self.all_model_parts.Get('FluidPart')

        if self.pp.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras":
            # apply the initial turbulent viscosity on all of the nodes
            turb_visc = self.pp.FluidSolverConfiguration.TurbulentViscosity
            for node in fluid_model_part.Nodes:
                node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, turb_visc)
                visc = node.GetSolutionStepValue(VISCOSITY)
                node.SetSolutionStepValue(MOLECULAR_VISCOSITY, 0, visc)
                if node.IsFixed(VELOCITY_X) and node.GetSolutionStepValue(VELOCITY_X, 0) != 0.0 or \
                   node.IsFixed(VELOCITY_Y) and node.GetSolutionStepValue(VELOCITY_Y, 0) != 0.0 or \
                   node.IsFixed(VELOCITY_Z) and node.GetSolutionStepValue(VELOCITY_Z, 0) != 0.0:
                    node.Fix(TURBULENT_VISCOSITY)

            # select nodes on the wall
            self.fluid_solver.wall_nodes = []
            for i in SolverSettings.SA_wall_group_ids:
                # get the nodes of the wall for SA.
                nodes = fluid_model_part.GetNodes(i)
                for node in nodes:
                    self.fluid_solver.wall_nodes.append(node)
                    node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, 0.0)
                    node.Fix(TURBULENT_VISCOSITY)

    def GetFieldUtility(self):
        field_utility = None

        if self.pp.CFD_DEM.ElementType == "SwimmingNanoParticle":
            flow_field = ConstantVelocityField(0.00001, 0, 0)
            space_time_set = SpaceTimeSet()
            field_utility = FluidFieldUtility(space_time_set, flow_field, 1000.0, 1e-6)

        return field_utility

    def GetResultsCreator(self):
        return None

    def ApplyForwardCoupling(self, alpha = 'None'):
        self.projection_module.ApplyForwardCoupling(alpha)

    def ApplyForwardCouplingOfVelocityOnly(self):
        self.projection_module.ApplyForwardCouplingOfVelocityOnly()

    def PerformFinalOperations(self, time = None):
        pass
