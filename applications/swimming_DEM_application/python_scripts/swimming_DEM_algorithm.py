from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#TODO: set self.SolverSettings only once! 
#TODO: the same for self.fluid_model_part
#TODO: test DEM bounding box

import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import dem_main_script_ready_for_coupling as DEM_algorithm
import swimming_DEM_procedures as SDP
import variables_management as vars_man
import time as timer


sys.path.insert(0,'')
import DEM_explicit_solver_var as DEM_parameters
BaseAlgorithm = DEM_algorithm.Solution

class Algorithm(BaseAlgorithm):
    def __init__(self, pp):
        self.StartTimer()
        self.pp = pp
        self.SetBetaParamters()
        self.SetDoSolveDEMVariable()
        # creating a basset_force tool to perform the operations associated with the calculation of this force along the path of each particle
        self.GetBassetForceTools()
        BaseAlgorithm.__init__(self)
        self.CreateParts()

    def CreateParts(self):
        # Order must be respected here
        # defining a fluid model
        self.all_model_parts.Add(ModelPart("FluidPart"))
        # defining a model part for the mixed part
        self.all_model_parts.Add(ModelPart("MixedPart"))

    def StartTimer(self):
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
        self.pp.CFD_DEM.fluid_fraction_grad_type = 0
        self.pp.CFD_DEM.store_full_gradient = 0
        self.pp.CFD_DEM.laplacian_calculation_type = 0
        self.pp.CFD_DEM.do_search_neighbours = True
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
        self.pp.CFD_DEM.time_window = 0.8
        self.pp.CFD_DEM.number_of_exponentials = 10
        self.pp.CFD_DEM.number_of_quadrature_steps_in_window = int(self.pp.CFD_DEM.time_window / self.pp.CFD_DEM.delta_time_quadrature)
        self.pp.CFD_DEM.print_steps_per_plot_step = 1
        self.pp.CFD_DEM.PostCationConcentration = False
        self.pp.CFD_DEM.do_impose_flow_from_field = False
        self.pp.CFD_DEM.print_MATERIAL_ACCELERATION_option = True
        self.pp.CFD_DEM.print_FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED_option = False
        self.pp.CFD_DEM.print_VELOCITY_GRADIENT_option = 1
        self.pp.CFD_DEM.print_FLUID_FRACTION_GRADIENT_option = 0
        self.pp.CFD_DEM.print_FLUID_FRACTION_GRADIENT_PROJECTED_option = 0
        self.pp.CFD_DEM.print_VORTICITY_option = 1
        self.pp.CFD_DEM.print_MATERIAL_ACCELERATION_option = True
        self.pp.CFD_DEM.calculate_diffusivity_option = False
        self.pp.CFD_DEM.print_CONDUCTIVITY_option = False
        self.pp.CFD_DEM.filter_velocity_option = False
        self.pp.CFD_DEM.print_PARTICLE_VEL_option = False
        # Making the fluid step an exact multiple of the DEM step
        self.pp.Dt = int(self.pp.Dt / self.pp.CFD_DEM.MaxTimeStep) * self.pp.CFD_DEM.MaxTimeStep
        self.pp.viscosity_modification_type = 0.0
        self.domain_size = 3
        self.pp.type_of_inlet = 'VelocityImposed' # 'VelocityImposed' or 'ForceImposed'
        self.pp.force = Vector(3)
        self.pp.force[0] = 0
        self.pp.force[1] = 0
        self.pp.force[2] = 1e-10

        # defining and adding imposed porosity fields        
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
        
        self.spheres_model_part = self.all_model_parts.Get('SpheresPart')
        self.fluid_model_part = self.all_model_parts.Get('FluidPart')

        self.SetFluidSolverModule()
        
        self.vars_man = vars_man
        self.vars_man.ConstructListsOfVariables(self.pp)
                                        
        self.FluidInitialize()
        self.DispersePhaseInitialize()
        
        self.vars_man.AddingExtraProcessInfoVariables(self.pp, self.fluid_model_part, self.spheres_model_part)

    def FluidInitialize(self):

        self.AddFluidVariables()   
        self.ReadFluidModelPart()
        self.SetFluidBufferSizeAndAddDofs()        
        self.SetFluidSolver()        
        self.fluid_solver.Initialize()        
        self.ActivateTurbulenceModel()
        
    def DispersePhaseInitialize(self):
        self.vars_man.AddNodalVariables(self.spheres_model_part, self.pp.dem_vars)
        self.vars_man.AddNodalVariables(self.rigid_face_model_part, self.pp.rigid_faces_vars)
        self.vars_man.AddNodalVariables(self.DEM_inlet_model_part, self.pp.inlet_vars)
        BaseAlgorithm.Initialize(self)                
        
    def SetFluidSolverModule(self):
        self.SolverSettings = self.pp.FluidSolverConfiguration
        self.solver_module = import_solver(self.SolverSettings)
    
    def SetFluidSolver(self):
        self.fluid_solver = self.solver_module.CreateSolver(self.fluid_model_part, self.SolverSettings)
        
    def AddFluidVariables(self):
        # caution with breaking up this block (memory allocation)! {
        self.AddFluidVariablesByFluidSolver()
        self.AddFluidVariablesBySwimmingDEMAlgorithm()
                
    def AddFluidVariablesByFluidSolver(self):
        
        if "REACTION" in self.pp.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(REACTION)
        if "DISTANCE" in self.pp.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)
            
        self.solver_module.AddVariables(self.fluid_model_part, self.SolverSettings)
        
    def AddFluidVariablesBySwimmingDEMAlgorithm(self):
        self.vars_man.AddNodalVariables(self.fluid_model_part, self.pp.fluid_vars) 

    def SetFluidBufferSizeAndAddDofs(self):
        self.SetFluidBufferSizeAndAddDofsByFluidSolver()
        self.SetFluidBufferSizeAndAddDofsBySwimmingDEMAlgorithm()  

    def SetFluidBufferSizeAndAddDofsByFluidSolver(self):
        spheres_model_part = self.all_model_parts.Get('SpheresPart')        
        self.SolverSettings = self.pp.FluidSolverConfiguration
        self.fluid_model_part.SetBufferSize(3)
        self.solver_module.AddDofs(self.fluid_model_part, self.SolverSettings)
        
    def SetFluidBufferSizeAndAddDofsBySwimmingDEMAlgorithm(self):
        SDP.AddExtraDofs(self.pp, self.fluid_model_part, self.spheres_model_part, self.cluster_model_part, self.DEM_inlet_model_part)
        
    def DEMSolve(self, time = 'None'): # time is passed in case it is needed
        self.solver.Solve()

    def FluidSolve(self, time = 'None'):
        pass
        #self.fluid_solver.Solve()

    def PerformZeroStepInitializations(self):
        pass

    def PerformInitialDEMStepOperations(self, time = None):
        pass

    def SetInlet(self):
        if DEM_parameters.dem_inlet_option:
            #Constructing the inlet and initializing it (must be done AFTER the self.spheres_model_part Initialize)
            # Note that right now only inlets of a single type are possible. This should be generalized.
            if self.pp.type_of_inlet == 'VelocityImposed':
                self.DEM_inlet = DEM_Inlet(self.DEM_inlet_model_part)
            elif self.pp.type_of_inlet == 'ForceImposed':
                self.DEM_inlet = DEM_Force_Based_Inlet(self.DEM_inlet_model_part, self.pp.force)

            self.DEM_inlet.InitializeDEM_Inlet(self.spheres_model_part, self.creator_destructor)

    def SetAnalyticFaceWatcher(self):
        from analytic_tools import analytic_data_procedures
        self.watcher = AnalyticFaceWatcher()
        self.watcher_analyser = analytic_data_procedures.WatcherAnalyzer(analytic_face_watcher = self.watcher, path = self.main_path)

    def SetInletWatcher(self):
        self.watcher_analyser.SetInlet(self.DEM_inlet)

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

    def GetRecoveryCounter(self):
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

    def FillHistoryForcePrecalculatedVectors(self):
        # Warning: this estimation is based on a constant time step for DEM. This is usually the case, but could not be so. A more robust implementation is needed!
        N_steps = int(self.pp.CFD_DEM.FinalTime / self.pp.CFD_DEM.MaxTimeStep) + 20
        spheres_model_part = self.all_model_parts.Get('SpheresPart')
        if self.pp.CFD_DEM.basset_force_type > 0:
            self.basset_force_tool.FillDaitcheVectors(N_steps, self.pp.CFD_DEM.quadrature_order, self.pp.CFD_DEM.time_steps_per_quadrature_step)
        if self.pp.CFD_DEM.basset_force_type >= 3 or self.pp.CFD_DEM.basset_force_type == 1:
            self.basset_force_tool.FillHinsbergVectors(spheres_model_part, self.pp.CFD_DEM.number_of_exponentials, self.pp.CFD_DEM.number_of_quadrature_steps_in_window)

    def AppendValuesForTheHistoryForce(self):
        spheres_model_part = self.all_model_parts.Get('SpheresPart')

        if self.pp.CFD_DEM.basset_force_type == 1 or self.pp.CFD_DEM.basset_force_type >= 3:
            self.basset_force_tool.AppendIntegrandsWindow(spheres_model_part)
        elif self.pp.CFD_DEM.basset_force_type == 2:
            self.basset_force_tool.AppendIntegrands(spheres_model_part)

    def GetBassetForceTools(self):
        self.basset_force_tool = SDP.BassetForceTools()

    def ActivateTurbulenceModel(self):
        self.fluid_model_part = self.all_model_parts.Get('FluidPart')

        if self.pp.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras":
            # apply the initial turbulent viscosity on all of the nodes
            turb_visc = self.pp.FluidSolverConfiguration.TurbulentViscosity
            for node in self.fluid_model_part.Nodes:
                node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, turb_visc)
                visc = node.GetSolutionStepValue(VISCOSITY)
                node.SetSolutionStepValue(MOLECULAR_VISCOSITY, 0, visc)
                if node.IsFixed(VELOCITY_X) and node.GetSolutionStepValue(VELOCITY_X, 0) != 0.0 or \
                   node.IsFixed(VELOCITY_Y) and node.GetSolutionStepValue(VELOCITY_Y, 0) != 0.0 or \
                   node.IsFixed(VELOCITY_Z) and node.GetSolutionStepValue(VELOCITY_Z, 0) != 0.0:
                    node.Fix(TURBULENT_VISCOSITY)

            # select nodes on the wall
            self.fluid_solver.wall_nodes = []
            for i in self.SolverSettings.SA_wall_group_ids:
                # get the nodes of the wall for SA.
                nodes = self.fluid_model_part.GetNodes(i)
                for node in nodes:
                    self.fluid_solver.wall_nodes.append(node)
                    node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, 0.0)
                    node.Fix(TURBULENT_VISCOSITY)

    def GetFieldUtility(self):
        return None

    def GetResultsCreator(self):
        return None

    def ApplyForwardCoupling(self, alpha = 'None'):
        self.projection_module.ApplyForwardCoupling(alpha)

    def ApplyForwardCouplingOfVelocityOnly(self, time = None):
        self.projection_module.ApplyForwardCouplingOfVelocityOnly()

    def PerformFinalOperations(self, time = None):
        pass
