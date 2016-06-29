from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import utilities
#~ import connectivity_mapper                  # Auxiliary matching meshes communicator
import residual_definitions                 # Residual definitions
#~ import mvqn_strategy                        # MultiVector Quasi-Newton method strategy
#~ import jfnk_strategy                        # Jacobian Free Newton-Krylov method strategy
#~ import relaxation_strategy                  # Relaxation strategy

# Import libraries
import time as timemodule                   # Import time library as timemodule (avoid interferences with "time" var)
import json                                 # Encoding library (for data exchange)

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(structure_main_model_part, fluid_main_model_part, project_parameters):
    return PartitionedFSISolver(structure_main_model_part, fluid_main_model_part, project_parameters)
    

class PartitionedFSISolver:
    def __init__(self, structure_main_model_part, fluid_main_model_part, project_parameters): 
        
        # Initial tests
        if project_parameters["structure_solver_settings"]["problem_data"]["time_step"].GetDouble() != project_parameters["fluid_solver_settings"]["problem_data"]["time_step"].GetDouble():
            raise("ERROR: Different time step among subdomains!")
        if project_parameters["structure_solver_settings"]["problem_data"]["start_time"].GetDouble() != project_parameters["fluid_solver_settings"]["problem_data"]["start_step"].GetDouble():
            raise("ERROR: Different number of time steps among subdomains!")
        if project_parameters["structure_solver_settings"]["problem_data"]["end_time"].GetDouble() != project_parameters["fluid_solver_settings"]["problem_data"]["end_time"].GetDouble():
            raise("ERROR: Different final time among subdomains!")
        
        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.structure_main_model_part = structure_main_model_part    
        self.fluid_main_model_part = fluid_main_model_part    
                
        # Settings string in JSON format
        default_settings = KratosMultiphysics.Parameters("""
        {
        "structure_solver_settings":
            {
            "solver_type": "solid_mechanics_implicit_dynamic_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "echo_level": 0,
            "time_integration_method": "Implicit",
            "analysis_type": "nonlinear",
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_iteration": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "solution_type": "Dynamic",
            "scheme_type": "Newmark",
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-4,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "Super LU",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "processes_sub_model_part_list": [""],
            "problem_domain_sub_model_part_list": ["solid_model_part"]
            },
        "fluid_solver_settings":
            {
            "solver_type": "navier_stokes_solver_vmsmonolithic",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "maximum_iterations": 10,
            "dynamic_tau": 0.0,
            "oss_switch": 0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": false,
            "divergence_clearance_steps": 0,
            "reform_dofs_at_each_iteration": true,
            "relative_velocity_tolerance": 1e-5,
            "absolute_velocity_tolerance": 1e-7,
            "relative_pressure_tolerance": 1e-5,
            "absolute_pressure_tolerance": 1e-7,
            "linear_solver_settings"        : {
                "solver_type" : "AMGCL_NS_Solver",
                "krylov_type" : "bicgstab",
                "velocity_block_preconditioner" : {
                    "tolerance" : 1e-3,
                    "precondioner_type" : "spai0"
                },
                "pressure_block_preconditioner" : {
                    "tolerance" : 1e-2,
                    "precondioner_type" : "spai0"
                },
                "tolerance" : 1e-6,
                "krylov_type": "bicgstab",
                "gmres_krylov_space_dimension": 50,
                "max_iteration": 50,
                "verbosity" : 0,
                "scaling": true,
                "coarse_enough" : 5000
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "alpha":-0.3,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "regularization_coef": 1000,
            "MoveMeshFlag": false,
            "use_slip_conditions": false,
            "turbulence_model": "None",
            "use_spalart_allmaras": false
            },
        "coupling_solver_settings":
            {
            "coupling_scheme"   : "DirichletNeumann",
            "solver_type"       : "partitioned_fsi_solver",
            "nl_tol"            : 1e-5,
            "nl_max_it"         : 50,
            "coupling_strategy" : {
                "solver_type"       : "relaxation_strategy",
                "acceleration_type" : "Aitken",
                "w_0"               : 0.825
                },
            "structure_interfaces_list" : [""],
            "fluid_interfaces_list" : [""]
            }
        }
        """)
        
        # Take the each one of the solvers settings from the ProjectParameters
        self.settings = KratosMultiphysics.Parameters("{}")
        self.settings.AddValue("structure_solver_settings",project_parameters["structure_solver_settings"]["solver_settings"])
        self.settings.AddValue("fluid_solver_settings",project_parameters["fluid_solver_settings"]["solver_settings"])
        self.settings.AddValue("coupling_solver_settings",project_parameters["coupling_solver_settings"]["solver_settings"])
        
        # Overwrite the default settings with user-provided parameters
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        print("*** Partitioned FSI solver construction starts...")
        # Construct the structure solver
        structure_solver_module = __import__(self.settings["structure_solver_settings"]["solver_type"].GetString())
        self.structure_solver = structure_solver_module.CreateSolver(self.structure_main_model_part, 
                                                                     self.settings["structure_solver_settings"])
        print("* Structure solver constructed.")
        
        # Construct the fluid solver
        fluid_solver_module = __import__(self.settings["fluid_solver_settings"]["solver_type"].GetString())
        self.fluid_solver = fluid_solver_module.CreateSolver(self.fluid_main_model_part, 
                                                             self.settings["fluid_solver_settings"])
        print("* Fluid solver constructed.")
        
        # Construct the coupling partitioned strategy
        self.max_nl_it = self.settings["coupling_solver_settings"]["nl_max_it"].GetInt()
        self.nl_tol = self.settings["coupling_solver_settings"]["nl_tol"].GetDouble()
        
        coupling_strategy_name = self.settings["coupling_solver_settings"]["coupling_strategy"]["solver_type"].GetString()
        interface_strategy_module = __import__(coupling_strategy_name)       
        self.coupling_strategy = interface_strategy_module.CreateStrategy(self.settings["coupling_solver_settings"]["coupling_strategy"])
        print("* Coupling strategy constructed.")
        print("*** FSI partitioned solver construction finished.")
        
        
    def GetMinimumBufferSize(self):
        # Get structure buffer size
        self.structure_solver.GetMinimumBufferSize()
        # Get fluid buffer size
        self.fluid_solver.GetMinimumBufferSize()


    def AddVariables(self):
        ## Structure variables addition
        # Standard CFD variables addition
        self.structure_solver.AddVariables()
        # FSI variables addition
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_INTERFACE)    # TODO: IS_INTERFACE is deprecated. Move to INTERFACE.
        # Mapper variables addition
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MAUX)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.AUX)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.VAUX)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.SCALAR_PROJECTED)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.VECTOR_PROJECTED)
        
        ## Fluid variables addition
        # Standard CFD variables addition
        self.fluid_solver.AddVariables()
        # FSI variables addition
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_INTERFACE)    # TODO: IS_INTERFACE is deprecated. Move to INTERFACE.
        # Mapper variables addition
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MAUX) 
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.AUX)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.VAUX)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.SCALAR_PROJECTED)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.VECTOR_PROJECTED)
        
                
    def ImportModelPart(self):
        # Import structure model part
        self.structure_solver.ImportModelPart()
        
        # Import fluid model part
        self.fluid_solver.ImportModelPart()
        
        
    def AddDofs(self):
        # Add DOFs structure
        self.structure_solver.AddDofs()     
           
        # Add DOFs fluid
        self.fluid_solver.AddDofs()
        
        
    def Initialize(self):
        # Initialize structure solver
        self.structure_solver.Initialize()
        
        # Initialize fluid solver
        self.fluid_solver.Initialize()

        # Get interface problem sizes
        interface_problem_sizes = self._GetInterfaceProblemSizes()
        self.solid_interface_problem_size = interface_problem_sizes[0]
        self.fluid_interface_problem_size = interface_problem_sizes[1]
        
        # Get the domain size
        self.domain_size = self._GetDomainSize()
                
        # Construct the interface mapper
        ### NOTE: The old mapper in FSI app is used --> Move to mortar mapper asap.
        ### NOTE2: The flag "INTERFACE" ("IS_INTERFACE" if the old mapper is used) has to be set to the interface nodes before the mapper construction.
        import NonConformant_OneSideMap
        self.interface_mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(self.fluid_solver.main_model_part,
                                                                                  self.structure_solver.main_model_part, 
                                                                                  1.0, 50, 1e-5)
        print("Non-conformant one side map constructed.")
                                                                                  
                                                                             
        # Construct the interface residual 
        self.interface_residual = residual_definitions.CreateResidual(self.settings["coupling_solver_settings"]["coupling_scheme"].GetString(),
                                                                      self.fluid_solver,  
                                                                      self.structure_solver, 
                                                                      self.interface_mapper) 
        print("Interface residual constructed.")
        
        # Interface solution initial guess (it might be velocity or fluxes depending on the type of coupling)
        # Note that the FSI problem is defined in terms of the fluid interface
        fluid_interface_residual_size = self.fluid_interface_problem_size*self.domain_size
        self.iteration_guess_value = KratosMultiphysics.Vector(fluid_interface_residual_size) 

        for i in range(0,fluid_interface_residual_size):
            self.iteration_guess_value[i] = 0.0001
                                                    
                
    def GetComputeModelPart(self):
        ##### THIS FUNCTION IS PROBABLY NOT NECESSARY IN THE COMBINED SOLVER
        pass
        
        
    def GetOutputVariables(self):
        pass
        
        
    def ComputeDeltaTime(self):
        pass
        
        
    def SaveRestart(self):
        pass #one should write the restart file here
        
        
    def Solve(self):
        
        ##### DEVELOPMENT FLAGS --> TO BE REMOVED AS SOON AS THEY ARE IMPLEMENTED IN THE JSON FILE.
        mesh_prediction = False
        fluid_prediction = False
        move_interface =  False
        
        # Compute mesh prediction # --> IMPLEMENT CORRECTLY THE MESH PREDICTION
        #~ if mesh_prediction == True:
            #~ print("Computing time step ",step," prediction...")
            #~ # Get the previous step fluid interface nodal fluxes
            #~ fluid_flux_prev_step = residual_definitions.GetInterfaceNodalResult(self.fluid_solver,"REACTION",1)
            #~ fluid_flux_prev_step_comm = wet_interface_comm.FluidToStructure_VectorMap(fluid_flux_prev_step,True)
                
            #~ # Solve the current step structure problem with the previous step fluid interface nodal fluxes
            #~ D2_problem.SolveNeumann(fluid_flux_prev_step_comm)
            #~ solid_interface_disp = residual_definitions.GetInterfaceNodalResult(D2_problem,"DISPLACEMENT",0)
            #~ solid_interface_disp_comm = wet_interface_comm.StructureToFluid_VectorMap(solid_interface_disp,False)
                
            #~ # Solve the fluid mesh problem with the obtained solid interface displacements
            #~ D1_problem.SolveMesh(solid_interface_disp_comm)
        
            #~ if fluid_prediction == True:
                #~ # Get the obtained solid interface velocity
                #~ solid_interface_velo = residual_definitions.GetInterfaceNodalResult(D2_problem,"VELOCITY",1)
                #~ solid_interface_velo_comm = wet_interface_comm.StructureToFluid_VectorMap(solid_interface_velo,False)
            
                #~ # Solve the fluid domain with the obtained solid interface velocity and the updated mesh
                #~ D1_problem.SolveDirichlet(solid_interface_velo_comm)
            
            #~ print("Time step ",step," prediction computed.")
        
        
        for nl_it in range(1,self.max_nl_it+1):
            
            print("    NL-ITERATION ",nl_it,"STARTS.")
            
            print("     Residual computation starts...")
            vel_residual = self.interface_residual.ComputeResidual(self.iteration_guess_value)            
            nl_res_norm = space.TwoNorm(vel_residual)
            print("     Residual computation finished. |res|=",nl_res_norm)      
            
            # Move interface nodes
            if move_interface == True:
                
                    solid_interface_disp = self.interface_residual.GetInterfaceNodalResult(self.structure_solver,"DISPLACEMENT",0)
                    #~ solid_interface_disp = D2_problem.GetInterfaceDisplacement()
                    solid_interface_disp_comm = wet_interface_comm.StructureToFluid_VectorMap(solid_interface_disp)
                    D1_problem.MoveInterface(solid_interface_disp_comm)
                    
            # Check convergence
            
            
                    
            
        
        
        ### ESTO DEBE HACERSE DENTRO DEL RESIDUAL
        # Initialize structure solver
        self.structure_solver.Solve()
        # Initialize fluid solver
        self.fluid_solver.Solve()
        # Initialize coupling solver
        
        # This solve must contain the current step resolution, that is to say the non-linear loop.
        
        #~ self.mechanical_solver.Solve()


    def SetEchoLevel(self, level):
        pass
        #~ self.solver.SetEchoLevel(level)


    def Clear(self):
        pass
        #~ self.solver.Clear()
        
        
    def Check(self):
        pass
        # How to do a check of the interface solver?¿?¿ Probably it is not necessary...
        #~ self.mechanical_solver.Check()
            
    def _GetInterfaceProblemSizes(self):
        
        # Get solid interface problem size
        structure_interface_pb_size = 0
        for i in range(self.settings["coupling_solver_settings"]["structure_interfaces_list"].size()):
            interface_submodelpart_name = self.settings["coupling_solver_settings"]["structure_interfaces_list"][i].GetString()
            structure_interface_pb_size += len(self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name).Nodes)
        
        # Get fluid interface problem size
        fluid_interface_pb_size = 0
        for i in range(self.settings["coupling_solver_settings"]["fluid_interfaces_list"].size()):
            interface_submodelpart_name = self.settings["coupling_solver_settings"]["fluid_interfaces_list"][i].GetString()
            fluid_interface_pb_size += len(self.fluid_solver.main_model_part.GetSubModelPart(interface_submodelpart_name).Nodes)
            
        return (structure_interface_pb_size, fluid_interface_pb_size)
        
    def _GetDomainSize(self):
        
        fluid_domain_size = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        structure_domain_size = self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        
        if fluid_domain_size !=structure_domain_size:
            raise("ERROR: Solid domain size and fluid domain size are not equal!")
        
        return fluid_domain_size
            
            
        
