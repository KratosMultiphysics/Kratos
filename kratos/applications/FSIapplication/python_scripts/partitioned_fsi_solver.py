from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import utilities
import residual_definitions                 # Residual definitions
import NonConformant_OneSideMap             # Import non-conformant mapper

# Import libraries
import time as timemodule                   # Import time library as timemodule (avoid interferences with "time" var)
import json                                 # Encoding library (for data exchange)

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as KratosALE
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
            "move_interface"    : true,
            "mesh_prediction"   : true,
            "coupling_strategy" : {
                "solver_type"       : "relaxation_strategy",
                "acceleration_type" : "Aitken",
                "w_0"               : 0.825
                },
            "mesh_solver"               : "mesh_solver_structural_similarity",
            "mesh_reform_dofs_each_step"     : false,
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
        
        # Auxiliar variables
        self.max_nl_it = self.settings["coupling_solver_settings"]["nl_max_it"].GetInt()
        self.nl_tol = self.settings["coupling_solver_settings"]["nl_tol"].GetDouble()
        self.move_interface = self.settings["coupling_solver_settings"]["move_interface"].GetBool()
        self.mesh_prediction = self.settings["coupling_solver_settings"]["mesh_prediction"].GetBool()
        
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
        coupling_strategy_name = self.settings["coupling_solver_settings"]["coupling_strategy"]["solver_type"].GetString()
        interface_strategy_module = __import__(coupling_strategy_name)       
        self.coupling_strategy = interface_strategy_module.CreateStrategy(self.settings["coupling_solver_settings"]["coupling_strategy"])
        print("* Coupling strategy constructed.")
        
        # Construct the ALE mesh solver
        mesh_solver_name = self.settings["coupling_solver_settings"]["mesh_solver"].GetString()
        self.mesh_solver_module = __import__(mesh_solver_name)
        self.mesh_solver = self.mesh_solver_module.CreateMeshSolver(self.fluid_main_model_part,
                                                                    self.settings["coupling_solver_settings"]["mesh_reform_dofs_each_step"].GetBool())
        print("* ALE mesh solver constructed.")
        
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
        
        ## Fluid variables addition
        # Standard CFD variables addition
        self.fluid_solver.AddVariables()
        # FSI variables addition
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_INTERFACE)    # TODO: IS_INTERFACE is deprecated. Move to INTERFACE.
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
        # Mesh solver variables addition
        self.mesh_solver_module.AddVariables(self.fluid_solver.main_model_part)
                
        ## Mapper variables addition
        NonConformant_OneSideMap.AddVariables(self.fluid_solver.main_model_part,self.structure_solver.main_model_part)
                
                
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
        self.mesh_solver_module.AddDofs(self.fluid_solver.main_model_part)
        
        
    def Initialize(self):
        # Initialize structure solver
        self.structure_solver.Initialize()
        
        # Initialize fluid solver
        self.fluid_solver.Initialize()
        
        # Mesh solver initialization
        self.mesh_solver.Initialize()

        # Get interface problem sizes
        interface_problem_sizes = self._GetInterfaceProblemSizes()
        self.solid_interface_problem_size = interface_problem_sizes[0]
        self.fluid_interface_problem_size = interface_problem_sizes[1]
        
        # Get the domain size
        self.domain_size = self._GetDomainSize()
                
        # Construct the interface mapper
        ### NOTE: The old mapper in FSI app is used --> Move to mortar mapper asap.
        ### NOTE2: The flag "INTERFACE" ("IS_INTERFACE" if the old mapper is used) has to be set to the interface nodes before the mapper construction.
        self.interface_mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(self.fluid_solver.main_model_part,
                                                                                  self.structure_solver.main_model_part, 
                                                                                  1.0, 50, 1e-5)
        print("Non-conformant one side map constructed.")
                                                                                  
        # Construct the interface residual 
        self.interface_residual = residual_definitions.CreateResidual(self.settings["coupling_solver_settings"],
                                                                      self.fluid_solver,  
                                                                      self.structure_solver, 
                                                                      self.interface_mapper) 
        print("Interface residual constructed.")
        
        # Set the Neumann B.C. in the structure interface
        self._SetStructureNeumannCondition()
        
        # Set the Neumann B.C. in the fluid interface
        if self.settings["coupling_solver_settings"]["coupling_scheme"].GetString() == "NeumannNeumann":
            self._SetFluidNeumannCondition()
            
        # Set ALE mesh boundary conditions
        self._SetFluidMeshBoundaryConditions()
        
        # Interface solution initial guess (it might be velocity or fluxes depending on the type of coupling)
        # Note that the FSI problem is defined in terms of the fluid interface
        fluid_interface_residual_size = self.fluid_interface_problem_size*self.domain_size
        self.iteration_guess_value = KratosMultiphysics.Vector(fluid_interface_residual_size) 

        for i in range(0,fluid_interface_residual_size):
            self.iteration_guess_value[i] = 0.0001

                
    def GetComputeModelPart(self):
        pass
        
        
    def GetOutputVariables(self):
        pass
        
        
    def ComputeDeltaTime(self):
        pass
        
        
    def SaveRestart(self):
        pass
        
        
    def Solve(self):
                
        ## Compute mesh prediction ##
        if self.mesh_prediction == True:
            print("Computing time step ",step," prediction...")
            # Get the previous step fluid interface nodal fluxes
            self.interface_mapper.FluidToStructure_VectorMap(KratosMultiphysics.REACTION,
                                                             KratosMultiphysics.POINT_LOAD,
                                                             False,False)                                              
                
            # Solve the current step structure problem with the previous step fluid interface nodal fluxes
            self.structure_solver.Solve()
            # Map the obtained structure displacement to the fluid interface
            self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                             KratosMultiphysics.DISPLACEMENT,
                                                             True,False)
            # Solve the mesh problem
            self.mesh_solver.Solve()
            self.mesh_solver.MoveNodes()
                                                             
            # Map the obtained structure velocity to the fluid interface
            self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.VELOCITY,
                                                             KratosMultiphysics.VELOCITY,
                                                             True,False)
        
            # Solve the fluid problem with the updated mesh and the solid interface velocity as Dirichlet B.C.
            self.fluid_solver.Solve()
            
        ## Non-Linear interface coupling iteration ##
        for nl_it in range(1,self.max_nl_it+1):
            
            self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] = nl_it
            self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] = nl_it
            
            print("    NL-ITERATION ",nl_it,"STARTS.")
            
            print("     Residual computation starts...")
            vel_residual = self.interface_residual.ComputeResidual(self.iteration_guess_value)            
            nl_res_norm = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.RESIDUAL_NORM]
            print("     Residual computation finished. |res|=",nl_res_norm)      
            
            # Move interface nodes
            if self.move_interface == True:
                    
                    self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                                     KratosMultiphysics.DISPLACEMENT,
                                                                     True,False)                      # Project the structure interface displacement onto the fluid interface
                    self.mesh_solver.MoveNodes()                                                      # Move interface nodes
                    
            # Check convergence
            if nl_res_norm < self.nl_tol:
            
                print("    NON-LINEAR ITERATION CONVERGENCE ACHIEVED")
                print("    Total non-linear iterations: ",nl_it," NL residual norm: ",nl_res_norm)
            
                break 
                
            else:
                
                self.iteration_guess_value = self.coupling_strategy.InterfaceSolutionUpdate(self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS],
                                                                                            nl_it,
                                                                                            self.iteration_guess_value,
                                                                                            vel_residual)
              
        ## Mesh update ##
        self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                         KratosMultiphysics.DISPLACEMENT,
                                                         True,False)                           # Project the structure interface displacement onto the fluid interface
        self.mesh_solver.Solve()                                                               # Solve the mesh problem
        self.mesh_solver.MoveNodes()                                                           # Move nodes
                
        # Compute mesh residual
        mesh_res_norm = self.interface_residual.ComputeMeshResidual()
        print("MESH RESIDUAL NORM: ",mesh_res_norm)
        

    def SetEchoLevel(self, level):
        pass
        #~ self.solver.SetEchoLevel(level)
        
        
    def SetTimeStep(self, step):
        
        self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS] = step
        self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS] = step


    def Clear(self):
        pass
        
        
    def Check(self):
        pass
            
            
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
        
        
    def _SetStructureNeumannCondition(self):
        
        structure_computational_submodelpart = self.structure_solver.main_model_part.GetSubModelPart("solid_computational_model_part")
        
        aux_count = len(self.structure_solver.main_model_part.Conditions)       # Get the last existing condition numbering
        aux_count += 1
        print("max aux_count",aux_count)
        aux_count = 0
        for cond in self.structure_solver.main_model_part.Conditions:
            if(cond.Id > aux_count):
                aux_count = cond.Id
        aux_count += 1
        print("max aux_count",aux_count)
        
        
        for i in range(self.settings["coupling_solver_settings"]["structure_interfaces_list"].size()):
            interface_submodelpart_name = self.settings["coupling_solver_settings"]["structure_interfaces_list"][i].GetString()
            interface_submodelpart_i = self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name)
            # NOTE: In this manner, two interface submodelparts cannot share a node (it would be repeated in the pointload conditions...)
        
            for node in interface_submodelpart_i.Nodes:
                
                # Create the point load condition
                if self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                    structure_computational_submodelpart.CreateNewCondition("PointLoadCondition2D1N",aux_count,[node.Id],self.structure_solver.main_model_part.Properties[0]) 
                elif self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
                    structure_computational_submodelpart.CreateNewCondition("PointLoadCondition3D1N",aux_count,[node.Id],self.structure_solver.main_model_part.Properties[0]) 
                
                aux_count+=1
        
        
    def _SetFluidNeumannCondition(self):
        
        fluid_computational_volume_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart("volume_model_part_name")
        
        aux_count = len(self.fluid_solver.main_model_part.Conditions)       # Get the last existing condition numbering
        aux_count += 1
        print("max aux_count",aux_count)
        aux_count = 0
        for cond in self.fluid_solver.main_model_part.Conditions:
            if(cond.Id > aux_count):
                aux_count = cond.Id
        aux_count += 1
        print("max aux_count",aux_count)
        
        
        for i in range(self.settings["coupling_solver_settings"]["fluid_interfaces_list"].size()):
            interface_submodelpart_name = self.settings["coupling_solver_settings"]["fluid_interfaces_list"][i].GetString()
            interface_submodelpart_i = self.fluid_solver.main_model_part.GetSubModelPart(interface_submodelpart_name)
            # NOTE: In this manner, two interface submodelparts cannot share a node (it would be repeated in the pointload conditions...)
        
            for node in interface_submodelpart_i.Nodes:
                
                # NOTE: THIS CONDITION IS LOCAL IN MY MACHINE, DECIDE WHAT TO DO.
                # Create the fluid load condition
                if self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                    fluid_computational_volume_submodelpart.CreateNewCondition("PointForce2Dfluid",aux_count,[node.Id],self.fluid_solver.main_model_part.Properties[0]) 
                elif self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
                    fluid_computational_volume_submodelpart.CreateNewCondition("PointForce3Dfluid",aux_count,[node.Id],self.fluid_solver.main_model_part.Properties[0]) 
                
                aux_count+=1
                
            
    def _SetFluidMeshBoundaryConditions(self):
        # Note that this function also initializes the interface displacement as 0 (the interface submodelpart is within the skin_parts group).
        
        for i in range(self.settings["fluid_solver_settings"]["skin_parts"].size()):
            skin_submodelpart_name = self.settings["fluid_solver_settings"]["skin_parts"][i].GetString()
            skin_submodelpart_i = self.fluid_solver.main_model_part.GetSubModelPart(skin_submodelpart_name)
            
            for node in skin_submodelpart_i.Nodes:
                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
                node.SetValue(KratosMultiphysics.DISPLACEMENT_X,0.0)
                node.SetValue(KratosMultiphysics.DISPLACEMENT_Y,0.0)
                node.SetValue(KratosMultiphysics.DISPLACEMENT_Z,0.0)
