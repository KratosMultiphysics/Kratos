from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import utilities
#~ import connectivity_mapper                  # Auxiliary matching meshes communicator
#~ import residual_definitions                 # Residual definitions
#~ import mvqn_strategy                        # MultiVector Quasi-Newton method strategy
#~ import jfnk_strategy                        # Jacobian Free Newton-Krylov method strategy
#~ import relaxation_strategy                  # Relaxation strategy

# Import libraries
import numpy
import scipy
import scipy.sparse
import scipy.sparse.linalg
import time as timemodule                   # Import time library as timemodule (avoid interferences with "time" var)
import json                                 # Encoding library (for data exchange)

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

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
            "solver_type": "mvqn_strategy"
            }
        }
        """)
        
        # Take the each one of the solvers settings from the ProjectParameters
        self.settings = KratosMultiphysics.Parameters("{}")
        self.settings.AddValue("structure_solver_settings",project_parameters["structure_solver_settings"]["solver_settings"])
        self.settings.AddValue("fluid_solver_settings",project_parameters["fluid_solver_settings"]["solver_settings"])
        #~ custom_settings.AddValue("coupling_solver_settings",project_parameters["coupling_solver_settings"]["solver_settings"])           # TO BE IMPLEMENTED           
        
        # Overwrite the default settings with user-provided parameters
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        # Construct the structure solver
        structure_solver_module = __import__(self.settings["structure_solver_settings"]["solver_type"].GetString())
        self.structure_solver = structure_solver_module.CreateSolver(self.structure_main_model_part, self.settings["structure_solver_settings"])
        
        # Construct the fluid solver
        fluid_solver_module = __import__(self.settings["fluid_solver_settings"]["solver_type"].GetString())
        self.fluid_solver = fluid_solver_module.CreateSolver(self.fluid_main_model_part, self.settings["fluid_solver_settings"])
        
        # Construct the coupling strategy
        ### TO BE DONE IN A SIMILAR FASHION AS THE ONES ABOVE
               
        print("Construction of PartitionedFSISolver finished.")
        
        
    def GetMinimumBufferSize(self):
        # Get structure buffer size
        self.structure_solver.GetMinimumBufferSize()
        # Get fluid buffer size
        self.fluid_solver.GetMinimumBufferSize()


    def AddVariables(self):
        # Add structure variables
        self.structure_solver.AddVariables()
        # Add fluid variables
        self.fluid_solver.AddVariables()
        # Add coupling variables
        
                
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
        # Initialize coupling solver
        
        
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
        
        ## INITIAL TEST SOLVING THE FLUID AND THE STRUCTURE SEPARATELY
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
