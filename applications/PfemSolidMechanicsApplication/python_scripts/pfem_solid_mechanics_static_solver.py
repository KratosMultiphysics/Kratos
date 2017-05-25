from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_static_solver

def CreateSolver(main_model_part, custom_settings):
    return PfemStaticMechanicalSolver(main_model_part, custom_settings)

class PfemStaticMechanicalSolver(solid_mechanics_static_solver.StaticMechanicalSolver):
    #derived class from StaticMechanicalSolver in order to add KratosPfemSolid variables    
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "solid_mechanics_static_solver",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Static",
            "analysis_type": "Non-Linear",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "water_pressure_dofs": false,
            "jacobian_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "implex": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": false,
            "clear_storage": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-4,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "bodies_list": [],
            "problem_domain_sub_model_part_list": ["solid_model_part"],
            "processes_sub_model_part_list": [""]
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        print("Construction of (PFEM) Static Mechanical Solver finished")

    def AddVariables(self):
        
        # Add displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        # Add reactions for the displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        # Add nodal force variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
        # Add specific variables for the problem conditions
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosSolid.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosSolid.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosSolid.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

        if self.settings["rotation_dofs"].GetBool():
            # Add specific variables for the problem (rotation dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
        if self.settings["pressure_dofs"].GetBool():
            # Add specific variables for the problem (pressure dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            self.main_model_part.AddNodalSolutionStepVariable(KratosSolid.PRESSURE_REACTION)
        if self.settings.Has("water_pressure_dofs"):
            if self.settings["water_pressure_dofs"].GetBool():
                # Add specific variables for the problem (pressure dofs)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.WATER_PRESSURE)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        if self.settings.Has("jacobian_dofs"):
            if self.settings["jacobian_dofs"].GetBool():
               self.main_model_part.AddNodalSolutionStepVariable(KratosPfemSolid.JACOBIAN)
               self.main_model_part.AddNodalSolutionStepVariable(KratosPfemSolid.REACTION_JACOBIAN)
        print("::[(PFEM)Mechanical Solver]:: Variables ADDED")

    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X);
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y);
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z);
            
        if(self.settings["solution_type"].GetString() == "Dynamic"):
            for node in self.main_model_part.Nodes:
                # adding first derivatives as dofs
                node.AddDof(KratosMultiphysics.VELOCITY_X);
                node.AddDof(KratosMultiphysics.VELOCITY_Y);
                node.AddDof(KratosMultiphysics.VELOCITY_Z);
                # adding second derivatives as dofs
                node.AddDof(KratosMultiphysics.ACCELERATION_X);
                node.AddDof(KratosMultiphysics.ACCELERATION_Y);
                node.AddDof(KratosMultiphysics.ACCELERATION_Z);                
            
        if self.settings["rotation_dofs"].GetBool():
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X);
                node.AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y);
                node.AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z);
                
        if(self.settings["solution_type"].GetString() == "Dynamic" and self.settings["rotation_dofs"].GetBool()):
            for node in self.main_model_part.Nodes:       
                # adding first derivatives as dofs
                node.AddDof(KratosMultiphysics.ANGULAR_VELOCITY_X);
                node.AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Y);
                node.AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Z);
                # adding second derivatives as dofs
                node.AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_X);
                node.AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Y);
                node.AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Z);
                    
        if self.settings["pressure_dofs"].GetBool():                
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.PRESSURE, KratosSolid.PRESSURE_REACTION);
            if not self.settings["stabilization_factor"].IsNull():
                self.main_model_part.ProcessInfo[KratosMultiphysics.STABILIZATION_FACTOR] = self.settings["stabilization_factor"].GetDouble()

        if self.settings.Has("water_pressure_dofs"):
            if self.settings["water_pressure_dofs"].GetBool():
                for node in self.main_model_part.Nodes:       
                    node.AddDof(KratosMultiphysics.WATER_PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE);

        if self.settings.Has("jacobian_dofs"):
            if self.settings["jacobian_dofs"].GetBool():
                for node in self.main_model_part.Nodes:       
                    node.AddDof(KratosPfemSolid.JACOBIAN, KratosPfemSolid.REACTION_JACOBIAN);
        print("::[(PFEM)Mechanical Solver]:: DOF's ADDED")


 
