from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()
KratosMultiphysics.CheckRegisteredApplications("ParticleMechanicsApplication")

# Import applications and dependencies
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Importing the base class
from python_solver import PythonSolver

def CreateSolver(model, custom_settings):

    domain_size = custom_settings["domain_size"].GetInt()

    ### In MPM three model parts are needed
    ## Material model part definition
    material_model_part_name = custom_settings["model_part_name"].GetString()
    material_model_part = KratosMultiphysics.ModelPart(material_model_part_name) # Equivalent to model_part3 in the old format
    material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    ## Initial material model part definition
    initial_material_model_part_name = "Initial_" + material_model_part_name
    initial_material_model_part = KratosMultiphysics.ModelPart(initial_material_model_part_name) #Equivalent to model_part2 in the old format
    initial_material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    ## Grid model part definition
    grid_model_part = KratosMultiphysics.ModelPart("Background_Grid") #Equivalent to model_part1 in the old format
    grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    ## Adding into one model
    model.AddModelPart(grid_model_part)
    model.AddModelPart(initial_material_model_part)
    model.AddModelPart(material_model_part)

    return ParticleMPMSolver(model, custom_settings)

class ParticleMPMSolver(PythonSolver):

    ### Solver constructor
    def __init__(self, model, custom_settings):
        super(ParticleMPMSolver, self).__init__(model, custom_settings)

        # Default settings
        material_model_part_name = custom_settings["model_part_name"].GetString()
        self.model_part1 = self.model.GetModelPart("Background_Grid")                        #grid_model_part
        self.model_part2 = self.model.GetModelPart("Initial_" + material_model_part_name)    #initial_model_part
        self.model_part3 = self.model.GetModelPart(material_model_part_name)                 #mpm_model_part
        self.min_buffer_size = 3

        # There is only a single rank in OpenMP, we always print
        self.is_printing_rank = True

        # Default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "MPM_Material",
            "domain_size"     : 2,
            "time_stepping" : {
                "time_step"       : 0.0005
            },
            "solver_type"                        : "StaticSolver",
            "echo_level"                         : 0,
            "time_integration_method"            : "Implicit",
            "analysis_type"                      : "linear",
            "scheme_type"                        : "Newmark",
            "grid_model_import_settings"              : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name_Grid"
            },
            "model_import_settings"              : {
                "input_type"        : "mdpa",
                "input_filename"    : "unknown_name_Body"
            },
            "material_import_settings"           : {
                "materials_filename" : "ParticleMaterials.json"
            },
            "explicit_integration_scheme"        : "CentralDifferences",
            "time_step_prediction_level"         : "Automatic",
            "rayleigh_damping"                   : false,
            "rotation_dofs"                      : false,
            "pressure_dofs"                      : false,
            "reform_dof_set_at_each_step"        : false,
            "line_search"                        : false,
            "implex"                             : false,
            "compute_reactions"                  : true,
            "compute_contact_forces"             : false,
            "convergence_criterion"              : "Residual_criteria",
            "displacement_relative_tolerance"    : 1.0E-4,
            "displacement_absolute_tolerance"    : 1.0E-9,
            "residual_relative_tolerance"        : 1.0E-4,
            "residual_absolute_tolerance"        : 1.0E-9,
            "max_iteration"                      : 10,
            "geometry_element"                   : "Triangle",
            "number_of_material"                 : 1,
            "particle_per_element"               : 3,
            "impenetrability_condition"          : true,
            "move_mesh_flag"                     : false,
            "problem_domain_sub_model_part_list" : [],
            "processes_sub_model_part_list"      : [],
            "linear_solver_settings": {
                "solver_type" : "AMGCL",
                "smoother_type":"damped_jacobi",
                "krylov_type": "cg",
                "coarsening_type": "aggregation",
                "max_iteration": 200,
                "provide_coordinates": false,
                "gmres_krylov_space_dimension": 100,
                "verbosity" : 0,
                "tolerance": 1e-7,
                "scaling": false,
                "block_size": 3,
                "use_block_matrices_if_possible" : true,
                "coarse_enough" : 50
            }
        }""")

        # Overwrite the default settings with user-provided parameters
        self.settings.ValidateAndAssignDefaults(default_settings)

        # Construct the linear solvers
        import linear_solver_factory
        if(self.settings["linear_solver_settings"]["solver_type"].GetString() == "AMGCL"):
            self.block_builder = True
        else:
            self.block_builder = False
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        self.print_on_rank_zero("::[ParticleMPMSolver]:: ", "Solver is constructed correctly.")


    ### Solver public functions
    def AddVariables(self):
        # Add variables to different model parts
        self._add_variables_to_model_part(self.model_part1)
        self._add_variables_to_model_part(self.model_part2)

        self.print_on_rank_zero("::[ParticleMPMSolver]:: ", "Variables are added.")
 
    def ImportModelPart(self):
        # Read model part
        self._model_part_reading()

        self.print_on_rank_zero("::[ParticleMPMSolver]:: ","Models are imported.")

    def PrepareModelPart(self):
        # Set buffer size
        self._set_buffer_size()

        # Executes the check and prepare model process
        self._execute_check_and_prepare()

        self.print_on_rank_zero("::[ParticleMPMSolver]:: ", "ModelPart prepared for Solver.")

    def GetComputingModelPart(self):
        if not self.model.HasModelPart(self.settings["model_part_name"].GetString()):
            raise Exception("The ComputingModelPart was not created yet!")
        return self.model.GetModelPart(self.settings["model_part_name"].GetString())

    def GetGridModelPart(self):
        if not self.model.HasModelPart("Background_Grid"):
            raise Exception("The GridModelPart was not created yet!")
        return self.model.GetModelPart("Background_Grid")

    def AddDofs(self):
        # Add dofs to different model parts
        self._add_dofs_to_model_part(self.model_part1)
        self._add_dofs_to_model_part(self.model_part2)
        
        self.print_on_rank_zero("::[ParticleMPMSolver]:: ","DOFs are added.")

    def Initialize(self):
        #TODO: implement solver_settings and change the input of the constructor in MPM_strategy.h
       
        # Set definition of the convergence criteria
        self.convergence_criterion_type = self.settings["convergence_criterion"].GetString()
        self.rel_disp_tol               = self.settings["displacement_relative_tolerance"].GetDouble()
        self.abs_disp_tol               = self.settings["displacement_absolute_tolerance"].GetDouble()
        self.rel_res_tol                = self.settings["residual_relative_tolerance"].GetDouble()
        self.abs_res_tol                = self.settings["residual_absolute_tolerance"].GetDouble()
        self.max_iters                  = self.settings["max_iteration"].GetInt()

        # Set definition of the global solver type
        self.solver_type                    = self.settings["solver_type"].GetString()
        self.time_integration_method        = self.settings["time_integration_method"].GetString()

        # Set definition of the solver parameters
        self.compute_reactions      = self.settings["compute_reactions"].GetBool()
        self.compute_contact_forces = self.settings["compute_contact_forces"].GetBool()
        self.rotation_dofs          = self.settings["rotation_dofs"].GetBool()
        self.pressure_dofs          = self.settings["pressure_dofs"].GetBool()
        self.line_search            = self.settings["line_search"].GetBool()
        self.implex                 = self.settings["implex"].GetBool()
        self.move_mesh_flag         = self.settings["move_mesh_flag"].GetBool()

        # Set default solver_settings parameters
        self.geometry_element   = self.settings["geometry_element"].GetString()
        self.number_particle    = self.settings["particle_per_element"].GetInt()
        if self.geometry_element == "Triangle":
            if (self.domain_size == 2):
                if (self.pressure_dofs):
                    self.new_element = KratosParticle.CreateUpdatedLagragianUP2D3N()
                else:
                    self.new_element = KratosParticle.CreateUpdatedLagragian2D3N()
            else:
                if (self.pressure_dofs):
                    raise Exception("Element for mixed U-P formulation in 3D for Tetrahedral Element is not yet implemented.")
                else:
                    self.new_element = KratosParticle.CreateUpdatedLagragian3D4N()
        elif self.geometry_element == "Quadrilateral":
            if (self.domain_size == 2):
                if (self.pressure_dofs):
                    raise Exception("Element for mixed U-P formulation in 2D for Quadrilateral Element is not yet implemented.")
                else:
                    self.new_element = KratosParticle.CreateUpdatedLagragian2D4N()
            else:
                if (self.pressure_dofs):
                    raise Exception("Element for mixed U-P formulation in 3D for Hexahedral Element is not yet implemented.")
                else:
                    self.new_element = KratosParticle.CreateUpdatedLagragian3D8N()

        # Initialize solver
        if(self.domain_size==2):
            self.solver = KratosParticle.MPM2D(self.model_part1, self.model_part2, self.model_part3, self.linear_solver, self.new_element, self.move_mesh_flag, self.solver_type, self.geometry_element, self.number_particle, self.block_builder, self.pressure_dofs)
        else:
            self.solver = KratosParticle.MPM3D(self.model_part1, self.model_part2, self.model_part3, self.linear_solver, self.new_element, self.move_mesh_flag, self.solver_type, self.geometry_element,  self.number_particle, self.block_builder, self.pressure_dofs)
      
        # Set echo level
        self._set_echo_level()
        
        # Check if everything is assigned correctly
        self._check()
      
        self.print_on_rank_zero("::[ParticleMPMSolver]:: ","Solver is initialized correctly.")

    def SolveSolutionStep(self):
        (self.solver).Solve()

    def Solve(self):
        warn_message  = 'Using "Solve" is deprecated and will be removed in the future!\n'
        warn_message += 'Use the "SolveSolutionStep" instead'
        self.print_warning_on_rank_zero("::[ParticleMPMSolver]:: ", warn_message)
        (self.solver).Solve()

    ### Solver private functions
    def _add_variables_to_model_part(self, model_part):
        # Add displacements
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        
        # Add dynamic variables
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        
        # Add reactions for the displacements
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        
        # Add nodal force variables
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FORCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCE)
        #model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
        
        # Add specific variables for the problem conditions
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MOMENTUM)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_INERTIA)
        model_part.AddNodalSolutionStepVariable(KratosParticle.DISPLACEMENT_AUX)
        model_part.AddNodalSolutionStepVariable(KratosParticle.AUX_VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosParticle.AUX_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(KratosParticle.AUX_R)
        model_part.AddNodalSolutionStepVariable(KratosParticle.AUX_T)
        model_part.AddNodalSolutionStepVariable(KratosParticle.AUX_R_VEL)
        model_part.AddNodalSolutionStepVariable(KratosParticle.AUX_T_VEL)
        model_part.AddNodalSolutionStepVariable(KratosParticle.AUX_R_ACC)
        model_part.AddNodalSolutionStepVariable(KratosParticle.AUX_T_ACC)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_LUMPED_MASS)
        
        # Add variables for arbitrary slope with slip
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        
        # Add variables for specific cases
        if self.settings["rotation_dofs"].GetBool():
            # add specific variables for the problem (rotation dofs)
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
        if self.settings["pressure_dofs"].GetBool():
            # add specific variables for the problem (pressure dofs)
            model_part.AddNodalSolutionStepVariable(KratosSolid.PRESSURE_REACTION)
            model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MPRESSURE)
            model_part.AddNodalSolutionStepVariable(KratosParticle.AUX_PRESSURE)
        if self.settings["time_integration_method"].GetString() == "explicit" :
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MIDDLE_VELOCITY)

    def _model_part_reading(self):
        # reading the model part of the background grid
        if(self.settings["grid_model_import_settings"]["input_type"].GetString() == "mdpa"):
            KratosMultiphysics.ModelPartIO(self.settings["grid_model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.model_part1)
        else:
            raise Exception("Other input options are not implemented yet.")
        
        # reading the model part of the material point
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.model_part2)
        else:
            raise Exception("Other input options are not implemented yet.")

    def _execute_check_and_prepare(self):
        # Specific active node and element check for particle MPM solver
        for node in self.model_part1.Nodes:
            if (node.Is(KratosMultiphysics.ACTIVE)):
                print(node.Id)
                
        for element in self.model_part2.Elements:
            element.Set(KratosMultiphysics.ACTIVE, True)

        # Specify domain size
        self.domain_size = self.model_part3.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

         # Read material property
        materials_imported = self._import_constitutive_laws()
        if materials_imported:
            self.print_on_rank_zero("::[ParticleMPMSolver]:: ","Constitutive law was successfully imported.")
        else:
            self.print_warning_on_rank_zero("::[ParticleMPMSolver]:: ","Constitutive law was not imported.")

        # Clone property of model_part2 to model_part3
        self.model_part3.Properties = self.model_part2.Properties

    def _import_constitutive_laws(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            import read_materials_process
            # Create a dictionary of model parts.
            Model = KratosMultiphysics.Model()
            Model.AddModelPart(self.model_part2)
            # Add constitutive laws and material properties from json file to model parts.
            read_materials_process.ReadMaterialsProcess(Model, self.settings["material_import_settings"])
                       
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    
    def _add_dofs_to_model_part(self, model_part):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, model_part)
        
        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X, model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y, model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z, model_part)
        
        if self.settings["pressure_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosSolid.PRESSURE_REACTION, model_part)

    def _set_buffer_size(self):
        current_buffer_size = self.model_part1.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.model_part1.SetBufferSize(self.min_buffer_size)
        else:
            self.model_part1.SetBufferSize(current_buffer_size)
        
        current_buffer_size = self.model_part2.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.model_part2.SetBufferSize(self.min_buffer_size)
        else:
            self.model_part2.SetBufferSize(current_buffer_size)

    def _set_echo_level(self):
        self.solver.SetEchoLevel(self.echo_level)

    def _check(self):
        self.solver.Check()

    def _is_printing_rank(self):
        return self.is_printing_rank