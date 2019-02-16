from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Importing the base class
from python_solver import PythonSolver

def CreateSolver(model, custom_settings):
    return ParticleMPMSolver(model, custom_settings)

class ParticleMPMSolver(PythonSolver):

    ### Solver constructor
    def __init__(self, model, custom_settings):
        super(ParticleMPMSolver, self).__init__(model, custom_settings)

        # Add model part containers
        self._add_model_part_containers()

        # Default settings
        self.min_buffer_size = 2

        # There is only a single rank in OpenMP, we always print
        self.is_printing_rank = True

        # Default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "MPM_Material",
            "domain_size"     : -1,
            "time_stepping"   : {},
            "solver_type"                        : "StaticSolver",
            "echo_level"                         : 0,
            "analysis_type"                      : "linear",
            "scheme_type"                        : "",
            "grid_model_import_settings"              : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name_Grid"
            },
            "model_import_settings"              : {
                "input_type"        : "mdpa",
                "input_filename"    : "unknown_name_Body"
            },
            "material_import_settings"           : {
                "materials_filename" : ""
            },
            "pressure_dofs"                      : false,
            "compute_reactions"                  : false,
            "convergence_criterion"              : "Residual_criteria",
            "displacement_relative_tolerance"    : 1.0E-4,
            "displacement_absolute_tolerance"    : 1.0E-9,
            "residual_relative_tolerance"        : 1.0E-4,
            "residual_absolute_tolerance"        : 1.0E-9,
            "max_iteration"                      : 20,
            "axis_symmetric_flag"                : false,
            "move_mesh_flag"                     : false,
            "problem_domain_sub_model_part_list" : [],
            "processes_sub_model_part_list"      : [],
            "element_search_settings"            : {
                "max_number_of_results"          : 1000,
                "searching_tolerance"            : 1.0E-5
            },
            "linear_solver_settings"             : {
                "solver_type" : "amgcl",
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

        # Temporary warnings, to be removed
        if custom_settings.Has("geometry_element"):
            custom_settings.RemoveValue("geometry_element")
            warning = '\n::[ParticleMPMSolver]:: W-A-R-N-I-N-G: You have specified "geometry_element", '
            warning += 'which is deprecated and will be removed soon. \nPlease remove it from the "solver settings"!\n'
            self.print_warning_on_rank_zero("Geometry element", warning)
        if custom_settings.Has("particle_per_element"):
            custom_settings.RemoveValue("particle_per_element")
            warning = '\n::[ParticleMPMSolver]:: W-A-R-N-I-N-G: You have specified "particle_per_element", '
            warning += 'which is deprecated and will be removed soon. \nPlease remove it from the "solver settings"!\n'
            self.print_warning_on_rank_zero("Particle per element", warning)
        if custom_settings.Has("line_search"):
            custom_settings.RemoveValue("line_search")
            warning = '\n::[ParticleMPMSolver]:: W-A-R-N-I-N-G: You have specified "line_search", '
            warning += 'which is deprecated and will be removed soon. \nPlease remove it from the "solver settings"!\n'
            self.print_warning_on_rank_zero("Geometry element", warning)

        # Overwrite the default settings with user-provided parameters
        self.settings.ValidateAndAssignDefaults(default_settings)

        # Construct the linear solvers
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        linear_solver_type = self.settings["linear_solver_settings"]["solver_type"].GetString()
        if(linear_solver_type == "amgcl" or linear_solver_type == "AMGCL"):
            self.block_builder = True
        else:
            self.block_builder = False
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        self.print_on_rank_zero("::[ParticleMPMSolver]:: ", "Solver is constructed correctly.")


    ### Solver public functions
    def AddVariables(self):
        # Add variables to different model parts
        self._add_variables_to_model_part(self.grid_model_part)
        self._add_variables_to_model_part(self.initial_material_model_part)

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
        self._add_dofs_to_model_part(self.grid_model_part)
        self._add_dofs_to_model_part(self.initial_material_model_part)

        self.print_on_rank_zero("::[ParticleMPMSolver]:: ","DOFs are added.")

    def Initialize(self):
        #TODO: implement solver_settings and change the input of the constructor in MPM_strategy.h

        # Set definition of the convergence criteria
        self.convergence_criterion_type = self.settings["convergence_criterion"].GetString()
        self.rel_disp_tol               = self.settings["displacement_relative_tolerance"].GetDouble()
        self.abs_disp_tol               = self.settings["displacement_absolute_tolerance"].GetDouble()
        self.rel_res_tol                = self.settings["residual_relative_tolerance"].GetDouble()
        self.abs_res_tol                = self.settings["residual_absolute_tolerance"].GetDouble()
        self.max_iteration              = self.settings["max_iteration"].GetInt()

        # Set definition of the global solver type
        self.solver_type                    = self.settings["solver_type"].GetString()

        # Set definition of the solver parameters
        self.compute_reactions      = self.settings["compute_reactions"].GetBool()
        self.pressure_dofs          = self.settings["pressure_dofs"].GetBool()
        self.axis_symmetric_flag    = self.settings["axis_symmetric_flag"].GetBool()
        self.move_mesh_flag         = self.settings["move_mesh_flag"].GetBool()

        # Set definition of search element
        self.max_number_of_search_results = self.settings["element_search_settings"]["max_number_of_results"].GetInt()
        self.searching_tolerance          = self.settings["element_search_settings"]["searching_tolerance"].GetDouble()

        # Identify geometry type
        self._identify_geometry_type()

        # Set default solver_settings parameters
        if self.geometry_element == "Triangle":
            if (self.domain_size == 2):
                if (self.pressure_dofs):
                    self.new_element = KratosParticle.CreateUpdatedLagragianUP2D3N()
                else:
                    if (self.axis_symmetric_flag):
                        self.new_element = KratosParticle.CreateUpdatedLagragianAxis2D3N()
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
                    if (self.axis_symmetric_flag):
                        self.new_element = KratosParticle.CreateUpdatedLagragianAxis2D4N()
                    else:
                        self.new_element = KratosParticle.CreateUpdatedLagragian2D4N()
            else:
                if (self.pressure_dofs):
                    raise Exception("Element for mixed U-P formulation in 3D for Hexahedral Element is not yet implemented.")
                else:
                    self.new_element = KratosParticle.CreateUpdatedLagragian3D8N()

        # Initialize solver
        if(self.domain_size==2):
            self.solver = KratosParticle.MPM2D(self.grid_model_part, self.initial_material_model_part, self.material_model_part,
                                self.linear_solver, self.new_element, self.solver_type, self.max_iteration, self.compute_reactions,
                                self.block_builder, self.pressure_dofs, self.move_mesh_flag)
        else:
            self.solver = KratosParticle.MPM3D(self.grid_model_part, self.initial_material_model_part, self.material_model_part,
                                self.linear_solver, self.new_element, self.solver_type, self.max_iteration, self.compute_reactions,
                                self.block_builder, self.pressure_dofs, self.move_mesh_flag)

        # Set echo level
        self._set_echo_level()

        # Check if everything is assigned correctly
        self._check()

        self.print_on_rank_zero("::[ParticleMPMSolver]:: ","Solver is initialized correctly.")

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt
        self.grid_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.grid_model_part.CloneTimeStep(new_time)

        return new_time

    def ComputeDeltaTime(self):
        return self.settings["time_stepping"]["time_step"].GetDouble()

    def SearchElement(self):
        self.solver.SearchElement(self.max_number_of_search_results, self.searching_tolerance)

    def InitializeSolutionStep(self):
        self.SearchElement()
        self.solver.Initialize()
        self.solver.InitializeSolutionStep()

    def Predict(self):
        self.solver.Predict()

    def SolveSolutionStep(self):
        is_converged = self.solver.SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.solver.FinalizeSolutionStep()

        self.solver.Clear()

    ### Solver private functions
    def _add_model_part_containers(self):

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size not in [2,3]:
            err_msg  = "The input \"domain_size\" is wrong!"
            err_msg += "Available options are: \"2\" or \"3\""
            raise Exception(err_msg)

        ### In MPM three model parts are needed
        ## Material model part definition
        material_model_part_name = self.settings["model_part_name"].GetString()
        if not self.model.HasModelPart(material_model_part_name):
            self.material_model_part = self.model.CreateModelPart(material_model_part_name) # Equivalent to model_part3 in the old format
            self.material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        ## Initial material model part definition
        initial_material_model_part_name = "Initial_" + material_model_part_name
        if not self.model.HasModelPart(initial_material_model_part_name):
            self.initial_material_model_part = self.model.CreateModelPart(initial_material_model_part_name) #Equivalent to model_part2 in the old format
            self.initial_material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        ## Grid model part definition
        if not self.model.HasModelPart("Background_Grid"):
            self.grid_model_part = self.model.CreateModelPart("Background_Grid") #Equivalent to model_part1 in the old format
            self.grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    def _add_variables_to_model_part(self, model_part):
        # Add displacements
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # Add dynamic variables
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        # Add reactions for the displacements
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # Add specific variables for the problem conditions
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MOMENTUM)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_INERTIA)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)

        # Add variables for arbitrary slope with slip
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        # Add variables for specific cases
        if self.settings["pressure_dofs"].GetBool():
            # add specific variables for the problem (pressure dofs)
            model_part.AddNodalSolutionStepVariable(KratosParticle.PRESSURE_REACTION)
            model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MPRESSURE)

    def _model_part_reading(self):
        # reading the model part of the background grid
        if(self.settings["grid_model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.grid_model_part, self.settings["grid_model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")

        # reading the model part of the material point
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.initial_material_model_part, self.settings["model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")

    def _execute_check_and_prepare(self):
        # Specific active node and element check for particle MPM solver
        for node in self.grid_model_part.Nodes:
            if (node.Is(KratosMultiphysics.ACTIVE)):
                self.print_on_rank_zero("::[ParticleMPMSolver]:: ","WARNING: This grid node have been set active: ", node.Id)

        # Setting active initial elements
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, self.initial_material_model_part.Elements)

        # Specify domain size
        self.domain_size = self.material_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # Read material property
        materials_imported = self._import_constitutive_laws()
        if materials_imported:
            self.print_on_rank_zero("::[ParticleMPMSolver]:: ","Constitutive law was successfully imported.")
        else:
            self.print_warning_on_rank_zero("::[ParticleMPMSolver]:: ","Constitutive law was not imported.")

        # Clone property of model_part2 to model_part3
        self.material_model_part.Properties = self.initial_material_model_part.Properties

    def _import_constitutive_laws(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    def _identify_geometry_type(self):
        for mpm in self.grid_model_part.Elements:
            num_nodes = len(mpm.GetNodes())
            break

        if (self.domain_size == 2 and num_nodes == 3) or (self.domain_size == 3 and num_nodes == 4):
            self.geometry_element = "Triangle"
        elif (self.domain_size == 2 and num_nodes == 4) or (self.domain_size == 3 and num_nodes == 8):
            self.geometry_element = "Quadrilateral"

    def _add_dofs_to_model_part(self, model_part):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, model_part)

        if self.settings["pressure_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosParticle.PRESSURE_REACTION, model_part)

    def _set_buffer_size(self):
        current_buffer_size = self.grid_model_part.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.grid_model_part.SetBufferSize(self.min_buffer_size)
        else:
            self.grid_model_part.SetBufferSize(current_buffer_size)

        current_buffer_size = self.initial_material_model_part.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.initial_material_model_part.SetBufferSize(self.min_buffer_size)
        else:
            self.initial_material_model_part.SetBufferSize(current_buffer_size)

    def _set_echo_level(self):
        self.solver.SetEchoLevel(self.echo_level)

    def _check(self):
        self.solver.Check()

    def _is_printing_rank(self):
        return self.is_printing_rank
