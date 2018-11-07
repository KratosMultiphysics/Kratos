from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from python_solver import PythonSolver

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("ShallowWaterApplication")

# Import applications
import KratosMultiphysics.ShallowWaterApplication as Shallow

def CreateSolver(model, custom_settings):
    return ShallowWaterBaseSolver(model, custom_settings)

class ShallowWaterBaseSolver(PythonSolver):
    def __init__(self, model, custom_settings):  # Constructor of the class
        settings = self._ValidateSettings(custom_settings)

        super(ShallowWaterBaseSolver, self).__init__(model, settings)

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        ## Set the element and condition names for the replace settings
        ## These should be defined in derived classes
        self.element_name = None
        self.condition_name = None
        self.min_buffer_size = 2

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        # Initialize shallow water variables utility
        self.ShallowVariableUtils = Shallow.ShallowWaterVariablesUtility(self.main_model_part, self.settings["dry_height"].GetDouble())

    def AddVariables(self):
        # Primitive variables
        self.main_model_part.AddNodalSolutionStepVariable(Shallow.HEIGHT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        # Physic problem parameters
        self.main_model_part.AddNodalSolutionStepVariable(Shallow.FREE_SURFACE_ELEVATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.GRAVITY)
        self.main_model_part.AddNodalSolutionStepVariable(Shallow.BATHYMETRY)
        self.main_model_part.AddNodalSolutionStepVariable(Shallow.MANNING)
        self.main_model_part.AddNodalSolutionStepVariable(Shallow.RAIN)
        # Auxiliary variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

    def AddDofs(self):
        raise Exception("Calling the base class instead of the derived one")

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        # Defining the variables
        gravity = self.settings["gravity"].GetDouble()
        time_scale = self.settings["time_scale"].GetString()
        water_height_scale = self.settings["water_height_scale"].GetString()

        # Time unit converter
        if   time_scale == "seconds":
            time_unit_converter =     1
        elif time_scale == "minutes":
            time_unit_converter =    60
        elif time_scale == "hours":
            time_unit_converter =  3600
        elif time_scale == "days":
            time_unit_converter = 86400
        else:
            raise Exception("unknown time scale")

        # Water height unit converter
        if   water_height_scale == "meters":
            water_height_unit_converter = 1.0
        elif water_height_scale == "millimeters":
            water_height_unit_converter = 0.001
        else:
            raise Exception("unknown water height scale")

        # Set ProcessInfo variables
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_Z, gravity * time_unit_converter**2)
        self.main_model_part.ProcessInfo.SetValue(Shallow.TIME_UNIT_CONVERTER, time_unit_converter)
        self.main_model_part.ProcessInfo.SetValue(Shallow.WATER_HEIGHT_UNIT_CONVERTER, water_height_unit_converter)

        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Replace default elements and conditions
            self._ReplaceElementsAndConditions()
            ## Executes the check and prepare model process (Create computing_model_part)
            self._ExecuteCheckAndPrepare()

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def GetComputingModelPart(self):
        return self.main_model_part

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            raise Exception("Estimation Dt utility not yet implemented")

        # Creating the solution strategy for the mesh stage
        self.conv_criteria = KratosMultiphysics.DisplacementCriteria(self.settings["relative_tolerance"].GetDouble(),
                                                                     self.settings["absolute_tolerance"].GetDouble())
        (self.conv_criteria).SetEchoLevel(self.echo_level)

        #~ self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(domain_size,   # DomainSize
                                                                                             domain_size+1) # BlockSize

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.main_model_part,
                                                                            self.time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())

        (self.solver).SetEchoLevel(max(0, self.echo_level-1))
        (self.solver).Check()

        (self.solver).Initialize()

        self.print_on_rank_zero("::[ShallowWaterBaseSolver]::", "Mesh stage solver initialization finished")

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.solver.InitializeSolutionStep()

    def Predict(self):
        if self._TimeBufferIsInitialized():
            self.solver.Predict()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self.solver.SolveSolutionStep()
            return is_converged

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.solver.FinalizeSolutionStep()

    def Check(self):
        self.solver.Check()

    def Clear(self):
        self.solver.Clear()

    #### Specific internal functions ####

    def _IsPrintingRank(self):
        return self._is_printing_rank

    def _TimeBufferIsInitialized(self):
        # We always have one extra old step (step 0, read from input)
        return self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1 >= self.GetMinimumBufferSize()

    def _ComputeDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            delta_time = self.EstimateDeltaTimeUtility.EstimateDt()
        # User-defined delta time
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble()

        # Move particles utility needs to read delta_time
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, delta_time)

        return delta_time

    def _ValidateSettings(self, settings):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type"              : "shallow_water_base_solver",
            "model_part_name"          : "main_model_part",
            "domain_size"              : 2,
            "gravity"                  : 9.81,
            "time_scale"               : "seconds",
            "water_height_scale"       : "meters",
            "model_import_settings"    : {
                "input_type"               : "mdpa",
                "input_filename"           : "unknown_name"
            },
            "echo_level"               : 0,
            "buffer_size"              : 2,
            "dynamic_tau"              : 0.005,
            "dry_height"               : 0.01,
            "relative_tolerance"       : 1e-6,
            "absolute_tolerance"       : 1e-9,
            "maximum_iterations"       : 20,
            "compute_reactions"        : false,
            "reform_dofs_at_each_step" : false,
            "calculate_norm_dx"        : true,
            "move_mesh_flag"           : false,
            "volume_model_part_name"   : "volume_model_part",
            "skin_parts"               : [""],
            "no_skin_parts"            : [""],
            "linear_solver_settings"   : {
                "solver_type"              : "AMGCL",
                "smoother_type"            : "ilu0",
                "krylov_type"              : "gmres",
                "coarsening_type"          : "aggregation",
                "max_iteration"            : 100,
                "tolerance"                : 1.0e-6,
                "scaling"                  : false
            },
            "time_stepping"            : {
                "automatic_time_step"      : false,
                "time_step"                : 0.01
            }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    def _ReplaceElementsAndConditions(self):
        ## Get number of nodes and domain size
        elem_num_nodes = self._GetElementNumNodes()
        cond_num_nodes = self._GetConditionNumNodes()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## Complete the element name
        if (self.element_name is not None):
            new_elem_name = self.element_name + str(int(domain_size)) + "D" + str(int(elem_num_nodes)) + "N"
        else:
            raise Exception("There is no element name. Define the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if (self.condition_name is not None):
            new_cond_name = self.condition_name + str(int(domain_size)) + "D" + str(int(cond_num_nodes)) + "N"
        else:
            raise Exception("There is no condition name. Define the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

    def _GetElementNumNodes(self):
        if self.main_model_part.NumberOfElements() != 0:
            element_num_nodes = len(self.main_model_part.Elements.__iter__().__next__().GetNodes()) # python3 syntax
        else:
            element_num_nodes = 0

        element_num_nodes = self.main_model_part.GetCommunicator().MaxAll(element_num_nodes)
        return element_num_nodes

    def _GetConditionNumNodes(self):
        if self.main_model_part.NumberOfConditions() != 0:
                condition_num_nodes = len(self.main_model_part.Conditions.__iter__().__next__().GetNodes()) # python3 syntax
        else:
            condition_num_nodes = 0

        condition_num_nodes = self.main_model_part.GetCommunicator().MaxAll(condition_num_nodes)
        return condition_num_nodes

    def _ExecuteCheckAndPrepare(self):
        pass