from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from importlib import import_module

# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.PFEM2Application as PFEM2

def CreateSolver(model, custom_settings):
    return PFEM2BaseSolver(model, custom_settings)

class PFEM2BaseSolver(PythonSolver):
    def __init__(self, model, settings):  # Constructor of the class
        self._validate_settings_in_baseclass = True
        super(PFEM2BaseSolver, self).__init__(model, settings)

        ## Set the element and condition names for the replace settings
        ## These should be defined in derived classes
        self.element_name = self.settings["formulation"]["element_type"].GetString()
        self.condition_name = "LineCondition2D3N"
        self.min_buffer_size = 2

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        if self.model.HasModelPart(model_part_name):
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            self.model_part = self.model.CreateModelPart(model_part_name)

        self.domain_size = self.settings["domain_size"].GetInt()
        self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, self.domain_size)

        self.get_particles_stage().ExecuteInitialize()

    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KM.DISTANCE)
        self.model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.model_part.AddNodalSolutionStepVariable(KM.PRESS_PROJ)
        self.model_part.AddNodalSolutionStepVariable(PFEM2.PRESS_PROJ_NO_RO)
        self.model_part.AddNodalSolutionStepVariable(PFEM2.PREVIOUS_ITERATION_PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(PFEM2.PROJECTED_VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(PFEM2.DELTA_VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(PFEM2.MEAN_SIZE)
        self.model_part.AddNodalSolutionStepVariable(KM.YP)
        self.model_part.AddNodalSolutionStepVariable(KM.RHS)
        self.model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        self.model_part.AddNodalSolutionStepVariable(KM.NODAL_AREA)
        self.model_part.AddNodalSolutionStepVariable(KM.NODAL_MASS)
        self.model_part.AddNodalSolutionStepVariable(KM.BODY_FORCE)
        self.model_part.AddNodalSolutionStepVariable(KM.MESH_VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KM.VISCOSITY_AIR)
        self.model_part.AddNodalSolutionStepVariable(KM.VISCOSITY_WATER)
        self.model_part.AddNodalSolutionStepVariable(KM.DENSITY_AIR)
        self.model_part.AddNodalSolutionStepVariable(KM.DENSITY_WATER)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.PRESSURE, self.model_part)
        KM.VariableUtils().AddDof(KM.VELOCITY_X, self.model_part)
        KM.VariableUtils().AddDof(KM.VELOCITY_Y, self.model_part)
        KM.VariableUtils().AddDof(KM.VELOCITY_Z, self.model_part)
        KM.VariableUtils().AddDof(KM.DISTANCE, self.model_part)

    def ImportModelPart(self):
        self._ImportModelPart(self.model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        # Set ProcessInfo variables
        self.model_part.ProcessInfo.SetValue(KM.STEP, 0)

        if not self.model_part.ProcessInfo[KM.IS_RESTARTED]:
            ## Replace default elements and conditions
            self._ReplaceElementsAndConditions()
            ## Set buffer size
            self.model_part.SetBufferSize(self.min_buffer_size)
            ## Import materials and send properties to the nodes and process info
            self._ImportMaterials()

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def GetComputingModelPart(self):
        return self.model_part

    def Initialize(self):
        KM.Logger.PrintInfo("::[PFEM2BaseSolver]:", "Initializing ...")

        self._ApplyBoundaryConditions()

        # The pfem2 mesh solution strategy is created here if it does not already exist.
        mesh_strategy = self.get_mesh_strategy()
        mesh_strategy.SetEchoLevel(max(0, self.echo_level) - 1)
        mesh_strategy.Initialize()

        # The cpp instance for the particles stage is created here
        self.get_particles_stage().ExecuteBeforeSolutionLoop()

        KM.Logger.PrintInfo("::[PFEM2BaseSolver]:", "Finished initialization.")

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt

        self.model_part.CloneTimeStep(new_time)
        self.model_part.ProcessInfo[KM.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.get_particles_stage().ExecuteInitializeSolutionStep()
            self.get_mesh_strategy().InitializeSolutionStep()

    def Predict(self):
        if self._TimeBufferIsInitialized():
            self.get_mesh_strategy().Predict()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self.get_mesh_strategy().SolveSolutionStep()
            return is_converged
        else:
            return True

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.get_mesh_strategy().FinalizeSolutionStep()
            self.get_particles_stage().ExecuteFinalizeSolutionStep()

    def Check(self):
        self.get_mesh_strategy().Check()

    def Clear(self):
        self.get_mesh_strategy().Clear()

    #### Specific internal functions ####

    def _TimeBufferIsInitialized(self):
        return self.model_part.ProcessInfo[KM.STEP] >= self.GetMinimumBufferSize()

    def _ComputeDeltaTime(self):
        delta_time = 0.0
        if self.settings["time_stepping"]["automatic_time_step"].GetBool():
            raise Exception("Automatic time step not implemented in PFEM2 Application")
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
        return delta_time

    @classmethod
    def GetDefaultSettings(cls):
        default_settings = KM.Parameters("""
        {
            "solver_type"              : "pfem_2_solver",
            "model_part_name"          : "model_part",
            "domain_size"              : 2,
            "model_import_settings"    : {
                "input_type"               : "mdpa",
                "input_filename"           : "unknown_name"
            },
            "material_import_settings": {
                "materials_filename"       : "ProjectParameters_materials.json"
            },
            "formulation"              : {
                "element_type"             : "MonolithicPFEM2"
            },
            "echo_level"               : 0,
            "buffer_size"              : 2,
            "relative_tolerance"       : 1e-6,
            "absolute_tolerance"       : 1e-9,
            "maximum_iterations"       : 20,
            "compute_reactions"        : false,
            "reform_dofs_at_each_step" : false,
            "calculate_norm_dx"        : true,
            "move_mesh_flag"           : false,
            "skin_parts"               : [],
            "no_skin_parts"            : [],
            "volume_model_part_name"   : "name",
            "linear_solver_settings"   : {
                "solver_type"              : "amgcl"
            },
            "time_stepping"            : {
                "automatic_time_step"      : false,
                "time_step"                : 0.01
            },
            "particles_stage_settings" : {
                "convection_type"          : "pfem_2"
            }
        }""")
        default_settings.AddMissingParameters(super(PFEM2BaseSolver,cls).GetDefaultSettings())
        return default_settings

    def _ReplaceElementsAndConditions(self):
        ## Complete the element name
        if self.element_name is None:
            raise Exception("There is no element name. Define the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if self.condition_name is None:
            raise Exception("There is no condition name. Define the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        self.settings.AddValue("element_replace_settings", KM.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(self.element_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(self.condition_name)

        ## Call the replace elements and conditions process
        KM.ReplaceElementsAndConditionsProcess(self.model_part, self.settings["element_replace_settings"]).Execute()

    def _ApplyBoundaryConditions(self):
        KM.BodyNormalCalculationUtils().CalculateBodyNormals(self.model_part, self.domain_size)
        if self.domain_size == 2:
            self.addBC = PFEM2.AddFixedPressureCondition2D(self.model_part)
        else:
            self.addBC = PFEM2.AddFixedPressureCondition3D(self.model_part)
        self.addBC.AddThem()

    def _ImportMaterials(self):
        print(self.settings)
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            material_settings = KM.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KM.ReadMaterialsUtility(material_settings, self.model)

            self._assign_nodally_properties()
            KM.Logger.PrintInfo("::[PFEM2BaseSolver]:", "Materials were successfully imported.")
        else:
            KM.Logger.PrintInfo("::[PFEM2BaseSolver]:", "Materials were not imported.")

    def _assign_nodally_properties(self):
        with open(self.settings["material_import_settings"]["materials_filename"].GetString(), 'r') as parameter_file:
            materials = KM.Parameters(parameter_file.read())

        for i in range(materials["properties"].size()):
            model_part = self.model.GetModelPart(materials["properties"][i]["model_part_name"].GetString())
            mat = materials["properties"][i]["Material"]

            for key, value in mat["Variables"].items():
                var = KM.KratosGlobals.GetVariable(key)
                if value.IsDouble():
                    KM.VariableUtils().SetScalarVar(var, value.GetDouble(), model_part.Nodes)
                    model_part.ProcessInfo.SetValue(var, value.GetDouble())
                elif value.IsVector():
                    KM.VariableUtils().SetVectorVar(var, value.GetVector(), model_part.Nodes)
                    model_part.ProcessInfo.SetValue(var, value.GetVector())
                else:
                    raise ValueError("Type of value is not available")

    def get_time_scheme(self):
        if not hasattr(self, '_time_scheme'):
            self._time_scheme = self._create_time_scheme()
        return self._time_scheme

    def get_convergence_criterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._create_convergence_criterion()
        return self._convergence_criterion

    def get_linear_solver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._create_linear_solver()
        return self._linear_solver

    def get_builder_and_solver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._create_builder_and_solver()
        return self._builder_and_solver

    def get_mesh_strategy(self):
        if not hasattr(self, '_mesh_strategy'):
            self._mesh_strategy = self._create_mesh_strategy()
        return self._mesh_strategy

    def get_particles_stage(self):
        if not hasattr(self, '_particles_stage'):
            self._particles_stage = self._create_particles_stage()
        return self._particles_stage

    def _create_convergence_criterion(self):
        convergence_criterion = KM.DisplacementCriteria(self.settings["relative_tolerance"].GetDouble(),
                                                        self.settings["absolute_tolerance"].GetDouble())
        convergence_criterion.SetEchoLevel(self.echo_level)
        return convergence_criterion

    def _create_linear_solver(self):
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        builder_and_solver = KM.ResidualBasedBlockBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _create_time_scheme(self):
        """Create the solution scheme for the Navier-Stokes problem.
        """
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _create_mesh_strategy(self):
        model_part = self.GetComputingModelPart()
        time_scheme = self.get_time_scheme()
        linear_solver = self.get_linear_solver()
        convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return KM.ResidualBasedNewtonRaphsonStrategy(model_part,
                                                     time_scheme,
                                                     linear_solver,
                                                     convergence_criterion,
                                                     builder_and_solver,
                                                     self.settings["maximum_iterations"].GetInt(),
                                                     self.settings["compute_reactions"].GetBool(),
                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                     self.settings["move_mesh_flag"].GetBool())

    def _create_particles_stage(self):
        convection_settings = KM.Parameters('''{"Parameters" : {}}''')
        convection_settings["Parameters"] = self.settings["particles_stage_settings"].Clone()
        convection_settings["Parameters"].RemoveValue("convection_type")
        process_name = "KratosMultiphysics."
        convection_type = self.settings["particles_stage_settings"]["convection_type"].GetString()
        if convection_type == "pfem_2":
            process_name += "PFEM2Application.pfem_2_process"
            if not convection_settings["Parameters"].Has("model_part_name"):
                convection_settings["Parameters"].AddValue("model_part_name", self.settings["model_part_name"])
        else:
            raise Exception("The requested particles stage type is not available: " + convection_type)
        python_module = import_module(process_name)
        print(convection_settings)
        return python_module.Factory(convection_settings, self.model)
