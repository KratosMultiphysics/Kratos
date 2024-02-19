import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
from KratosMultiphysics.python_solver import PythonSolver


def CreateSolver(model, parameters):
    return PfemFluidSolver(model, parameters)

class PfemFluidSolver(PythonSolver):

    def __init__(self, model, parameters):

        self._validate_settings_in_baseclass=True # To be removed eventually

        super(PfemFluidSolver,self).__init__(model, parameters)

        #construct the linear solver
        from KratosMultiphysics import python_linear_solver_factory
        self.pressure_linear_solver = python_linear_solver_factory.ConstructSolver(self.settings["pressure_linear_solver_settings"])
        self.velocity_linear_solver = python_linear_solver_factory.ConstructSolver(self.settings["velocity_linear_solver_settings"])

        self.compute_reactions = self.settings["compute_reactions"].GetBool()

        super(PfemFluidSolver, self).__init__(model, parameters)

        model_part_name = self.settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)


    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
             "solver_type": "pfem_fluid_solver",
            "model_part_name": "PfemFluidModelPart",
            "physics_type"   : "fluid",
            "domain_size": 2,
            "time_stepping"               : {
                "automatic_time_step" : false,
                "time_step"           : 0.001
            },
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings"           : {
                "materials_filename" : "unknown_name"
            },
            "buffer_size": 3,
            "echo_level": 1,
            "reform_dofs_at_each_step": false,
            "clear_storage": false,
            "compute_reactions": true,
            "move_mesh_flag": true,
            "dofs"                : [],
            "stabilization_factor": 1.0,
            "line_search": false,
            "compute_contact_forces": false,
            "block_builder": false,
            "component_wise": false,
            "predictor_corrector": true,
            "time_order": 2,
            "maximum_velocity_iterations": 1,
            "maximum_pressure_iterations": 7,
            "velocity_tolerance": 1e-5,
            "pressure_tolerance": 1e-5,
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "amgcl",
                "max_iteration"                  : 5000,
                "tolerance"                      : 1e-9,
                "provide_coordinates"            : false,
                "scaling"                        : false,
                "smoother_type"                  : "damped_jacobi",
                "krylov_type"                    : "cg",
                "coarsening_type"                : "aggregation",
                "verbosity"                      : 0
            },
            "velocity_linear_solver_settings": {
                "solver_type"                    : "bicgstab",
                "max_iteration"                  : 5000,
                "tolerance"                      : 1e-9,
                "preconditioner_type"            : "none",
                "scaling"                        : false
            },
            "solving_strategy_settings":{
               "time_step_prediction_level": 0,
               "max_delta_time": 1.0e-5,
               "fraction_delta_time": 0.9,
               "rayleigh_damping": false,
               "rayleigh_alpha": 0.0,
               "rayleigh_beta" : 0.0
            },
        "bodies_list": [],
        "problem_domain_sub_model_part_list": [],
        "constitutive_laws_list": [],
        "processes_sub_model_part_list": [],
        "constraints_process_list": [],
        "loads_process_list"       : [],
        "output_process_list"      : [],
        "output_configuration"     : {},
        "problem_process_list"     : [],
        "processes"                : {},
        "output_processes"         : {},
        "check_process_list": [],
            "penalty_coefficient" : 0.0
        }""")
        this_defaults.AddMissingParameters(super(PfemFluidSolver, cls).GetDefaultParameters())
        return this_defaults


    def GetMinimumBufferSize(self):
        return 3

    def Initialize(self):

        # Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        self.fluid_solver = KratosPfemFluid.TwoStepVPStrategy(self.computing_model_part,
                                                               self.velocity_linear_solver,
                                                               self.pressure_linear_solver,
                                                               self.settings["reform_dofs_at_each_step"].GetBool(),
                                                               self.settings["velocity_tolerance"].GetDouble(),
                                                               self.settings["pressure_tolerance"].GetDouble(),
                                                               self.settings["maximum_pressure_iterations"].GetInt(),
                                                               self.settings["time_order"].GetInt(),
                                                               self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION])

        # Set echo_level
        echo_level = self.settings["echo_level"].GetInt()
        self.fluid_solver.SetEchoLevel(echo_level)

        # Self initialize strategy
        self.fluid_solver.Initialize()

        # Check if everything is assigned correctly
        self.fluid_solver.Check()
        
        # Set penalty coefficient for Cut-PFEM
        #TODO: Create a Cut-PFEM solver deriving from this one and do this in there
        self.computing_model_part.ProcessInfo[KratosPfemFluid.PENALTY_COEFFICIENT] = self.settings["penalty_coefficient"].GetDouble()

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.PRESSURE)
            node.AddDof(KratosMultiphysics.VELOCITY_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z)

    def ImportModelPart(self):

        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):

        self.computing_model_part_name = "fluid_computing_domain"

        # Import materials and relative tables
        materials_imported = self.ImportMaterials()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[PfemFluidSolver]:: ", "Materials were successfully imported.")
        else:
            KratosMultiphysics.Logger.PrintInfo("::[PfemFluidSolver]:: ", "Materials were not imported.")

        # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("constitutive_laws_list",self.settings["constitutive_laws_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
        params.AddValue("material_import_settings",self.settings["material_import_settings"])
        if( self.settings.Has("bodies_list") ):
            params.AddValue("bodies_list",self.settings["bodies_list"])

        self.CheckAndPrepareModelProcess(params)
        self.main_model_part.SetBufferSize( self.settings["buffer_size"].GetInt() )

        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            current_buffer_size = self.GetMinimumBufferSize()

        self.main_model_part.SetBufferSize( current_buffer_size )

        # Fill buffer
        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        time = time - delta_time * (current_buffer_size)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for size in range(0, current_buffer_size):
            step = size - (current_buffer_size -1)
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            self.main_model_part.CloneTimeStep(time)

        self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

        if (abs(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) < 1e-5 * delta_time):
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)


    def CheckAndPrepareModelProcess(self, params):
        # CheckAndPrepareModelProcess creates the fluid_computational model part
        from KratosMultiphysics.PfemFluidDynamicsApplication import pfem_check_and_prepare_fluid_model_process
        pfem_check_and_prepare_fluid_model_process.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

    def _ComputeDeltaTime(self):

        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        return delta_time

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        self.fluid_solver.InitializeSolutionStep()
        # if self._TimeBufferIsInitialized():
        #     self.fluid_solver.InitializeSolutionStep()

        ## Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            adaptive_time_interval = KratosPfemFluid.AdaptiveTimeIntervalProcess(self.main_model_part,self.settings["echo_level"].GetInt())
            adaptive_time_interval.Execute()

    def Predict(self):
        pass

    def SolveSolutionStep(self):
        converged = self.fluid_solver.SolveSolutionStep()
        return converged

    def FinalizeSolutionStep(self):
        #pass
        self.fluid_solver.FinalizeSolutionStep()

        unactive_peak_elements = False
        unactive_sliver_elements = False
        if(unactive_peak_elements == True or unactive_sliver_elements == True):
            set_active_flag = KratosPfemFluid.SetActiveFlagProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.settings["echo_level"].GetInt())
            set_active_flag.Execute()

    def SetEchoLevel(self, level):
        self.fluid_solver.SetEchoLevel(level)

    def Clear(self):
        self.fluid_solver.Clear()

    def Check(self):
        self.fluid_solver.Check()

    def ImportMaterials(self):
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

    def _TimeBufferIsInitialized(self):
        # We always have one extra old step (step 0, read from input)
        return self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1 >= self.GetMinimumBufferSize()


