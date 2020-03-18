from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.DelaunayMeshingApplication  as KratosDelaunay
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
        print("Construction of 2-step Pfem Fluid Solver finished.")
        super(PfemFluidSolver, self).__init__(model, parameters)

        model_part_name = self.settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)


    @classmethod
    def GetDefaultSettings(cls):
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
        "processes_sub_model_part_list": [],
        "constraints_process_list": [],
        "loads_process_list"       : [],
        "output_process_list"      : [],
        "output_configuration"     : {},
        "problem_process_list"     : [],
        "processes"                : {},
        "output_processes"         : {},
        "check_process_list": []
        }""")
        this_defaults.AddMissingParameters(super(PfemFluidSolver, cls).GetDefaultSettings())
        return this_defaults


    def GetMinimumBufferSize(self):
        return 3

    def Initialize(self):

        print("::[Pfem Fluid Solver]:: -START-")

        print(self.main_model_part.SetBufferSize(self.settings["buffer_size"].GetInt()))

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

        echo_level = self.settings["echo_level"].GetInt()

        # Set echo_level
        self.fluid_solver.SetEchoLevel(echo_level)

        # Set initialize flag
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.mechanical_solver.SetInitializePerformedFlag(True)

        # Check if everything is assigned correctly
        self.fluid_solver.Check()
        print("::[Pfem Fluid Solver]:: -END- ")


    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POISSON_RATIO)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YOUNG_MODULUS)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_ERROR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)


        #VARIABLES FOR PAPANASTASIOU MODEL
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FLOW_INDEX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELD_SHEAR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT)

        #VARIABLES FOR MU-I RHEOLOGY MODEL
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.STATIC_FRICTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.DYNAMIC_FRICTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ZERO)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DIAMETER)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INFINITE_FRICTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ONE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ALPHA_PARAMETER)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION_OLD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION_RATE)

        # PFEM fluid variables
        # self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NORMVELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELDED)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FREESURFACE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_ACCELERATION)

        # Pfem Extra Vars
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_NORMAL)

        self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.OFFSET)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.SHRINK_FACTOR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.MEAN_ERROR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.RIGID_WALL)

        self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.PROPERTY_ID)

        print("::[PfemFluidSolver]:: Variables ADDED")


    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.PRESSURE)
            node.AddDof(KratosMultiphysics.DENSITY)
            node.AddDof(KratosMultiphysics.VELOCITY_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)
        print("::[Pfem Fluid Solver]:: DOF's ADDED")

    def ImportModelPart(self):

        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):

        print("::[PfemFluidSolver]:: Model preparing started.")

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

        print ("::[Pfem Fluid Solver]:: Model preparing finished.")


    def CheckAndPrepareModelProcess(self, params):
        # CheckAndPrepareModelProcess creates the fluid_computational model part
        from KratosMultiphysics.PfemFluidDynamicsApplication import pfem_check_and_prepare_model_process_fluid
        pfem_check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

    def _ComputeDeltaTime(self):

        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        return delta_time

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.fluid_solver.Solve()

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeStrategy(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.fluid_solver.Initialize()

    def InitializeSolutionStep(self):
        #self.fluid_solver.InitializeSolutionStep()
        if self._TimeBufferIsInitialized():
            self.fluid_solver.InitializeSolutionStep()

        ## Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            adaptive_time_interval = KratosPfemFluid.AdaptiveTimeIntervalProcess(self.main_model_part,self.settings["echo_level"].GetInt())
            adaptive_time_interval.Execute()

    def Predict(self):
        pass

    def SolveSolutionStep(self):
        is_converged = True
        self.fluid_solver.Solve()
        return is_converged

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


