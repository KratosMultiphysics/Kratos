from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return MechanicalSolver(main_model_part, custom_settings)

#Base class to develop other solvers
class MechanicalSolver(object):

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
            "solver_type": "solid_mechanics_solver",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Dynamic",
            "analysis_type": "Non-Linear",
            "time_integration_method": "Implicit",
            "scheme_type": "Newmark",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": null,
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
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "bodies_list": [],
            "problem_domain_sub_model_part_list": ["solid"],
            "processes_sub_model_part_list": [""]
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Warning: Construction of Base Mechanical Solver finished")


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

        print("::[Mechanical Solver]:: Variables ADDED")


    def GetMinimumBufferSize(self):
        return 2;

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

        print("::[Mechanical Solver]:: DOF's ADDED")


    def ImportModelPart(self):

        print("::[Mechanical Solver]:: Model reading starts.")

        self.computing_model_part_name = "computing_domain" #this submodelpart will be labeled with KratosMultiphysics.ACTIVE flag, you can recover it checking the flag.

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):

            # Model part reading
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            print("    Import input model part.")

            # Check and prepare model process and construct constitutive law
            self._ExecuteAfterReading()

            # Set and fill buffer
            self._SetAndFillBuffer()

        elif(self.settings["model_import_settings"]["input_type"].GetString() == "rest"):

            problem_path = os.getcwd()
            restart_path = os.path.join(problem_path, self.settings["model_import_settings"]["input_filename"].GetString() + "__" + self.settings["model_import_settings"]["input_file_label"].GetString() )

            if(os.path.exists(restart_path+".rest") == False):
                print("    rest file does not exist , check the restart step selected ")

            print("    Load Restart file: ", self.settings["model_import_settings"]["input_filename"].GetString() + "__" + self.settings["model_import_settings"]["input_file_label"].GetString())
            # set serializer flag
            self.serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE      # binary
            # self.serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR # ascii
            # self.serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL   # ascii

            serializer = KratosMultiphysics.Serializer(restart_path, self.serializer_flag)

            serializer.Load(self.main_model_part.Name, self.main_model_part)

            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True

            print(self.main_model_part)

        else:
            raise Exception("Other input options are not yet implemented.")

        print ("::[Mechanical Solver]:: Model reading finished.")

    def Initialize(self):
        raise Exception("please implement the Custom Initialization of your solver")

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.mechanical_solver.Solve()

    # solve :: sequencial calls

    def InitializeStrategy(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False ):
            self.mechanical_solver.Initialize()

    def InitializeSolutionStep(self):
        self.mechanical_solver.InitializeSolutionStep()

    def Predict(self):
        self.mechanical_solver.Predict()

    def SolveSolutionStep(self):
        self.mechanical_solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.mechanical_solver.FinalizeSolutionStep()

    # solve :: sequencial calls

    def SetEchoLevel(self, level):
        self.mechanical_solver.SetEchoLevel(level)

    def Clear(self):
        self.mechanical_solver.Clear()

    def Check(self):
        self.mechanical_solver.Check()

    #### Specific internal functions ####

    def _ExecuteAfterReading(self):
        self.computing_model_part_name = "computing_domain" #this submodelpart will be labeled with KratosMultiphysics.ACTIVE flag, you can recover it checking the flag.

        # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])

        if( self.settings.Has("bodies_list") ):
            params.AddValue("bodies_list",self.settings["bodies_list"])

        # CheckAndPrepareModelProcess creates the solid_computational model part
        import check_and_prepare_model_process_solid
        check_and_prepare_model_process_solid.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # Constitutive law import
        import constitutive_law_python_utility as constitutive_law_utils
        constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(self.main_model_part,
                                                                         self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]);
        constitutive_law.Initialize();
        print("    Constitutive law initialized.")

    def _SetAndFillBuffer(self):
        # Set buffer size
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
            #delta_time is computed from previous time in process_info
            self.main_model_part.CloneTimeStep(time)

        self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    def _GetSolutionScheme(self, scheme_type, component_wise, compute_contact_forces):
        raise Exception("please implement the Custom Choice of your Scheme (_GetSolutionScheme) in your solver")

    def _GetConvergenceCriterion(self):
        # Creation of an auxiliar Kratos parameters object to store the convergence settings
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("rotation_dofs",self.settings["rotation_dofs"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("component_wise",self.settings["component_wise"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])

        # Construction of the class convergence_criterion
        import convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.convergence_criterion(conv_params)

        return convergence_criterion.mechanical_convergence_criterion

    def _GetBuilderAndSolver(self, component_wise, block_builder):
        # Creating the builder and solver
        if(component_wise):
            builder_and_solver = KratosSolid.ComponentWiseBuilderAndSolver(self.linear_solver)
        else:
            if(block_builder):
                # To keep matrix blocks in builder
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(self.linear_solver)

        return builder_and_solver

    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search, implex):
        raise Exception("please implement the Custom Choice of your Mechanical Solver (_GetMechanicalSolver) in your solver")
