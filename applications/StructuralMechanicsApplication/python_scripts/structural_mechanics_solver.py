from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return MechanicalSolver(main_model_part, custom_settings)


class MechanicalSolver(object):
    """The base class for structural mechanics solvers.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes must override the function _create_solution_scheme which
    constructs and returns a solution scheme. Depending on the type of
    solver, derived classes may also need to override the following functions:

    _create_solution_scheme
    _create_convergence_criterion
    _create_linear_solver
    _create_builder_and_solver
    _create_mechanical_solver

    The mechanical_solver, builder_and_solver, etc. should alway be retrieved
    using the getter functions get_mechanical_solver, get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    settings -- Kratos parameters containing solver settings.
    main_model_part -- the model part used to construct the solver.
    """
    def __init__(self, main_model_part, custom_settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "Static",
            "echo_level": 0,
            "buffer_size": 2,
            "analysis_type": "non_linear",
            "time_integration_method": "implicit",
            "scheme_type": "newmark",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "computing_model_part_name" : "computing_domain",
            "material_import_settings" :{
                "materials_filename": ""
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": true,
            "clear_storage": false,
            "move_mesh_flag": true,
            "multi_point_constraints_used": false,
            "convergence_criterion": "residual_criterion",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "component_wise" : false,
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

        # Overwrite the default settings with user-provided parameters.
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part
        print("::[MechanicalSolver]:: Construction finished")

    def AddVariables(self):
        # Add displacements.
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        # Add specific variables for the problem conditions.
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        if self.settings["rotation_dofs"].GetBool():
            # Add specific variables for the problem (rotation dofs).
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LOCAL_POINT_MOMENT)
            # TODO: Can we combine POINT_TORQUE and POINT_MOMENT???
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_TORQUE)
        if self.settings["pressure_dofs"].GetBool(): # TODO: The creation of UP and USigma elements is pending
            # Add specific variables for the problem (pressure dofs).
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.PRESSURE_REACTION)
        print("::[MechanicalSolver]:: Variables ADDED")

    def GetMinimumBufferSize(self):
        return 2

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z,self.main_model_part)
        if self.settings["pressure_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE_REACTION,self.main_model_part)
        print("::[MechanicalSolver]:: DOF's ADDED")

    def ImportModelPart(self):
        print("::[MechanicalSolver]:: Importing model part.")
        problem_path = os.getcwd()
        input_filename = self.settings["model_import_settings"]["input_filename"].GetString()
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            # Import model part from mdpa file.
            print("    Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa")
            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)
            print("    Finished reading model part from mdpa file.")
            # Check and prepare computing model part and import constitutive laws.
            self._execute_after_reading()
            self._set_and_fill_buffer()
        elif(self.settings["model_import_settings"]["input_type"].GetString() == "rest"):
            # Import model part from restart file.
            problem_path = os.getcwd()
            restart_path = os.path.join(problem_path, input_filename + "__" + self.settings["model_import_settings"]["input_file_label"].GetString())
            if(os.path.exists(restart_path+".rest") == False):
                raise Exception("Restart file not found: " + restart_path + ".rest")
            print("    Loading Restart file: ", restart_path + ".rest")
            serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE      # binary
            # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR # ascii
            # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL   # ascii
            serializer = KratosMultiphysics.Serializer(restart_path, serializer_flag)
            serializer.Load(self.main_model_part.Name, self.main_model_part)
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
            #I use it to rebuild the contact conditions.
            load_step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1
            self.main_model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step
            print("    Finished loading model part from restart file.")
        else:
            raise Exception("Other model part input options are not yet implemented.")
        print(self.main_model_part)
        print("::[MechanicalSolver]:: Finished importing model part.")

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        print("::[MechanicalSolver]:: Initializing ...")
        # The mechanical solver is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        mechanical_solver = self.get_mechanical_solver()
        mechanical_solver.SetEchoLevel(self.settings["echo_level"].GetInt())
        if (self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
            mechanical_solver.Initialize()
        else:
            # SetInitializePerformedFlag is not a member of SolvingStrategy but
            # is used by ResidualBasedNewtonRaphsonStrategy.
            if hasattr(mechanical_solver, SetInitializePerformedFlag):
                mechanical_solver.SetInitializePerformedFlag(True)
        self.Check()
        print("::[MechanicalSolver]:: Finished initialization.")

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        mechanical_solver = self.get_mechanical_solver()
        mechanical_solver.Solve()

    def InitializeSolutionStep(self):
        self.get_mechanical_solver().InitializeSolutionStep()

    def Predict(self):
        self.get_mechanical_solver().Predict()

    def SolveSolutionStep(self):
        is_converged = self.get_mechanical_solver().SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.get_mechanical_solver().FinalizeSolutionStep()

    def SetEchoLevel(self, level):
        self.get_mechanical_solver().SetEchoLevel(level)

    def Clear(self):
        self.get_mechanical_solver().Clear()

    def Check(self):
        self.get_mechanical_solver().Check()

    #### Specific internal functions ####

    def get_solution_scheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._create_solution_scheme()
        return self._solution_scheme

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

    def get_mechanical_solver(self):
        if not hasattr(self, '_mechanical_solver'):
            self._mechanical_solver = self._create_mechanical_solver()
        return self._mechanical_solver

    def import_constitutive_laws(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            import read_materials_process
            # Create a dictionary of model parts.
            Model = {self.main_model_part.Name : self.main_model_part}
            for i in range(self.settings["problem_domain_sub_model_part_list"].size()):
                part_name = self.settings["problem_domain_sub_model_part_list"][i].GetString()
                Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})
            for i in range(self.settings["processes_sub_model_part_list"].size()):
                part_name = self.settings["processes_sub_model_part_list"][i].GetString()
                Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})
            # Add constitutive laws and material properties from json file to model parts.
            read_materials_process.ReadMaterialsProcess(Model, self.settings["material_import_settings"])
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    def validate_and_transfer_matching_settings(self, origin_settings, destination_settings):
        """Transfer matching settings from origin to destination.

        If a name in origin matches a name in destination, then the setting is
        validated against the destination.

        The typical use is for validating and extracting settings in derived classes:

        class A:
            def __init__(self, model_part, a_settings):
                default_a_settings = Parameters('''{
                    ...
                }''')
                a_settings.ValidateAndAssignDefaults(default_a_settings)
        class B(A):
            def __init__(self, model_part, custom_settings):
                b_settings = Parameters('''{
                    ...
                }''') # Here the settings contain default values.
                self.validate_and_transfer_matching_settings(custom_settings, b_settings)
                super().__init__(model_part, custom_settings)
        """
        for name, dest_value in destination_settings.items():
            if origin_settings.Has(name): # Validate and transfer value.
                orig_value = origin_settings[name]
                if dest_value.IsDouble() and orig_value.IsDouble():
                    destination_settings[name].SetDouble(origin_settings[name].GetDouble())
                elif dest_value.IsInt() and orig_value.IsInt():
                    destination_settings[name].SetInt(origin_settings[name].GetInt())
                elif dest_value.IsBool() and orig_value.IsBool():
                    destination_settings[name].SetBool(origin_settings[name].GetBool())
                elif dest_value.IsString() and orig_value.IsString():
                    destination_settings[name].SetString(origin_settings[name].GetString())
                elif dest_value.IsArray() and orig_value.IsArray():
                    if dest_value.size() != orig_value.size():
                        raise Exception('len("' + name + '") != ' + str(dest_value.size()))
                    for i in range(dest_value.size()):
                        if dest_value[i].IsDouble() and orig_value[i].IsDouble():
                            dest_value[i].SetDouble(orig_value[i].GetDouble())
                        elif dest_value[i].IsInt() and orig_value[i].IsInt():
                            dest_value[i].SetInt(orig_value[i].GetInt())
                        elif dest_value[i].IsBool() and orig_value[i].IsBool():
                            dest_value[i].SetBool(orig_value[i].GetBool())
                        elif dest_value[i].IsString() and orig_value[i].IsString():
                            dest_value[i].SetString(orig_value[i].GetString())
                        elif dest_value[i].IsSubParameter() and orig_value[i].IsSubParameter():
                            self.validate_and_transfer_matching_settings(orig_value[i], dest_value[i])
                            if len(orig_value[i].items()) != 0:
                                raise Exception('Json settings not found in default settings: ' + orig_value[i].PrettyPrintJsonString())
                        else:
                            raise Exception('Unsupported parameter type.')
                elif dest_value.IsSubParameter() and orig_value.IsSubParameter():
                    self.validate_and_transfer_matching_settings(orig_value, dest_value)
                    if len(orig_value.items()) != 0:
                        raise Exception('Json settings not found in default settings: ' + orig_value.PrettyPrintJsonString())
                else:
                    raise Exception('Unsupported parameter type.')
                origin_settings.RemoveValue(name)

    #### Private functions ####

    def _execute_after_reading(self):
        """Prepare computing model part and import constitutive laws. """
        # Auxiliary parameters object for the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("computing_model_part_name",self.settings["computing_model_part_name"])
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
        if( self.settings.Has("bodies_list") ):
            params.AddValue("bodies_list",self.settings["bodies_list"])
        # Assign mesh entities from domain and process sub model parts to the computing model part.
        import check_and_prepare_model_process_structural
        check_and_prepare_model_process_structural.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # Import constitutive laws.
        materials_imported = self.import_constitutive_laws()
        if materials_imported:
            print("    Constitutive law was successfully imported.")
        else:
            print("    Constitutive law was not imported.")

    def _set_and_fill_buffer(self):
        """Prepare nodal solution step data containers and time step information. """
        # Set the buffer size for the nodal solution steps data. Existing nodal
        # solution step data may be lost.
        buffer_size = self.settings["buffer_size"].GetInt()
        if buffer_size < self.GetMinimumBufferSize():
            buffer_size = self.GetMinimumBufferSize()
        self.main_model_part.SetBufferSize(buffer_size)
        # Cycle the buffer. This sets all historical nodal solution step data to
        # the current value and initializes the time stepping in the process info.
        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        step =-buffer_size
        time = time - delta_time * buffer_size
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for i in range(0, buffer_size):
            step = step + 1
            time = time + delta_time
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            self.main_model_part.CloneTimeStep(time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    def _add_dynamic_variables(self):
        # For being consistent for Serial and Trilinos
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)

    def _add_dynamic_dofs(self):
        # For being consistent for Serial and Trilinos
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_Z,self.main_model_part)
        if(self.settings["rotation_dofs"].GetBool()):
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_VELOCITY_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Z,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Z,self.main_model_part)

    def _create_convergence_criterion(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("rotation_dofs",self.settings["rotation_dofs"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])
        import convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.convergence_criterion(conv_params)
        return convergence_criterion.mechanical_convergence_criterion

    def _create_linear_solver(self):
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        if self.settings["block_builder"].GetBool():
            if self.settings["multi_point_constraints_used"].GetBool():
                builder_and_solver = KratosMultiphysics.StructuralMechanicsApplication.ResidualBasedBlockBuilderAndSolverWithMpc(linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        else:
            if self.settings["multi_point_constraints_used"].GetBool():
                raise Exception("To use MPCs you also have to set \"block_builder\" to \"true\"")
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _create_solution_scheme(self):
        """Create the solution scheme for the structural problem.
        """
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _create_mechanical_solver(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            mechanical_solver = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            mechanical_solver = self._create_newton_raphson_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return mechanical_solver

    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part, 
                                                              mechanical_scheme, 
                                                              linear_solver, 
                                                              builder_and_solver, 
                                                              self.settings["compute_reactions"].GetBool(), 
                                                              self.settings["reform_dofs_at_each_step"].GetBool(), 
                                                              False, 
                                                              self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part, 
                                                                     mechanical_scheme, 
                                                                     linear_solver, 
                                                                     mechanical_convergence_criterion, 
                                                                     builder_and_solver,
                                                                     self.settings["max_iteration"].GetInt(),
                                                                     self.settings["compute_reactions"].GetBool(),
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.settings["move_mesh_flag"].GetBool())
