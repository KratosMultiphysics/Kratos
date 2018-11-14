from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return FemDemMechanicalSolver(main_model_part, custom_settings)

#Base class to develop other solvers
class FemDemMechanicalSolver(object):
    """The base class for solid mechanics solvers.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes must override the function _create_mechanical_solver which
    constructs and returns a valid solving strategy. Depending on the type of
    solver, derived classes may also need to override the following functions:

    _create_solution_scheme
    _create_convergence_criterion
    _create_linear_solver
    _create_builder_and_solver
    _create_mechanical_solver

    The mechanical_solver, builder_and_solver, etc. should alway be retrieved
    using the getter functions _get_mechanical_solver, _get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    settings -- Kratos parameters containing solver settings.
    main_model_part -- the model part used to construct the solver.
    """
    def __init__(self, main_model_part, custom_settings):         
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "solid_mechanics_solver",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Dynamic",
            "time_integration_method": "Implicit",
            "scheme_type": "Newmark",
	    "analysis_type": "Non-Linear",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "computing_model_part_name" : "computing_domain",
            "dofs": [],
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "implex": false,
            "stabilization_factor": null,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "axisymmetric": false,
            "block_builder": true,
            "move_mesh_flag": true,
            "clear_storage": false,
            "component_wise": false,
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "extrapolation_required" : false,
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


        #trick to allow null value in a stabilization_factor variable
        if(custom_settings.Has("stabilization_factor")):
            if(custom_settings["stabilization_factor"].IsDouble()):
                default_settings["stabilization_factor"].SetDouble(0.0)

        
        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        print("::[Solid_Mechanical_Solver]:: Constructed")

        
    def GetMinimumBufferSize(self):
        return 2;

    def SetVariables(self):
        
        self.nodal_variables = []
        self.dof_variables   = []
        self.dof_reactions   = [] 
        
        # Add displacements
        self.dof_variables = self.dof_variables + ['DISPLACEMENT']
        self.dof_reactions = self.dof_reactions + ['REACTION'] 
        
        # Add dynamic variables
        if(self.settings["solution_type"].GetString() == "Dynamic" or (self.settings["scheme_type"].GetString() != "Linear")):
            self.dof_variables = self.dof_variables + ['VELOCITY','ACCELERATION']
            self.dof_reactions = self.dof_reactions + ['NOT_DEFINED','NOT_DEFINED']
        
        # Add specific variables for the problem conditions
        self.nodal_variables = self.nodal_variables + ['VOLUME_ACCELERATION','POSITIVE_FACE_PRESSURE','NEGATIVE_FACE_PRESSURE','POINT_LOAD','LINE_LOAD','SURFACE_LOAD']
        
        # Add nodal force variables for component wise calculation
        if( self.settings.Has("component_wise") ):
            if self.settings["component_wise"].GetBool():
                self.nodal_variables = self.nodal_variables + ['INTERNAL_FORCE','EXTERNAL_FORCE']
 
        # Add rotational variables
        if self._check_input_dof("ROTATION"):                    
            # Add specific variables for the problem (rotation dofs)
            self.dof_variables = self.dof_variables + ['ROTATION']
            self.dof_reactions = self.dof_reactions + ['TORQUE']
            if(self.settings["solution_type"].GetString() == "Dynamic" or (self.settings["scheme_type"].GetString() != "Linear")):
                self.dof_variables = self.dof_variables + ['ANGULAR_VELOCITY','ANGULAR_ACCELERATION']
                self.dof_reactions = self.dof_reactions + ['NOT_DEFINED','NOT_DEFINED']
            # Add specific variables for the problem conditions
            self.nodal_variables = self.nodal_variables + ['POINT_MOMENT']
            # Add large rotation variables
            self.nodal_variables = self.nodal_variables + ['STEP_DISPLACEMENT','STEP_ROTATION','DELTA_ROTATION']

            
        # Add pressure variables
        if self._check_input_dof("PRESSURE"):
            # Add specific variables for the problem (pressure dofs)
            self.dof_variables = self.dof_variables + ['PRESSURE']
            self.dof_reactions = self.dof_reactions + ['PRESSURE_REACTION']
            if not self.settings["stabilization_factor"].IsNull():
                self.main_model_part.ProcessInfo[KratosMultiphysics.STABILIZATION_FACTOR] = self.settings["stabilization_factor"].GetDouble()

        # Add contat variables
        if self._check_input_dof("LAGRANGE_MULTIPLIER"):
            # Add specific variables for the problem (contact dofs)
            self.dof_variables = self.dof_variables + ['LAGRANGE_MULTIPLIER_NORMAL']
            self.dof_reactions = self.dof_reactions + ['LAGRANGE_MULTIPLIER_NORMAL_REACTION']
        
        
        
    def AddVariables(self):

        self.SetVariables()
        
        self.nodal_variables = self.nodal_variables + self.dof_variables + self.dof_reactions 

        self.nodal_variables = [self.nodal_variables[i] for i in range(0,len(self.nodal_variables)) if self.nodal_variables[i] != 'NOT_DEFINED']

        for variable in self.nodal_variables:            
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable))
            #print(" Added variable ", KratosMultiphysics.KratosGlobals.GetVariable(variable),"(",variable,")")
            
        print("::[Mechanical_Solver]:: General Variables ADDED")
                                                              
        
    def AddDofs(self):
        AddDofsProcess = KratosSolid.AddDofsProcess(self.main_model_part, self.dof_variables, self.dof_reactions)
        AddDofsProcess.Execute()
                
        print("::[Mechanical_Solver]:: DOF's ADDED")
        
    def ImportModelPart(self):

        print("::[Mechanical_Solver]:: Importing model part.")
        problem_path = os.getcwd()
        input_filename = self.settings["model_import_settings"]["input_filename"].GetString()
        
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):            
            # Import model part from mdpa file.
            print("   Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa ")
            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)
            # print("   Finished reading model part from mdpa file ")
            # Check and prepare computing model part and import constitutive laws.
            self._execute_after_reading()
            self._set_and_fill_buffer()
            
        elif(self.settings["model_import_settings"]["input_type"].GetString() == "rest"):
            # Import model part from restart file.
            restart_path = os.path.join(problem_path, self.settings["model_import_settings"]["input_filename"].GetString() + "__" + self.settings["model_import_settings"]["input_file_label"].GetString() )
            if(os.path.exists(restart_path+".rest") == False):
                raise Exception("Restart file not found: " + restart_path + ".rest")
            print("   Loading Restart file: ", restart_path + ".rest ")
            # set serializer flag
            serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE      # binary
            # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR # ascii
            # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL   # ascii

            serializer = KratosMultiphysics.Serializer(restart_path, serializer_flag)
            serializer.Load(self.main_model_part.Name, self.main_model_part)

            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
            #I use it to rebuild the contact conditions.
            load_step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] +1;
            self.main_model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step
            # print("   Finished loading model part from restart file ")            

        else:
            raise Exception("Other input options are not yet implemented.")


        print(self.main_model_part)
        print ("::[Mechanical_Solver]:: Finished importing model part.")
            
    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        # Model part writing
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)
    
    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        print("::[Mechanical_Solver]:: -START-")
        
        # The mechanical solver is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        mechanical_solver = self._get_mechanical_solver()
        mechanical_solver.SetEchoLevel(self.settings["echo_level"].GetInt())
        if (self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
            mechanical_solver.Initialize()
        else:
            # SetInitializePerformedFlag is not a member of SolvingStrategy but
            # is used by ResidualBasedNewtonRaphsonStrategy.
            if hasattr(mechanical_solver, SetInitializePerformedFlag):
                mechanical_solver.SetInitializePerformedFlag(True)
        self.Check()
        print("::[Mechanical_Solver]:: -END-")        
        
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
        self._get_mechanical_solver().Solve()
    
    def InitializeSolutionStep(self):
        self._get_mechanical_solver().InitializeSolutionStep()

    def Predict(self):
        self._get_mechanical_solver().Predict()                                            
    def SolveSolutionStep(self):
        is_converged = self._get_mechanical_solver().SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self._get_mechanical_solver().FinalizeSolutionStep()
                                                    
    def SetEchoLevel(self, level):
        self._get_mechanical_solver().SetEchoLevel(level)

    def Clear(self):
        self._get_mechanical_solver().Clear()

    def Check(self):
        self._get_mechanical_solver().Check()
                                                    
    #### Solver internal methods ####

    def _check_input_dof(self, variable):
        dofs_list = self.settings["dofs"]
        for i in range(0, dofs_list.size() ):
            if dofs_list[i].GetString() == variable:
                return True
        return False
            
    def _get_solution_scheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._create_solution_scheme()
        return self._solution_scheme

    def _get_convergence_criterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._create_convergence_criterion()
        return self._convergence_criterion

    def _get_linear_solver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._create_linear_solver()
        return self._linear_solver

    def _get_builder_and_solver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._create_builder_and_solver()
        return self._builder_and_solver

    def _get_mechanical_solver(self):
        if not hasattr(self, '_mechanical_solver'):
            self._mechanical_solver = self._create_mechanical_solver()
        return self._mechanical_solver
                                                    
    def _execute_after_reading(self):
        # The computing_model_part is labeled 'KratosMultiphysics.ACTIVE' flag (in order to recover it)
        self.computing_model_part_name = "computing_domain" 

        # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
       
        if( self.settings.Has("bodies_list") ):
            params.AddValue("bodies_list",self.settings["bodies_list"])

        # CheckAndPrepareModelProcess creates the computating_model_part
        import check_and_prepare_model_process
        check_and_prepare_model_process.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # Import constitutive laws
        materials_imported = self._import_constitutive_laws()
        if materials_imported:
            print("   Constitutive law was successfully imported.")
        else:
            print("   Constitutive law was not imported.")

    def _import_constitutive_laws(self):
        
        if os.path.isfile("materials.py"):
            # Constitutive law import
            import constitutive_law_python_utility as constitutive_law_utils
            constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(self.main_model_part,
                                                                             self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]);
            constitutive_law.Initialize();
            
            return True
        else:
            return False        

    def _validate_and_transfer_matching_settings(self, origin_settings, destination_settings):
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
                            self._validate_and_transfer_matching_settings(orig_value[i], dest_value[i])
                            if len(orig_value[i].items()) != 0:
                                raise Exception('Json settings not found in default settings: ' + orig_value[i].PrettyPrintJsonString())
                        else:
                            raise Exception('Unsupported parameter type.')
                elif dest_value.IsSubParameter() and orig_value.IsSubParameter():
                    self._validate_and_transfer_matching_settings(orig_value, dest_value)
                    if len(orig_value.items()) != 0:
                        raise Exception('Json settings not found in default settings: ' + orig_value.PrettyPrintJsonString())
                else:
                    raise Exception('Unsupported parameter type.')
                origin_settings.RemoveValue(name)
        
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
        
    def _create_solution_scheme(self):
        raise Exception("please implement the Custom Choice of your Scheme (_create_solution_scheme) in your solver")
    
    def _create_convergence_criterion(self):
        # Creation of an auxiliar Kratos parameters object to store the convergence settings
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddEmptyValue("rotation_dofs").SetBool(self._check_input_dof("ROTATION"))
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("component_wise",self.settings["component_wise"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])
        
        # Construction of the class convergence_criterion
        import convergence_criteria_factory_fem_dem
        convergence_criteria_factory_fem_dem = convergence_criteria_factory_fem_dem.convergence_criterion(conv_params)
        
        return convergence_criteria_factory_fem_dem.mechanical_convergence_criterion

    def _create_linear_solver(self):
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver
    
    def _create_builder_and_solver(self):
        linear_solver = self._get_linear_solver()
        if(self.settings["component_wise"].GetBool() == True):
            builder_and_solver = KratosSolid.ComponentWiseBuilderAndSolver(linear_solver)
        else:
            if(self.settings["block_builder"].GetBool() == True):
                # To keep matrix blocks in builder
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
                
        return builder_and_solver

    def _create_mechanical_solver(self):
        raise Exception("please implement the Custom Choice of your Mechanical Solver (_create_mechanical_solver) in your solver")
