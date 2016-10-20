from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PfemFluidDynamicsApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return QuasiIncompressibleFluidSolver(main_model_part, custom_settings)

class QuasiIncompressibleFluidSolver:

    #def __init__(self, model_part, domain_size, periodic=False):
    def __init__(self, main_model_part, custom_settings):
        
        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = Parameters("""
        {  
            "solver_type": "two_step_v_p_solver",
            "model_import_settings"              : {
                "input_type"       : "mdpa",  
                "input_filename"   : "WaveCoarseFluid",
                "input_file_label" : "0"
            },
            "echo_level": 1,
            "buffer_size": 3,
            "predictor_corrector": true,
            "maximum_velocity_iterations": 1,
            "maximum_pressure_iterations": 7,
            "velocity_tolerance": 1e-5,
            "pressure_tolerance": 1e-5,
            "time_order": 2,
            "dynamic_tau": 0.01,
            "compute_reactions": true,
            "divergence_clearance_steps": 0,
            "reform_dofs_at_each_step"           : true,
            "pressure_linear_solver_config":  {
                "solver_type"                    : "AMGCL",
                "max_iteration"                  : 1000,
                "tolerance"                      : 1e-12,
                "scaling"                        : false,
                "smoother_type"                  : "damped_jacobi",
                "krylov_type"                    : "cg"
            },
            "velocity_linear_solver_config": {
                "solver_type"                    : "bicgstab",
                "max_iteration"                  : 10000,
                "tolerance"                      : 1e-12,
                "preconditioner_type"            : "ILU0Preconditioner",
                "scaling"                        : false
            },
            "problem_domain_sub_model_part_list": ["fluid_model_part"],
            "bodies_list": [
                {"body_name":"body1",
                "parts_list":["Part1"]
                },
                {"body_name":"body2",
                "parts_list":["Part2","Part3"]
                }
            ],
            "neighbour_search": true,
            "processes_sub_model_part_list": [""],
            "solution_type": "Dynamic",
            "time_integration_method": "Implicit",
            "scheme_type": "Newmark",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "line_search": false,
            "compute_contact_forces": false,
            "block_builder": false,
            "clear_storage": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "convergence_criterion": "Residual_criteria",
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            }
        } 
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
       
        #construct the linear solver
        import linear_solver_factory
        self.pressure_linear_solver = linear_solver_factory.ConstructSolver(self.settings["pressure_linear_solver_settings"])
        self.velocity_linear_solver = linear_solver_factory.ConstructSolver(self.settings["velocity_linear_solver_settings"])

        self.compute_reactions = self.settings["compute_reactions"].GetBool()

        print("Construction of QuasiIncompressibleFluidSolver finished.")


    def GetMinimumBufferSize(self):
        return 2;

    def Initialize(self):

        print("Initialize two_step_v_p_solver !!")


        compute_model_part = self.GetComputingModelPart()
        
        MoveMeshFlag = True

    # default settings
        self.settings.velocity_tolerance = self.settings["velocity_tolerance"].GetDouble()
        if(hasattr(self.settings, "velocity_tolerance")):
            self.settings.velocity_tolerance = self.settings["velocity_tolerance"].GetDouble()
        if(hasattr(self.settings, "pressure_tolerance")):
            self.settings.pressure_tolerance = self.settings["pressure_tolerance"].GetDouble()
        self.settings.pressure_tolerance = self.settings["pressure_tolerance"].GetDouble()
        if(hasattr(self.settings, "maximum_velocity_iterations")):
            self.settings.maximum_velocity_iterations = self.settings["maximum_velocity_iterations"].GetInt()
        if(hasattr(self.settings, "maximum_pressure_iterations")):
            self.settings.maximum_pressure_iterations = self.settings["maximum_pressure_iterations"].GetInt()
        if(hasattr(self.settings, "time_order")):
            self.settings.time_order = self.settings["time_order"].GetInt()
        if(hasattr(self.settings, "compute_reactions")):
            self.settings.compute_reactions = self.settings["compute_reactions"].GetBool()
        if(hasattr(self.settings, "reform_dofs_at_each_step")):
            self.settings.ReformDofAtEachStep = self.settings["reform_dofs_at_each_step"].GetBool()
        self.settings.ReformDofAtEachStep = self.settings["reform_dofs_at_each_step"].GetBool()
        if(hasattr(self.settings, "predictor_corrector")):
            self.settings.predictor_corrector = self.settings["predictor_corrector"].GetBool()
        if(hasattr(self.settings, "echo_level")):
            self.settings.echo_level = self.settings["echo_level"].GetInt()
        if(hasattr(self.settings, "dynamic_tau")):
            self.settings.dynamic_tau = self.settings["dynamic_tau"].GetDouble()

        # linear solver settings
        import linear_solver_factory
        if(hasattr(self.settings, "pressure_linear_solver_custom_settings")):
            self.settings.pressure_linear_solver = linear_solver_factory.ConstructSolver(
                self.settings.pressure_linear_solver_self.settings)
        if(hasattr(self.settings, "velocity_linear_solver_custom_settings")):
            self.settings.velocity_linear_solver = linear_solver_factory.ConstructSolver(
                self.settings.velocity_linear_solver_self.settings)
        if(hasattr(self.settings, "divergence_cleareance_step")):
            self.settings.divergence_clearance_steps = self.settings.divergence_cleareance_step
 
        """
        if self.periodic == True:
            self.solver_settings = TwoStepVPSettingsPeriodic(self.model_part,
                                                             self.domain_size,
                                                             self.time_order,
                                                             self.use_slip_conditions,
                                                             MoveMeshFlag,
                                                             self.ReformDofAtEachIteration,PATCH_INDEX)
        else:
            self.solver_settings = TwoStepVPSettings(self.model_part,
                                                     self.domain_size,
                                                     self.time_order,
                                                     self.use_slip_conditions,
                                                     MoveMeshFlag,
                                                     self.ReformDofAtEachIteration)
        self.solver_settings.SetEchoLevel(self.echo_level)

        self.solver_settings.SetStrategy(TwoStepVPStrategyLabel.Velocity,
                                         self.velocity_linear_solver,
                                         self.velocity_tolerance,
                                         self.maximum_velocity_iterations)

        self.solver_settings.SetStrategy(TwoStepVPStrategyLabel.Pressure,
                                         self.pressure_linear_solver,
                                         self.pressure_tolerance,
                                         self.maximum_pressure_iterations)
        """


        """
        if self.periodic == True:
            self.solver = TwoStepVPStrategy(self.model_part,
                                            self.solver_settings, 
                                            self.predictor_corrector, 
                                            PATCH_INDEX)
        else:
            self.solver = TwoStepVPStrategy(self.model_part,
                                            self.solver_settings, 
                                            self.predictor_corrector)
        """

        self.solver = TwoStepVPStrategy(self.main_model_part,
                                        self.velocity_linear_solver,
                                        self.pressure_linear_solver,
                                        MoveMeshFlag,
                                        self.ReformDofAtEachStep,
                                        self.velocity_tolerance,
                                        self.pressure_tolerance,
                                        self.maximum_velocity_iterations,
                                        self.maximum_pressure_iterations,
                                        self.time_order,
                                        self.main_model_part.ProcessInfo[DOMAIN_SIZE],
                                        self.predictor_corrector)

        self.solver.Check()

# self.solver.ApplyFractionalVelocityFixity()
        # generating the slip conditions
# self.create_slip_conditions.Execute()
# (self.solver).SetSlipProcess(self.create_slip_conditions);
# self.slip_conditions_initialized = True
# (self.solver).SetEchoLevel(self.echo_level)
        print("finished initialization of the fluid strategy")


    def InitializeStressStrain(self):

        self.solver.InitializeStressStrain()

        print("Initialize Stress Strain finished ")


    def Solve(self):
        #if(self.ReformDofAtEachStep):
            # (self.neighbour_search).Execute()
# self.slip_conditions_initialized = False
            #if self.use_slip_conditions:
               # self.normal_util.CalculateOnSimplex(
                #    self.model_part, self.domain_size, IS_STRUCTURE)
                    
        (self.solver).Solve()

        self.solver.CalculateAccelerations()  # ACCELERATION
        #self.solver.CalculateDisplacements()  # DISPLACEMENTS
        self.solver.CalculateHistoricalVariables()  # STRESS-STRAIN

         #if(self.compute_reactions):
         #    self.solver.CalculateReactions()  # REACTION)

     
    def Clear(self):
        (self.solver).Clear()
        self.slip_conditions_initialized = True


    def WriteRestartFile(self, FileName):
        restart_file = open(FileName + ".mdpa", 'w')
        import new_restart_utilities
        new_restart_utilities.PrintProperties(restart_file)
        new_restart_utilities.PrintNodes(self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintElements(
            "Fluid3D", self.model_part.Elements, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VELOCITY_X, "VELOCITY_X", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VELOCITY_Y, "VELOCITY_Y", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VELOCITY_Z, "VELOCITY_Z", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            ACCELERATION_X, "ACCELERATION_X", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            ACCELERATION_Y, "ACCELERATION_Y", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            ACCELERATION_Z, "ACCELERATION_Z", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            ACCELERATION_X, "DISPLACEMENT_X", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            ACCELERATION_Y, "DISPLACEMENT_Y", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            ACCELERATION_Z, "DISPLACEMENT_Z", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            PRESSURE, "PRESSURE", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VISCOSITY, "VISCOSITY", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            BULK_MODULUS, "BULK_MODULUS", self.model_part.Nodes, restart_file)
        restart_file.close()

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

        
    def ImportModelPart(self):
        
        print("::[PFEM Fluid Mechanics Solver]:: Model reading starts.")

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            
            print("    Importing input model part...")
            
            ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            print("    Imported input model part.")
            
            self.computing_model_part_name = "fluid_computing_domain"
            # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
            params = Parameters("{}")
            params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)
            params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
            params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])         
            params.AddValue("bodies_list",self.settings["bodies_list"])         

            # CheckAndPrepareModelProcess creates the solid_computational model part
            import check_and_prepare_model_process_fluid
            check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()
            
            for el in self.main_model_part.Elements:
                density = el.Properties.GetValue(DENSITY)
                viscosity = el.Properties.GetValue(VISCOSITY)
                bulk_modulus = el.Properties.GetValue(BULK_MODULUS)
                break
            print ("density: ",density)
            print ("viscosity: ",viscosity)
            print ("bulk_modulus: ",bulk_modulus)
            VariableUtils().SetScalarVar(DENSITY, density, self.main_model_part.Nodes)              # Set density
            VariableUtils().SetScalarVar(VISCOSITY, viscosity, self.main_model_part.Nodes)  # Set kinematic viscosity
            VariableUtils().SetScalarVar(BULK_MODULUS, bulk_modulus, self.main_model_part.Nodes)  # Set kinematic viscosity
            self.main_model_part.SetBufferSize( self.settings["buffer_size"].GetInt() )
        
            current_buffer_size = self.main_model_part.GetBufferSize()
            if(self.GetMinimumBufferSize() > current_buffer_size):
                current_buffer_size = self.GetMinimumBufferSize()

            self.main_model_part.SetBufferSize( current_buffer_size )

            #fill buffer
            delta_time = self.main_model_part.ProcessInfo[DELTA_TIME]
            time = self.main_model_part.ProcessInfo[TIME]
            time = time - delta_time * (current_buffer_size)
            self.main_model_part.ProcessInfo.SetValue(TIME, time)            
            for size in range(0, current_buffer_size):
                step = size - (current_buffer_size -1)
                self.main_model_part.ProcessInfo.SetValue(STEP, step)
                time = time + delta_time
                #delta_time is computed from previous time in process_info
                self.main_model_part.CloneTimeStep(time)

           #  self.main_model_part.ProcessInfo[IS_RESTARTED] = False
            

        elif(self.settings["model_import_settings"]["input_type"].GetString() == "rest"):

            problem_path = os.getcwd()
            restart_path = os.path.join(problem_path, self.settings["model_import_settings"]["input_filename"].GetString() + "__" + self.settings["model_import_settings"]["input_file_label"].GetString() )

            if(os.path.exists(restart_path+".rest") == False):
                print("    rest file does not exist , check the restart step selected ")

            # set serializer flag
            self.serializer_flag = SerializerTraceType.SERIALIZER_NO_TRACE      # binary
            # self.serializer_flag = SerializerTraceType.SERIALIZER_TRACE_ERROR # ascii
            # self.serializer_flag = SerializerTraceType.SERIALIZER_TRACE_ALL   # ascii

            serializer = Serializer(restart_path, self.serializer_flag)

            serializer.Load(self.main_model_part.Name, self.main_model_part)
            print("    Load input restart file.")

            self.main_model_part.ProcessInfo[IS_RESTARTED] = True

            print(self.main_model_part)

        else:
            raise Exception("Other input options are not yet implemented.")
        
        
        print ("::[Mechanical Solver]:: Model reading finished.")

#
#
def CreateSolver(model_part, config, custom_settings):
    #fluid_solver = QuasiIncompressibleFluidSolver(model_part, config.domain_size, periodic)
    fluid_solver = QuasiIncompressibleFluidSolver(model_part, custom_settings)

    # default settings
    #fluid_solver.velocity_tolerance = config.velocity_tolerance
    #if(hasattr(config, "velocity_tolerance")):
        #fluid_solver.velocity_tolerance = config.velocity_tolerance
    #if(hasattr(config, "pressure_tolerance")):
        #fluid_solver.pressure_tolerance = config.pressure_tolerance
    #if(hasattr(config, "maximum_velocity_iterations")):
        #fluid_solver.maximum_velocity_iterations = config.maximum_velocity_iterations
    #if(hasattr(config, "maximum_pressure_iterations")):
        #fluid_solver.maximum_pressure_iterations = config.maximum_pressure_iterations
    #if(hasattr(config, "time_order")):
        #fluid_solver.time_order = config.time_order
    #if(hasattr(config, "compute_reactions")):
        #fluid_solver.compute_reactions = config.compute_reactions
    #if(hasattr(config, "ReformDofAtEachIteration")):
        #fluid_solver.ReformDofAtEachIteration = config.ReformDofAtEachIteration
    #if(hasattr(config, "predictor_corrector")):
        #fluid_solver.predictor_corrector = config.predictor_corrector
    #if(hasattr(config, "echo_level")):
        #fluid_solver.echo_level = config.echo_level
    #if(hasattr(config, "dynamic_tau")):
        #fluid_solver.dynamic_tau = config.dynamic_tau

    fluid_solver.velocity_tolerance = config["velocity_tolerance"].GetDouble()
    if(hasattr(config, "velocity_tolerance")):
        fluid_solver.velocity_tolerance = config["velocity_tolerance"].GetDouble()
    if(hasattr(config, "pressure_tolerance")):
        fluid_solver.pressure_tolerance = config["pressure_tolerance"].GetDouble()
    if(hasattr(config, "maximum_velocity_iterations")):
        fluid_solver.maximum_velocity_iterations = config["maximum_velocity_iterations"].GetInt()
    if(hasattr(config, "maximum_pressure_iterations")):
        fluid_solver.maximum_pressure_iterations = config["maximum_pressure_iterations"].GetInt()
    if(hasattr(config, "time_order")):
        fluid_solver.time_order = config["time_order"].GetInt()
    if(hasattr(config, "compute_reactions")):
        fluid_solver.compute_reactions = config["compute_reactions"].GetBool()
    if(hasattr(config, "reform_dofs_at_each_step")):
        fluid_solver.ReformDofAtEachStep = config["reform_dofs_at_each_step"].GetBool()
    if(hasattr(config, "predictor_corrector")):
        fluid_solver.predictor_corrector = config["predictor_corrector"].GetBool()
    if(hasattr(config, "echo_level")):
        fluid_solver.echo_level = config["echo_level"].GetInt()
    if(hasattr(config, "dynamic_tau")):
        fluid_solver.dynamic_tau = config["dynamic_tau"].GetDouble()

    # linear solver settings
    import linear_solver_factory
    if(hasattr(config, "pressure_linear_solver_config")):
        fluid_solver.pressure_linear_solver = linear_solver_factory.ConstructSolver(
            config.pressure_linear_solver_config)
    if(hasattr(config, "velocity_linear_solver_config")):
        fluid_solver.velocity_linear_solver = linear_solver_factory.ConstructSolver(
            config.velocity_linear_solver_config)
    if(hasattr(config, "divergence_cleareance_step")):
        fluid_solver.divergence_clearance_steps = config.divergence_cleareance_step

    return fluid_solver

def CreateSolver(model_part, custom_settings):
    #fluid_solver = QuasiIncompressibleFluidSolver(model_part, config.domain_size, periodic)
    fluid_solver = QuasiIncompressibleFluidSolver(model_part, custom_settings)

    # default settings
    fluid_solver.velocity_tolerance = custom_settings["velocity_tolerance"].GetDouble()
    if(hasattr(custom_settings, "velocity_tolerance")):
        fluid_solver.velocity_tolerance = custom_settings["velocity_tolerance"].GetDouble()
    if(hasattr(custom_settings, "pressure_tolerance")):
        fluid_solver.pressure_tolerance = custom_settings["pressure_tolerance"].GetDouble()
    if(hasattr(custom_settings, "maximum_velocity_iterations")):
        fluid_solver.maximum_velocity_iterations = custom_settings["maximum_velocity_iterations"].GetInt()
    if(hasattr(custom_settings, "maximum_pressure_iterations")):
        fluid_solver.maximum_pressure_iterations = custom_settings["maximum_pressure_iterations"].GetInt()
    if(hasattr(custom_settings, "time_order")):
        fluid_solver.time_order = custom_settings["time_order"].GetInt()
    if(hasattr(custom_settings, "compute_reactions")):
        fluid_solver.compute_reactions = custom_settings["compute_reactions"].GetBool()
    if(hasattr(custom_settings, "reform_dofs_at_each_step")):
        fluid_solver.ReformDofAtEachStep = custom_settings["reform_dofs_at_each_step"].GetBool()
    

    if(hasattr(custom_settings, "predictor_corrector")):
        fluid_solver.predictor_corrector = custom_settings["predictor_corrector"].GetBool()
    if(hasattr(custom_settings, "echo_level")):
        fluid_solver.echo_level = custom_settings["echo_level"].GetInt()
    if(hasattr(custom_settings, "dynamic_tau")):
        fluid_solver.dynamic_tau = custom_settings["dynamic_tau"].GetDouble()

    # linear solver settings
    import linear_solver_factory
    if(hasattr(custom_settings, "pressure_linear_solver_custom_settings")):
        fluid_solver.pressure_linear_solver = linear_solver_factory.ConstructSolver(
            custom_settings.pressure_linear_solver_custom_settings)
    if(hasattr(custom_settings, "velocity_linear_solver_custom_settings")):
        fluid_solver.velocity_linear_solver = linear_solver_factory.ConstructSolver(
            custom_settings.velocity_linear_solver_custom_settings)
    if(hasattr(custom_settings, "divergence_cleareance_step")):
        fluid_solver.divergence_clearance_steps = custom_settings.divergence_cleareance_step


    fluid_solver.pressure_tolerance = custom_settings["pressure_tolerance"].GetDouble()
    fluid_solver.maximum_velocity_iterations = custom_settings["maximum_velocity_iterations"].GetInt()
    fluid_solver.maximum_pressure_iterations = custom_settings["maximum_pressure_iterations"].GetInt()
    fluid_solver.time_order = custom_settings["time_order"].GetInt()
    fluid_solver.compute_reactions = custom_settings["compute_reactions"].GetBool()
    fluid_solver.ReformDofAtEachStep = custom_settings["reform_dofs_at_each_step"].GetBool()
    fluid_solver.predictor_corrector = custom_settings["predictor_corrector"].GetBool()
    fluid_solver.echo_level = custom_settings["echo_level"].GetInt()
    fluid_solver.dynamic_tau = custom_settings["dynamic_tau"].GetDouble()
    fluid_solver.divergence_clearance_steps = custom_settings["divergence_clearance_steps"].GetInt()


    return fluid_solver



def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(FRACT_VEL)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT)
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ)
    model_part.AddNodalSolutionStepVariable(CONV_PROJ)
# model_part.AddNodalSolutionStepVariable(ADVPROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(BULK_MODULUS)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(FREESURFACE)
    model_part.AddNodalSolutionStepVariable(INTERF)
    model_part.AddNodalSolutionStepVariable(Y_WALL)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    # Stokes needs it (in case periodic conditions are required)
    # model_part.AddNodalSolutionStepVariable(PATCH_INDEX)


    print("variables for the two step v p solver added correctly")



def AddDofs(model_part, config=None):

    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(PRESSURE)
        node.AddDof(VELOCITY_X)
        node.AddDof(VELOCITY_Y)
        node.AddDof(VELOCITY_Z)

    print("dofs for the two step v p solver added correctly")



def NodalChecksAndAssignations(model_part):

    numFluid=0
    numRigid=0
    numRigidFluid=0
    numRigidNotFluid=0
    numBoundary=0
    numIsolated=0
    numFreeSurface=0
    numBlocked=0
    mean_nodal_h=0
    for node in model_part.Nodes:
        # adding dofs
        if (node.Is(FLUID)):
            numFluid+=1

            nodal_h=node.GetSolutionStepValue(NODAL_H)
            mean_nodal_h+=nodal_h

            #density=node.GetSolutionStepValue(DENSITY)
            #print("density ",density)

            #viscosity=node.GetSolutionStepValue(VISCOSITY)
            #print("VISCOSITY ",viscosity)

            #bulk_modulus=node.GetSolutionStepValue(BULK_MODULUS)
            #print("bulk_modulus ",bulk_modulus)

        if (node.Is(RIGID)):
            numRigid+=1
            if (node.Is(FLUID)):
                numRigidFluid+=1
            else:
                numRigidNotFluid+=1   
                node.SetSolutionStepValue(PRESSURE,0.0)
        if (node.Is(BOUNDARY)):
            numBoundary+=1
        if (node.Is(ISOLATED)):
            numIsolated+=1
            node.SetSolutionStepValue(PRESSURE,0.0)

        if (node.Is(FREE_SURFACE)):
            numFreeSurface+=1
        if (node.Is(BLOCKED)):
            numBlocked+=1
 
    mean_nodal_h*=1.0/numFluid;
    print("nodal_h is  ",nodal_h)
    print("numFluid ",numFluid)
    print("numRigid ",numRigid)
    print("numRigidFluid ",numRigidFluid)
    print("numRigidNotFluid ",numRigidNotFluid)
    print("numBoundary ",numBoundary)
    print("numIsolated ",numIsolated)
    print("numFreeSurface ",numFreeSurface)
    print("numBlocked ",numBlocked)

