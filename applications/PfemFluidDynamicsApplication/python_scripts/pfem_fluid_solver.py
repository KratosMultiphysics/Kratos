from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return PfemFluidSolver(main_model_part, custom_settings)

class PfemFluidSolver:

    def __init__(self, main_model_part, custom_settings):
        
        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {  
            "echo_level": 1,
            "buffer_size": 3,
            "solver_type": "pfem_fluid_solver",
             "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": false,
            "clear_storage": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "predictor_corrector": true,
            "time_order": 2,
            "maximum_velocity_iterations": 1,
            "maximum_pressure_iterations": 7,
            "velocity_tolerance": 1e-5,
            "pressure_tolerance": 1e-5,
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "AMGCL",
                "max_iteration"                  : 1000,
                "tolerance"                      : 1e-12,
                "provide_coordinates"            : true,
                "scaling"                        : false,
                "smoother_type"                  : "damped_jacobi",
                "krylov_type"                    : "cg",
                "coarsening_type"                : "aggregation",
                "verbosity"                      : 0
            },
            "velocity_linear_solver_settings": {
                "solver_type"                    : "bicgstab",
                "max_iteration"                  : 10000,
                "tolerance"                      : 1e-12,
                "preconditioner_type"            : "ILU0Preconditioner",
                "scaling"                        : false
            },
            "bodies_list": [
                {"body_name":"body1",
                "parts_list":["Part1"]
                },
                {"body_name":"body2",
                "parts_list":["Part2","Part3"]
                }
            ],
            "problem_domain_sub_model_part_list": ["fluid_model_part"],
            "processes_sub_model_part_list": [""]
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

        print("Construction of 2-step Pfem Fluid Solver finished.")


    def GetMinimumBufferSize(self):
        return 2;

    def Initialize(self):

        print("::[Pfem Fluid Solver]:: -START-")

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
                                                              self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Set echo_level
        self.fluid_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

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
        
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
                        
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POISSON_RATIO)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YOUNG_MODULUS)
 
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)


        # PFEM fluid variables
        # self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NORMVELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FREESURFACE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_ACCELERATION)

        print("::[Pfem Fluid Solver]:: Variables ADDED")
                

    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.PRESSURE)
            node.AddDof(KratosMultiphysics.VELOCITY_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)
        print("::[Pfem Fluid Solver]:: DOF's ADDED")
                      
        
    def ImportModelPart(self):
        
        print("::[Pfem Fluid Solver]:: Model reading starts.")

        self.computing_model_part_name = "fluid_computing_domain"
        
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            
            print("    Importing input model part...")
            
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            print("    Import input model part.")
            
            
            # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
            params = KratosMultiphysics.Parameters("{}")
            params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)
            params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
            params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
            if( self.settings.Has("bodies_list") ):
                params.AddValue("bodies_list",self.settings["bodies_list"])         

            # CheckAndPrepareModelProcess creates the fluid_computational model part
            import pfem_check_and_prepare_model_process_fluid
            pfem_check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

            # Set Properties to nodes : Deprecated
            #self.SetProperties()

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
            

        elif(self.settings["model_import_settings"]["input_type"].GetString() == "rest"):

            problem_path = os.getcwd()
            restart_path = os.path.join(problem_path, self.settings["model_import_settings"]["input_filename"].GetString() + "__" + self.settings["model_import_settings"]["input_file_label"].GetString() )

            if(os.path.exists(restart_path+".rest") == False):
                print("    rest file does not exist , check the restart step selected ")

            print("    Load Restart file: ", self.settings["model_import_settings"]["input_filename"].GetString() + "__" + self.settings["model_import_settings"]["input_file_label"].GetString())
            # set serializer flag
            self.serializer_flag = SerializerTraceType.SERIALIZER_NO_TRACE      # binary
            # self.serializer_flag = SerializerTraceType.SERIALIZER_TRACE_ERROR # ascii
            # self.serializer_flag = SerializerTraceType.SERIALIZER_TRACE_ALL   # ascii

            serializer = Serializer(restart_path, self.serializer_flag)

            serializer.Load(self.main_model_part.Name, self.main_model_part)
            print("    Load input restart file.")

            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True

            print(self.main_model_part)

        else:
            raise Exception("Other input options are not yet implemented.")
        
        
        print ("::[Pfem Fluid Solver]:: Model reading finished.")

        
    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()
            
        self.fluid_solver.Solve()
        
        #self.fluid_solver.CalculateAccelerations()  # ACCELERATION
        #self.fluid_solver.CalculateDisplacements()  # DISPLACEMENTS

    # solve :: sequencial calls
    
    def InitializeStrategy(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        #self.fluid_solver.Initialize()

    def InitializeSolutionStep(self):
        #self.fluid_solver.InitializeSolutionStep()

        adaptive_time_interval = KratosPfemFluid.AdaptiveTimeIntervalProcess(self.main_model_part,self.settings["echo_level"].GetInt())
        adaptive_time_interval.Execute()

        unactive_peak_elements = True
        unactive_sliver_elements = True
        set_active_flag = KratosPfemFluid.SetActiveFlagProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.settings["echo_level"].GetInt())
        set_active_flag.Execute()

        #split_elements = KratosPfemFluid.SplitElementsProcess(self.main_model_part,self.settings["echo_level"].GetInt())
        #split_elements.ExecuteInitialize()

    def Predict(self):
        pass
        #self.fluid_solver.Predict()

    def SolveSolutionStep(self):
        #self.fluid_solver.SolveSolutionStep()
        self.fluid_solver.Solve()

    def FinalizeSolutionStep(self):
        #pass
        self.fluid_solver.FinalizeSolutionStep()  

        unactive_peak_elements = True
        unactive_sliver_elements = True
        set_active_flag = KratosPfemFluid.SetActiveFlagProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.settings["echo_level"].GetInt())
        set_active_flag.ExecuteFinalize()

        #split_elements = KratosPfemFluid.SplitElementsProcess(self.main_model_part,self.settings["echo_level"].GetInt())
        #split_elements.ExecuteFinalize()
        

    # solve :: sequencial calls

    
    def SetEchoLevel(self, level):
        self.fluid_solver.SetEchoLevel(level)

    def Clear(self):
        self.fluid_solver.Clear()
        
    def Check(self):
        self.fluid_solver.Check()
#

#   Extra methods:: custom AFranci...
#
    def SetProperties(self):
        for el in self.main_model_part.Elements:
            density = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
            bulk_modulus = el.Properties.GetValue(KratosMultiphysics.BULK_MODULUS)
            young_modulus = el.Properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)
            poisson_ratio = el.Properties.GetValue(KratosMultiphysics.POISSON_RATIO)
            break
            
        print ("density: ",density)
        print ("viscosity: ",viscosity)
        print ("bulk_modulus: ",bulk_modulus)
        print ("young_modulus: ",young_modulus)
        print ("poisson_ratio: ",poisson_ratio)

#


#


