from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.ExternalSolversApplication as KratosExternal
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
#check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return DamMechanicalSolver(main_model_part, custom_settings)


class DamMechanicalSolver:

    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
         {
            "general_data"                        : {
                "problem_name"    : "pruebas",
                "model_part_name" : "MainModelPart",
                "domain_size"     : 2,
                "NumberofThreads" : 1,
                "type_of_problem" : "Thermo-Mechanical",
                "time_scale"      : "Seconds",
                "delta_time"      : 1,
                "ending_time"     : 10
            },
            "diffusion_settings"                  : {
                "variables"             : {
                    "unknown_variable"       : "KratosMultiphysics.TEMPERATURE",
                    "difussion_variable"     : "KratosMultiphysics.CONDUCTIVITY",
                    "specific_heat_variable" : "KratosMultiphysics.SPECIFIC_HEAT",
                    "density_variable"       : "KratosMultiphysics.DENSITY"
                },
                "temporal_scheme"       : "Backward-Euler",
                "reference_temperature" : "Reservoir_Information"
            },
            "mechanical_settings"                 : {
                "solver_type"                     : "dam_new_mechanical_solver",
                "model_import_settings"           : {
                    "input_type"     : "mdpa",
                    "input_filename" : "pruebas"
                },
                "solution_type"                   : "Quasi-Static",
                "analysis_type"                   : "Linear",
                "strategy_type"                   : "Newton-Raphson",
                "scheme_type"                     : "Newmark",
                "convergence_criterion"           : "Residual_criterion",
                "displacement_relative_tolerance" : 0.0001,
                "displacement_absolute_tolerance" : 1e-9,
                "residual_relative_tolerance"     : 0.0001,
                "residual_absolute_tolerance"     : 1e-9,
                "max_iteration"                   : 10,
                "max_radius_factor"               : 5.0,
                "min_radius_factor"               : 0.5,
                "max_iteration"                   : 10,
                "echo_level"                      : 0,
                "buffer_size"                     : 2,
                "compute_reactions"               : true,
                "reform_step_dofs"                : false,
                "move_mesh_flag"                  : true,
                "type_of_builder"                 : "Elimination",
                "type_of_solver"                  : "Iterative",
                "solver_class"                    : "AMGCL"
            },
            "problem_domain_sub_model_part_list"  : [""],
            "problem_domain_body_sub_model_part_list"  : [""],
            "problem_domain_joint_sub_model_part_list" : [""],
            "processes_sub_model_part_list"       : [""],
            "nodal_processes_sub_model_part_list" : [""],
            "load_processes_sub_model_part_list"  : [""],
            "loads_sub_model_part_list": [""],
            "loads_variable_list": [""],
            "output_configuration"                : {
                "result_file_configuration" : {
                    "gidpost_flags"       : {
                        "GiDPostMode"           : "GiD_PostBinary",
                        "WriteDeformedMeshFlag" : "WriteDeformed",
                        "WriteConditionsFlag"   : "WriteConditions",
                        "MultiFileFlag"         : "SingleFile"
                    },
                    "output_frequency"    : 1.0,
                    "nodal_results"       : [""],
                    "gauss_point_results" : [""]
                }
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        # Get the linear solver
        self.linear_solver = self.LinearSolver(self.settings["mechanical_settings"]["type_of_solver"].GetString(),
                                    self.settings["mechanical_settings"]["solver_class"].GetString())
    
        print("Construction of Dam Mechanical Solver finished")
    
    def GetMinimumBufferSize(self):
        return 2;

    def AddVariables(self):
        #add displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        #add reactions for the displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        #add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        #add variables for the solid conditions
        self.main_model_part.AddNodalSolutionStepVariable(KratosSolid.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosSolid.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosSolid.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        #add volume acceleration
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.NODAL_CAUCHY_STRESS_TENSOR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.Vi_POSITIVE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.Viii_POSITIVE)
        
        print("Mechanical variables correctly added")

    def ImportModelPart(self):
                        
        if(self.settings["mechanical_settings"]["model_import_settings"]["input_type"].GetString() == "mdpa"):
            
            # Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["mechanical_settings"]["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            print("    Import input model part.")
            
            self.computing_model_part_name = "solid_computing_domain"
            ## Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
            aux_params = KratosMultiphysics.Parameters("{}")
            aux_params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)
            aux_params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
            aux_params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
            
            ## CheckAndPrepareModelProcess creates the solid_computational_model_part
            import check_and_prepare_model_process_poro
            check_and_prepare_model_process_poro.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()
            
            # Constitutive law import
            import constitutive_law_utility
            constitutive_law_utility.SetConstitutiveLaw(self.main_model_part)
            
        else:
            raise Exception("Other input options are not yet implemented.")
        
        self.main_model_part.SetBufferSize( self.settings["mechanical_settings"]["buffer_size"].GetInt() )
        
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )
                
        print ("Model reading finished.")


    def AddDofs(self):
        for node in self.main_model_part.Nodes:
            ##Solid dofs
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X,KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y,KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z,KratosMultiphysics.REACTION_Z)

        print(" Mechanical DOFs correctly added")


    def Initialize(self):
    
        # Get the computational_model_part 
        compute_model_part = self.GetComputingModelPart()
        
        # Builder and solver creation
        builder_and_solver = self.BuilderAndSolverCreator(self.settings["mechanical_settings"]["type_of_builder"].GetString(), self.linear_solver)
        
        # Solution scheme creation
        scheme = self.SchemeCreator(self.settings["mechanical_settings"]["solution_type"].GetString(), self.settings["mechanical_settings"]["scheme_type"].GetString())
        
        # Convergence criterion
        convergence_criterion = self.ConvergenceCriterion(self.settings["mechanical_settings"]["convergence_criterion"].GetString())

        # Solver creation
        self.solver = self.ConstructTheSolver(self.main_model_part,
                                            scheme,
                                            convergence_criterion,
                                            builder_and_solver,
                                            self.settings["mechanical_settings"]["strategy_type"].GetString()) 
        
        # Set echo_level
        self.solver.SetEchoLevel(self.settings["mechanical_settings"]["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.solver.Check();

        print ("Initialization DamSolver finished")
        
    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart("solid_computing_domain")
    
    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass
        
    def SaveRestart(self):
        pass #one should write the restart file here
        
    def Solve(self):
        
        self.solver.Solve()

    def SetEchoLevel(self, level):
        
        self.solver.SetEchoLevel(level)

    def Clear(self):
        
        self.solver.Clear()
        
    def Check(self):
        
        self.solver.Check()

    #### Specific internal functions ####
    
    def LinearSolver(self, damTypeofSolver, solver_class):
    
        if(damTypeofSolver == "Direct"):
            if(solver_class == "SuperLU"):
                linear_solver = KratosMultiphysics.ExternalSolversApplication.SuperLUSolver()
            elif(solver_class == "SkylineLUFactorization"):
                linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        elif(damTypeofSolver == "Iterative"):
            if(solver_class == "BICGSTAB"):
                tolerance = 1e-5
                max_iterations = 1000
                precond = KratosMultiphysics.ILU0Preconditioner()
                linear_solver = KratosMultiphysics.BICGSTABSolver(tolerance,max_iterations,precond)
            elif(solver_class == "AMGCL"):
                amgcl_smoother = KratosMultiphysics.AMGCLSmoother.ILU0
                amgcl_krylov_type = KratosMultiphysics.AMGCLIterativeSolverType.BICGSTAB
                tolerance = 1e-5
                max_iterations = 1000
                verbosity = 0 #0->shows no information, 1->some information, 2->all the information
                gmres_size = 50
                linear_solver =  KratosMultiphysics.AMGCLSolver(amgcl_smoother,amgcl_krylov_type,tolerance,max_iterations,verbosity,gmres_size)
                  
        return linear_solver
    
    
    def BuilderAndSolverCreator(self, typebuilder, linear_solver):
        
        if(typebuilder == "Elimination"):
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        
        return builder_and_solver
        
    def SchemeCreator(self, solution_type, scheme_type):
        
        smoothing = "Yes"
        
        for i in range(self.settings["output_configuration"]["result_file_configuration"]["nodal_results"].size()):
            outputs = self.settings["output_configuration"]["result_file_configuration"]["nodal_results"][i].GetString()
            if (outputs == "NODAL_CAUCHY_STRESS_TENSOR"):
                smoothing = "Yes"
        

            if (solution_type == "Quasi-Static"):       
                scheme =  KratosDam.IncrementalUpdateStaticSmoothingScheme()
            else:
                if(scheme_type == "Newmark"):
                    damp_factor_m = 0.0
                else:
                    damp_factor_m = -0.01
                    
                scheme = KratosDam.BossakDisplacementSmoothingScheme(damp_factor_m)

        return scheme
         
         
    def ConvergenceCriterion(self, convergence_criterion):
        dis_rel_tol = self.settings["mechanical_settings"]["displacement_relative_tolerance"].GetDouble()
        dis_abs_tol = self.settings["mechanical_settings"]["displacement_absolute_tolerance"].GetDouble()
        res_rel_tol = self.settings["mechanical_settings"]["residual_relative_tolerance"].GetDouble()
        res_abs_tol = self.settings["mechanical_settings"]["residual_absolute_tolerance"].GetDouble()
        echo_level = self.settings["mechanical_settings"]["echo_level"].GetInt()
        
        if(convergence_criterion == "Displacement_criterion"):
            convergence_criterion = KratosSolid.DisplacementConvergenceCriterion(dis_rel_tol, dis_abs_tol)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion == "Residual_criterion"):
            convergence_criterion = KratosMultiphysics.ResidualCriteria(res_rel_tol, res_abs_tol)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion == "And_criterion"):
            Displacement = KratosSolid.DisplacementConvergenceCriterion(dis_rel_tol, dis_abs_tol)
            Displacement.SetEchoLevel(echo_level)
            Residual = KratosMultiphysics.ResidualCriteria(res_rel_tol, res_abs_tol)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
        elif(convergence_criterion == "Or_criterion"):
            Displacement = KratosSolid.DisplacementConvergenceCriterion(dis_rel_tol, dis_abs_tol)
            Displacement.SetEchoLevel(echo_level)
            Residual = KratosMultiphysics.ResidualCriteria(res_rel_tol, res_abs_tol)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
        
        return convergence_criterion
    
    def ConstructTheSolver(self, compute_model_part, scheme, convergence_criterion, builder_and_solver, strategy_type):
    
        max_iters = self.settings["mechanical_settings"]["max_iteration"].GetInt()
        compute_reactions = self.settings["mechanical_settings"]["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["mechanical_settings"]["reform_step_dofs"].GetBool()
        move_mesh_flag = self.settings["mechanical_settings"]["move_mesh_flag"].GetBool()
        
        if(strategy_type == "Newton-Raphson"):
            self.main_model_part.ProcessInfo[KratosPoro.IS_CONVERGED]=True
            self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.main_model_part, 
                                                                            scheme, 
                                                                            self.linear_solver, 
                                                                            convergence_criterion, 
                                                                            builder_and_solver, 
                                                                            max_iters, 
                                                                            compute_reactions, 
                                                                            reform_step_dofs, 
                                                                            move_mesh_flag)
                                                                            
        else:
            
            # Arc-Length strategy
            arc_length_params = KratosMultiphysics.Parameters("{}")
            arc_length_params.AddValue("desired_iterations",self.settings["mechanical_settings"]["desired_iterations"])
            arc_length_params.AddValue("max_radius_factor",self.settings["mechanical_settings"]["max_radius_factor"])
            arc_length_params.AddValue("min_radius_factor",self.settings["mechanical_settings"]["min_radius_factor"])
            arc_length_params.AddValue("loads_sub_model_part_list",self.settings["loads_sub_model_part_list"])
            arc_length_params.AddValue("loads_variable_list",self.settings["loads_variable_list"])
            
            self.solver = KratosPoro.PoromechanicsRammArcLengthStrategy(self.main_model_part,
                                                                   scheme,
                                                                   self.linear_solver,
                                                                   convergence_criterion,
                                                                   builder_and_solver,
                                                                   arc_length_params,
                                                                   max_iters,
                                                                   compute_reactions,
                                                                   reform_step_dofs,
                                                                   move_mesh_flag)
            
        return self.solver
                                        
            
            
