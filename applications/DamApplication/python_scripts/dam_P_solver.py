from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.ExternalSolversApplication as KratosExternal
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
#check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return DamUPSolver(main_model_part, custom_settings)

class DamUPSolver:
    def __init__(self, model_part, custom_settings):

        self.main_model_part = model_part    
        
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
        }""")
        
               
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
               
        amgcl_smoother = KratosMultiphysics.AMGCLSmoother.ILU0
        amgcl_krylov_type = KratosMultiphysics.AMGCLIterativeSolverType.BICGSTAB
        tolerance = 1e-5
        max_iterations = 1000
        verbosity = 0 #0->shows no information, 1->some information, 2->all the information
        gmres_size = 50
        self.linear_solver =  KratosMultiphysics.AMGCLSolver(amgcl_smoother,amgcl_krylov_type,tolerance,max_iterations,verbosity,gmres_size)
 

    def AddVariables(self):
        ## Add pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        ## Add Dynamic pressure Variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.Dt_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.Dt2_PRESSURE)
        
        print("P variables correctly added")
        
        
    def AddDofs(self):
        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.PRESSURE)
            
        print("P DOFs correctly added")
        
    def Initialize(self):
        beta=0.25
        gamma=0.5
        time_scheme = KratosDam.DamPScheme(beta,gamma)
        move_mesh_flag = False #USER SHOULD NOT CHANGE THIS
        
        conv_criteria = KratosMultiphysics.ResidualCriteria(1e-5,1e-9)
        max_iterations = 10
        
        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
            self.main_model_part, 
            time_scheme, 
            self.linear_solver,
            conv_criteria,
            max_iterations,
            False, 
            False, 
            move_mesh_flag)
                                                                                   

        (self.solver).SetEchoLevel(self.settings["mechanical_settings"]["echo_level"].GetInt())
        
        print ("Initialization P Solver finished")
        
    def ImportModelPart(self):
        
        if(self.settings["mechanical_settings"]["model_import_settings"]["input_type"].GetString() == "mdpa"):
            
            # Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["mechanical_settings"]["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            print("    Import input model part.")
            
            self.computing_model_part_name = "acoustic_computing_domain"
            ## Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
            aux_params = KratosMultiphysics.Parameters("{}")
            aux_params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)
            aux_params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
            aux_params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
            
            ## CheckAndPrepareModelProcess creates the solid_computational_model_part
            import check_and_prepare_model_process_poro
            check_and_prepare_model_process_poro.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()
            
            # Constitutive law import
            import dam_constitutive_law_utility
            dam_constitutive_law_utility.SetConstitutiveLaw(self.main_model_part)       
        else:
            raise Exception("other input options are not yet implemented")
        
        self.main_model_part.SetBufferSize( self.settings["mechanical_settings"]["buffer_size"].GetInt() )
        current_buffer_size = self.main_model_part.GetBufferSize()
        
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )
                
        print ("model reading finished")     
   
        
    def GetMinimumBufferSize(self):
        return 2;
        
    def GetComputingModelPart(self):
        return self.main_model_part
                        
    def Solve(self):
        (self.solver).Solve()

    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)
        
    def Clear(self):
        (self.solver).Clear()


