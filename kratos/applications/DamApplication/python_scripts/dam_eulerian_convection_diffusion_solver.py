from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.DamApplication
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return DamThermalSolver(main_model_part, custom_settings)

class DamThermalSolver:
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
                "evolution_type"  : "Exact",
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
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solvers
        pDiagPrecond = KratosMultiphysics.DiagonalPreconditioner()
        self.linear_solver = KratosMultiphysics.BICGSTABSolver(1e-9, 5000, pDiagPrecond)

    def AddVariables(self):
        
        # We must copy the settings to a C++ object.  
        thermal_settings = KratosMultiphysics.ConvectionDiffusionSettings()
        thermal_settings.SetDiffusionVariable(KratosMultiphysics.CONDUCTIVITY)
        thermal_settings.SetUnknownVariable(KratosMultiphysics.TEMPERATURE)
        thermal_settings.SetSpecificHeatVariable(KratosMultiphysics.SPECIFIC_HEAT)
        thermal_settings.SetDensityVariable(KratosMultiphysics.DENSITY)
        # We save it in the model part
        (self.main_model_part.ProcessInfo).SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, thermal_settings)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
    

    def AddDofs(self):
        thermal_settings = (self.main_model_part.ProcessInfo).GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.TEMPERATURE)
        
        print(" Thermal DOFs correctly added")
    
        
    def Initialize(self):
        self.ReformDofAtEachIteration = False
        domain_size = self.settings["general_data"]["domain_size"].GetInt()
        
        parameters = self.settings
        
        self.solver = KratosMultiphysics.DamApplication.DamEulerianConvectionDiffusionStrategy(self.main_model_part,
                                                                                self.linear_solver,
                                                                                parameters,
                                                                                self.ReformDofAtEachIteration,
                                                                                domain_size)
                                                                                       
        #~ self.solver = KratosMultiphysics.ConvectionDiffusionApplication.ResidualBasedEulerianConvectionDiffusionStrategy(self.main_model_part,
                                                                                                                        #~ self.linear_solver,
                                                                                                                        #~ self.ReformDofAtEachIteration,
                                                                                                                        #~ domain_size)
                                                                                        
                                                                                        

        (self.solver).SetEchoLevel(self.settings["mechanical_settings"]["echo_level"].GetInt())
        
        print ("Initialization Dam Thermal Solver finished")
        
    def GetMinimumBufferSize(self):
        return 2;
                        
    def Solve(self):
        (self.solver).Solve()

    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

