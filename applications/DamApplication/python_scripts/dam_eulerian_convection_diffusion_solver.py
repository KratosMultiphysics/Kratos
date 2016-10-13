from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.DamApplication as KratosDam
#check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return DamThermalSolver(main_model_part, custom_settings)


class DamThermalSolver:

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

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
            
        print("Construction of Dam Thermal Solver finished")
    
    def GetMinimumBufferSize(self):
        return 2;

    def AddVariables(self, py_settings=None):
        #if python(string) settings are given, we copy all the settings to the model part.
        #otherwise, we assume the user has already defined the settings in the model part.
        
        # TODO: This information must be obtained from json file.
        py_settings = type("settings", (object,), {
            "unknown_variable"       : "KratosMultiphysics.TEMPERATURE",
            "diffusion_variable"     : "KratosMultiphysics.CONDUCTIVITY",
            "specific_heat_variable" : "KratosMultiphysics.SPECIFIC_HEAT",
            "density_variable"       : "KratosMultiphysics.DENSITY"
        })
                
        if py_settings is not None:
            
            #we must copy the settings to a c++ object, from now on we will only use the model part, input (py) settings will be ignored in the adddofs, initialize, etc
            thermal_settings = KratosMultiphysics.ConvectionDiffusionSettings() #c++ settings            
            if hasattr(py_settings, "velocity_variable"):
                thermal_settings.SetVelocityVariable(eval(py_settings.velocity_variable))
            if hasattr(py_settings, "mesh_velocity_variable"):
                thermal_settings.SetMeshVelocityVariable(eval(py_settings.mesh_velocity_variable))
            if hasattr(py_settings, "diffusion_variable"):                
                thermal_settings.SetDiffusionVariable(eval(py_settings.diffusion_variable))
            if hasattr(py_settings, "unknown_variable"):
                thermal_settings.SetUnknownVariable(eval(py_settings.unknown_variable))
            if hasattr(py_settings, "specific_heat_variable"):
                thermal_settings.SetSpecificHeatVariable(eval(py_settings.specific_heat_variable))
            if hasattr(py_settings, "density_variable"):
                thermal_settings.SetDensityVariable(eval(py_settings.density_variable))
            #and now we save it in the model part.            
            (self.main_model_part.ProcessInfo).SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, thermal_settings)
            
        #now we can add the variable, as the name suggest
        if (self.main_model_part.ProcessInfo).Has(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS):
            thermal_settings  = (self.main_model_part.ProcessInfo).GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
            if thermal_settings.IsDefinedVelocityVariable():
                self.main_model_part.AddNodalSolutionStepVariable(thermal_settings.GetVelocityVariable())
            if thermal_settings.IsDefinedMeshVelocityVariable():
                self.main_model_part.AddNodalSolutionStepVariable(thermal_settings.GetMeshVelocityVariable())
            if thermal_settings.IsDefinedDiffusionVariable():
                self.main_model_part.AddNodalSolutionStepVariable(thermal_settings.GetDiffusionVariable())
            if thermal_settings.IsDefinedUnknownVariable():
                self.main_model_part.AddNodalSolutionStepVariable(thermal_settings.GetUnknownVariable())
            if thermal_settings.IsDefinedSpecificHeatVariable():
                self.main_model_part.AddNodalSolutionStepVariable(thermal_settings.GetSpecificHeatVariable())
            if thermal_settings.IsDefinedDensityVariable():
                self.main_model_part.AddNodalSolutionStepVariable(thermal_settings.GetDensityVariable())
        else:
            raise ValueError('CONVECTION_DIFFUSION_SETTINGS not defined in the model part!')
            
        print(" Thermal variables correctly added")

    def AddDofs(self, settings=None):
                
        thermal_settings  = (self.main_model_part.ProcessInfo).GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)
        for node in self.main_model_part.Nodes:
            node.AddDof(KratosMultiphysics.TEMPERATURE);
        print("user should include the variables as needed by nodes")
        
        print(" Thermal DOFs correctly added")


    def Initialize(self):
    
        # Solver creation
        self.solver = self.ConstructTheSolver() 
        
        # Set echo_level
        self.solver.SetEchoLevel(self.settings["mechanical_settings"]["echo_level"].GetInt())

        print ("Initialization DamThermalSolver finished")
        
    def Solve(self):
        
        self.solver.Solve()

    def SetEchoLevel(self, level):
        
        self.solver.SetEchoLevel(level)


    #### Specific internal functions ####    
    
    def ConstructTheSolver(self):
        
        domain_size = self.settings["general_data"]["domain_size"].GetInt()
        # assignation of parameters to be used
        self.ReformDofAtEachIteration = False;
        
        # definition of the solvers
        pDiagPrecond = KratosMultiphysics.DiagonalPreconditioner()
        self.linear_solver = KratosMultiphysics.BICGSTABSolver(1e-9, 5000, pDiagPrecond)
        
        parameters = self.settings        
           
        self.solver = KratosDam.DamEulerianConvectionDiffusionStrategy(self.main_model_part,
                                                                                self.linear_solver,
                                                                                parameters,
                                                                                self.ReformDofAtEachIteration,
                                                                                domain_size)
        return self.solver

