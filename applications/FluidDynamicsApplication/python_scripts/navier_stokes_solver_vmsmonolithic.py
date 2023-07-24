# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

import math

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

class StabilizedFormulation(object):
    """Helper class to define stabilization-dependent parameters."""
    def __init__(self,settings):
        self.element_name = None
        self.element_integrates_in_time = False
        self.element_has_nodal_properties = False
        self.historical_nodal_properties_variables_list = []
        self.non_historical_nodal_properties_variables_list = []
        self.process_data = {}

        #TODO: Keep this until the MonolithicWallCondition is removed to ensure backwards compatibility in solvers with no defined condition_name
        self.condition_name = "MonolithicWallCondition"

        if settings.Has("element_type"):
            formulation = settings["element_type"].GetString()
            if formulation == "vms":
                self._SetUpClassicVMS(settings)
            elif formulation == "qsvms":
                self._SetUpQSVMS(settings)
            elif formulation == "dvms":
                self._SetUpDVMS(settings)
            elif formulation == "fic":
                self._SetUpFIC(settings)
            elif formulation == "symbolic":
                warn_msg  = 'Provided \'element_name\' is \'symbolic\'. This is has been renamed to \'weakly_compressible\'. Use this instead.'
                KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)
                self._SetUpWeaklyCompressible(settings)
            elif formulation == "weakly_compressible":
                self._SetUpWeaklyCompressible(settings)
        else:
            print(settings)
            raise RuntimeError("Argument \'element_type\' not found in stabilization settings.")

    def SetProcessInfo(self,model_part):
        for variable,value in self.process_data.items():
            model_part.ProcessInfo[variable] = value

    def _SetUpClassicVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "vms",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.01
        }""")

        default_non_newtonian_settings = KratosMultiphysics.Parameters(r"""{
            "power_law_k": 1e-6,
            "power_law_n": 1.0,
            "yield_stress": 0.0,
            "regularization_coefficient" : 100.0
        }""")

        # if non-newtonian, there are some extra options
        if settings.Has("non_newtonian_fluid_parameters"):
            self.non_newtonian_option = True
            default_settings.AddValue("non_newtonian_fluid_parameters", default_non_newtonian_settings)
            self.element_name = 'HerschelBulkleyVMS'
        else:
            self.non_newtonian_option = False
            self.element_name = 'VMS'
        self.condition_name = "MonolithicWallCondition"

        settings.ValidateAndAssignDefaults(default_settings)

        # set the nodal material properties flag
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = [KratosMultiphysics.DENSITY, KratosMultiphysics.VISCOSITY]

        # validate the non-newtonian parameters if necessary
        if self.non_newtonian_option:
            settings["non_newtonian_fluid_parameters"].ValidateAndAssignDefaults(default_non_newtonian_settings)
            self.process_data[KratosMultiphysics.POWER_LAW_K] = settings["non_newtonian_fluid_parameters"]["power_law_k"].GetDouble()
            self.process_data[KratosMultiphysics.POWER_LAW_N] = settings["non_newtonian_fluid_parameters"]["power_law_n"].GetDouble()
            self.process_data[KratosMultiphysics.YIELD_STRESS] = settings["non_newtonian_fluid_parameters"]["yield_stress"].GetDouble()
            self.process_data[KratosCFD.REGULARIZATION_COEFFICIENT] = settings["non_newtonian_fluid_parameters"]["regularization_coefficient"].GetDouble()

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpQSVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "qsvms",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.0,
            "element_manages_time_integration": false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        if settings["element_manages_time_integration"].GetBool() == False:
            self.element_name = "QSVMS"
            self.element_integrates_in_time = False
        else:
            self.element_name = "TimeIntegratedQSVMS"
            self.element_integrates_in_time = True
        self.condition_name = "NavierStokesWallCondition"

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpDVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "dvms",
            "dynamic_tau": 0.0,
            "use_orthogonal_subscales": false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "DVMS"
        self.condition_name = "NavierStokesWallCondition"

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpFIC(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "fic",
            "beta": 0.8,
            "adjust_beta_dynamically": false,
            "dynamic_tau": 0.0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        dynamic_beta = settings["adjust_beta_dynamically"].GetBool()
        if dynamic_beta:
            KratosMultiphysics.Logger.PrintWarning("NavierStokesSolverVMSMonolithic","FIC with dynamic beta not yet implemented, using provided beta as a constant value")
        else:
            self.element_name = "FIC"
        self.condition_name = "NavierStokesWallCondition"

        self.process_data[KratosCFD.FIC_BETA] = settings["beta"].GetDouble()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = 0
        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()

    def _SetUpWeaklyCompressible(self,settings):
        #TODO: Remove this after deprecation period is over
        if (settings["element_type"].GetString() == "symbolic"):
            settings["element_type"].SetString("weakly_compressible")

        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "weakly_compressible",
            "dynamic_tau": 1.0,
            "sound_velocity": 1.0e+12
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "WeaklyCompressibleNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.element_integrates_in_time = True

        # set the nodal material properties flag
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = [KratosMultiphysics.DENSITY]
        self.non_historical_nodal_properties_variables_list = [KratosMultiphysics.SOUND_VELOCITY]

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        #TODO: Remove SOUND_VELOCITY from ProcessInfo. Should be obtained from the properties.
        self.process_data[KratosMultiphysics.SOUND_VELOCITY] = settings["sound_velocity"].GetDouble()

def CreateSolver(model, custom_settings):
    return NavierStokesSolverMonolithic(model, custom_settings)

class NavierStokesSolverMonolithic(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "monolithic",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "formulation": {
                "element_type": "vms"
            },
            "maximum_iterations": 10,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": true,
            "assign_neighbour_elements_to_conditions": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"        : {
                "solver_type" : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01,
                "time_step"           : 0.0
            },
            "time_scheme":"bossak",
            "alpha":-0.3,
            "velocity_relaxation":0.9,
            "pressure_relaxation":0.9,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "move_mesh_flag": false
        }""")

        default_settings.AddMissingParameters(super(NavierStokesSolverMonolithic, cls).GetDefaultParameters())
        return default_settings

    def _BackwardsCompatibilityHelper(self,settings):
        ## Backwards compatibility -- deprecation warnings
        if settings.Has("stabilization"):
            msg  = "Input JSON data contains deprecated setting \'stabilization\'.\n"
            msg += "Please rename it to \'formulation\' (and rename \'stabilization/formulation\' to \'formulation/element_type\' if it exists).\n"
            KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            settings.AddValue("formulation", settings["stabilization"])
            settings.RemoveValue("stabilization")
            settings["formulation"].AddValue("element_type", settings["formulation"]["formulation"])
            settings["formulation"].RemoveValue("formulation")

        if settings.Has("oss_switch"):
            msg  = "Input JSON data contains deprecated setting \'oss_switch\' (int).\n"
            msg += "Please define \'formulation/element_type\' (set it to \'vms\')\n"
            msg += "and set \'formulation/use_orthogonal_subscales\' (bool) instead."
            KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            if not settings.Has("formulation"):
                settings.AddValue("formulation",KratosMultiphysics.Parameters(r'{"element_type":"vms"}'))
            settings["formulation"].AddEmptyValue("use_orthogonal_subscales")
            settings["formulation"]["use_orthogonal_subscales"].SetBool(bool(settings["oss_switch"].GetInt()))
            settings.RemoveValue("oss_switch")

        if settings.Has("dynamic_tau"):
            msg  = "Input JSON data contains deprecated setting \'dynamic_tau\' (float).\n"
            msg += "Please define \'formulation/element_type\' (set it to \'vms\') and \n"
            msg += "set \'formulation/dynamic_tau\' (float) instead."
            KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            if not settings.Has("formulation"):
                settings.AddValue("formulation",KratosMultiphysics.Parameters(r'{"element_type":"vms"}'))
            settings["formulation"].AddEmptyValue("dynamic_tau")
            settings["formulation"]["dynamic_tau"].SetDouble(settings["dynamic_tau"].GetDouble())
            settings.RemoveValue("dynamic_tau")

        if settings.Has("turbulence_model") and settings["turbulence_model"].IsString():
            if settings["turbulence_model"].GetString().lower()!="none":
                msg = "Ignoring deprecated \"turbulence_model\" (string) setting."
                KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            settings.RemoveValue("turbulence_model")

        return settings


    def __init__(self, model, custom_settings):
        custom_settings = self._BackwardsCompatibilityHelper(custom_settings)
        super(NavierStokesSolverMonolithic,self).__init__(model,custom_settings)

        # Set up the auxiliary class with the formulation settings
        self._SetFormulation()

        # Update the default buffer size according to the selected time scheme
        self._SetTimeSchemeBufferSize()

        ## Get time-step size
        self.timestep=self.settings["time_stepping"]["time_step"].GetDouble()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesSolverMonolithic finished.")

    def AddVariables(self):
        ## Add base class variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)

        # Adding variables required for the nodal material properties
        if self.element_has_nodal_properties:
            for variable in self.historical_nodal_properties_variables_list:
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        # Adding variables required for the periodic conditions
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def Initialize(self):
        # If the solver requires an instance of the stabilized formulation class, set the process info variables
        if hasattr(self, 'formulation'):
            self.formulation.SetProcessInfo(self.GetComputingModelPart())

        # Construct and initialize the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def CalculateNodalArea(self):    
        nodal_area_process = KratosMultiphysics.CalculateNodalAreaProcess(self.main_model_part,2)
        nodal_area_process.Execute()

    def CalculateTheError(self):
         
        time=self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]  
        print("dt=", self.timestep)
        #if (time==self.timestep or time>4.998):  
        if (True):
        #if (time>1.27888999):       
         print("CalculateTheError")  
         self.CalculateNodalArea()   
         #mu=0.001
         mu=1.0
         rho=1.0
         totalerrorvelocitywithnodalarea=0.0
         totalerrorpressurewithnodalarea=0.0
         totalerrorvelocity=0.0
         totalerrorpressure=0.0
         nodalareasum=0.0
         totalarea=0.0 
         for node in self.main_model_part.Nodes:
          vel_x=-1.0*math.sin(node.X)*math.cos(node.Y)*math.exp(-2.0*(mu/rho)*time)  
          vel_y=math.cos(node.X)*math.sin(node.Y)*math.exp(-2.0*(mu/rho)*time)
          vel_z=0.0
          pressure=(rho/4.0)*(math.cos(2.0*node.X)+math.cos(2.0*node.Y))*math.exp(-4.0*(mu/rho)*time) 
          numerical_vel_x=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
          numerical_vel_y=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
          numerical_pressure=node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
          nodalarea=node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
          nodalareasum+=nodalarea
          errorvel_x = vel_x-numerical_vel_x
          errorvel_y = vel_y-numerical_vel_y
          errorpressure = pressure-numerical_pressure
          totalerrorvelocitywithnodalarea +=(errorvel_x*errorvel_x+errorvel_y*errorvel_y)*nodalarea
          totalerrorpressurewithnodalarea +=(errorpressure*errorpressure)*nodalarea
          totalerrorvelocity +=(errorvel_x*errorvel_x+errorvel_y*errorvel_y)
          totalerrorpressure +=(errorpressure*errorpressure)
         totalerrorvelocitywithnodalarea=totalerrorvelocitywithnodalarea/nodalareasum
         totalerrorpressurewithnodalarea=totalerrorpressurewithnodalarea/nodalareasum
         totalerrorvelocitywithnodalarea  = math.sqrt(totalerrorvelocitywithnodalarea)  
         totalerrorpressurewithnodalarea  = math.sqrt(totalerrorpressurewithnodalarea)
         totalerrorvelocity  = math.sqrt(totalerrorvelocity)  
         totalerrorpressure  = math.sqrt(totalerrorpressure)
         #error calculation with 7 gauss points https://kratos-wiki.cimne.upc.edu/index.php/Numerical_Integration
         totalerror=0.0
         error_vel=0.0
         error_pressure=0.0
         totalarea=0.0
         for element in self.main_model_part.Elements:
          x10 = element.GetNode(1).X - element.GetNode(0).X
          y10 = element.GetNode(1).Y - element.GetNode(0).Y
          x20 = element.GetNode(2).X - element.GetNode(0).X
          y20 = element.GetNode(2).Y - element.GetNode(0).Y
          detJ = x10 * y20-y10 * x20
          Area = 0.5*detJ
          totalarea=totalarea+Area
          #node_a
          xa_g=0.0*(element.GetNode(0).X)+0.0*(element.GetNode(1).X)+1.0*(element.GetNode(2).X)
          ya_g=0.0*(element.GetNode(0).Y)+0.0*(element.GetNode(1).Y)+1.0*(element.GetNode(2).Y)
          weight_a=1.0/40.0
          vela_x=-1.0*math.sin(xa_g)*math.cos(ya_g)*math.exp(-2.0*(mu/rho)*time)  
          vela_y=math.cos(xa_g)*math.sin(ya_g)*math.exp(-2.0*(mu/rho)*time)
          pressure_a=(rho/4.0)*(math.cos(2.0*xa_g)+math.cos(2.0*ya_g))*math.exp(-4.0*(mu/rho)*time) 
          numerical_vela_x=0.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+1.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
          numerical_vela_y=0.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+1.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
          numerical_pressure_a=0.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+1.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.PRESSURE))
          error_vela_x=vela_x-numerical_vela_x
          error_vela_y=vela_y-numerical_vela_y
          error_pressurea=pressure_a-numerical_pressure_a
          error_vel=error_vel+(error_vela_x*error_vela_x+error_vela_y*error_vela_y)*weight_a*detJ
          error_pressure=error_pressure+(error_pressurea*error_pressurea)*weight_a*detJ
          #node_b
          xb_g=0.5*(element.GetNode(0).X)+0.0*(element.GetNode(1).X)+0.5*(element.GetNode(2).X)
          yb_g=0.5*(element.GetNode(0).Y)+0.0*(element.GetNode(1).Y)+0.5*(element.GetNode(2).Y)
          weight_b=1.0/15.0
          velb_x=-1.0*math.sin(xb_g)*math.cos(yb_g)*math.exp(-2.0*(mu/rho)*time)  
          velb_y=math.cos(xb_g)*math.sin(yb_g)*math.exp(-2.0*(mu/rho)*time)
          pressure_b=(rho/4.0)*(math.cos(2.0*xb_g)+math.cos(2.0*yb_g))*math.exp(-4.0*(mu/rho)*time) 
          numerical_velb_x=0.5*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.5*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
          numerical_velb_y=0.5*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.5*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
          numerical_pressure_b=0.5*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.5*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.PRESSURE))
          errorb_x=velb_x-numerical_velb_x
          errorb_y=velb_y-numerical_velb_y
          error_pressureb=pressure_b-numerical_pressure_b
          error_vel=error_vel+(errorb_x*errorb_x+errorb_y*errorb_y)*weight_b*detJ
          error_pressure=error_pressure+(error_pressureb*error_pressureb)*weight_b*detJ
          #node_c
          xc_g=1.0*(element.GetNode(0).X)+0.0*(element.GetNode(1).X)+0.0*(element.GetNode(2).X)
          yc_g=1.0*(element.GetNode(0).Y)+0.0*(element.GetNode(1).Y)+0.0*(element.GetNode(2).Y)
          weight_c=1.0/40.0
          velc_x=-1.0*math.sin(xc_g)*math.cos(yc_g)*math.exp(-2.0*(mu/rho)*time)  
          velc_y=math.cos(xc_g)*math.sin(yc_g)*math.exp(-2.0*(mu/rho)*time)
          pressure_c=(rho/4.0)*(math.cos(2.0*xc_g)+math.cos(2.0*yc_g))*math.exp(-4.0*(mu/rho)*time) 
          numerical_velc_x=1.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
          numerical_velc_y=1.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
          numerical_pressure_c=1.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.PRESSURE))
          errorc_x=velc_x-numerical_velc_x
          errorc_y=velc_y-numerical_velc_y
          error_pressurec=pressure_c-numerical_pressure_c
          error_vel=error_vel+(errorc_x*errorc_x+errorc_y*errorc_y)*weight_c*detJ
          error_pressure=error_pressure+(error_pressurec*error_pressurec)*weight_c*detJ
          #node_d
          xd_g=0.5*(element.GetNode(0).X)+0.5*(element.GetNode(1).X)+0.0*(element.GetNode(2).X)
          yd_g=0.5*(element.GetNode(0).Y)+0.5*(element.GetNode(1).Y)+0.0*(element.GetNode(2).Y)
          weight_d=1.0/15.0
          veld_x=-1.0*math.sin(xd_g)*math.cos(yd_g)*math.exp(-2.0*(mu/rho)*time)  
          veld_y=math.cos(xd_g)*math.sin(yd_g)*math.exp(-2.0*(mu/rho)*time)
          pressure_d=(rho/4.0)*(math.cos(2.0*xd_g)+math.cos(2.0*yd_g))*math.exp(-4.0*(mu/rho)*time) 
          numerical_veld_x=0.5*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.5*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
          numerical_veld_y=0.5*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.5*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
          numerical_pressure_d=0.5*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.5*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.PRESSURE))
          errord_x=veld_x-numerical_veld_x
          errord_y=veld_y-numerical_veld_y
          error_pressured=pressure_d-numerical_pressure_d
          error_vel=error_vel+(errord_x*errord_x+errord_y*errord_y)*weight_d*detJ
          error_pressure=error_pressure+(error_pressured*error_pressured)*weight_d*detJ
          #node_e
          xe_g=0.0*(element.GetNode(0).X)+1.0*(element.GetNode(1).X)+0.0*(element.GetNode(2).X)
          ye_g=0.0*(element.GetNode(0).Y)+1.0*(element.GetNode(1).Y)+0.0*(element.GetNode(2).Y)
          weight_e=1.0/40.0
          vele_x=-1.0*math.sin(xe_g)*math.cos(ye_g)*math.exp(-2.0*(mu/rho)*time)  
          vele_y=math.cos(xe_g)*math.sin(ye_g)*math.exp(-2.0*(mu/rho)*time)
          pressure_e=(rho/4.0)*(math.cos(2.0*xe_g)+math.cos(2.0*ye_g))*math.exp(-4.0*(mu/rho)*time) 
          numerical_vele_x=0.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+1.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
          numerical_vele_y=0.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+1.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
          numerical_pressure_e=0.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+1.0*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.0*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.PRESSURE))
          errore_x=vele_x-numerical_vele_x
          errore_y=vele_y-numerical_vele_y
          error_pressuree=pressure_e-numerical_pressure_e
          error_vel=error_vel+(errore_x*errore_x+errore_y*errore_y)*weight_e*detJ
          error_pressure=error_pressure+(error_pressuree*error_pressuree)*weight_e*detJ
          #node_f
          xf_g=0.0*(element.GetNode(0).X)+0.5*(element.GetNode(1).X)+0.5*(element.GetNode(2).X)
          yf_g=0.0*(element.GetNode(0).Y)+0.5*(element.GetNode(1).Y)+0.5*(element.GetNode(2).Y)
          weight_f=1.0/15.0
          velf_x=-1.0*math.sin(xf_g)*math.cos(yf_g)*math.exp(-2.0*(mu/rho)*time)  
          velf_y=math.cos(xf_g)*math.sin(yf_g)*math.exp(-2.0*(mu/rho)*time)
          pressure_f=(rho/4.0)*(math.cos(2.0*xf_g)+math.cos(2.0*yf_g))*math.exp(-4.0*(mu/rho)*time) 
          numerical_velf_x=0.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.5*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.5*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
          numerical_velf_y=0.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.5*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.5*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
          numerical_pressure_f=0.0*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.5*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.5*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.PRESSURE))
          errorf_x=velf_x-numerical_velf_x
          errorf_y=velf_y-numerical_velf_y
          error_pressuref=pressure_f-numerical_pressure_f
          error_vel=error_vel+(errorf_x*errorf_x+errorf_y*errorf_y)*weight_f*detJ
          error_pressure=error_pressure+(error_pressuref*error_pressuref)*weight_f*detJ
          #node_g
          xg_g=0.333333333333333*(element.GetNode(0).X)+0.333333333333333*(element.GetNode(1).X)+0.333333333333333*(element.GetNode(2).X)
          yg_g=0.333333333333333*(element.GetNode(0).Y)+0.333333333333333*(element.GetNode(1).Y)+0.333333333333333*(element.GetNode(2).Y)
          weight_g=9.0/40.0
          velg_x=-1.0*math.sin(xg_g)*math.cos(yg_g)*math.exp(-2.0*(mu/rho)*time)  
          velg_y=math.cos(xg_g)*math.sin(yg_g)*math.exp(-2.0*(mu/rho)*time)
          pressure_g=(rho/4.0)*(math.cos(2.0*xg_g)+math.cos(2.0*yg_g))*math.exp(-4.0*(mu/rho)*time) 
          numerical_velg_x=0.333333333333333*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.333333333333333*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))+0.333333333333333*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
          numerical_velg_y=0.333333333333333*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.333333333333333*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))+0.333333333333333*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
          numerical_pressure_g=0.333333333333333*(element.GetNode(0).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.333333333333333*(element.GetNode(1).GetSolutionStepValue(KratosMultiphysics.PRESSURE))+0.333333333333333*(element.GetNode(2).GetSolutionStepValue(KratosMultiphysics.PRESSURE))
          errorg_x=velg_x-numerical_velg_x
          errorg_y=velg_y-numerical_velg_y
          error_pressureg=pressure_g-numerical_pressure_g
          error_vel=error_vel+(errorg_x*errorg_x+errorg_y*errorg_y)*weight_g*detJ
          error_pressure=error_pressure+(error_pressureg*error_pressureg)*weight_g*detJ

         error_vel=error_vel/totalarea
         error_vel=math.sqrt(error_vel)
         error_pressure=error_pressure/totalarea
         error_pressure=math.sqrt(error_pressure)
         print("dt=", self.timestep)
         print("totalerrorvelocitywithnodalarea=", totalerrorvelocitywithnodalarea) 
         print("totalerrorpressurewithnodalarea=", totalerrorpressurewithnodalarea)
         print("totalerrorvelocity=", totalerrorvelocity) 
         print("totalerrorpressure=", totalerrorpressure)   
         print("totalerrorvelocitywith7gp=", error_vel)  
         print("totalerrorpressurewith7gp=", error_pressure) 
    
    def InitializeSolutionStep(self):
        # If required, compute the BDF coefficients
        if hasattr(self, 'time_discretization'):
            (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)
        # Perform the solver InitializeSolutionStep
        self._GetSolutionStrategy().InitializeSolutionStep()

    def _SetFormulation(self):
        self.formulation = StabilizedFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.element_integrates_in_time = self.formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties
        self.historical_nodal_properties_variables_list = self.formulation.historical_nodal_properties_variables_list
        self.non_historical_nodal_properties_variables_list = self.formulation.non_historical_nodal_properties_variables_list

    def _SetTimeSchemeBufferSize(self):
        scheme_type = self.settings["time_scheme"].GetString()
        if scheme_type == "bossak":
            self.min_buffer_size = 2
        elif scheme_type == "bdf2":
            self.min_buffer_size = 3
        elif scheme_type == "steady":
            self.min_buffer_size = 1
            self._SetUpSteadySimulation()
        else:
            msg  = "Unknown time_scheme option found in project parameters:\n"
            msg += "\"" + scheme_type + "\"\n"
            msg += "Accepted values are \"bossak\", \"bdf2\" or \"steady\".\n"
            raise Exception(msg)

    def _SetUpSteadySimulation(self):
        '''Overwrite time stepping parameters so that they do not interfere with steady state simulations.'''
        self.settings["time_stepping"]["automatic_time_step"].SetBool(False)
        if self.settings["formulation"].Has("dynamic_tau"):
            self.settings["formulation"]["dynamic_tau"].SetDouble(0.0)

    def _SetNodalProperties(self):
        set_density = KratosMultiphysics.DENSITY in self.historical_nodal_properties_variables_list
        set_viscosity = KratosMultiphysics.VISCOSITY in self.historical_nodal_properties_variables_list
        set_sound_velocity = KratosMultiphysics.SOUND_VELOCITY in self.non_historical_nodal_properties_variables_list

        # Get density and dynamic viscostity from the properties of the first element
        for el in self.main_model_part.Elements:
            # Get DENSITY from properties
            if set_density:
                rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                if rho <= 0.0:
                    raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            # Get DYNAMIC_VISCOSITY from properties and calculate the kinematic one (VISCOSITY)
            if set_viscosity:
                dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
                if dyn_viscosity <= 0.0:
                    raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
                kin_viscosity = dyn_viscosity / rho
            # Get SOUND_VELOCITY
            if set_sound_velocity:
                if el.Properties.Has(KratosMultiphysics.SOUND_VELOCITY):
                    sound_velocity = el.Properties.GetValue(KratosMultiphysics.SOUND_VELOCITY)
                else:
                    sound_velocity = 1.0e+12 # Default sound velocity value
                    KratosMultiphysics.Logger.PrintWarning('No \'SOUND_VELOCITY\' value found in Properties {0}. Setting default value {1}'.format(el.Properties.Id, sound_velocity))
                if sound_velocity <= 0.0:
                    raise Exception("SOUND_VELOCITY set to {0} in Properties {1}, positive number expected.".format(sound_velocity, el.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")

        # Transfer the obtained properties to the nodes
        if set_density:
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        if set_viscosity:
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)
        if set_sound_velocity:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.SOUND_VELOCITY, sound_velocity, self.main_model_part.Nodes)
