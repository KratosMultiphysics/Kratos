import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector, Logger

import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.PlasmaDynamicsApplication as PlasmaDynamics
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from KratosMultiphysics.SwimmingDEMApplication.variables_management import VariablesManager as SDEMVariablesManager
import KratosMultiphysics.SwimmingDEMApplication.parameters_tools as PT

def Say(*args):
    Logger.PrintInfo("PlasmaDynamics", *args)
    Logger.Flush()
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)


def GetGlobalVariableByName(variable_name):
    modules = [Kratos, DEM, PlasmaDynamics]
    for mod in modules:
        try:
            return Kratos.KratosGlobals.GetVariable(variable_name)
        except Exception:
            pass
    names = [mod.__name__ for mod in modules]
    error_message = ('No variable with name \''
        + variable_name + '\' exists in either of the modules:\n')
    for name in names[:-1]:
        error_message += name + ', '
    error_message += 'or ' + names[-1] + '.'
    raise AttributeError(error_message)



#TODO: complete the class
class VariablesManager(SDEMVariablesManager):
    
    def __init__(self, parameters):
        super().__init__(parameters)
        
    def ConstructListsOfVariables(self, parameters):
        # PRINTING VARIABLES
        # constructing lists of variables to be printed

        self.nodal_results, self.gauss_points_results = [], []
        #self.fluid_parameters = parameters['fluid_parameters']
        self.fluid_phase_parameters = parameters['fluid_phase_parameters'] 
        self.fluid_parameters = {
            "problem_data"        : self.fluid_phase_parameters["problem_data"],
            "output_processes"    : self.fluid_phase_parameters["output_processes"],               
            "processes"           : self.fluid_phase_parameters["processes"],                
            "solver_settings"     : self.fluid_phase_parameters["solver_settings"]["fluid_solver_settings"]
                                  }  
        
        if 'sdem_output_processes' in self.fluid_parameters:
            gid_output_options = self.fluid_parameters["sdem_output_processes"]["gid_output"][0]["Parameters"]
            result_file_configuration = gid_output_options["postprocess_parameters"]["result_file_configuration"]
            gauss_point_results = result_file_configuration["gauss_point_results"]
            nodal_variables = self.fluid_parameters["sdem_output_processes"]["gid_output"][0]["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_results"]
            self.nodal_results = [nodal_variables[i].GetString() for i in range(nodal_variables.size())]
            self.gauss_points_results = [gauss_point_results[i].GetString() for i in range(gauss_point_results.size())]

        self.ConstructListsOfResultsToPrint(parameters)


        # COUPLING VARIABLES
        # listing the variables involved in the fluid-particles coupling

        if parameters["coupling"]["coupling_level_type"].GetInt():
            self.ConstructListsOfVariablesForCoupling(parameters)

        # VARIABLES TO ADD
        # listing nodal variables to be added to the model parts (memory will be allocated for them)

        # fluid variables
        self.fluid_vars = []
        self.fluid_vars += [Kratos.TORQUE]
        # self.fluid_vars += [Kratos.TEMPERATURE]
        self.fluid_vars += self.fluid_printing_vars
        self.fluid_vars += self.coupling_fluid_vars

        if parameters["pressure_grad_recovery_type"].GetInt() > 0:
            self.fluid_vars += [Fluid.RECOVERED_PRESSURE_GRADIENT]

        if (parameters["gradient_calculation_type"].GetInt() > 1
            or parameters["pressure_grad_recovery_type"].GetInt() > 1
            or parameters["material_acceleration_calculation_type"].GetInt() == 7
            or parameters["laplacian_calculation_type"].GetInt() > 1):
            self.fluid_vars += [Kratos.NODAL_WEIGHTS]

        if parameters["material_acceleration_calculation_type"].GetInt():
            self.fluid_vars += [Kratos.MATERIAL_ACCELERATION]
            self.fluid_vars += [Kratos.VELOCITY_COMPONENT_GRADIENT]

            if (parameters["material_acceleration_calculation_type"].GetInt() == 5
                or parameters["material_acceleration_calculation_type"].GetInt() == 6):
                if parameters["store_full_gradient_option"].GetBool():
                    self.fluid_vars += [Kratos.VELOCITY_X_GRADIENT]
                    self.fluid_vars += [Kratos.VELOCITY_Y_GRADIENT]
                    self.fluid_vars += [Kratos.VELOCITY_Z_GRADIENT]

        if (parameters["vorticity_calculation_type"].GetInt() > 0
            or PT.RecursiveFindParametersWithCondition(parameters["properties"], 'vorticity_induced_lift_parameters')):
            self.fluid_vars += [Kratos.VORTICITY]

        if parameters["laplacian_calculation_type"].GetInt():
            self.fluid_vars += [Kratos.VELOCITY_LAPLACIAN]

        if PT.RecursiveFindTrueBoolInParameters(parameters["properties"], 'do_apply_faxen_corrections'):
            self.fluid_vars += [Kratos.VELOCITY_LAPLACIAN_RATE]

        if parameters["coupling"]["backward_coupling"]["calculate_diffusivity_option"].GetBool():
            self.fluid_vars += [Kratos.CONDUCTIVITY]

        # dem variables
        self.dem_vars = []
        self.dem_vars += self.dem_printing_vars
        self.dem_vars += self.coupling_dem_vars
        self.dem_vars += [Kratos.BUOYANCY]
        self.dem_vars += [Kratos.VELOCITY_OLD]

        if self.do_include_history_force:
            self.dem_vars += [Kratos.BASSET_FORCE]

        if parameters["frame_of_reference"]["frame_type"].GetInt() and self.do_include_history_force > 0:
            self.dem_vars += [SDEM.DISPLACEMENT_OLD]
            self.dem_vars += [Kratos.VELOCITY_OLD_OLD]

        if (parameters["custom_dem"]["translational_integration_scheme"].GetString()
            in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}
            or self.do_include_history_force):
            self.dem_vars += [Kratos.VELOCITY_OLD]
            self.dem_vars += [Kratos.ADDITIONAL_FORCE_OLD]
            self.dem_vars += [Kratos.AUX_VEL]

        if parameters["add_each_hydro_force_option"].GetBool():
            self.dem_vars += [Kratos.DRAG_FORCE]

        self.dem_vars += [Kratos.PARTICLE_SPHERICITY] # TODO: add only when needed

        if (PT.RecursiveFindParametersWithCondition(parameters["properties"], 'vorticity_induced_lift_parameters')
            and parameters["add_each_hydro_force_option"].GetBool()):
            self.dem_vars += [Kratos.LIFT_FORCE]

        if parameters["add_each_hydro_force_option"].GetBool():
            self.dem_vars += [Kratos.VIRTUAL_MASS_FORCE]

        # clusters variables
        self.clusters_vars = []

        # rigid faces variables
        self.rigid_faces_vars = [Kratos.VELOCITY,
                                  Kratos.ANGULAR_VELOCITY,
                                  Kratos.DISPLACEMENT,
                                  DEM.DELTA_DISPLACEMENT,
                                  Kratos.DELTA_ROTATION,
                                  DEM.CONTACT_FORCES,
                                  DEM.DEM_PRESSURE,
                                  DEM.ELASTIC_FORCES,
                                  Kratos.PRESSURE,
                                  DEM.TANGENTIAL_ELASTIC_FORCES,
                                  DEM.SHEAR_STRESS,
                                  Kratos.NODAL_AREA,
                                  Kratos.VELOCITY_OLD]

        if parameters["custom_fluid"]["embedded_option"].GetBool():
            self.rigid_faces_vars += [Kratos.FORCE]
            self.rigid_faces_vars += [Kratos.POSITIVE_FACE_PRESSURE]
            self.rigid_faces_vars += [Kratos.NEGATIVE_FACE_PRESSURE]

        self.fluid_vars += self.rigid_faces_vars

        # inlet variables
        self.inlet_vars = self.dem_vars

    

        
