import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM

from KratosMultiphysics.DEMApplication.sphere_strategy import ExplicitStrategy as SphereStrategy


class PlasmaStrategy(SphereStrategy):

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, parameters, procedures):
        self.project_parameters = parameters
        super().__init__(all_model_parts, creator_destructor, dem_fem_search, parameters['dem_parameters'], procedures)
        
    def TranslationalIntegrationSchemeTranslator(self, name):
        return SphereStrategy.TranslationalIntegrationSchemeTranslator(self, name)
    
    def RotationalIntegrationSchemeTranslator(self, name_translational, name_rotational):
        return SphereStrategy.RotationalIntegrationSchemeTranslator(self, name_translational, name_rotational)  
        
    def CreateCPlusPlusStrategy(self):
        self.SetVariablesAndOptions()

        if (self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Velocity_Verlet'):
            self.cplusplus_strategy = DEM.IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                                  self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                                  self.search_strategy, self.solver_settings)
        else:
            self.cplusplus_strategy = DEM.ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                                 self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                                 self.search_strategy, self.solver_settings)
    
    def GetTranslationalSchemeInstance(self, class_name):
        try:
            translational_scheme = super().GetTranslationalSchemeInstance(class_name)
        except Exception:
            translational_scheme = PlasmaStrategy.class_name()
        return translational_scheme
    
    def GetRotationalSchemeInstance(self, class_name):
        try:
            rotational_scheme = super().GetRotationalSchemeInstance(class_name)
        except Exception:
            rotational_scheme = PlasmaStrategy.class_name()
        return rotational_scheme
    
    def GetPlasmaConstitutiveLawParametersIfItExists(self, properties):
        if self.project_parameters.Has('properties'):
            for p in self.project_parameters["properties"]:
                if p['properties_id'].GetInt() == int(properties.Id) and p.Has('plasma_dynamics_law_parameters'):
                    return p['plasma_dynamics_law_parameters']
        return None
    
    #TODO
    @staticmethod
    def CreatePlasmaDynamicsLaw(properties, plasma_dynamics_law_parameters):
        pass
    
    def ModifyProperties(self, properties, param = 0):
        super().ModifyProperties(properties, param)
        
        plasma_dynamics_law_parameters = self.GetPlasmaConstitutiveLawParametersIfItExists(properties)
        if plasma_dynamics_law_parameters:
            PlasmaStrategy.CreatePlasmaDynamicsLaw(properties, plasma_dynamics_law_parameters)

        if not param:
            if not properties.Has(Kratos.PARTICLE_SPHERICITY):
                properties[Kratos.PARTICLE_SPHERICITY] = 1.0
        
