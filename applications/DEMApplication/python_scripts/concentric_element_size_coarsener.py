import math
import scipy.stats as stats

import KratosMultiphysics
import KratosMultiphysics.DEMApplication as DEM
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage

def ComputeMeanRadiusOfThisParticle(x, y, z, fine_radius):

    distance_to_origin = math.sqrt(x*x + y*y)

    specimen_type = 1 # 1: CTW16, # 2: CTW10
    if specimen_type == 1:
        max_distance_for_fine_radius = 0.00762
    else:
        max_distance_for_fine_radius = 0.01651

    slope = 0.01
    if distance_to_origin < max_distance_for_fine_radius:
        radius = fine_radius
    else:
        radius = fine_radius + slope * (distance_to_origin - max_distance_for_fine_radius)

    return radius

class ElementSizeModifier(DEMAnalysisStage):

    def __init__(self, model, project_parameters, size_modifier_parameters, flush_frequency=10.0):


        self.size_modifier_parameters = size_modifier_parameters

        default_input_parameters = KratosMultiphysics.Parameters("""
        {
            "mean_diameter_of_particles": 0.1,
            "standard_deviation": 0.5,
            "max_diameter_of_particles": 0.2,
            "min_diameter_of_particles": 0.05,
            "initiation_time": 0.0,
            "process_duration": 1.0,
            "global_damping": 0.8,
            "time_step": 5e-7,
            "geometry_settings": {
                "geometry_type": "Radial",
                "inner_radius": 1.0,
                "outer_radius": 5.0,
                "tolerance": 0.001
            },
            "material_settings":{
            "static_friction": 0.0,
            "dynamic_friciton": 0.0,
            "young_modulus": 5.0e7,
            "density_for_artificial_concentric_weight": 2500
            }
        }
        """)
        self.size_modifier_parameters.ValidateAndAssignDefaults(default_input_parameters)
        self.eraser_counter = 0
        self.list_of_static_frictions_at_start = []
        self.list_of_dynamic_frictions_at_start = []
        self.list_of_young_modulus_at_start = []
        self.list_of_coefficients_of_restittution_at_start = []
        self.list_of_damping_gammas_at_start = []
        project_parameters["GlobalDamping"].SetDouble(self.size_modifier_parameters["global_damping"].GetDouble())
        project_parameters["GravityX"].SetDouble(0.0)
        project_parameters["GravityY"].SetDouble(0.0)
        project_parameters["GravityZ"].SetDouble(0.0)
        project_parameters["RollingFrictionOption"].SetBool(False)
        project_parameters["OutputFileType"].SetString("Ascii")
        project_parameters["MaxTimeStep"].SetDouble(self.size_modifier_parameters["time_step"].GetDouble())
        total_needed_time = self.size_modifier_parameters["initiation_time"].GetDouble() + self.size_modifier_parameters["process_duration"].GetDouble()
        if project_parameters["FinalTime"].GetDouble() < total_needed_time:
            project_parameters["FinalTime"].SetDouble(total_needed_time)
        super().__init__(model, project_parameters)

    def Initialize(self):
        super().Initialize()
        self._SaveInitialRadiusOfAllParticles()
        self._GetDeviationFromMeanSizeOfAllParticles()
        for props in self.spheres_model_part.Properties:
            self.list_of_static_frictions_at_start.append(props[DEM.STATIC_FRICTION])
            props.SetValue(DEM.STATIC_FRICTION, self.size_modifier_parameters["material_settings"]["static_friction"].GetDouble())
            self.list_of_dynamic_frictions_at_start.append(props[DEM.DYNAMIC_FRICTION])
            props.SetValue(DEM.DYNAMIC_FRICTION, self.size_modifier_parameters["material_settings"]["dynamic_friciton"].GetDouble())
            self.list_of_young_modulus_at_start.append(props[KratosMultiphysics.YOUNG_MODULUS])
            props.SetValue(KratosMultiphysics.YOUNG_MODULUS, self.size_modifier_parameters["material_settings"]["young_modulus"].GetDouble())
            self.list_of_coefficients_of_restittution_at_start.append(props[DEM.COEFFICIENT_OF_RESTITUTION])
            props.SetValue(DEM.COEFFICIENT_OF_RESTITUTION, self.size_modifier_parameters["material_settings"]["coefficient_of_restitution"].GetDouble())

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        center = KratosMultiphysics.Array3()
        center[0] = center[1] = center[2] = 0.0 #TODO: input
        density = self.size_modifier_parameters["material_settings"]["density_for_artificial_concentric_weight"].GetDouble()
        self.PreUtilities.ApplyConcentricForceOnParticles(self.spheres_model_part, center, density)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self._ModifyRadiusOfAllParticles()
        self._EraseElementsOutsideDomainAndMarkSkin()

    def Finalize(self):
        counter = 0
        for props in self.spheres_model_part.Properties:
            props.SetValue(DEM.STATIC_FRICTION, self.list_of_static_frictions_at_start[counter])
            props.SetValue(DEM.DYNAMIC_FRICTION, self.list_of_dynamic_frictions_at_start[counter])
            props.SetValue(KratosMultiphysics.YOUNG_MODULUS, self.list_of_young_modulus_at_start[counter])
            props.SetValue(DEM.COEFFICIENT_OF_RESTITUTION, self.list_of_young_modulus_at_start[counter])
            counter += 1
        super().Finalize()

    def _GetDeviationFromMeanSizeOfAllParticles(self):
        min_radius = (self.size_modifier_parameters["max_diameter_of_particles"].GetDouble())/ 2.0
        max_radius = (self.size_modifier_parameters["min_diameter_of_particles"].GetDouble())/ 2.0
        mean_radius_of_particles = 0.5 * self.size_modifier_parameters["mean_diameter_of_particles"].GetDouble()
        radius_standard_deviation = self.size_modifier_parameters["standard_deviation"].GetDouble() / 2.0

        distribution = stats.truncnorm((min_radius - mean_radius_of_particles)/ radius_standard_deviation, (max_radius - mean_radius_of_particles) / radius_standard_deviation, loc=mean_radius_of_particles, scale=radius_standard_deviation)

        for node in self.spheres_model_part.Nodes:
            random_radius = distribution.rvs(1)
            deviation = random_radius[0] - mean_radius_of_particles
            node.SetValue(DEM.DEVIATION, deviation)

    def _SaveInitialRadiusOfAllParticles(self):
        for node in self.spheres_model_part.Nodes:
            node.SetValue(DEM.INITIAL_RADIUS, node.GetSolutionStepValue(KratosMultiphysics.RADIUS))

    def _ModifyRadiusOfAllParticles(self):
        if self.time > self.size_modifier_parameters["initiation_time"].GetDouble() + self.size_modifier_parameters["process_duration"].GetDouble():
            return
        mean_diameter_of_particles = self.size_modifier_parameters["mean_diameter_of_particles"].GetDouble()
        for node in self.spheres_model_part.Nodes:
            x = node.X
            y = node.Y
            z = node.Z
            radius_defined_by_function = ComputeMeanRadiusOfThisParticle(x, y, z, mean_diameter_of_particles/2.0)
            deviation_for_mean_radius = node.GetValue(DEM.DEVIATION)
            initial_radius = node.GetValue(DEM.INITIAL_RADIUS)
            actual_deviation_for_current_radius = deviation_for_mean_radius * radius_defined_by_function / initial_radius
            intended_radius_at_end = radius_defined_by_function + actual_deviation_for_current_radius
            portion_of_process = (self.time - self.size_modifier_parameters["initiation_time"].GetDouble()) / self.size_modifier_parameters["process_duration"].GetDouble()
            intended_radius_at_current_time = initial_radius + portion_of_process * (intended_radius_at_end - initial_radius)
            node.SetSolutionStepValue(KratosMultiphysics.RADIUS, intended_radius_at_current_time)

    def _EraseElementsOutsideDomainAndMarkSkin(self):
        if self.eraser_counter == 1:
            self.eraser_counter = 0
            mean_diameter_of_particles = self.size_modifier_parameters["mean_diameter_of_particles"].GetDouble()
            if self.size_modifier_parameters["geometry_settings"]["geometry_type"].GetString() == "Radial":
                max_radius = self.size_modifier_parameters["geometry_settings"]["outer_radius"].GetDouble()
                tolerance = self.size_modifier_parameters["geometry_settings"]["tolerance"].GetDouble()
                center = KratosMultiphysics.Array3()
                center[0] = center[1] = center[2] = 0.0
                self.PreUtilities.ResetSkinParticles(self.spheres_model_part)
                self.PreUtilities.MarkToEraseParticlesOutsideRadius(self.spheres_model_part, max_radius, center, tolerance)
                inner_radius = self.size_modifier_parameters["geometry_settings"]["inner_radius"].GetDouble()
                radius_at_inner_boundary = mean_diameter_of_particles/2.0
                self.PreUtilities.SetSkinParticlesInnerBoundary(self.spheres_model_part, inner_radius, 2.0 * radius_at_inner_boundary)
                radius_at_outer_boundary = ComputeMeanRadiusOfThisParticle(max_radius, 0.0, 0.0, mean_diameter_of_particles/2.0)
                portion_of_process = (self.time - self.size_modifier_parameters["initiation_time"].GetDouble()) / self.size_modifier_parameters["process_duration"].GetDouble()
                radius_at_outer_boundary_at_current_time = mean_diameter_of_particles/2.0 + portion_of_process * (radius_at_outer_boundary - mean_diameter_of_particles/2.0)
                self.PreUtilities.SetSkinParticlesOuterBoundary(self.spheres_model_part, max_radius, 1.4 * radius_at_outer_boundary_at_current_time)
        else:
            self.eraser_counter += 1

if __name__ == "__main__":
    from KratosMultiphysics import Logger
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.INFO)
    with open("ProjectParametersDEM.json",'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    with open("SizeModifierParametersDEM.json",'r') as size_modifier_parameter_file:
        size_modifier_parameters = KratosMultiphysics.Parameters(size_modifier_parameter_file.read())

    model = KratosMultiphysics.Model()
    ElementSizeModifier(model, project_parameters, size_modifier_parameters).Run()
