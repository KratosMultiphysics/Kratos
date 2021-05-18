import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import os
import math
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis

class BackwardCouplingTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super().__init__(model, varying_parameters)
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'backward_coupling_tests/')

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = 1)

    def Finalize(self):
        super().Finalize()
        self._CalculateVolumes()

    def _CalculateVolumes(self):
        spheres_mp = self.model.GetModelPart('SpheresPart')
        fluid_mp = self.model.GetModelPart('FluidModelPart')

        # Adding up DEM particles' volume and making sure it is not zero
        total_particles_volume = 0.0

        for node in spheres_mp.Nodes:
            total_particles_volume += node.GetSolutionStepValue(Kratos.RADIUS) ** 3
        total_particles_volume *= 4/3 * math.pi

        # Calculating the integral of the particles' volume, which should be the
        # same amount (up to finite-precision errors)
        particles_component_volume = 0.0
        for node in fluid_mp.Nodes:
            nodal_area = node.GetSolutionStepValue(Kratos.NODAL_AREA)
            particles_volume_fraction = 1 - node.GetSolutionStepValue(Kratos.FLUID_FRACTION)
            particles_component_volume += nodal_area * particles_volume_fraction

        results = Kratos.Parameters("{}")
        Add = results.AddEmptyValue
        Add("particles_component_volume").SetDouble(particles_component_volume)
        Add("total_particles_volume").SetDouble(total_particles_volume)
        self.project_parameters.AddValue('results', results)