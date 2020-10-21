import KratosMultiphysics as Kratos
from KratosMultiphysics import Model, Parameters
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import os
import math
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis

from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import Say as Say

class BackwardCouplingTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super(BackwardCouplingTestAnalysis, self).__init__(model, varying_parameters)
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'backward_coupling_tests/')

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = 1)

    def FinalizeSolutionStep(self):
        super(BackwardCouplingTestAnalysis, self).FinalizeSolutionStep()

    def Finalize(self):
        super(BackwardCouplingTestAnalysis, self).Finalize()
        spheres_mp = self.model.GetModelPart('SpheresPart')
        fluid_mp = self.model.GetModelPart('FluidModelPart')

        # Adding DEM particles' volume and making sure it is not zero
        total_particles_volume = 0.0

        for node in spheres_mp.Nodes:
            total_particles_volume += node.GetSolutionStepValue(Kratos.RADIUS) ** 3
        total_particles_volume *= 4/3 * math.pi

        assert(not math.isclose(total_particles_volume, 0.0))

        # Checking that the integral of the particles' volume fraction gives the same amount
        particles_component_volume = 0.0
        for node in fluid_mp.Nodes:
            nodal_area = node.GetSolutionStepValue(Kratos.NODAL_AREA)
            particles_volume_fraction = 1 - node.GetSolutionStepValue(Kratos.FLUID_FRACTION)
            particles_component_volume += nodal_area * particles_volume_fraction

        assert(math.isclose(particles_component_volume, total_particles_volume))