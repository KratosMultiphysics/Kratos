from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import main_script

class HIRMAProjectSolution(main_script.Solution):

    def BeforeSolveOperations(self, time):
        super(HIRMAProjectSolution, self).BeforeSolveOperations(time)
        self.ApplyWaterPressureUpstream(time, self.spheres_model_part)
        
    def ApplyWaterPressureUpstream(self, time, spheres_model_part):

        for element in spheres_model_part.Elements:

            node = element.GetNode(0)                
            is_skin = node.GetSolutionStepValue(SKIN_SPHERE)
            x_coord = node.X
            y_coord = node.Y
            z_coord = node.Z
            external_force = Array3()
            
            water_density = 1000.0
            gravity = 9.81
            # Multiplier in order to accelerate the results
            K = 30.0
            
            # TODO: Hardcoded data which depends on each particular HIRMA geometry and mesh. Getting this data may be further improved
            # This data will be used to define the set of spheres that will be subjected to water pressure
            max_X_to_apply_pressure = 0.0
            dam_Z_max = 2.728
            spheres_diam = 0.2
            dam_Y_max = 3.2
            dam_Y_min = -0.8
            dam_Z_min = -1.4

            if is_skin:
                external_force[0] = external_force[1] = external_force[2] = 0.0
            else:
                continue

            if x_coord < max_X_to_apply_pressure and y_coord > (dam_Y_min + spheres_diam) and y_coord < (dam_Y_max - spheres_diam) and z_coord > (dam_Z_min + spheres_diam):
                # Hydrostatic pressure (converted to forces in spheres)
                # The external force is multiplied by time in order to apply the load progressively. This time has no units
                external_force[0] = max(0.0, K * water_density * gravity * spheres_diam * spheres_diam * time * (dam_Z_max - z_coord))

            node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, external_force)

if __name__=="__main__":
    model = Model()
    HIRMAProjectSolution(model).Run()
