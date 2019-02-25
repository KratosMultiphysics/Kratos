from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
import copy
import random

def RandomPositive(supremum = 1.0):
    value = random.random() * supremum
    
    if value < 0.00000001:
        value = 0.00000001
    
    return value

def CalculateInertiaOfBall(r, sqrt_of_m):
    return 0.4 * math.pi * r ** 2 * sqrt_of_m ** 2

def PadWithSpaces(lines):   
    width = GetColumMaxWidth(lines, 2)
    aux = copy.deepcopy(lines)

    del lines[:]
    
    for line_aux in aux:        
        lines += [line_aux.ljust(width)]
    
def GetColumMaxWidth(lines, margin):
    width = 0

    for line in lines:
        width = max(width, len(line))
        
    width += margin
    
    return width

class Parameters:
    def __init__(self):
        
        self.delta_time = RandomPositive()        
        self.radius = RandomPositive()
        self.sphericity = RandomPositive()
        self.sqrt_of_mass = RandomPositive()
        self.inertia = CalculateInertiaOfBall(self.radius, self.sqrt_of_mass)
        self.fluid_density = RandomPositive()
        self.kinematic_viscosity = RandomPositive()
        self.fluid_fraction = RandomPositive()
        self.gel_strength = RandomPositive()
        self.power_law_n = RandomPositive(3.0)
        self.power_law_k = RandomPositive(1000)
        self.yield_stress = RandomPositive(1000)
        self.initial_drag_force = RandomPositive(1000)
        self.drag_law_slope = RandomPositive(1000)
        self.power_law_tol = RandomPositive(0.0001)
        
        self.coor_x = 0.0
        self.coor_y = 0.0
        self.coor_z = 0.0
        
        self.displ_x = 0.0
        self.displ_y = 0.0
        self.displ_z = 0.0

        self.velocity_x = 0.0
        self.velocity_y = 0.0
        self.velocity_z = 0.0
        
        self.velocity_old_x = 0.0
        self.velocity_old_y = 0.0
        self.velocity_old_z = 0.0
        
        self.accel_x = 0.0
        self.accel_y = 0.0
        self.accel_z = 0.0                

        self.rot_angle_x = 0.0
        self.rot_angle_y = 0.0
        self.rot_angle_z = 0.0

        self.euler_angle_x = 0.0
        self.euler_angle_y = 0.0
        self.euler_angle_z = 0.0        

        self.ang_vel_x = 0.0
        self.ang_vel_y = 0.0
        self.ang_vel_z = 0.0                
        
        self.gravity_x = 0.0
        self.gravity_y = 0.0
        self.gravity_z = -9.81
        
        self.tota_force_x = 0.0
        self.tota_force_y = 0.0
        self.tota_force_z = 0.0

        self.fluid_velocity_x = 0.0
        self.fluid_velocity_y = 0.0
        self.fluid_velocity_z = 0.0

        self.fluid_acceleration_x = 0.0
        self.fluid_acceleration_y = 0.0
        self.fluid_acceleration_z = 0.0

        self.fluid_vorticity_x = 0.0
        self.fluid_vorticity_y = 0.0
        self.fluid_vorticity_z = 0.0

        self.shear_rate_projected = 0.2

        self.pressure_gradient_x = 0.0
        self.pressure_gradient_y = 0.0
        self.pressure_gradient_z = 0.0
        
        self.non_newtonian_option = 1.0
        self.manually_imposed_drag_law_option = 0
        self.fluid_model_type = 1       
        self.buoyancy_force_type = 0
        self.drag_force_type = 0
        self.virtual_mass_force_type = 0
        self.lift_force_type = 0
        self.magnus_force_type = 0        
        self.hydro_torque_type = 0                
        self.drag_porosity_correction_type = 0
        self.drag_modifier_type = 1

        self.nodal_mass_coeff = 1

        self.problem_name="suspended_particles"
        self.kratos_path="D:\Kratos"
   
    def GetParametersString(self):
        my_string = ""
        keys = []
        values = []
        
        for var in vars(self).keys():
            keys += [str(var)]
        
        PadWithSpaces(keys)
        
        for var in vars(self).values():
            values += [str(var)]
         
        i = 0
        
        for var in keys:
            my_string += "\n" + var + "= " + values[i]
            i += 1
            
        return my_string