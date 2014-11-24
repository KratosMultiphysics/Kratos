import math



class Parameters:
    def __init__(self):
        
        self.delta_time = 0.00000001
        
        self.radius = 1.2
        self.sphericity = 1.8
        self.sqrt_of_mass = 1.0
        self.inertia = CalculateInertiaOfBall(self.radius, self.sqrt_of_mass)
        self.fluid_density = 1.0
        self.kinematic_viscosity = 1.0
        self.sphericity = 1.0
        
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

        self.fluid_fraction = 0.8
        self.fluid_model_type = 1
        
        self.buoyancy_force_type = 1
        self.drag_force_type = 1
        self.virtual_mass_force_type = 1
        self.lift_force_type = 1
        self.magnus_force_type = 1        
        self.hydro_torque_type = 1                
        self.drag_porosity_correction_type = 0

        self.gel_strength = 1.0
        self.power_law_n = 1.0
        self.power_law_k = 1.0
        self.yield_stress = 1.0
        self.non_newtonian_option = 1.0
        self.initial_drag_force = 1.0
        self.drag_law_slope = 1.0
        self.power_law_tol = 0.0001
        self.manually_imposed_drag_law_option = 0
        self.drag_modifier_type = 1
        self.nodal_mass_coeff = 1

        self.problem_name="suspended_particles"
        self.kratos_path="D:\Kratos"

def CalculateInertiaOfBall(r, sqrt_of_m):
    return 0.4 * math.pi * r ** 2 * sqrt_of_m ** 2
    