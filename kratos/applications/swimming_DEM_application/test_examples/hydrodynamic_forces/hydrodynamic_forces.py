from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
import copy
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.MeshingApplication import *
import swimming_DEM_gid_output
import ProjectParameters
import random

def AddVariables(model_part, pp):
    AddNodalVariables(model_part)
    AddAndInitializeProcessInfoVariables(model_part, pp)
        
def InitializeVariables(model_part, pp):
    InitializeNodalVariables(model_part, pp)
    AddAndInitializeProcessInfoVariables(model_part, pp)       

def InitializeNodalVariables(model_part, pp):
    for node in model_part.Nodes:
        node.X = pp.coor_x
        node.Y = pp.coor_y
        node.Z = pp.coor_z
        
        node.SetSolutionStepValue(DISPLACEMENT_X, pp.displ_x)
        node.SetSolutionStepValue(DISPLACEMENT_Y, pp.displ_y)
        node.SetSolutionStepValue(DISPLACEMENT_Z, pp.displ_z)
        
        node.SetSolutionStepValue(VELOCITY_X, pp.velocity_x)
        node.SetSolutionStepValue(VELOCITY_Y, pp.velocity_y)
        node.SetSolutionStepValue(VELOCITY_Z, pp.velocity_z)
              
        node.SetSolutionStepValue(VELOCITY_X, 1, pp.velocity_old_x)
        node.SetSolutionStepValue(VELOCITY_Y, 1, pp.velocity_old_y)
        node.SetSolutionStepValue(VELOCITY_Z, 1, pp.velocity_old_z)
        
        node.SetSolutionStepValue(PARTICLE_ROTATION_ANGLE_X, pp.rot_angle_x)
        node.SetSolutionStepValue(PARTICLE_ROTATION_ANGLE_Y, pp.rot_angle_y)
        node.SetSolutionStepValue(PARTICLE_ROTATION_ANGLE_Z, pp.rot_angle_z)     
        
        node.SetSolutionStepValue(EULER_ANGLES_X, pp.euler_angle_x)
        node.SetSolutionStepValue(EULER_ANGLES_Y, pp.euler_angle_y)
        node.SetSolutionStepValue(EULER_ANGLES_Z, pp.euler_angle_z)          
        
        node.SetSolutionStepValue(ANGULAR_VELOCITY_X, pp.ang_vel_x)
        node.SetSolutionStepValue(ANGULAR_VELOCITY_Y, pp.ang_vel_y)
        node.SetSolutionStepValue(ANGULAR_VELOCITY_Z, pp.ang_vel_z)
               
        node.SetSolutionStepValue(TOTAL_FORCES_X, pp.tota_force_x)
        node.SetSolutionStepValue(TOTAL_FORCES_Y, pp.tota_force_y)
        node.SetSolutionStepValue(TOTAL_FORCES_Z, pp.tota_force_z)       
        
        node.SetSolutionStepValue(FLUID_VEL_PROJECTED_X, pp.fluid_velocity_x)
        node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Y, pp.fluid_velocity_y)
        node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Z, pp.fluid_velocity_z)        
        
        node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_X, pp.fluid_acceleration_x)
        node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_Y, pp.fluid_acceleration_y)
        node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_Z, pp.fluid_acceleration_z)        
        
        node.SetSolutionStepValue(FLUID_VORTICITY_PROJECTED_X, pp.fluid_vorticity_x)
        node.SetSolutionStepValue(FLUID_VORTICITY_PROJECTED_Y, pp.fluid_vorticity_y)
        node.SetSolutionStepValue(FLUID_VORTICITY_PROJECTED_Z, pp.fluid_vorticity_z)
        
        node.SetSolutionStepValue(PRESSURE_GRAD_PROJECTED_X, pp.pressure_gradient_x)
        node.SetSolutionStepValue(PRESSURE_GRAD_PROJECTED_Y, pp.pressure_gradient_y)
        node.SetSolutionStepValue(PRESSURE_GRAD_PROJECTED_Z, pp.pressure_gradient_z)                    
        
        node.SetSolutionStepValue(RADIUS, pp.radius)
        node.SetSolutionStepValue(PARTICLE_SPHERICITY, pp.sphericity)
        node.SetSolutionStepValue(SQRT_OF_MASS, pp.sqrt_of_mass)
        node.SetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA, pp.inertia)
        node.SetSolutionStepValue(FLUID_FRACTION_PROJECTED, pp.fluid_fraction)
        node.SetSolutionStepValue(FLUID_DENSITY_PROJECTED, pp.fluid_density)
        node.SetSolutionStepValue(FLUID_VISCOSITY_PROJECTED, pp.kinematic_viscosity)
        node.SetSolutionStepValue(POWER_LAW_N, pp.power_law_n)
        node.SetSolutionStepValue(POWER_LAW_K, pp.power_law_k)    
        node.SetSolutionStepValue(YIELD_STRESS, pp.yield_stress)   
        node.SetSolutionStepValue(SHEAR_RATE_PROJECTED, pp.shear_rate_projected)
           
    
def AddVariables(model_part):
    # COMMON
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
    model_part.AddNodalSolutionStepVariable(GROUP_ID)   

    # KINEMATIC
    model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
    
    # FORCES
    model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
    model_part.AddNodalSolutionStepVariable(DAMP_FORCES)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_FORCE)

    # BASIC PARTICLE PROPERTIES
    model_part.AddNodalSolutionStepVariable(RADIUS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_SPHERICITY)
    model_part.AddNodalSolutionStepVariable(SQRT_OF_MASS)
    model_part.AddNodalSolutionStepVariable(PRINCIPAL_MOMENTS_OF_INERTIA)
    model_part.AddNodalSolutionStepVariable(CHARACTERISTIC_LENGTH)
    model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)

    # ROTATION RELATED PROPERTIES
    model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA)
    model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_DAMP_RATIO)
    model_part.AddNodalSolutionStepVariable(ROLLING_FRICTION)

    # OTHER PROPERTIES
    model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD

    # LOCAL AXIS
    model_part.AddNodalSolutionStepVariable(EULER_ANGLES)

    # FLAGS
    model_part.AddNodalSolutionStepVariable(GROUP_ID)            # Differencied groups for plotting, etc..
    
    # ONLY VISUALIZATION
    model_part.AddNodalSolutionStepVariable(EXPORT_ID)
    #model_part.AddNodalSolutionStepVariable(EXPORT_GROUP_ID)
    
    # SWIMMING
    model_part.AddNodalSolutionStepVariable(REYNOLDS_NUMBER)
    model_part.AddNodalSolutionStepVariable(PRESSURE_GRAD_PROJECTED)
    model_part.AddNodalSolutionStepVariable(HYDRODYNAMIC_FORCE)
    model_part.AddNodalSolutionStepVariable(FLUID_VEL_PROJECTED)
    model_part.AddNodalSolutionStepVariable(FLUID_ACCEL_PROJECTED)
    model_part.AddNodalSolutionStepVariable(FLUID_FRACTION_PROJECTED)
    model_part.AddNodalSolutionStepVariable(FLUID_DENSITY_PROJECTED)
    model_part.AddNodalSolutionStepVariable(FLUID_VISCOSITY_PROJECTED)
    model_part.AddNodalSolutionStepVariable(BUOYANCY)
    model_part.AddNodalSolutionStepVariable(DRAG_FORCE)
    model_part.AddNodalSolutionStepVariable(VIRTUAL_MASS_FORCE)
    model_part.AddNodalSolutionStepVariable(LIFT_FORCE)    
    model_part.AddNodalSolutionStepVariable(FLUID_VORTICITY_PROJECTED)
    model_part.AddNodalSolutionStepVariable(SHEAR_RATE_PROJECTED)
    model_part.AddNodalSolutionStepVariable(FLUID_ACCEL_PROJECTED)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(POWER_LAW_N)
    model_part.AddNodalSolutionStepVariable(POWER_LAW_K)
    model_part.AddNodalSolutionStepVariable(POWER_LAW_K)
    model_part.AddNodalSolutionStepVariable(YIELD_STRESS)   
    
def AddAndInitializeProcessInfoVariables(model_part, pp):    
    # SIMULATION FLAGS
    
    model_part.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, 0)
    model_part.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, 0)
    model_part.ProcessInfo.SetValue(ROTATION_OPTION, 1)       

    # GLOBAL MATERIAL PROPERTIES
    model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, pp.nodal_mass_coeff)

    # PRINTING VARIABLES

    model_part.ProcessInfo.SetValue(FORCE_CALCULATION_TYPE, 0)
    model_part.ProcessInfo.SetValue(DAMP_TYPE, 0)
    model_part.ProcessInfo.SetValue(ROLLING_FRICTION_OPTION, 0)
    model_part.ProcessInfo.SetValue(PRINT_EXPORT_ID, 1)

    # TIME RELATED PARAMETERS
    model_part.ProcessInfo.SetValue(DELTA_TIME, pp.delta_time)
    
    # SWIMMING
    model_part.ProcessInfo.SetValue(BUOYANCY_FORCE_TYPE, pp.buoyancy_force_type)
    model_part.ProcessInfo.SetValue(DRAG_FORCE_TYPE, pp.drag_force_type)
    model_part.ProcessInfo.SetValue(VIRTUAL_MASS_FORCE_TYPE, pp.virtual_mass_force_type)
    model_part.ProcessInfo.SetValue(LIFT_FORCE_TYPE, pp.lift_force_type)
    model_part.ProcessInfo.SetValue(MAGNUS_FORCE_TYPE, pp.magnus_force_type)
    model_part.ProcessInfo.SetValue(DRAG_POROSITY_CORRECTION_TYPE, pp.drag_porosity_correction_type)
    model_part.ProcessInfo.SetValue(FLUID_MODEL_TYPE, pp.fluid_model_type)
    model_part.ProcessInfo.SetValue(MANUALLY_IMPOSED_DRAG_LAW_OPTION, pp.manually_imposed_drag_law_option)
    model_part.ProcessInfo.SetValue(DRAG_MODIFIER_TYPE, pp.drag_modifier_type)
    model_part.ProcessInfo.SetValue(INIT_DRAG_FORCE, pp.initial_drag_force)
    model_part.ProcessInfo.SetValue(DRAG_LAW_SLOPE, pp.drag_law_slope)
    model_part.ProcessInfo.SetValue(POWER_LAW_TOLERANCE, pp.power_law_tol)    
    model_part.ProcessInfo.SetValue(GRAVITY_X, pp.gravity_x)
    model_part.ProcessInfo.SetValue(GRAVITY_Y, pp.gravity_y)
    model_part.ProcessInfo.SetValue(GRAVITY_Z, pp.gravity_z)
    
           
    model_part.ProcessInfo.SetValue(CASE_OPTION, 1)
    model_part.ProcessInfo.SetValue(TRIHEDRON_OPTION, 0)
    model_part.ProcessInfo.SetValue(BOUNDING_BOX_OPTION, 0)
    model_part.ProcessInfo.SetValue(FIX_VELOCITIES_FLAG, 0)
    model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0);
    model_part.ProcessInfo.SetValue(TOTAL_CONTACTS, 0);
    model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, 0)
    model_part.ProcessInfo.SetValue(ACTIVATE_SEARCH, 1)  # needed in the basic for the continuum.
    
def AddDofs(model_part):

    for node in model_part.Nodes:
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(ANGULAR_VELOCITY_X, REACTION_X)
        node.AddDof(ANGULAR_VELOCITY_Y, REACTION_Y)
        node.AddDof(ANGULAR_VELOCITY_Z, REACTION_Z)

def PadWithSpacesa(lines, lines_to_print):
    widths = GetColumsMaxWidths(lines, 2)

    for line in lines:
        i = 0
        
        for i_entry in range(0, len(line)):
            line[i_entry] = line[i_entry].ljust(widths[i_entry])
    
    i = 0
    
    for line in lines:
        line_to_print = ""
        
        for entry in line:
            line_to_print += entry
                
        lines_to_print += [line_to_print]
    
    return sum(widths)    

def PadWithSpaces(lines):    
    widths = GetColumsMaxWidths(lines, 2)

    for line in lines:
        
        for i_entry in range(0, len(line)):
            line[i_entry] = line[i_entry].ljust(widths[i_entry])     
   
    aux = copy.deepcopy(lines)
    
    del lines[:]
    
    for line_aux in aux:
        line = ""
        
        for entry in line_aux:
            line += entry
         
        lines += [line]
   
    return sum(widths)       
    
def GetColumsMaxWidths(lines, margin):
    widths = []
    
    if (not len(lines)):
        return widths    
        
    for entry in lines[0]:
        widths += [1]
    
    for line in lines:
        i = 0
        
        for entry in line:
            
            if (len(str(entry)) > widths[i]):
                widths[i] = len(str(entry))              
                
            i += 1
        
    for i in range(0, len(widths)):
        widths[i] += margin

    return widths
   
def Norm(v):
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    
def Distance(v1, v2):
    return math.sqrt((v1[0] - v2 [0]) ** 2 + (v1[1] - v2 [1]) ** 2 + (v1[2] - v2 [2]) ** 2)

def Dot(v1, v2):
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    
def Normalize(v, modulus = 1.0):
    
    if (Norm(v) == 0.0):
        return
    
    else:
        coeff = modulus / Norm(v)
        v[0] *= coeff
        v[1] *= coeff
        v[2] *= coeff
        
def RandomVector(modulus = 1.0):
    v = Array3()
    v[0] = random.choice([-1, 1]) * random.random()
    v[1] = random.choice([-1, 1]) * random.random()
    v[2] = random.choice([-1, 1]) * random.random()
    Normalize(v, modulus)
    
    return v
    
def PrintResults(title, tests):
    lines = []
    lines += [["Test Number", "Target", "Calculated", "Error", "Veredict", "Description"]]
    
    i_test = 0
    
    for test in tests:
        lines += [[str(i_test)] + test.string_results]
        i_test += 1

    total_width = PadWithSpaces(lines)    
    print()
    print ("=" * total_width)
    print (title)
    print ("-" * total_width)
    print (lines[0])
    print ("-" * total_width)
    print (lines[1])
    print (lines[2])
    print (lines[3])
    print ("=" * total_width)
    print()        

class Benchmark:
    def __init__(self):
        pass
    
    tests = []
    
    @staticmethod
    def PrintResults(title, tests):
        lines = []
        lines += [["Test id", "Target", "Calculated", "Error", "Veredict", "Description"]]
        
        i_test = 0
        
        for test in tests:
            if (test.has_results):
                lines += [[str(i_test)] + test.string_results]            
                i_test += 1

        total_width = PadWithSpaces(lines)    
        print()
        print ("=" * total_width)
        print (title)
        print ("-" * total_width)
        print (lines[0])
        print ("-" * total_width)
        
        for i in range(1, i_test + 1):
            print (lines[i])

        print ("=" * total_width)
        print()
     
    @staticmethod  
    def ErrorMetric(v1, v2):
        if (Norm(v1) == 0.0 and Norm(v2) == 0.0):
            return 0.0
            
        else:
            return Distance(v1, v2) / (Norm(v1) + Norm(v2))
    
    
class BuoyancyBenchmark(Benchmark):    
    title = "Buoyancy force test results"
    tests = []
    
    @staticmethod
    def PrintResults():
        Benchmark.PrintResults(BuoyancyBenchmark.title, BuoyancyBenchmark.tests)
        
    @staticmethod
    def ConvertToStrings(results):
         
        if results[3]:
            word = "OK"
        else:
            word = "Fail"
        print()
        return [str(results[0]), str(results[1]), "{:.2e}".format(results[2]), word, str(results[4])]
    
    def __init__(self, pp, buoyancy_force_type, drag_force_type, pressure_gradient, description):
        self.pp = copy.deepcopy(pp) 
        self.buoyancy_tol = 10e-6

        self.pp.virtual_mass_force_type = 0
        self.pp.lift_force_type = 0
        self.magnus_force_type = 0             
        self.pp.pressure_gradient_x = pressure_gradient[0]
        self.pp.pressure_gradient_y = pressure_gradient[1]
        self.pp.pressure_gradient_z = pressure_gradient[2]
        self.pp.drag_force_type = drag_force_type
        self.pp.buoyancy_force_type = buoyancy_force_type
        
        self.description = description
        self.has_results = False
        
        BuoyancyBenchmark.tests += [self]
        Benchmark.tests += [self]
     
    def Test(self, model_part, benchmark_utils, target_buoyancy):
        self.target_buoyancy = target_buoyancy
        InitializeVariables(model_part, self.pp)
        benchmark_utils.ComputeHydrodynamicForces(model_part)

        for node in model_part.Nodes:
            buoyancy = node.GetSolutionStepValue(BUOYANCY)
                  
        error = Benchmark.ErrorMetric(buoyancy, self.target_buoyancy)
        
        if (error < self.buoyancy_tol):
            veredict = True
            
        else:
            veredict = False
            
        self.results = [self.target_buoyancy, buoyancy, error, veredict, self.description] 
        self.string_results = BuoyancyBenchmark.ConvertToStrings(self.results)
        self.has_results = True

class DragBenchmark(Benchmark):
    title = "Drag force test results"
    tests = []
    
    @staticmethod
    def PrintResults():
        Benchmark.PrintResults(DragBenchmark.title, DragBenchmark.tests)
        
    @staticmethod
    def ConvertToStrings(results):
         
        if results[3]:
            word = "OK"
        else:
            word = "Fail"
            
        return [str(results[0]), str(results[1]), "{:.2e}".format(results[2]), word, str(results[4])]
        
    def __init__(self, pp, drag_force_type, particle_reynolds, fluid_fraction, sphericity, description):
        self.pp = copy.deepcopy(pp) 
        self.pp.drag_force_type = drag_force_type 
        self.pp.fluid_fraction = fluid_fraction
        self.pp.sphericity = sphericity
        self.description = description

        self.drag_tol = 10e-6
        self.pp.buoyancy_force_type = 0
        self.pp.virtual_mass_force_type = 0
        self.pp.lift_force_type = 0  
        self.magnus_force_type = 0             
        self.CalculateFlowVariables(reynolds)

        self.has_results = False        
        DragBenchmark.tests += [self]
        Benchmark.tests += [self]
    
    def CalculateFlowVariables(self, reynolds):
        
        if (reynolds < 0):
            raise ValueError("The particle's Reynolds number should be non-negative")
        
        elif (reynolds == 0.0):
            pass
            
        else:
            velocity = RandomVector()
            fluid_velocity = RandomVector()
            self.pp.velocity_x = velocity[0]
            self.pp.velocity_y = velocity[1]
            self.pp.velocity_z = velocity[2]
            self.pp.fluid_velocity_x = fluid_velocity[0]
            self.pp.fluid_velocity_y = fluid_velocity[1]
            self.pp.fluid_velocity_z = fluid_velocity[2]
            slip_vel = fluid_velocity - velocity
            r = self.pp.radius
            dens = self.pp.fluid_density
           
            self.pp.kinematic_viscosity = 2 * r * dens * Norm(slip_vel) / reynolds
            
    def Test(self, model_part, benchmark_utils, target_drag):
        self.target_drag = target_drag  
        InitializeVariables(model_part, self.pp)
        benchmark_utils.ComputeHydrodynamicForces(model_part)

        for node in model_part.Nodes:
            drag = node.GetSolutionStepValue(DRAG_FORCE)
                  
        error = Benchmark.ErrorMetric(drag, self.target_drag)
        
        if (error < self.drag_tol):
            veredict = True
            
        else:
            veredict = False
            
        self.results = [self.target_drag, drag, error, veredict, self.description] 
        self.string_results = DragBenchmark.ConvertToStrings(self.results)
        self.has_results = True

    
class VirtualMassBenchmark(Benchmark):    
    title = "Virtual mass force test results"
    tests = []
    
    @staticmethod
    def PrintResults():
        Benchmark.PrintResults(VirtualMassBenchmark.title, VirtualMassBenchmark.tests)
        
    @staticmethod
    def ConvertToStrings(results):
         
        if results[3]:
            word = "OK"
        else:
            word = "Fail"
            
        return [str(results[0]), str(results[1]), "{:.2e}".format(results[2]), word, str(results[4])]
    
    def __init__(self, pp, virtual_mass_force_type, acceleration_number, fluid_fraction, description):
        self.pp = copy.deepcopy(pp) 
        self.virtual_mass_tol = 10e-6
        
        self.pp.buoyancy_force_type = 0
        self.pp.lift_force_type = 0
        self.magnus_force_type = 0             
        self.pp.drag_force_type = 0
        self.pp.virtual_mass_force_type = virtual_mass_force_type
        self.pp.fluid_fraction = fluid_fraction
        self.description = description
        self.has_results = False
        self.CalculateFlowVariables(acceleration_number)
        
        VirtualMassBenchmark.tests += [self]
        Benchmark.tests += [self]
        
    def CalculateFlowVariables(self, acceleration_number):
        
        if (acceleration_number < 0):
            raise ValueError("The particle's Reynolds number should be non-negative")
        
        elif (acceleration_number == 0.0):
            pass
            
        else:
            velocity = RandomVector()
            fluid_velocity = RandomVector()
            acceleration = RandomVector()
            velocity_old = velocity - self.pp.delta_time * acceleration
            fluid_acceleration = RandomVector()
            self.pp.velocity_x = velocity[0]
            self.pp.velocity_y = velocity[1]
            self.pp.velocity_z = velocity[2]
            self.pp.velocity_old_x = velocity_old[0]
            self.pp.velocity_old_y = velocity_old[1]
            self.pp.velocity_old_z = velocity_old[2]
            self.pp.accel_x = acceleration[0]
            self.pp.accel_y = acceleration[1]
            self.pp.accel_z = acceleration[2]                        
            self.pp.fluid_velocity_x = fluid_velocity[0]
            self.pp.fluid_velocity_y = fluid_velocity[1]
            self.pp.fluid_velocity_z = fluid_velocity[2]
            self.pp.fluid_acceleration_x = fluid_acceleration[0]
            self.pp.fluid_acceleration_y = fluid_acceleration[1]
            self.pp.fluid_acceleration_z = fluid_acceleration[2]            
            
            slip_vel = fluid_velocity - velocity
            slip_accel = fluid_acceleration - acceleration
            self.pp.radius = Dot(slip_vel, slip_vel) / abs(4 * Dot(slip_vel, slip_accel) * acceleration_number)

    def Test(self, model_part, benchmark_utils, target_virtual_mass):
        self.target_virtual_mass = target_virtual_mass
        InitializeVariables(model_part, self.pp)
        benchmark_utils.ComputeHydrodynamicForces(model_part)

        for node in model_part.Nodes:
            virtual_mass = node.GetSolutionStepValue(VIRTUAL_MASS_FORCE)
                  
        error = Benchmark.ErrorMetric(virtual_mass, self.target_virtual_mass)
        
        if (error < self.virtual_mass_tol):
            veredict = True
            
        else:
            veredict = False
            
        self.results = [self.target_virtual_mass, virtual_mass, error, veredict, self.description] 
        self.string_results = VirtualMassBenchmark.ConvertToStrings(self.results)
        self.has_results = True   
 
    
class SaffmanBenchmark(Benchmark):    
    title = "Saffman force test results"
    tests = []
    
    @staticmethod
    def PrintResults():
        Benchmark.PrintResults(SaffmanBenchmark.title, SaffmanBenchmark.tests)
        
    @staticmethod
    def ConvertToStrings(results):
         
        if results[3]:
            word = "OK"
        else:
            word = "Fail"
            
        return [str(results[0]), str(results[1]), "{:.2e}".format(results[2]), word, str(results[4])]
    
    def __init__(self, pp, saffman_force_type, reynolds, reynolds_shear, description):
        self.pp = copy.deepcopy(pp) 
        self.saffman_tol = 10e-6
        
        self.pp.buoyancy_force_type = 0
        self.pp.lift_force_type = saffman_force_type
        self.magnus_force_type = 0     
        self.pp.drag_force_type = 0
        self.pp.fluid_fraction = fluid_fraction
        self.description = description
        
        self.has_results = False
        self.CalculateFlowVariables(reynolds, reynolds_shear)
        
        SaffmanBenchmark.tests += [self]
        Benchmark.tests += [self]
        
    def CalculateFlowVariables(self, reynolds, reynolds_shear):
        
        if (reynolds < 0 or reynolds_shear < 0):
            raise ValueError("The particle's Reynold's number and shear Reynold's numbers should be non-negative")
        
        elif (reynolds_shear * reynolds == 0.0):
            pass
            
        else:
            velocity = RandomVector()
            fluid_velocity = RandomVector()
            ang_velocity = RandomVector()
            vorticity = RandomVector()
            
            self.pp.velocity_x = velocity[0]
            self.pp.velocity_y = velocity[1]
            self.pp.velocity_z = velocity[2]
            self.ang_vel_x = ang_velocity[0]
            self.ang_vel_y = ang_velocity[1]
            self.ang_vel_z = ang_velocity[2]                                   
            self.pp.fluid_velocity_x = fluid_velocity[0]
            self.pp.fluid_velocity_y = fluid_velocity[1]
            self.pp.fluid_velocity_z = fluid_velocity[2]    
            
            slip_vel = fluid_velocity - velocity
            slip_vort = vorticity - ang_velocity
            r = self.pp.radius
            dens = self.pp.fluid_density
            self.pp.kinematic_viscosity = 2 * r * dens * Norm(slip_vel) / reynolds
            vorticity_norm = self.pp.kinematic_viscosity * reynolds_shear / (dens * 4 * r ** 2 )
            Normalize(vorticity, vorticity_norm)
            self.pp.fluid_vorticity_x = vorticity[0]
            self.pp.fluid_vorticity_y = vorticity[1]
            self.pp.fluid_vorticity_z = vorticity[2]  

    def Test(self, model_part, benchmark_utils, target_saffman):
        self.target_saffman = target_saffman
        InitializeVariables(model_part, self.pp)
        benchmark_utils.ComputeHydrodynamicForces(model_part)

        for node in model_part.Nodes:
            saffman = node.GetSolutionStepValue(LIFT_FORCE)
                  
        error = Benchmark.ErrorMetric(saffman, self.target_saffman)
        
        if (error < self.saffman_tol):
            veredict = True
            
        else:
            veredict = False
            
        self.results = [self.target_saffman, saffman, error, veredict, self.description] 
        self.string_results = SaffmanBenchmark.ConvertToStrings(self.results)
        self.has_results = True  
        
#***************************************************************************************************************************
   
# BENCHMARKS

#***************************************************************************************************************************  

pp = ProjectParameters.Parameters()
model_part = ModelPart("OneBallModelPart")
AddVariables(model_part)

model_part_io_solid = ModelPartIO("hydrodynamic_forces", True)
model_part_io_solid.ReadModelPart(model_part)
model_part.SetBufferSize(2)
AddDofs(model_part)

benchmark_utils = BenchmarkUtils()

# Buoyancy
#***************************************************************************************************************************
#***************************************************************************************************************************

pressure_gradient = RandomVector()

#***************************************************************************************************************************
buoyancy_test_0 = BuoyancyBenchmark(pp, 0, 1, pressure_gradient, "Inactive")

buoyancy_target_0 = RandomVector(0)

buoyancy_test_0.Test(model_part, benchmark_utils, buoyancy_target_0)

#***************************************************************************************************************************
buoyancy_test_1 = BuoyancyBenchmark(pp, 1, 1, pressure_gradient, "Standard")

buoyancy_target_1 = Array3()
buoyancy_target_1[0] = - 4/3 * math.pi * pressure_gradient[0]
buoyancy_target_1[1] = - 4/3 * math.pi * pressure_gradient[1]
buoyancy_target_1[2] = - 4/3 * math.pi * pressure_gradient[2]

buoyancy_test_1.Test(model_part, benchmark_utils, buoyancy_target_1)

#***************************************************************************************************************************
buoyancy_test_2 = BuoyancyBenchmark(pp, 1, 2, pressure_gradient, "Weatherford: hydrostatic buoyancy")

buoyancy_target_2 = Array3()
buoyancy_target_2[0] = 0.0
buoyancy_target_2[1] = 0.0
buoyancy_target_2[2] = 9.81 * 4/3 * math.pi

buoyancy_test_2.Test(model_part, benchmark_utils, buoyancy_target_2)

BuoyancyBenchmark.PrintResults()

# Drag
#***************************************************************************************************************************
#***************************************************************************************************************************

viscosity = 10e-4
fluid_fraction = 0.5
reynolds = 0.0
sphericity = 1.0

#***************************************************************************************************************************

drag_test_0 = DragBenchmark(pp, 0, reynolds, fluid_fraction, sphericity, "Inactive")

drag_target_0 = RandomVector(0)

drag_test_0.Test(model_part, benchmark_utils, drag_target_0)

#***************************************************************************************************************************
reynolds = 1.0
drag_test_1 = DragBenchmark(pp, 1, reynolds, fluid_fraction, sphericity, "Stokes regime")

viscosity = drag_test_1.pp.kinematic_viscosity
density = drag_test_1.pp.fluid_density
radius = drag_test_1.pp.radius
slip_vel = Array3()
slip_vel[0] = drag_test_1.pp.fluid_velocity_x - drag_test_1.pp.velocity_x
slip_vel[1] = drag_test_1.pp.fluid_velocity_y - drag_test_1.pp.velocity_y
slip_vel[2] = drag_test_1.pp.fluid_velocity_z - drag_test_1.pp.velocity_z
drag_target_1 = 6 * math.pi * viscosity * radius * slip_vel

drag_test_1.Test(model_part, benchmark_utils, drag_target_1)

#***************************************************************************************************************************
drag_test_2 = DragBenchmark(pp, 5, reynolds, fluid_fraction, sphericity, "Newtonian regime")

# calculating target
slip_vel[0] = drag_test_2.pp.fluid_velocity_x - drag_test_2.pp.velocity_x
slip_vel[1] = drag_test_2.pp.fluid_velocity_y - drag_test_2.pp.velocity_y
slip_vel[2] = drag_test_2.pp.fluid_velocity_z - drag_test_2.pp.velocity_z
density = drag_test_2.pp.fluid_density
radius = drag_test_2.pp.radius
drag_target_2 = 0.5 * math.pi * radius ** 2 * density * Norm(slip_vel) * slip_vel
drag_target_2 *= 0.44 

drag_test_2.Test(model_part, benchmark_utils, drag_target_2)

#***************************************************************************************************************************
reynolds = 1.0
drag_test_3 = DragBenchmark(pp, 6, reynolds, fluid_fraction, sphericity, "Intermediate regime, Re = " + str(reynolds))

# calculating target
slip_vel[0] = drag_test_3.pp.fluid_velocity_x - drag_test_3.pp.velocity_x
slip_vel[1] = drag_test_3.pp.fluid_velocity_y - drag_test_3.pp.velocity_y
slip_vel[2] = drag_test_3.pp.fluid_velocity_z - drag_test_3.pp.velocity_z
density = drag_test_3.pp.fluid_density
radius = drag_test_3.pp.radius
drag_target_3 = 0.5 * math.pi * radius ** 2 * density * Norm(slip_vel) * slip_vel
drag_target_3 *= 24.0 / reynolds * (1.0 + 0.15 * math.pow(reynolds, 0.687))

drag_test_3.Test(model_part, benchmark_utils, drag_target_3)

#***************************************************************************************************************************
reynolds = 250
drag_test_31 = DragBenchmark(pp, 6, reynolds, fluid_fraction, sphericity, "Intermediate regime, Re = " + str(reynolds))

# calculating target
slip_vel[0] = drag_test_31.pp.fluid_velocity_x - drag_test_31.pp.velocity_x
slip_vel[1] = drag_test_31.pp.fluid_velocity_y - drag_test_31.pp.velocity_y
slip_vel[2] = drag_test_31.pp.fluid_velocity_z - drag_test_31.pp.velocity_z
density = drag_test_31.pp.fluid_density
radius = drag_test_31.pp.radius
drag_target_31 = 0.5 * math.pi * radius ** 2 * density * Norm(slip_vel) * slip_vel
drag_target_31 *= 24.0 / reynolds * (1.0 + 0.15 * math.pow(reynolds, 0.687))

drag_test_31.Test(model_part, benchmark_utils, drag_target_31)

#***************************************************************************************************************************
reynolds = 250
sphericity = 1.0
drag_test_4 = DragBenchmark(pp, 7, reynolds, fluid_fraction, sphericity, "Haider and Levenspiel, Re = " + str(reynolds) + ", sphericity = " + str(sphericity))

# calculating target
slip_vel[0] = drag_test_4.pp.fluid_velocity_x - drag_test_4.pp.velocity_x
slip_vel[1] = drag_test_4.pp.fluid_velocity_y - drag_test_4.pp.velocity_y
slip_vel[2] = drag_test_4.pp.fluid_velocity_z - drag_test_4.pp.velocity_z
density = drag_test_4.pp.fluid_density
radius = drag_test_4.pp.radius
drag_target_4 = 0.5 * math.pi * radius ** 2 * density * Norm(slip_vel) * slip_vel
A = math.exp(2.3288 - 6.4581 * sphericity + 2.4486 * sphericity ** 2)
B = 0.0964 + 0.5565 * sphericity
C = math.exp(4.905 - 13.8944 * sphericity + 18.4222 * sphericity ** 2 - 10.2599 * sphericity ** 3)
D = math.exp(1.4681 + 12.2584 * sphericity - 20.7322 * sphericity ** 2 + 15.8855 * sphericity ** 3)
drag_target_4 *= (24.0 / reynolds * (1.0 + A * math.pow(reynolds, B)) + C / (1 + D / reynolds))

drag_test_4.Test(model_part, benchmark_utils, drag_target_4)

#***************************************************************************************************************************
reynolds = 250
sphericity = 0.3
drag_test_41 = DragBenchmark(pp, 7, reynolds, fluid_fraction, sphericity, "Haider and Levenspiel, Re = " + str(reynolds) + ", sphericity = " + str(sphericity))

# calculating target
slip_vel[0] = drag_test_41.pp.fluid_velocity_x - drag_test_41.pp.velocity_x
slip_vel[1] = drag_test_41.pp.fluid_velocity_y - drag_test_41.pp.velocity_y
slip_vel[2] = drag_test_41.pp.fluid_velocity_z - drag_test_41.pp.velocity_z
density = drag_test_41.pp.fluid_density
radius = drag_test_41.pp.radius
drag_target_41 = 0.5 * math.pi * radius ** 2 * density * Norm(slip_vel) * slip_vel
A = math.exp(2.3288 - 6.4581 * sphericity + 2.4486 * sphericity ** 2)
B = 0.0964 + 0.5565 * sphericity
C = math.exp(4.905 - 13.8944 * sphericity + 18.4222 * sphericity ** 2 - 10.2599 * sphericity ** 3)
D = math.exp(1.4681 + 12.2584 * sphericity - 20.7322 * sphericity ** 2 + 15.8855 * sphericity ** 3)
drag_target_41 *= (24.0 / reynolds * (1.0 + A * math.pow(reynolds, B)) + C / (1 + D / reynolds))

drag_test_41.Test(model_part, benchmark_utils, drag_target_41)

DragBenchmark.PrintResults()

# Virtual Mass
#***************************************************************************************************************************
#***************************************************************************************************************************

acceleration_number = 1.0
virtual_mass_test_0 = VirtualMassBenchmark(pp, 0, acceleration_number, fluid_fraction, "Inactive")

virtual_mass_target_0 = RandomVector(0)

virtual_mass_test_0.Test(model_part, benchmark_utils, virtual_mass_target_0)

#***************************************************************************************************************************

virtual_mass_test_1 = VirtualMassBenchmark(pp, 1, acceleration_number, fluid_fraction, "Stokes regime")
slip_accel = Array3()
slip_accel[0] = virtual_mass_test_1.pp.fluid_acceleration_x - virtual_mass_test_1.pp.accel_x
slip_accel[1] = virtual_mass_test_1.pp.fluid_acceleration_y - virtual_mass_test_1.pp.accel_y
slip_accel[2] = virtual_mass_test_1.pp.fluid_acceleration_z - virtual_mass_test_1.pp.accel_z
volume = 4.0 / 3.0  * math.pi * virtual_mass_test_1.pp.radius ** 3
virtual_mass_coeff = 0.5;
virtual_mass_target_1 = virtual_mass_coeff * virtual_mass_test_1.pp.fluid_density * volume * slip_accel

virtual_mass_test_1.Test(model_part, benchmark_utils, virtual_mass_target_1)

#***************************************************************************************************************************

fluid_fraction = 0.6
virtual_mass_test_2 = VirtualMassBenchmark(pp, 2, acceleration_number, fluid_fraction, "Zuber, fluid fraction = " + str(fluid_fraction))
slip_accel = Array3()
slip_accel[0] = virtual_mass_test_2.pp.fluid_acceleration_x - virtual_mass_test_2.pp.accel_x
slip_accel[1] = virtual_mass_test_2.pp.fluid_acceleration_y - virtual_mass_test_2.pp.accel_y
slip_accel[2] = virtual_mass_test_2.pp.fluid_acceleration_z - virtual_mass_test_2.pp.accel_z
volume = 4.0 / 3.0  * math.pi * virtual_mass_test_2.pp.radius ** 3
virtual_mass_coeff = 0.5 + 1.5 * (1 - fluid_fraction);
virtual_mass_target_2 = virtual_mass_coeff * virtual_mass_test_2.pp.fluid_density * volume * slip_accel

virtual_mass_test_2.Test(model_part, benchmark_utils, virtual_mass_target_2)

#***************************************************************************************************************************

acceleration_number = 3.6
virtual_mass_test_3 = VirtualMassBenchmark(pp, 3, acceleration_number, fluid_fraction, "Odar and Hamilton, acceleration number = " + str(acceleration_number))
slip_accel = Array3()
slip_accel[0] = virtual_mass_test_3.pp.fluid_acceleration_x - virtual_mass_test_3.pp.accel_x
slip_accel[1] = virtual_mass_test_3.pp.fluid_acceleration_y - virtual_mass_test_3.pp.accel_y
slip_accel[2] = virtual_mass_test_3.pp.fluid_acceleration_z - virtual_mass_test_3.pp.accel_z
volume = 4.0 / 3.0  * math.pi * virtual_mass_test_3.pp.radius ** 3
virtual_mass_coeff = 0.5
virtual_mass_coeff *= 2.1 - 0.132 / (acceleration_number ** 2 + 0.12);
virtual_mass_target_3 = virtual_mass_coeff * virtual_mass_test_3.pp.fluid_density * volume * slip_accel

virtual_mass_test_3.Test(model_part, benchmark_utils, virtual_mass_target_3)

#***************************************************************************************************************************

fluid_fraction = 0.6
acceleration_number = 3.6
virtual_mass_test_4 = VirtualMassBenchmark(pp, 4, acceleration_number, fluid_fraction, "Odar + Zuber; fluid frac. = " + str(fluid_fraction) + ", acc. num. = " + str(acceleration_number))
slip_accel = Array3()
slip_accel[0] = virtual_mass_test_4.pp.fluid_acceleration_x - virtual_mass_test_4.pp.accel_x
slip_accel[1] = virtual_mass_test_4.pp.fluid_acceleration_y - virtual_mass_test_4.pp.accel_y
slip_accel[2] = virtual_mass_test_4.pp.fluid_acceleration_z - virtual_mass_test_4.pp.accel_z
volume = 4.0 / 3.0  * math.pi * virtual_mass_test_4.pp.radius ** 3
virtual_mass_coeff = 0.5 + 1.5 * (1 - fluid_fraction)
virtual_mass_coeff *= 2.1 - 0.132 / (acceleration_number ** 2 + 0.12);
virtual_mass_target_4 = virtual_mass_coeff * virtual_mass_test_4.pp.fluid_density * volume * slip_accel

virtual_mass_test_4.Test(model_part, benchmark_utils, virtual_mass_target_4)

VirtualMassBenchmark.PrintResults()

# Saffman
#***************************************************************************************************************************
#***************************************************************************************************************************

reynolds = 1.0
reynolds_shear = 1.0
saffman_test_0 = SaffmanBenchmark(pp, 0, reynolds, reynolds_shear, "Inactive")

saffman_target_0 = RandomVector(0)

saffman_test_0.Test(model_part, benchmark_utils, saffman_target_0)

SaffmanBenchmark.PrintResults()
