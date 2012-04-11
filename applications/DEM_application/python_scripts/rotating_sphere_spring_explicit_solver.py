# -*- coding: utf-8 -*-
#importing the Kratos Library
import time as timer

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(NUMBER_OF_NEIGHBOURS)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(ROTATION)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(MOMENT)
    model_part.AddNodalSolutionStepVariable(RADIUS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)
    model_part.AddNodalSolutionStepVariable(PARTICLE_STATIC_FRICTION_COEF)
    model_part.AddNodalSolutionStepVariable(PARTICLE_DYNAMIC_FRICTION_COEF)
    model_part.AddNodalSolutionStepVariable(PARTICLE_STIFFNESS)
    model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS)
    model_part.AddNodalSolutionStepVariable(POISSON_RATIO)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MASS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_COEF_RESTITUTION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_ZETA)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)

    print "variables for the explicit solver added correctly"

def AddDofs(model_part):
    print "Dofs added"
 
class ExplicitSolver:
    
    def __init__(self, solver_id, Dt, search_radius, n_step_search, prox_tol, model_part):
 
        #data of the problem
        self.particle_destructor_and_constructor = Rotating_Spheres_Creator_Destructor()
        self.solver_id = solver_id
	self.step = 0;
        self.time_step = Dt;
        self.use_parallel_distance_calculation = False
        self.n_step_search = n_step_search
        self.acc_time_calc_forces = 0.0
        self.acc_time_motion = 0.0
        self.acc_time_search = 0.0

        #first neighbours search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.explicit_solver_object = Rotating_Spheres_Explicit_Solver(solver_id, search_radius, prox_tol, model_part)
        self.explicit_solver_object.Search_Neighbours()
        self.explicit_solver_object.Update_Contacting_Neighbours()
  
    def EvolveSystem(self, Dt, gravity):
 	#calculate forces
        before_forces_time = timer.time()
        self.explicit_solver_object.Calculate_Forces(Dt, gravity)
        elapsed_time_forces = timer.time() - before_forces_time
        self.acc_time_calc_forces += elapsed_time_forces
#        print 'acumulated computing time for the forces: ', self.acc_time_calc_forces

	#evolve motion
        before_motion_time = timer.time()
        self.explicit_solver_object.Evolve_Motion(Dt, gravity)
        elapsed_time_motion = timer.time() - before_motion_time
        self.acc_time_motion += elapsed_time_motion

	#neighbours search
        before_search_time = timer.time()
        if (self.step % self.n_step_search == 0):
            self.explicit_solver_object.Search_Neighbours()
            self.explicit_solver_object.Update_Contacting_Neighbours()
        elapsed_time_search = timer.time() - before_search_time
        self.acc_time_search += elapsed_time_search
#        print 'acumulated computing time for the search: ', self.acc_time_search

    def Solve(self, Dt, gravity):
	self.EvolveSystem(Dt, gravity)
	self.step += 1

    def EstimateTimeStep(self, risk_factor, Dt, max_Dt, min_Dt):
        dt = self.explicit_solver_object.Estimate_Time_Step(risk_factor);
        print 'the computed maximum time step is '
        print dt
        if(dt < min_Dt):
            dt = min_Dt
            print 'Estimated Dt smaller than min_Dt = ', min_Dt
        
        if(dt > max_Dt):
            dt = max_Dt
	    print 'Estimated Dt greater than min_Dt = ', max_Dt

        return dt

    def GetParticlePointers(self):
        return self.explicit_solver_object.Get_List_Of_Particle_Pointers()

    def Calculate_Model_Surrounding_Bounding_Box(self, particles_pointers, model_part, enlargement_factor):
        self.particle_destructor_and_constructor.Calculate_Surrounding_Bounding_Box(particles_pointers, model_part, enlargement_factor)

    def Destroy_Particles(self, particles_pointers, model_part):
        self.particle_destructor_and_constructor.Destroy_Distant_Particles(particles_pointers, model_part)

    


        
        

