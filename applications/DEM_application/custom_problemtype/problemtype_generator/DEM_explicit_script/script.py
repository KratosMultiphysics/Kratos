#import cProfile

import DEM_explicit_solver_var
import time as timer
#import matplotlib
#from numpy import *
#from pylab import *

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

#defining a model part for the solid part
solid_model_part = ModelPart("SolidPart");  
#############################################

##importing the solvers needed & defining rotation cases
rotation =False #default vaulue, changed on rotation cases.
contact_type = DEM_explicit_solver_var.ContactType
if (contact_type == 'spring_circle'):
    import circle_spring_explicit_solver as explicit_solver
elif (contact_type == 'hertzian_circle'):
    import circle_hertzian_explicit_solver as explicit_solver
elif (contact_type == 'spring_sphere'):
    import sphere_spring_explicit_solver as explicit_solver
elif (contact_type == 'hertzian_sphere'):
    import sphere_hertzian_explicit_solver as explicit_solver
elif (contact_type == 'rotating_spring_sphere'):
    import rotating_sphere_spring_explicit_solver as explicit_solver
    rotation = True
elif (contact_type == 'rotating_hertzian_sphere'):
    rotation = True
    import rotating_sphere_hertzian_explicit_solver as explicit_solver

#cProfile.run(explicit_solver.AddVariables(solid_model_part))
explicit_solver.AddVariables(solid_model_part)
#introducing input file name
input_file_name = DEM_explicit_solver_var.problem_name
#reading the solid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
#print "aqui hi arribem1"
##selecting output format
#cProfile.run(gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions))
gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
#print "aqui hi arribem2"
#gid_io.ReadModelPart(solid_model_part) 
model_part_io_solid = ModelPartIO(input_file_name)
#print "aqui hi arribem3"
model_part_io_solid.ReadModelPart(solid_model_part)
#print "aqui hi arribem4"
#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
solid_model_part.SetBufferSize(2)
#print "aqui hi arribem4"

##adding dofs
explicit_solver.AddDofs(solid_model_part)
    
#constructing the solver
#constants
ruta = './resultados.py'
p = open(ruta,'w')

solver_strategy = DEM_explicit_solver_var.Solver
if (solver_strategy == 'forward_euler'):
    solver_id = 1
elif (solver_strategy == 'runge_kutta_4 '):
    solver_id = 2
else :
    solver_id = 3
#print "aqui hi arribem5"
gravity = Vector(3)
gravity[0] = DEM_explicit_solver_var.gravity_x
gravity[1] = DEM_explicit_solver_var.gravity_y
gravity[2] = DEM_explicit_solver_var.gravity_z
final_time = DEM_explicit_solver_var.max_time
output_dt = DEM_explicit_solver_var.output_dt
safety_factor = DEM_explicit_solver_var.dt_safety_factor
number_of_inital_steps = DEM_explicit_solver_var.number_of_inital_steps
initial_time_step = DEM_explicit_solver_var.initial_time_step
introduced_max_dt =  DEM_explicit_solver_var.max_time_step
min_dt = DEM_explicit_solver_var.min_time_step

#settings to be changed
n_step_estimation = 1
n_step_destroy_distant = 1      # how many steps between elimination of distant particles?
n_step_search = 1
bounding_box_enlargement_factor = 2.0    # by what factor do we want to enlarge the strict bounding box
min_risk_factor = 1.0 / safety_factor

#mesh to be printed
if(DEM_explicit_solver_var.print_layers == False): # true by default
    mesh_name = 0.0
    gid_io.InitializeMesh( mesh_name)
    gid_io.WriteSphereMesh( solid_model_part.GetMesh() )
    gid_io.FinalizeMesh()
    gid_io.Flush()
    gid_io.InitializeResults(mesh_name, (solid_model_part).GetMesh());

#initializations
time = 0.0
time_old_print = 0.0
step = 0
dt = initial_time_step
risk_factor = min_risk_factor
search_radius = 0.0
max_radius = 0.0
min_radius = 0.0
first_it = True
#calculation of search radius
for node in solid_model_part.Nodes:
    rad = node.GetSolutionStepValue(RADIUS)
    if rad > max_radius:  
        max_radius = rad
    if first_it == True:
        min_radius = rad
        first_it = False
    if rad < min_radius:  
        min_radius = rad

search_radius = 2.5 * max_radius
prox_tol = 0.000001 * min_radius
bounding_box_enlargement_factor = max(1.0 + search_radius, bounding_box_enlargement_factor)

print ('The search radius is: ')+str(search_radius)
#print "aqui hi arribem 6"
solver = explicit_solver.ExplicitSolver(solver_id, dt, search_radius, n_step_search, prox_tol, solid_model_part)
#print "aqui hi arribem 7"
dt = solver.EstimateTimeStep(risk_factor, dt, introduced_max_dt, min_dt)
#print "aqui hi arribem 8"
list_of_particles_pointers = solver.GetParticlePointers()
#print "aqui hi arribem 9"
solver.Calculate_Model_Surrounding_Bounding_Box(list_of_particles_pointers, solid_model_part, bounding_box_enlargement_factor)
print ('The initial estimation of the time step is: Dt = ')+str(dt)
#print "explicit solver created"
solver.EvolveSystem(dt, gravity) #also creates particles
#print 'finished initializing solver and creating particles'
current_pr_time = timer.clock()
current_real_time = timer.time()
print ('Calculation starts at instant: ') + str(current_pr_time)
while(time < final_time):

    if(step < number_of_inital_steps):
        max_dt = initial_time_step
    else:
        max_dt = introduced_max_dt
        #progressively increment the safety factor
        #in the steps that follow a reduction of it
        risk_factor = risk_factor * 0.8
        if(risk_factor < min_risk_factor):
            risk_factor = min_risk_factor

    if ((step + 1) % n_step_destroy_distant == 0):
        solver.Destroy_Particles(list_of_particles_pointers, solid_model_part)

    time = time + dt
    solid_model_part.CloneTimeStep(time)
    solver.Solve(dt, gravity)
##############     GiD IO        ################################################################################
    time_to_print = time - time_old_print
    print str(time)
    if(time_to_print >= DEM_explicit_solver_var.output_dt):
        if(DEM_explicit_solver_var.print_layers == True):
            gid_io.InitializeMesh(time);
            gid_io.WriteSphereMesh(solid_model_part.GetMesh());
            gid_io.FinalizeMesh();
        print(time)
        gid_io.InitializeResults(time, solid_model_part.GetMesh());   
        gid_io.WriteNodalResults(VELOCITY, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(FORCE, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(RADIUS, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(NUMBER_OF_NEIGHBOURS, solid_model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(IS_STRUCTURE, solid_model_part.Nodes, time, 0)

        if (rotation == True):
            gid_io.WriteNodalResults(ANGULAR_VELOCITY, solid_model_part.Nodes, time, 0)
            gid_io.WriteNodalResults(MOMENT, solid_model_part.Nodes, time, 0)
        gid_io.Flush()
        if(DEM_explicit_solver_var.print_layers == True):
            gid_io.FinalizeResults()    
	time_old_print = time
    step += 1
    
p.close()
print 'Calculation ends at instant: ' + str(timer.time())
elapsed_pr_time = timer.clock() - current_pr_time
elapsed_real_time = timer.time() - current_real_time
print 'Calculation ends at processing time instant: ' + str(timer.clock())
print 'Elapsed processing time: ' + str(elapsed_pr_time)
print 'Elapsed real time: ' + str(elapsed_real_time)
if(DEM_explicit_solver_var.print_layers == False):
    gid_io.FinalizeResults()    
        

