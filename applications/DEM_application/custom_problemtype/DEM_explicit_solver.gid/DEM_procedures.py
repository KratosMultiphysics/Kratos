from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
#from KratosMultiphysics.MetisApplication import *

from DEM_explicit_solver_var import *
from pressure_script import *

from numpy import *

#from KratosMultiphysics.mpi import *

def AddMpiVariables(model_part):
    
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    model_part.AddNodalSolutionStepVariable(INTERNAL_ENERGY)
    model_part.AddNodalSolutionStepVariable(OSS_SWITCH)
    
def PerformInitialPartition(model_part,model_part_io_solid,input_file_name):
    
    domain_size = 3
    
    print "("+str(mpi.rank)+","+str(mpi.size)+")"+"before performing the division"
    number_of_partitions = mpi.size #we set it equal to the number of processors
    if mpi.rank == 0 :
        print "("+str(mpi.rank)+","+str(mpi.size)+")"+"start partition process"
        partitioner = MortonDivideInputToPartitionsProcess(model_part_io_solid, number_of_partitions, domain_size);
        partitioner.Execute()

    print "("+str(mpi.rank)+","+str(mpi.size)+")"+"division performed"
    mpi.world.barrier()

    MPICommSetup = SetMPICommunicatorProcess(model_part)
    MPICommSetup.Execute()

    print "("+str(mpi.rank)+","+str(mpi.size)+")"+"Comunicator Set"

    print "("+str(mpi.rank)+","+str(mpi.size)+")"+"Reading: "+input_file_name+"_"+str(mpi.rank)

    my_input_filename = input_file_name+"_"+str(mpi.rank)
    model_part_io_solid = ModelPartIO(my_input_filename)
    
    return model_part_io_solid
    


def ProcModelData(solid_model_part,solver):

  # Previous Calculations.

  Model_Data = open('Model_Data.txt','w')
  
  #mean radius, and standard deviation:

  i = 0
  sum_radi = 0
  sum_squared = 0
  for node in solid_model_part.Nodes:
  
    sum_radi += node.GetSolutionStepValue(RADIUS)
    sum_squared += node.GetSolutionStepValue(RADIUS)**2
    i+=1

  mean=sum_radi/i
  var =sum_squared/i-mean**2
  if(abs(var)<1e-05):
    var=0
  std_dev=var**0.5
  
  Model_Data.write("Radius Mean: "+str(mean)+'\n')
  Model_Data.write("Std Deviation: "+str(std_dev)+'\n')
  Model_Data.write('\n')
  
  Total_Particles     = len(solid_model_part.Nodes)
  Total_Contacts      = solver.model_part.ProcessInfo.GetValue(TOTAL_CONTACTS)/2
  Coordination_Number    = 1.0*(Total_Contacts*2)/Total_Particles
  
  Model_Data.write("Total Number of Particles: "+str(Total_Particles)+'\n')
  Model_Data.write("Total Number of Contacts: "+str(Total_Contacts)+'\n')
  Model_Data.write("Coordination Number NC: "+str(Coordination_Number)+'\n')
  Model_Data.write('\n')
  
  Model_Data.write("Volume Elements: "+str(mass_elements)+'\n')
  
  Model_Data.close()    

sup_layer = list()
inf_layer = list()
fix_particles = list()
force_measurement = list()
special_selection = list()
others = list()
   
def ProcListDefinition(model_part,solver):
  
  # Defining lists (FOR COMPRESSION TESTS)
  
  for node in model_part.Nodes:
    if (node.GetSolutionStepValue(GROUP_ID)==1):
      sup_layer.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==2):
      inf_layer.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==3):
      fix_particles.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==4):
      force_measurement.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==5):
      special_selection.append(node)
    else:
      others.append(node)

    
def ProcGiDSolverTransfer(model_part,solver):
    
    if (Integration_Scheme == 'forward_euler'):
        time_scheme = FowardEulerScheme()
    elif (Integration_Scheme == 'mid_point_rule'):
        time_scheme = MidPointScheme()
    elif (Integration_Scheme == 'const_average_acc'):
        time_scheme = ConstAverageAccelerationScheme()
    else:
        print('scheme not defined')
        
    if (NormalForceCalculation == "Linear"):
        force_calculation_type_id = 0
    elif (NormalForceCalculation == "Hertz"):
        force_calculation_type_id = 1

    if(DampId == "ViscDamp"):
        damp_id = 1
    else:
        damp_id = 0
        
    solver.damp_id=damp_id

    if(RotaDampId == "LocalDamp"):
        rota_damp_id = 1
    elif(RotaDampId == "RollingFric"):
        rota_damp_id = 2
    else:
        rota_damp_id = 0
        
    solver.rota_damp_id=rota_damp_id
    
    solver.magic_factor = MagicFactor


    gravity = Vector(3)
    gravity[0] = gravity_x
    gravity[1] = gravity_y
    gravity[2] = gravity_z

    solver.gravity=gravity
    
    m_search_radius_extension = search_radius_extension

    #options for the solver
    
    if(VirtualMassOption == "ON"):
        solver.virtual_mass_OPTION=1 #xapuza
        
    solver.nodal_mass_coeff=VirtualMassCoefficient    

    solver.final_time = final_time
    
    if(DeltaOption=="OFF"):
        m_search_radius_extension = 0.0;

    solver.time_scheme=time_scheme
    solver.force_calculation_type_id=force_calculation_type_id

    if (CriticalTimeOption =="ON"):
        solver.critical_time_OPTION=1; #xapuza

    if(DeltaOption =="ON"):
        solver.delta_OPTION=True

    if(ContinuumOption =="ON"):
        solver.continuum_simulating_OPTION=True
      
        if(ContactMeshOption =="ON"):
            solver.contact_mesh_OPTION=1  #xapuza
               
        if(ConcreteTestOption =="ON"):
           solver.concrete_test_OPTION=1  #xapuza 
           
           if(TriaxialOption =="ON"):
           
             solver.initial_pressure_time = InitialTime
             solver.time_increasing_ratio = IncreasingTemporaily
           
      
        if(FailureCriterionOption =="Mohr-Coulomb"):
            solver.failure_criterion_OPTION=1 
        elif(FailureCriterionOption =="Uncoupled"):
            solver.failure_criterion_OPTION=2
        
        solver.tau_zero       = TauZero
        solver.sigma_max      = SigmaMax
        solver.sigma_min      = SigmaMin
        solver.internal_fricc = InternalFricc
      
    solver.search_radius_extension = m_search_radius_extension

    if(RotationOption =="ON"):
        solver.rotation_OPTION=1  #xapuza
    if(TrihedronOption =="ON"):
        solver.trihedron_OPTION=1  #xapuza 
    if(RotationalSpringOption =="ON"):
        solver.rotation_spring_OPTION=1  #xapuza
      
    solver.safety_factor = dt_safety_factor #for critical time step calculation 

    #Prints
    if (print_export_id=="1"):
      solver.print_export_id =1
      
    if (print_export_skin_sphere=="1"): 
      solver.print_export_skin_sphere = 1
      
    
    
    # global variable settings

    if(GlobalVariablesOption =="ON"):
        solver.global_variables_OPTION = 1  #xapuza

    solver.global_kn    = global_kn
    solver.global_kt    = global_kt
    solver.global_kr    = global_kr
    solver.global_rn    = global_rn
    solver.global_rt    = global_rt
    solver.global_rr    = global_rr
    
    solver.global_fri_ang = global_fri_ang

    # time settings
    
    solver.delta_time = max_time_step

    # bounding box

    n_step_destroy_distant = search_step      # how many steps between elimination of distant particles?
    solver.n_step_search   = search_step

    extra_radius = 0.0
    max_radius = 0.0
    min_radius = 0.0
    first_it = True

    #calculation of search radius
    for node in model_part.Nodes:
          
      rad = node.GetSolutionStepValue(RADIUS)
      if rad > max_radius:  
          max_radius = rad
      if first_it == True:
          min_radius = rad
          first_it = False
      if rad < min_radius:  
          min_radius = rad
          
    if(BoundingBoxOption =="ON"):
      solver.bounding_box_OPTION=1  #xapuza
      
    extra_radius = 2.5 * max_radius
    prox_tol = 0.000001 * min_radius  #currently not in use.
    m_bounding_box_enlargement_factor = max(1.0 + extra_radius, bounding_box_enlargement_factor)

    solver.enlargement_factor = m_bounding_box_enlargement_factor
    
    Pressure = ConfinementPressure*1e6 #Mpa
    
    if(Pressure!=0):
      
      solver.external_pressure = 1

def ProcSkinAndPressure(model_part,solver):
    
    #Defining list of skin particles (For a test tube of height 30 cm and diameter 15 cm)
    
    Pressure = ConfinementPressure*1e6 #Mpa
    
    print(" ")
    print(solver.external_pressure)
    print("")
    SKIN = list()  
    LAT = list()
    BOT = list()
    TOP = list()
    XLAT = list()  #only lat, not the corner ones
    XTOP = list()  #only top, not corner ones...
    XBOT = list()
    XTOPCORNER = list()
    XBOTCORNER = list()
 
    total_cross_section = 0.0
    
    #Cylinder dimensions
    
    h   = 0.3
    d   = 0.15
    eps = 2
    
    surface = 2*(3.141592*d*d*0.25)+(3.141592*d*h)
    
    top_pressure = 0.0
    bot_pressure = 0.0
      
    #SKIN DETERMINATION
    
    for element in model_part.Elements:
    
      element.SetValue(SKIN_SPHERE,0)
      if (element.GetValue(PREDEFINED_SKIN)==1):
		
        element.SetValue(SKIN_SPHERE,1)
      
      node = element.GetNode(0)
      r = node.GetSolutionStepValue(RADIUS,0)
      x = node.X
      y = node.Y
      z = node.Z

      cross_section = 3.141592*r*r

      if ( (x*x+z*z)>=((d/2-eps*r)*(d/2-eps*r)) ): 
      
        element.SetValue(SKIN_SPHERE,1)     
        total_cross_section = total_cross_section + cross_section
      
        LAT.append(node)
            
        if ( (y>eps*r ) and (y<(h-eps*r)) ) :
        
          SKIN.append(element)
        
          XLAT.append(node)
    
      if ( (y<=eps*r ) or (y>=(h-eps*r)) ): 

          element.SetValue(SKIN_SPHERE,1)
        
          SKIN.append(element)
        
          if ( y<=eps*r ):

              BOT.append(node)

          elif ( y>=(h-eps*r) ):

              TOP.append(node)

          if ( (x*x+z*z) >= ((d/2-eps*r)*(d/2-eps*r) ) ) :
         
              if ( y>h/2 ):

                  XTOPCORNER.append(node)
                
              else:

                  XBOTCORNER.append(node)
          else:

              if ( y<=eps*r ):
                
                  XBOT.append(node)
                
              elif ( y>=(h-eps*r) ):
                    
                  XTOP.append(node)
                                    
    if ( (TriaxialOption == "ON") and (Pressure != 0.0) ):
 
      ApplyPressure(Pressure,model_part,solver,SKIN,BOT,TOP,LAT,XLAT,XBOT,XBOTCORNER,XTOP,XTOPCORNER) 
      print("End Applying Imposed Forces")
     
def ProcPrintingVariables(gid_io,solid_model_part,contact_model_part,time):
  
	if (print_displacement=="1"):
	  gid_io.WriteNodalResults(DISPLACEMENT, contact_model_part.Nodes, time, 0)       
	if (print_velocity=="1"):
	  gid_io.WriteNodalResults(VELOCITY, contact_model_part.Nodes, time, 0)	  
	if (print_rhs=="1"):
	  gid_io.WriteNodalResults(RHS, contact_model_part.Nodes, time, 0)       
	if (print_applied_forces=="1"):
	  gid_io.WriteNodalResults(APPLIED_FORCE, contact_model_part.Nodes, time, 0)       
	if (print_total_forces=="1"):	  
	  gid_io.WriteNodalResults(TOTAL_FORCES, contact_model_part.Nodes, time, 0)	  
	if (print_damp_forces=="1"):
	  gid_io.WriteNodalResults(DAMP_FORCES, contact_model_part.Nodes, time, 0)        
	if (print_radius=="1"):
	  gid_io.WriteNodalResults(RADIUS, contact_model_part.Nodes, time, 0)       
	if (print_particle_cohesion=="1"):
	  gid_io.WriteNodalResults(PARTICLE_COHESION, contact_model_part.Nodes, time, 0)       
	if (print_particle_tension=="1"):
	  gid_io.WriteNodalResults(PARTICLE_TENSION, contact_model_part.Nodes, time, 0)
	if (print_group_id=="1"):
	  gid_io.WriteNodalResults(GROUP_ID, contact_model_part.Nodes, time, 0)
	if (print_export_id=="1"):
	  gid_io.WriteNodalResults(EXPORT_ID, contact_model_part.Nodes, time, 0)
	if (print_export_particle_failure_id=="1"):
	  gid_io.WriteNodalResults(EXPORT_PARTICLE_FAILURE_ID, contact_model_part.Nodes, time, 0)
	if (print_export_skin_sphere=="1"):
	  gid_io.WriteNodalResults(EXPORT_SKIN_SPHERE, contact_model_part.Nodes, time, 0)
	gid_io.WriteNodalResults(DEM_STRESS_XX, contact_model_part.Nodes, time, 0)
	gid_io.WriteNodalResults(DEM_STRESS_XY, contact_model_part.Nodes, time, 0)
	gid_io.WriteNodalResults(DEM_STRESS_XZ, contact_model_part.Nodes, time, 0)
	gid_io.WriteNodalResults(DEM_STRESS_YX, contact_model_part.Nodes, time, 0)
	gid_io.WriteNodalResults(DEM_STRESS_YY, contact_model_part.Nodes, time, 0)
	gid_io.WriteNodalResults(DEM_STRESS_YZ, contact_model_part.Nodes, time, 0)
	gid_io.WriteNodalResults(DEM_STRESS_ZX, contact_model_part.Nodes, time, 0)
	gid_io.WriteNodalResults(DEM_STRESS_ZY, contact_model_part.Nodes, time, 0)
	gid_io.WriteNodalResults(DEM_STRESS_ZZ, contact_model_part.Nodes, time, 0)
    
  #Aixo sempre per que si no hi ha manera de debugar
  #gid_io.WriteNodalResults(PARTITION_INDEX, contact_model_part.Nodes, time, 0)
  #gid_io.WriteNodalResults(INTERNAL_ENERGY, contact_model_part.Nodes, time, 0)

	if (ContactMeshOption == "ON"): ##xapuza
	  if (print_local_contact_force_low=="1"):
		  gid_io.PrintOnGaussPoints(LOCAL_CONTACT_FORCE_LOW,contact_model_part,time)
	  if (print_local_contact_force_high=="1"):
		  gid_io.PrintOnGaussPoints(LOCAL_CONTACT_FORCE_HIGH,contact_model_part,time)
	  if (print_contact_failure=="1"): 
		  gid_io.PrintOnGaussPoints(CONTACT_FAILURE,contact_model_part,time)	 
	  if (print_failure_criterion_state=="1"):
		  gid_io.PrintOnGaussPoints(FAILURE_CRITERION_STATE,contact_model_part,time)  	    
	  if (print_contact_tau=="1"):
		  gid_io.PrintOnGaussPoints(CONTACT_TAU,contact_model_part,time)
	  if (print_contact_sigma=="1"):
		  gid_io.PrintOnGaussPoints(CONTACT_SIGMA,contact_model_part,time)

	if (RotationOption == "ON"): ##xapuza
	  if (print_angular_velocity=="1"):
		  gid_io.WriteNodalResults(ANGULAR_VELOCITY, contact_model_part.Nodes, time, 0)
	  if (print_particle_moment=="1"):
		  gid_io.WriteNodalResults(PARTICLE_MOMENT, contact_model_part.Nodes, time, 0)
	  if (print_euler_angles=="1"):
		  gid_io.WriteLocalAxesOnNodes(EULER_ANGLES, contact_model_part.Nodes, time, 0)

	gid_io.Flush()
	sys.stdout.flush()


 
    
