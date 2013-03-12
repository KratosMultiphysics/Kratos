from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
#from KratosMultiphysics.MetisApplication import * #CARLOS

from DEM_explicit_solver_var import *
from pressure_script import *

from numpy import *

#from KratosMultiphysics.mpi import * #CARLOS


# GLOBAL VARIABLES OF THE SCRIPT
#Defining list of skin particles (For a test tube of height 30 cm and diameter 15 cm)
    
sup_layer_fm = list()
inf_layer_fm = list()
sup_plate_fm = list()
inf_plate_fm = list()
special_selection = list()
others = list()    
SKIN = list()  
LAT = list()
BOT = list()
TOP = list()
XLAT = list()  #only lat, not the corner ones
XTOP = list()  #only top, not corner ones...
XBOT = list()
XTOPCORNER = list()
XBOTCORNER = list()

   
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


def ProcListDefinition(model_part,solver):
  
  # Defining lists (FOR COMPRESSION TESTS)
  
  for node in model_part.Nodes:
    if (node.GetSolutionStepValue(GROUP_ID)==1):      #reserved for speciment particles with imposed displacement and strain-stress measurement (superior). Doesn't recive pressure
      sup_layer_fm.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==2):    #reserved for speciment particles with imposed displacement and strain-stress measurement (superior). Doesn't recive pressure
      inf_layer_fm.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==3):    #reserved for auxiliar strain-stress measurement plate (superior)
      sup_plate_fm.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==4):    #reserved for auxiliar strain-stress measurement plate (inferior)
      inf_plate_fm.append(node)
    elif (node.GetSolutionStepValue(GROUP_ID)==5):
      special_selection.append(node)
    else:
      others.append(node)
    
  return (sup_layer_fm, inf_layer_fm, sup_plate_fm, inf_plate_fm)

    
def ProcGiDSolverTransfer(model_part,solver):
    
    if (Integration_Scheme == 'forward_euler'):
        time_scheme = ForwardEulerScheme()
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

    if(NormalDampId == "ViscDamp"):
      if(TangentialDampId == "ViscDamp"):
        damp_id = 11
      else:
        damp_id = 10
    else:
      if(TangentialDampId == "ViscDamp"):
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
    
    m_search_radius_extension = 0.0
    m_amplified_continuum_search_radius_extension = 0.0

    #options for the solver
    
    if(VirtualMassOption == "ON"):
        solver.virtual_mass_OPTION=1 #xapuza
        
    solver.nodal_mass_coeff=VirtualMassCoefficient    

    solver.final_time = final_time
    
    if(DeltaOption=="ON"):
        m_search_radius_extension = search_radius_extension;
        if(ContinuumOption=="ON"):
           m_amplified_continuum_search_radius_extension = amplified_continuum_search_radius_extension;

    solver.time_scheme=time_scheme
    solver.force_calculation_type_id=force_calculation_type_id

    if (CriticalTimeOption =="ON"):
        solver.critical_time_OPTION=1; #xapuza

    if(DeltaOption =="ON"):
        solver.delta_OPTION=True

    if(ContinuumOption =="ON"):
        solver.continuum_simulating_OPTION=True
        
        if(NonLinearOption =="ON"):
           solver.Non_Linear_Option=1
           solver.C1 = C1
           solver.C2 = C2
           solver.N1 = N1
           solver.N2 = N2
        
        if(StressStrainOperations =="ON"): #xapuza
          solver.stress_strain_operations = 1#xapuza
 
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
    solver.amplified_continuum_search_radius_extension = m_amplified_continuum_search_radius_extension

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
      
    if (print_radial_displacement=="1"): 
      solver.print_radial_displacement = 1


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
    
    # global variable settings

    if(LimitSurfaceOption =="ON"):
        solver.limit_surface_OPTION = 1  #xapuza

    surface_normal_dir = Vector(3)
    surface_normal_dir[0] = surface_normal_dir_x
    surface_normal_dir[1] = surface_normal_dir_y    
    surface_normal_dir[2] = surface_normal_dir_z
    solver.surface_normal_dir = surface_normal_dir    
    surface_point_coor = Vector(3)
    surface_point_coor[0] = surface_point_coor_x
    surface_point_coor[1] = surface_point_coor_y
    surface_point_coor[2] = surface_point_coor_z
    solver.surface_point_coor = surface_point_coor
    solver.surface_friction_angle = surface_friction_angle

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
    
    if ( TriaxialOption =="ON" ):
       Pressure = ConfinementPressure*1e6 #Mpa
    else:
       Pressure = 0.0
    
    if(Pressure!=0):
      
      solver.external_pressure = 1
   
    if(FixVelocities =="ON"):
        solver.fix_velocities = 1  #xapuza
    solver.time_step_percentage_fix_velocities = TimePercentageFixVelocities   
    
    return Pressure
    
def ProcSkinAndPressure(model_part,solver):
    
    #SKIN DETERMINATION

    Pressure = ConfinementPressure*1e6 #Mpa
    total_cross_section = 0.0

    #Cylinder dimensions

    h   = 0.3
    d   = 0.15
    eps = 2.0

    surface = 2*(3.141592*d*d*0.25)+(3.141592*d*h)

    top_pressure = 0.0
    bot_pressure = 0.0

    xlat_area = 0.0
    xbot_area = 0.0
    xtop_area = 0.0
    xbotcorner_area = 0.0
    xtopcorner_area = 0.0
    
    for element in model_part.Elements:
      
      element.SetValue(SKIN_SPHERE,0)
   
      if ( predefined_skin_option == "OFF" ):
      
        node = element.GetNode(0)
        r = node.GetSolutionStepValue(RADIUS,0)
        x = node.X
        y = node.Y
        z = node.Z
        node_group = node.GetSolutionStepValue(GROUP_ID,0)
        cross_section = 3.141592*r*r

        #if( (node_group!=2) and (node_group!=4) ):
      
        if ( (x*x+z*z)>=((d/2-eps*r)*(d/2-eps*r)) ): 
      
             element.SetValue(SKIN_SPHERE,1)     

             LAT.append(node)
            
             if ( (y>eps*r ) and (y<(h-eps*r)) ) :
        
               SKIN.append(element)
        
               XLAT.append(node)
          
               xlat_area = xlat_area + cross_section
    
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
                       
                       xtopcorner_area = xtopcorner_area + cross_section
                
                   else:

                       XBOTCORNER.append(node)
                       xbotcorner_area = xbotcorner_area + cross_section
               else:

                   if ( y<=eps*r ):
                
                       XBOT.append(node)
                       xbot_area = xbot_area + cross_section
                
                   elif ( y>=(h-eps*r) ):
                    
                       XTOP.append(node)
                       xtop_area = xtop_area + cross_section

    print "End CLASSIC TEST SKIN DETERMINATION", "\n"
              
    return (xtop_area,xbot_area,xlat_area,xtopcorner_area,xbotcorner_area) 
     
def ProcApplyPressure(Pressure,model_part,solver,alpha_top,alpha_bot,alpha_lat):
     
    if ( predefined_skin_option == "ON"):
      print "\n", "Predefined Skin by the user, In this case is not correct to apply pressure yet"  ,"\n" 
    else:
      ApplyPressure(Pressure,model_part,solver,SKIN,BOT,TOP,LAT,XLAT,XBOT,XTOP,XBOTCORNER,XTOPCORNER,alpha_top,alpha_bot,alpha_lat)

    
def ProcMeasureBOT(BOT,solver):

    tol = 2.0
    y_mean = 0.0
    counter = 0.0
    for node in BOT:
      r = node.GetSolutionStepValue(RADIUS,0)
      y = node.Y
      
      y_mean += (y-r)*r
      counter += r

    return (y_mean,counter)      

def ProcMeasureTOP(TOP,solver):

    tol = 2.0
    y_mean = 0.0
    counter = 0.0
    
    for node in TOP:
      r = node.GetSolutionStepValue(RADIUS,0)
      y = node.Y
      
      y_mean += (y+r)*r
      counter += r

    return (y_mean,counter)
    
def ProcPrintingVariables(gid_io,export_model_part,time):
  
    if (print_displacement=="1"):
      gid_io.WriteNodalResults(DISPLACEMENT, export_model_part.Nodes, time, 0)       
    if (print_radial_displacement=="1"):
      gid_io.WriteNodalResults(RADIAL_DISPLACEMENT, export_model_part.Nodes, time, 0)       
    if (print_velocity=="1"):
      gid_io.WriteNodalResults(VELOCITY, export_model_part.Nodes, time, 0)
    if (print_rhs=="1"):
      gid_io.WriteNodalResults(RHS, export_model_part.Nodes, time, 0)       
    if (print_applied_forces=="1"):
      gid_io.WriteNodalResults(APPLIED_FORCE, export_model_part.Nodes, time, 0)       
    if (print_total_forces=="1"):     
      gid_io.WriteNodalResults(TOTAL_FORCES, export_model_part.Nodes, time, 0)    
    if (print_damp_forces=="1"):
      gid_io.WriteNodalResults(DAMP_FORCES, export_model_part.Nodes, time, 0)        
    if (print_radius=="1"):
      gid_io.WriteNodalResults(RADIUS, export_model_part.Nodes, time, 0)       
    if (print_particle_cohesion=="1"):
      gid_io.WriteNodalResults(PARTICLE_COHESION, export_model_part.Nodes, time, 0)       
    if (print_particle_tension=="1"):
      gid_io.WriteNodalResults(PARTICLE_TENSION, export_model_part.Nodes, time, 0)
    if (print_group_id=="1"):
      gid_io.WriteNodalResults(GROUP_ID, export_model_part.Nodes, time, 0)
    if (print_export_id=="1"):
      gid_io.WriteNodalResults(EXPORT_ID, export_model_part.Nodes, time, 0)
    if (print_export_particle_failure_id=="1"):
      gid_io.WriteNodalResults(EXPORT_PARTICLE_FAILURE_ID, export_model_part.Nodes, time, 0)
    if (print_export_skin_sphere=="1"):
      gid_io.WriteNodalResults(EXPORT_SKIN_SPHERE, export_model_part.Nodes, time, 0)
    if (print_stress_tensor == "1"):
      gid_io.WriteNodalResults(DEM_STRESS_XX, export_model_part.Nodes, time, 0)
      gid_io.WriteNodalResults(DEM_STRESS_XY, export_model_part.Nodes, time, 0)
      gid_io.WriteNodalResults(DEM_STRESS_XZ, export_model_part.Nodes, time, 0)
      gid_io.WriteNodalResults(DEM_STRESS_YX, export_model_part.Nodes, time, 0)
      gid_io.WriteNodalResults(DEM_STRESS_YY, export_model_part.Nodes, time, 0)
      gid_io.WriteNodalResults(DEM_STRESS_YZ, export_model_part.Nodes, time, 0)
      gid_io.WriteNodalResults(DEM_STRESS_ZX, export_model_part.Nodes, time, 0)
      gid_io.WriteNodalResults(DEM_STRESS_ZY, export_model_part.Nodes, time, 0)
      gid_io.WriteNodalResults(DEM_STRESS_ZZ, export_model_part.Nodes, time, 0)
    if (print_representative_volume == "1"):
      gid_io.WriteNodalResults(REPRESENTATIVE_VOLUME, export_model_part.Nodes, time, 0)
    
    #Aixo sempre per que si no hi ha manera de debugar
    #gid_io.WriteNodalResults(PARTITION_INDEX, export_model_part.Nodes, time, 0)
    #gid_io.WriteNodalResults(INTERNAL_ENERGY, export_model_part.Nodes, time, 0)

    if (ContactMeshOption == "ON"): ##xapuza
      if (print_local_contact_force_low=="1"):
          gid_io.PrintOnGaussPoints(LOCAL_CONTACT_FORCE_LOW,export_model_part,time)
      if (print_local_contact_force_high=="1"):
          gid_io.PrintOnGaussPoints(LOCAL_CONTACT_FORCE_HIGH,export_model_part,time)
      if (print_mean_contact_area=="1"): 
          gid_io.PrintOnGaussPoints(MEAN_CONTACT_AREA,export_model_part,time)
      if (print_contact_failure=="1"): 
          gid_io.PrintOnGaussPoints(CONTACT_FAILURE,export_model_part,time)  
      if (print_failure_criterion_state=="1"):
          gid_io.PrintOnGaussPoints(FAILURE_CRITERION_STATE,export_model_part,time)         
      if (print_contact_tau=="1"):
          gid_io.PrintOnGaussPoints(CONTACT_TAU,export_model_part,time)
      if (print_contact_sigma=="1"):
          gid_io.PrintOnGaussPoints(CONTACT_SIGMA,export_model_part,time)
          gid_io.PrintOnGaussPoints(LOCAL_CONTACT_AREA_HIGH,export_model_part,time)
          gid_io.PrintOnGaussPoints(LOCAL_CONTACT_AREA_LOW,export_model_part,time)
      gid_io.PrintOnGaussPoints(NON_ELASTIC_STAGE,export_model_part,time)    

    if (RotationOption == "ON"): ##xapuza
      if (print_angular_velocity=="1"):
          gid_io.WriteNodalResults(ANGULAR_VELOCITY, export_model_part.Nodes, time, 0)
      if (print_particle_moment=="1"):
          gid_io.WriteNodalResults(PARTICLE_MOMENT, export_model_part.Nodes, time, 0)
      if (print_euler_angles=="1"):
          gid_io.WriteLocalAxesOnNodes(EULER_ANGLES, export_model_part.Nodes, time, 0)

    gid_io.Flush()
    sys.stdout.flush()