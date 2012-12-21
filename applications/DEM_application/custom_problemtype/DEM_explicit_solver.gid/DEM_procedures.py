from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
#from KratosMultiphysics.MetisApplication import *
#from KratosMultiphysics.mpi import *

from numpy import *

from DEM_explicit_solver_var import *

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
    
def InitializeContactModelPart(solver):
        
    contact_model_part = solver.contact_model_part
    
    #MPICommSetup = SetMPICommunicatorProcess(contact_model_part)
    #MPICommSetup.Execute()
    
def InitializeSolver(model_part,solver):
    
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
    else:
        rota_damp_id = 0
        
    solver.rota_damp_id=rota_damp_id


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
    
    #Defining list of skin particles (For a test tube of height 30 cm and diameter 15 cm)

    Pressure = ConfinementPressure*1e6 #Mpa
    skin_list = list()

    if(ConcreteTestOption =="ON"): 
      
      for element in model_part.Elements:
        
        element.SetValue(SKIN_SPHERE,0)
        
        node = element.GetNode(0)
        r = node.GetSolutionStepValue(RADIUS,0)
        x = node.X
        y = node.Y
        z = node.Z
        
        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0
          
        Cross_section = 3.141592*r*r
        h=0.3
        d=0.15
        eps=2
        vect = zeros(3, double) 

        if ( (x*x+z*z)>=((d/2-eps*r)*(d/2-eps*r)) ): 
          
          element.SetValue(SKIN_SPHERE,1)
          skin_list.append(element)
          
          #vector normal al centre:
          vect_moduli = sqrt(x*x+z*z)
          #print(vect_moduli)
          if(vect_moduli>0.0):
            vect[0]=-x/vect_moduli
            vect[1]=0
            vect[2]=-z/vect_moduli
            
          #radius_elem = element.Node[0].GetSolutionStepValue(RADIUS)
          
          values[0]=-Cross_section*Pressure*vect[0]
          values[1]= 0.0
          values[2]=-Cross_section*Pressure*vect[2]
          
          
        if ( (y<=eps*r ) or (y>=(h-eps*r)) ): 

            element.SetValue(SKIN_SPHERE,1)
            #vector normal al centre:      
            values[0]=0.0
            values[2]=0.0
            if ( y>h/2 ):
                values[1]=-Cross_section*Pressure
            else:
                values[1]= Cross_section*Pressure

            if ( (x*x+z*z) >= ((d/2-eps*r)*(d/2-eps*r) ) ) :
                #vector normal al centre:
                vect_moduli = sqrt(x*x+z*z)
                
                if ( vect_moduli>0.0 ) :
                    vect[0]=-x/vect_moduli
                    vect[1]=0.0
                    vect[2]=-z/vect_moduli

                values[0]=-Cross_section*Pressure*vect[0]*0.70710678
                values[2]=-Cross_section*Pressure*vect[2]*0.70710678
                if ( y>h/2 ):
                    values[1]=-Cross_section*Pressure*0.70710678
                else:
                    values[1]=Cross_section*Pressure*0.70710678 
            else:
                skin_list.append(element)  

        node.SetSolutionStepValue(APPLIED_FORCE,values)
        
    print("End Applying Imposed Forces")
      
         
    #Initialize ContactModelPart
    
    InitializeContactModelPart(solver)   
          
    #Initialize the problem.

    solver.Initialize()