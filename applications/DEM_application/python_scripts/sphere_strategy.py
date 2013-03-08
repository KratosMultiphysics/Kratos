
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

# Check that KratosMultiphysics was imported in the main script
#CheckForPreviousImport(

from DEM_explicit_solver_var import *

def AddVariables(model_part):
  
   
# BASIQUES

    model_part.AddNodalSolutionStepVariable(RADIUS)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VELOCITY)    
    model_part.AddNodalSolutionStepVariable(ORIENTATION_REAL)
    model_part.AddNodalSolutionStepVariable(ORIENTATION_IMAG)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
    
    #may be grouped
    model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)
    model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS)
    model_part.AddNodalSolutionStepVariable(POISSON_RATIO)
    model_part.AddNodalSolutionStepVariable(RESTITUTION_COEFF)  # pot posarse mes ifs.....
    model_part.AddNodalSolutionStepVariable(PARTICLE_COHESION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_FRICTION) 
    model_part.AddNodalSolutionStepVariable(PARTICLE_TENSION)
    
    
    #es podrien eliminar
    model_part.AddNodalSolutionStepVariable( NODAL_MASS )
    
    model_part.AddNodalSolutionStepVariable( NUM_OF_NEIGH ) #temporal    
    model_part.AddNodalSolutionStepVariable( DEM_STRESS_XX )    
    model_part.AddNodalSolutionStepVariable( DEM_STRESS_XY )
    model_part.AddNodalSolutionStepVariable( DEM_STRESS_XZ )
    model_part.AddNodalSolutionStepVariable( DEM_STRESS_YX )
    model_part.AddNodalSolutionStepVariable( DEM_STRESS_YY )
    model_part.AddNodalSolutionStepVariable( DEM_STRESS_YZ )
    model_part.AddNodalSolutionStepVariable( DEM_STRESS_ZX )
    model_part.AddNodalSolutionStepVariable( DEM_STRESS_ZY )
    model_part.AddNodalSolutionStepVariable( DEM_STRESS_ZZ )
    
    
    
# ADVANCED

    model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   #colour defined in GiD
    model_part.AddNodalSolutionStepVariable(PARTICLE_CONTINUUM)  #Continuum group
    model_part.AddNodalSolutionStepVariable(GROUP_ID)            #differencied groups for plotting, etc..
    model_part.AddNodalSolutionStepVariable(REPRESENTATIVE_VOLUME)
    
    model_part.AddNodalSolutionStepVariable(RHS)
    model_part.AddNodalSolutionStepVariable(DAMP_FORCES)
    model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
    model_part.AddNodalSolutionStepVariable(APPLIED_FORCE)

    #if ( (ConfinementPressure != 0.0) and (TriaxialOption == "ON" ) ):
    model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_FORCE) 

#ROTATION

    if(RotationOption =="ON"):

      model_part.AddNodalSolutionStepVariable(PARTICLE_INERTIA)
      model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)     
      model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA)
      model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
      model_part.AddNodalSolutionStepVariable(DELTA_ROTA_DISPLACEMENT)
      
      #if(TrihedronOption =="ON"):
      model_part.AddNodalSolutionStepVariable(EULER_ANGLES)
         
      #if(RotaDampId =="LocalDamp"):
      model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_DAMP_RATIO)
         
      #if(RotaDampId =="RollingFric"): 
      model_part.AddNodalSolutionStepVariable(ROLLING_FRICTION) 

#ONLY VISUALITZATION
 
    if(print_export_skin_sphere =="1"):
      model_part.AddNodalSolutionStepVariable(EXPORT_SKIN_SPHERE)
      
    if(print_export_id =="1"):      
      model_part.AddNodalSolutionStepVariable(EXPORT_ID) 
      
    if(print_radial_displacement =="1"):
      model_part.AddNodalSolutionStepVariable(RADIAL_DISPLACEMENT) 
         
    
#Temporarily Unused

    #model_part.AddNodalSolutionStepVariable(EXPORT_PARTICLE_FAILURE_ID) //tocar el EXPORT_IDalizestep del sphere.cpp quan l'activis
 

    print "variables for the explicit solver added correctly"

def AddDofs(model_part):
    
    for node in model_part.Nodes:
    #adding dofs
        node.AddDof(DISPLACEMENT_X,REACTION_X);
        node.AddDof(DISPLACEMENT_Y,REACTION_Y);
        node.AddDof(DISPLACEMENT_Z,REACTION_Z);
        node.AddDof(VELOCITY_X,REACTION_X);
        node.AddDof(VELOCITY_Y,REACTION_Y);
        node.AddDof(VELOCITY_Z,REACTION_Z);
        node.AddDof(ANGULAR_VELOCITY_X,REACTION_X);
        node.AddDof(ANGULAR_VELOCITY_Y,REACTION_Y);
        node.AddDof(ANGULAR_VELOCITY_Z,REACTION_Z);        

    print "dofs for the DEM solution added correctly"

 
class ExplicitStrategy:
    
    def __init__(self,model_part,domain_size):

        self.model_part                     	= model_part    
        self.contact_model_part             	= ModelPart("ContactModelPart")
        self.domain_size                    	= domain_size
        self.damping_ratio                  	= 0.00;
        self.penalty_factor                 	= 10.00
        self.max_delta_time                 	= 0.05;
        self.final_time							= 3.0;
        self.fraction_delta_time            	= 0.90;
        self.MoveMeshFlag                   	= True;
        self.gravity                        	= Vector(3)#(0.0,-9.81,0.0)
        self.gravity[0] = 0.0
        self.gravity[1] = -9.81
        self.gravity[2] = 0.0
        self.delta_time                     	= 0.00001;
        self.virtual_mass_OPTION            	= 0; #its 1/0 xapuza
        self.nodal_mass_coeff               	= 0.0;
      
        #type of problem:

        self.critical_time_OPTION           	= 0 #its 1/0 xapuza
        self.final_time                         = 3.0
        
        self.magic_factor						=1.0

        self.delta_OPTION                   	= False
        self.continuum_simulating_OPTION    	= False
        self.case_OPTION                    	= 0  #aixo es una xapuza fins que pooyan permeti bools a pyton o tinguis flags.
        self.trihedron_OPTION               	= 0

        self.rotation_OPTION                	= 0  #its 1/0 xapuza
        self.rotation_spring_OPTION         	= 0  #its 1/0 xapuza
        self.bounding_box_OPTION            	= 0  #its 1/0 xapuza
        self.activate_search					= 1  #its 1/0 xapuza
        self.concrete_test_OPTION				= 0  #its 1/0 xapuza
        
        self.stress_strain_operations          = 0
        
        self.external_pressure					= 0
        self.time_increasing_ratio				= 15 # Percentage%
        self.initial_pressure_time				= 0.0
        
        self.contact_mesh_OPTION               = 0 #its 1/0 xapuza
        self.failure_criterion_OPTION          = 1 #its 1/0 xapuza
        self.tau_zero                          = 0.0
        self.sigma_max                         = 0.0
        self.sigma_min                         = 0.0
        self.internal_fricc	                    = 0.0
        
        self.Non_Linear_Option                 = 0
        self.N1                                = 1.0
        self.N2                                = 1.0
        self.C1                                = 1.0
        self.C2                                = 1.0
        
        self.fix_velocities                      = 0  
        self.time_step_percentage_fix_velocities = 0 #int(final_time/delta_time)*10;

        #global parameters
        self.global_variables_OPTION          = 0 #its 1/0 xapuza
        self.global_kn                        = 1000.0
        self.global_kt                        = 1000.0
        self.global_kr                        = 1000.0
        self.global_rn                        = 1000.0
        self.global_rt                        = 1000.0
        self.global_rr                        = 1000.0
        self.global_fri_ang                   = 40
        
        #prints
        
        self.print_export_id                   = 0
        self.print_export_skin_sphere          = 0
        self.print_radial_displacement         = 0

        #problem specific parameters

        self.force_calculation_type_id      	= 1
        self.damp_id                        	= 1
        self.rota_damp_id                   	= 0
        self.search_radius_extension        	= 0.0

        self.dummy_switch                   	= 0

        #problem utilities
        self.enlargement_factor             	= 1;
        self.n_step_search                  	= 1;
        self.safety_factor                  	= 1.0; #for critical time step

        self.create_and_destroy             	= particle_destructor_and_constructor();
        
        self.use_mpi                        	= 0; #MPI Carlos change to 1
        
        if(self.use_mpi):
            self.time_scheme                 	= MpiForwardEulerScheme();
        else:
            self.time_scheme                  	= ForwardEulerScheme();

    ######################################################################

    def Initialize(self):

        #definir les variables del ProcessInfo:
      
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)

        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)
        self.model_part.ProcessInfo.SetValue(FINAL_SIMULATION_TIME, self.final_time)

        #LA CONTINUUM OPTION SHA DE LLEGIR MIRANT SI TOTS ELS CONTINUUM GROUPS SON 0 LLAVORS LA OPTION ES FALSE
        #SI HI HA ALGUN KE ES DIFERENT DE 0 LLAVORS LA OPTION ES TRUE.
        
        if(self.delta_OPTION==True):
            if(self.continuum_simulating_OPTION==True): self.case_OPTION = 2
            else: self.case_OPTION = 1
        elif(self.delta_OPTION==False):
            if(self.continuum_simulating_OPTION==False): self.case_OPTION = 0
            else: self.case_OPTION = 3

        self.model_part.ProcessInfo.SetValue(AREA_VERTICAL_TAPA,0.0);
        self.model_part.ProcessInfo.SetValue(AREA_VERTICAL_CENTRE,0.0);
        self.model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED,0);
        self.model_part.ProcessInfo.SetValue(TOTAL_CONTACTS,0);
        self.model_part.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, self.critical_time_OPTION)
        self.model_part.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, self.virtual_mass_OPTION)
        self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)
        self.model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_OPTION)
        self.model_part.ProcessInfo.SetValue(ACTIVATE_SEARCH, self.activate_search)
        self.model_part.ProcessInfo.SetValue(ROTATION_OPTION, self.rotation_OPTION)
        self.model_part.ProcessInfo.SetValue(ROTATION_SPRING_OPTION, self.rotation_spring_OPTION)
        self.model_part.ProcessInfo.SetValue(BOUNDING_BOX_OPTION, self.bounding_box_OPTION)
        self.model_part.ProcessInfo.SetValue(TRIHEDRON_OPTION, self.trihedron_OPTION)
        self.model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, self.contact_mesh_OPTION)
        self.model_part.ProcessInfo.SetValue(CONCRETE_TEST_OPTION, self.concrete_test_OPTION)
        
        
        self.model_part.ProcessInfo.SetValue(FAILURE_CRITERION_OPTION, self.failure_criterion_OPTION)
        self.model_part.ProcessInfo.SetValue(CONTACT_SIGMA_MAX, self.sigma_max)
        self.model_part.ProcessInfo.SetValue(CONTACT_SIGMA_MIN, self.sigma_min)
        self.model_part.ProcessInfo.SetValue(CONTACT_TAU_ZERO, self.tau_zero)
        self.model_part.ProcessInfo.SetValue(CONTACT_INTERNAL_FRICC, self.internal_fricc)
        
        self.model_part.ProcessInfo.SetValue(NON_LINEAR_OPTION, self.Non_Linear_Option)
        
        self.model_part.ProcessInfo.SetValue(SLOPE_FRACTION_N1, self.N1)
        self.model_part.ProcessInfo.SetValue(SLOPE_FRACTION_N2, self.N2)
        
        self.model_part.ProcessInfo.SetValue(SLOPE_LIMIT_COEFF_C1, self.C1)
    
        self.model_part.ProcessInfo.SetValue(SLOPE_LIMIT_COEFF_C2, self.C2)
        self.model_part.ProcessInfo.SetValue(INITIAL_PRESSURE_TIME, self.initial_pressure_time)
        self.model_part.ProcessInfo.SetValue(TIME_INCREASING_RATIO, self.time_increasing_ratio)
       
        #####

        self.model_part.ProcessInfo.SetValue(FORCE_CALCULATION_TYPE, self.force_calculation_type_id)    
        self.model_part.ProcessInfo.SetValue(DAMP_TYPE, self.damp_id)
        self.model_part.ProcessInfo.SetValue(ROTA_DAMP_TYPE, self.rota_damp_id)
        self.model_part.ProcessInfo.SetValue(SEARCH_RADIUS_EXTENSION, self.search_radius_extension)

        self.model_part.ProcessInfo.SetValue(GLOBAL_VARIABLES_OPTION, self.global_variables_OPTION)
        self.model_part.ProcessInfo.SetValue(GLOBAL_KN, self.global_kn)
        self.model_part.ProcessInfo.SetValue(GLOBAL_KT, self.global_kt)
        self.model_part.ProcessInfo.SetValue(GLOBAL_KR, self.global_kr)
        self.model_part.ProcessInfo.SetValue(GLOBAL_RN, self.global_rn)
        self.model_part.ProcessInfo.SetValue(GLOBAL_RT, self.global_rt)
        self.model_part.ProcessInfo.SetValue(GLOBAL_RR, self.global_rr)
        self.model_part.ProcessInfo.SetValue(GLOBAL_FRI_ANG, self.global_fri_ang)
        
        self.model_part.ProcessInfo.SetValue(DEM_MAGIC_FACTOR,self.magic_factor)

        self.model_part.ProcessInfo.SetValue(DUMMY_SWITCH, self.dummy_switch)
        
        #self.model_part.ProcessInfo.SetValue(FIXED_VEL_BOT, 0)
        
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_1, 0) #Reserved for: message when confinement ends.
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_2, self.external_pressure) #Reserved for: External Applied force is acting
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_3, self.print_export_id) #reserved for: Export Print Skin sphere
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_4, self.print_export_skin_sphere) #reserved for print_export_skin_sphere
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_5, 0) 
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_6, self.fix_velocities) #reserved for fix_velocities
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_7, 0)#int( self.time_step_percentage_fix_velocities * ( self.final_time/self.delta_time) ) ) #reserved for timestep fix_velocities
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_8, self.print_radial_displacement)#reserved for ON OFF print RADIAL_DISPLACEMENT
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_9, self.stress_strain_operations)#reserved for ON_OFF stress_strain_operations
        
        self.model_part.ProcessInfo.SetValue(DOUBLE_DUMMY_1, 0.0) 
        self.model_part.ProcessInfo.SetValue(DOUBLE_DUMMY_2, 0.0) 
        self.model_part.ProcessInfo.SetValue(DOUBLE_DUMMY_3, self.time_step_percentage_fix_velocities)# reserved for percentage when start the fixing of velocities
        self.model_part.ProcessInfo.SetValue(DOUBLE_DUMMY_4, 0.0)
        self.model_part.ProcessInfo.SetValue(DOUBLE_DUMMY_5, 0.0)
        self.model_part.ProcessInfo.SetValue(DOUBLE_DUMMY_6, 0.0)
        self.model_part.ProcessInfo.SetValue(DOUBLE_DUMMY_7, 0.0)
        self.model_part.ProcessInfo.SetValue(DOUBLE_DUMMY_8, 0.0)
        self.model_part.ProcessInfo.SetValue(DOUBLE_DUMMY_9, 0.0)
               
        
        #creating the solution strategy
        if(self.use_mpi):
            self.solver = MpiExplicitSolverStrategy(self.model_part, self.contact_model_part, self.domain_size, self.enlargement_factor, self.damping_ratio, self.fraction_delta_time, self.delta_time, self.n_step_search, self.safety_factor,
                                            self.MoveMeshFlag, self.delta_OPTION, self.continuum_simulating_OPTION, self.time_scheme)
        else:
            self.solver = ExplicitSolverStrategy(self.model_part, self.contact_model_part, self.domain_size, self.enlargement_factor, self.damping_ratio, self.fraction_delta_time, self.delta_time, self.n_step_search, self.safety_factor,
                                            self.MoveMeshFlag, self.delta_OPTION, self.continuum_simulating_OPTION, self.time_scheme)
        #self.solver.Check() #es sa fer sempre un check despres de montar una estrategia.
        self.solver.Initialize() #aqui definirem el initialize dels elements pero tamb funcions que vulguem fer en el primer pras.
        
        #copying the nodes to the new model part, also the properties (void)
        self.contact_model_part.Nodes = self.model_part.Nodes;
     

    #######################################################################
    def Initial_Critical_Time(self):
        (self.solver).InitialCriticalTime()   

    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    
     #######################################################################

    #def Calculate_Model_Surrounding_Bounding_Box(self, model_part, enlargement_factor):
        #self.create_and_destroy.calculate_surrounding_bounding_box( model_part, enlargement_factor)

    #def Destroy_Particles(self, model_part):
        #self.create_and_destroy.destroy_distant_particles( model_part)

        
        
        #self.explicit_solver_object.Search_Neighbours() #after destructing we need to recalculate the neighbours instantly
        ##aixo te el problema de que despres de destruir hem de buscar pero tornarema  buscar en la estrategia??? organitzar-ho be.

          ##M:  be sure that after this we search for neighbours again.
