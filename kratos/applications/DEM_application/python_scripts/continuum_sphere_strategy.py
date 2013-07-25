from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

def AddVariables(model_part, Param):

    # KINEMATIC
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
    model_part.AddNodalSolutionStepVariable(DELTA_ROTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(ORIENTATION_REAL)
    model_part.AddNodalSolutionStepVariable(ORIENTATION_IMAG)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)

    # FORCES
    model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
    model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
    model_part.AddNodalSolutionStepVariable(DAMP_FORCES)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_FORCE)

    # BASIC PARTICLE PROPERTIES
    model_part.AddNodalSolutionStepVariable(RADIUS)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(SQRT_OF_MASS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)
    model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS)
    model_part.AddNodalSolutionStepVariable(POISSON_RATIO)
    model_part.AddNodalSolutionStepVariable(LN_OF_RESTITUTION_COEFF)
    model_part.AddNodalSolutionStepVariable(PARTICLE_FRICTION)

    # ROTATION RELATED PROPERTIES
    if (Var_Translator(Param.RotationOption)):
        model_part.AddNodalSolutionStepVariable(PARTICLE_INERTIA)
        model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA)
        model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_DAMP_RATIO)
        model_part.AddNodalSolutionStepVariable(ROLLING_FRICTION)

    # OTHER PROPERTIES
    model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD
    model_part.AddNodalSolutionStepVariable(PARTICLE_CONTINUUM)  # Continuum group
    model_part.AddNodalSolutionStepVariable(REPRESENTATIVE_VOLUME)
    model_part.AddNodalSolutionStepVariable(MAX_INDENTATION)
    
    # LOCAL AXIS
    model_part.AddNodalSolutionStepVariable(EULER_ANGLES)

    # BOUNDARY SURFACE
    if (Var_Translator(Param.LimitSurfaceOption)):
        model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_CONTACT_FORCES)
        model_part.AddNodalSolutionStepVariable(PARTICLE_SURFACE_ROTATE_SPRING_MOMENT)

    # FLAGS
    model_part.AddNodalSolutionStepVariable(GROUP_ID)            # Differencied groups for plotting, etc..
    model_part.AddNodalSolutionStepVariable(ERASE_FLAG)
    # ONLY VISUALITZATION

    if (Var_Translator(Param.PostExportId)):
      model_part.AddNodalSolutionStepVariable(EXPORT_ID)
    if (Var_Translator(Param.PredefinedSkinOption) | Var_Translator(Param.ConcreteTestOption) ):
      model_part.AddNodalSolutionStepVariable(EXPORT_SKIN_SPHERE)
    if (Var_Translator(Param.PostGroupId)):
      model_part.AddNodalSolutionStepVariable(EXPORT_GROUP_ID)

    print "Variables for the explicit solver added correctly"

def AddDofs(model_part):

    for node in model_part.Nodes:
        node.AddDof(DISPLACEMENT_X, REACTION_X);
        node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        node.AddDof(VELOCITY_X, REACTION_X);
        node.AddDof(VELOCITY_Y, REACTION_Y);
        node.AddDof(VELOCITY_Z, REACTION_Z);
        node.AddDof(ANGULAR_VELOCITY_X, REACTION_X);
        node.AddDof(ANGULAR_VELOCITY_Y, REACTION_Y);
        node.AddDof(ANGULAR_VELOCITY_Z, REACTION_Z);

    print "DOFs for the DEM solution added correctly"

def Var_Translator(variable):

    if (variable == "OFF" or variable == "0" or variable == 0):
        variable = 0
    else:
        variable = 1

    return variable

class ExplicitStrategy:

    def __init__(self, model_part, creator_destructor, Param):

        # Initialization of member variables

        # SIMULATION FLAGS        
        self.virtual_mass_option            = Var_Translator(Param.VirtualMassOption)  #its 1/0 xapuza
        self.critical_time_option           = Var_Translator(Param.AutoReductionOfTimeStepOption)  #its 1/0 xapuza
        self.case_option                    = 3  #aixo es una xapuza fins que pooyan permeti bools a pyton o tinguis flags.
        self.trihedron_option               = Var_Translator(Param.TrihedronOption)
        self.rotation_option                = Var_Translator(Param.RotationOption)
        self.rotation_spring_option         = Var_Translator(Param.RotationalSpringOption)  #its 1/0 xapuza
        self.bounding_box_option            = Var_Translator(Param.BoundingBoxOption)  #its 1/0 xapuza
        self.activate_search                = 1  #its 1/0 xapuza
        self.fix_velocities                 = Var_Translator(Param.FixVelocitiesOption)
        self.limit_surface_option           = Var_Translator(Param.LimitSurfaceOption)  #its 1/0 xapuza
        self.clean_init_indentation_option  = Var_Translator(Param.CleanIndentationsOption)
        self.homogeneous_material_option    = Var_Translator(Param.HomogeneousMaterialOption)
        self.global_variables_option        = Var_Translator(Param.GlobalVariablesOption)
        self.Non_Linear_Option              = Var_Translator(Param.NonLinearNormalElasticOption)
        self.stress_strain_operations       = Var_Translator(Param.StressStrainOperationsOption)
        self.MoveMeshFlag                   = True
        self.delta_option                   = Var_Translator(Param.DeltaOption)
        self.continuum_simulating_option    = Var_Translator(Param.ContinuumOption)
        self.contact_mesh_option            = Var_Translator( Var_Translator(Param.ContactMeshOption) & Var_Translator(Param.ContinuumOption) ) 
        self.concrete_test_option           = Var_Translator( Var_Translator(Param.ConcreteTestOption) & Var_Translator(Param.ContinuumOption) ) 
        self.triaxial_option                = Var_Translator( Var_Translator(Param.TriaxialOption) & self.concrete_test_option )
       
        self.search_radius_extension        = 1e-6 #needed for the tangential contacts. Charlie will modify the search. 
        self.amplified_continuum_search_radius_extension    = 1.0;
        self.automatic_bounding_box_option  = Var_Translator(Param.AutomaticBoundingBoxOption)              

        if (self.delta_option ):
            self.delta_option               = True
            self.search_radius_extension    = Param.SearchRadiusExtension

        if (Var_Translator(Param.ContinuumOption)):
            self.amplified_continuum_search_radius_extension = Param.AmplifiedSearchRadiusExtension;
                 
        if(self.delta_option==True):
          if(self.continuum_simulating_option): self.case_option = 2
          else: self.case_option = 1
        elif(self.delta_option==False):
          if(self.continuum_simulating_option == False): self.case_option = 0
          else: self.case_option = 3     
            
        # MODEL
        self.model_part                     = model_part
        self.contact_model_part             = ModelPart("ContactModelPart") #funcio kratos
        self.contact_model_part.Nodes       = self.model_part.Nodes;
        self.domain_size                    = Param.Dimension

        # BOUNDARY
        self.surface_normal_dir             = Vector(3)
        self.surface_normal_dir[0]          = Param.SurfaceNormalDirX
        self.surface_normal_dir[1]          = Param.SurfaceNormalDirY
        self.surface_normal_dir[2]          = Param.SurfaceNormalDirZ
        self.surface_point_coor             = Vector(3)
        self.surface_point_coor[0]          = Param.SurfacePointCoorX
        self.surface_point_coor[1]          = Param.SurfacePointCoorY
        self.surface_point_coor[2]          = Param.SurfacePointCoorZ
        self.surface_friction_angle         = Param.SurfaceFrictionAngle


        # GLOBAL PHISICAL ASPECTS
        self.gravity                        = Vector(3)
        self.gravity[0]                     = Param.GravityX
        self.gravity[1]                     = Param.GravityY
        self.gravity[2]                     = Param.GravityZ

        # GLOBAL MATERIAL PROPERTIES
        self.nodal_mass_coeff               = Param.VirtualMassCoefficient
        self.magic_factor                   = Param.MagicFactor

        if (self.global_variables_option):
            self.global_kn                  = Param.GlobalKn
            self.global_kt                  = Param.GlobalKt
            self.global_kr                  = Param.GlobalKr
            self.global_rn                  = Param.GlobalRn
            self.global_rt                  = Param.GlobalRT
            self.global_rr                  = Param.GlobalRr
            self.global_fri_ang             = Param.GlobalFrictionAngle

        if (Param.NormalForceCalculationType == "Linear"):
            self.force_calculation_type_id  = 0
        elif (Param.NormalForceCalculationType == "Hertz"):
            self.force_calculation_type_id  = 1

        if (self.Non_Linear_Option):
            self.C1                         = Param.C1
            self.C2                         = Param.C2
            self.N1                         = Param.N1
            self.N2                         = Param.N2

        if (Param.NormalDampingType == "ViscDamp"):

            if (Param.TangentialDampingType == "ViscDamp"):
                self.damp_id                = 11

            else:
                self.damp_id                = 10
        else:

            if (Param.TangentialDampingType == "ViscDamp"):
                self.damp_id                = 1

            else:
                self.damp_id                = 0

        if (Param.RotaDampingType == "LocalDamp"):
            self.rota_damp_id               = 1

        elif (Param.RotaDampingType == "RollingFric"):
            self.rota_damp_id               = 2

        else:
            self.rota_damp_id               = 0

        if (Param.FailureCriterionType == "Mohr-Coulomb"):
            self.failure_criterion_option   = 1

        elif (Param.FailureCriterionType == "Uncoupled"):
            self.failure_criterion_option   = 2

        self.tau_zero                       = Param.TauZero
        self.sigma_max                      = Param.SigmaMax
        self.sigma_min                      = Param.SigmaMin
        self.internal_fricc                 = Param.InternalFriction
        
        # CONCRETE TEST
        if (self.triaxial_option):
            self.initial_pressure_time        = Param.InitialPressureAplicationTime
            self.time_increasing_ratio        = Param.TotalTimePercentAsForceAplTime # (%)
        
        
        # PRINTING VARIABLES
        self.print_export_id                = Var_Translator(Param.PostExportId)
        self.print_export_skin_sphere       = Var_Translator(Param.PostExportSkinSphere)
        self.print_radial_displacement      = Var_Translator(Param.PostRadialDisplacement)
        self.print_group_id                 = Var_Translator(Param.PostGroupId)        
        
        self.dummy_switch                   = 0

        # TIME RELATED PARAMETERS

        self.delta_time                     = Param.MaxTimeStep
        self.max_delta_time                 = Param.MaxTimeStep
        self.final_time                     = Param.FinalTime

        # RESOLUTION METHODS AND PARAMETERS

        self.n_step_search                  = int(Param.TimeStepsPerSearchStep)
        self.safety_factor                  = Param.DeltaTimeSafetyFactor # For critical time step

        # CREATOR-DESTRUCTOR

        self.creator_destructor             = creator_destructor

        b_box_low     = Array3()
        b_box_high    = Array3()
        b_box_low[0]  = Param.BoundingBoxMaxX
        b_box_low[1]  = Param.BoundingBoxMaxY
        b_box_low[2]  = Param.BoundingBoxMaxZ
        b_box_high[0] = Param.BoundingBoxMinX
        b_box_high[1] = Param.BoundingBoxMinY
        b_box_high[2] = Param.BoundingBoxMinZ

        self.creator_destructor.SetLowNode(b_box_low)
        self.creator_destructor.SetHighNode(b_box_high)

        #if (self.automatic_bounding_box_option):
            #self.creator_destructor.CalculateSurroundingBoundingBox()


        # STRATEGIES

        self.search_strategy                = OMP_DEMSearch()

        if (Param.IntegrationScheme == 'forward_euler'):
            self.time_scheme = ForwardEulerScheme()
        elif (Param.IntegrationScheme == 'mid_point_rule'):
            self.time_scheme = MidPointScheme()
        elif (Param.IntegrationScheme == 'const_average_acc'):
            self.time_scheme = ConstAverageAccelerationScheme()
        else:
            print('scheme not defined')

    ######################################################################

    def Initialize(self):

        print("CONTINUUM PYTHON STRATEGY!")

        # Setting ProcessInfo variables
        
        # SIMULATION FLAGS
        self.model_part.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, self.virtual_mass_option)
        self.model_part.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, self.critical_time_option)
        self.model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_option)
        self.model_part.ProcessInfo.SetValue(TRIHEDRON_OPTION, self.trihedron_option)
        self.model_part.ProcessInfo.SetValue(ROTATION_OPTION, self.rotation_option)
        self.model_part.ProcessInfo.SetValue(BOUNDING_BOX_OPTION, self.bounding_box_option)
        self.model_part.ProcessInfo.SetValue(ACTIVATE_SEARCH, self.activate_search)
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_6, self.fix_velocities) #reserved for fix_velocities
        self.model_part.ProcessInfo.SetValue(GLOBAL_VARIABLES_OPTION, self.global_variables_option)
        self.model_part.ProcessInfo.SetValue(UNIFORM_MATERIAL_OPTION, self.homogeneous_material_option)
        self.model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0);
        self.model_part.ProcessInfo.SetValue(TOTAL_CONTACTS, 0);
        self.model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_option);
        self.model_part.ProcessInfo.SetValue(ROTATION_SPRING_OPTION, self.rotation_spring_option);
    
        # TOLERANCES
        self.model_part.ProcessInfo.SetValue(DISTANCE_TOLERANCE, 0);
        
        # BOUNDARY
        self.model_part.ProcessInfo.SetValue(LIMIT_SURFACE_OPTION, self.limit_surface_option)
        self.model_part.ProcessInfo.SetValue(SURFACE_NORMAL_DIR, self.surface_normal_dir)
        self.model_part.ProcessInfo.SetValue(SURFACE_POINT_COOR, self.surface_point_coor)
        self.model_part.ProcessInfo.SetValue(SURFACE_FRICC, self.surface_friction_angle)

        # GLOBAL PHISICAL ASPECTS
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)
        self.model_part.ProcessInfo.SetValue(DEM_MAGIC_FACTOR, self.magic_factor)

        # GLOBAL MATERIAL PROPERTIES

        if(self.homogeneous_material_option):
            self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)
     
        if (self.global_variables_option):
            self.model_part.ProcessInfo.SetValue(GLOBAL_KN, self.global_kn)
            self.model_part.ProcessInfo.SetValue(GLOBAL_KT, self.global_kt)

        # PRINTING VARIABLES
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_10, self.print_radial_displacement)#reserved for ON OFF print RADIAL_DISPLACEMENT
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_8, self.print_group_id) # Reserved for: Export Print Group ID
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_3, self.print_export_id) # Reserved for: Export Id
        self.model_part.ProcessInfo.SetValue(INT_DUMMY_4, self.print_export_skin_sphere) # Reserved for: Export Print Skin sphere
        self.model_part.ProcessInfo.SetValue(FORCE_CALCULATION_TYPE, self.force_calculation_type_id)
        self.model_part.ProcessInfo.SetValue(DAMP_TYPE, self.damp_id)
        self.model_part.ProcessInfo.SetValue(ROTA_DAMP_TYPE, self.rota_damp_id)
        #self.model_part.ProcessInfo.SetValue(INT_DUMMY_1, 0) #currently unused.

        # TIME RELATED PARAMETERS
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)
        self.model_part.ProcessInfo.SetValue(FINAL_SIMULATION_TIME, self.final_time)
        
        #CONTINUUM
        
        self.model_part.ProcessInfo.SetValue(SEARCH_RADIUS_EXTENSION, self.search_radius_extension)
        self.model_part.ProcessInfo.SetValue(AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION, self.amplified_continuum_search_radius_extension)
        
        self.model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, self.contact_mesh_option)
                
        self.model_part.ProcessInfo.SetValue(FAILURE_CRITERION_OPTION, self.failure_criterion_option)
        self.model_part.ProcessInfo.SetValue(CONTACT_SIGMA_MAX, self.sigma_max)
        self.model_part.ProcessInfo.SetValue(CONTACT_SIGMA_MIN, self.sigma_min)
        self.model_part.ProcessInfo.SetValue(CONTACT_TAU_ZERO, self.tau_zero)
        self.model_part.ProcessInfo.SetValue(CONTACT_INTERNAL_FRICC, self.internal_fricc)
        
        self.model_part.ProcessInfo.SetValue(NON_LINEAR_OPTION, self.Non_Linear_Option)
        
        if (self.Non_Linear_Option):
            self.model_part.ProcessInfo.SetValue(SLOPE_FRACTION_N1, self.N1)
            self.model_part.ProcessInfo.SetValue(SLOPE_FRACTION_N2, self.N2)
            self.model_part.ProcessInfo.SetValue(SLOPE_LIMIT_COEFF_C1, self.C1)
            self.model_part.ProcessInfo.SetValue(SLOPE_LIMIT_COEFF_C2, self.C2)
        
        
        print(self.triaxial_option)
        if (self.triaxial_option):
            self.model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 1)
            self.model_part.ProcessInfo.SetValue(INITIAL_PRESSURE_TIME, self.initial_pressure_time)
            self.model_part.ProcessInfo.SetValue(TIME_INCREASING_RATIO, self.time_increasing_ratio)

        #OTHERS
        
        self.model_part.ProcessInfo.SetValue(DUMMY_SWITCH, self.dummy_switch)

        # RESOLUTION METHODS AND PARAMETERS
        # Creating the solution strategy

        self.solver = ContinuumExplicitSolverStrategy(self.model_part, self.contact_model_part, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                      self.MoveMeshFlag, self.delta_option, self.continuum_simulating_option, self.creator_destructor, self.time_scheme, self.search_strategy)
  
                                  
        self.solver.Initialize() # Calls the solver Initialized function (initializes all elements and performs other necessary tasks before iterating)

    #######################################################################

    def Initial_Critical_Time(self):
        (self.solver).InitialTimeStepCalculation()

    #######################################################################

    def Solve(self):
        (self.solver).Solve()

    #######################################################################
