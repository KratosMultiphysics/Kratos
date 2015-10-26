from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import sphere_strategy as SolverStrategy
BaseExplicitStrategy = SolverStrategy.ExplicitStrategy


def Var_Translator(variable):

    if (variable == "OFF" or variable == "0" or variable == 0):
        variable = 0
    else:
        variable = 1

    return variable


class ExplicitStrategy(BaseExplicitStrategy):   

    def __init__(self, model_part, fem_model_part, cluster_model_part, inlet_model_part, creator_destructor, Param):

        # Initialization of member variables

        # SIMULATION FLAGS        
  
        self.self_strain_option             = Var_Translator(Param.StressStrainOption)
        self.critical_time_option           = Var_Translator(Param.AutoReductionOfTimeStepOption)   
        self.case_option                    = 3  
        self.trihedron_option               = Var_Translator(Param.PostEulerAngles)
        self.rotation_option                = Var_Translator(Param.RotationOption)
        self.bounding_box_option            = Var_Translator(Param.BoundingBoxOption)
        self.search_control                 = 1
        if(Var_Translator(Param.DontSearchUntilFailure)):
          print ("Search is not active until a bond is broken.")
          self.search_control                = 0
          if (len(fem_model_part.Nodes)>0 or Param.TestType== "BTS" ):   #MSI. This activates the search since there are fem contact elements. however only the particle - fem search should be active.
            print ("WARNING: Search should be activated since there might contact with FEM.")

        self.fix_velocities_flag                 = 0       

        self.clean_init_indentation_option = Var_Translator(Param.CleanIndentationsOption)
        self.MoveMeshFlag = True
        
        self.virtual_mass_option            = 0
        self.nodal_mass_coeff = Param.VirtualMassCoefficient
        if(self.nodal_mass_coeff != 1.00):
           self.virtual_mass_option            = 1
        
        self.delta_option = Var_Translator(Param.DeltaOption)
        self.contact_mesh_option = Var_Translator(Param.ContactMeshOption)
        self.test_type = Param.TestType
        self.automatic_bounding_box_option = Var_Translator(Param.AutomaticBoundingBoxOption)

        self.search_tolerance = 0.0
        self.coordination_number = 10.0
        if (hasattr(Param, "LocalResolutionMethod")):
          if(Param.LocalResolutionMethod == "hierarchical"):
            self.local_resolution_method = 1
          elif(Param.LocalResolutionMethod == "area_distribution"):
            self.local_resolution_method = 2
          else:
            self.local_resolution_method = 1
        else:
          self.local_resolution_method = 1
          
        self.amplified_continuum_search_radius_extension = Param.AmplifiedSearchRadiusExtension

        if (Param.DeltaOption == "None"):
            self.delta_option = 0

        elif (Param.DeltaOption == "Absolute"):
            self.delta_option = 1
            self.search_tolerance = Param.SearchTolerance

        elif (Param.DeltaOption == "Coordination_Number"):
            self.delta_option = 2
            self.coordination_number = Param.CoordinationNumber
            self.search_tolerance = 0.01 * Param.MeanRadius

       
        if(self.delta_option > 0):
           self.case_option = 2     #MSIMSI. only 2 cases, with delta or without but continuum always.

                
        self.fixed_vel_top = Param.LoadingVelocityTop
        self.fixed_vel_bot = Param.LoadingVelocityBot

        # BOUNDING_BOX
        self.enlargement_factor = Param.BoundingBoxEnlargementFactor
        self.top_corner = Array3()
        self.bottom_corner = Array3()
        self.top_corner[0] = Param.BoundingBoxMaxX
        self.top_corner[0] = Param.BoundingBoxMaxY
        self.top_corner[0] = Param.BoundingBoxMaxZ
        self.bottom_corner[0] = Param.BoundingBoxMinX
        self.bottom_corner[0] = Param.BoundingBoxMinY
        self.bottom_corner[0] = Param.BoundingBoxMinZ

        # MODEL
        self.model_part = model_part
        self.fem_model_part = fem_model_part
        self.contact_model_part = ModelPart("ContactModelPart")  # funcio kratos
        self.cluster_model_part = cluster_model_part
        self.inlet_model_part = inlet_model_part        
        self.domain_size = Param.Dimension


        # GLOBAL PHYSICAL ASPECTS
        self.gravity = Vector(3)
        self.gravity[0] = Param.GravityX
        self.gravity[1] = Param.GravityY
        self.gravity[2] = Param.GravityZ
        
        if (Param.MaterialModel == "Linear"):
            self.force_calculation_type_id = 0
        elif (Param.MaterialModel == "Hertz"):
            self.force_calculation_type_id = 1
        elif (Param.MaterialModel == "1DPlasticity"):
            self.force_calculation_type_id = 2
        elif (Param.MaterialModel == "ExpHard"):
            self.force_calculation_type_id = 3


        if (self.force_calculation_type_id == 3):
            self.donze_g1 = Param.G1
            self.donze_g2 = Param.G2
            self.donze_g3 = Param.G3
            self.donze_max_def = Param.MaxDef

        if (Param.LocalContactDamping == "Both"):         
            self.damp_id = 11
              
        elif (Param.LocalContactDamping == "Normal"):
            self.damp_id = 10
            
        elif (Param.LocalContactDamping == "Tangential"):
            self.damp_id = 1
        else:
            self.damp_id = 0
            
        self.dempack_option = Var_Translator(Param.Dempack)
                         
        self.rolling_friction_option = 0
        if (Var_Translator(Param.RollingFrictionOption)):
            self.rolling_friction_option = 1

        if (Param.FailureCriterionType == "Mohr-Coulomb"):
            self.failure_criterion_option = 1

        elif (Param.FailureCriterionType == "Uncoupled"):
            self.failure_criterion_option = 2


        # PRINTING VARIABLES
        self.print_export_id = Var_Translator(Param.PostExportId)
        self.print_export_skin_sphere = Var_Translator(Param.PostExportSkinSphere)

        self.dummy_switch = 0

        # TIME RELATED PARAMETERS

        self.delta_time = Param.MaxTimeStep
        self.max_delta_time = Param.MaxTimeStep
        self.final_time = Param.FinalTime

        # RESOLUTION METHODS AND PARAMETERS
        self.n_step_search = int(Param.NeighbourSearchFrequency)
        self.safety_factor = Param.DeltaTimeSafetyFactor  # For critical time step

        # CREATOR-DESTRUCTOR
        self.creator_destructor = creator_destructor       

        # STRATEGIES

        self.search_strategy = OMP_DEMSearch()

        if (Param.IntegrationScheme == 'Forward_Euler'):
            self.time_integration_scheme = ForwardEulerScheme()
        elif (Param.IntegrationScheme == 'Mid_Point_Rule'):
            self.time_integration_scheme = MidPointScheme()
        else:
            print('scheme not defined')

    #

    def Initialize(self):

        # Setting ProcessInfo variables

        # SIMULATION FLAGS
        self.model_part.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, self.virtual_mass_option)
        self.model_part.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, self.critical_time_option)
        self.model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_option)
        self.model_part.ProcessInfo.SetValue(TRIHEDRON_OPTION, self.trihedron_option)
        self.model_part.ProcessInfo.SetValue(ROTATION_OPTION, self.rotation_option)
        self.model_part.ProcessInfo.SetValue(BOUNDING_BOX_OPTION, self.bounding_box_option)
        self.model_part.ProcessInfo.SetValue(DEMPACK_OPTION, self.dempack_option)
        self.model_part.ProcessInfo.SetValue(SEARCH_CONTROL, self.search_control)
        self.model_part.ProcessInfo.SetValue(FIX_VELOCITIES_FLAG, self.fix_velocities_flag)
        self.model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0)
        self.model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_option)
        self.model_part.ProcessInfo.SetValue(STRESS_STRAIN_OPTION, self.self_strain_option)

        # GLOBAL PHISICAL ASPECTS
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)

        # GLOBAL MATERIAL PROPERTIES

        self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)

        # PRINTING VARIABLES
        self.model_part.ProcessInfo.SetValue(PRINT_EXPORT_ID, self.print_export_id)
        self.model_part.ProcessInfo.SetValue(PRINT_SKIN_SPHERE, self.print_export_skin_sphere)
        self.model_part.ProcessInfo.SetValue(DAMP_TYPE, self.damp_id)
        self.model_part.ProcessInfo.SetValue(ROLLING_FRICTION_OPTION, self.rolling_friction_option)

        # TIME RELATED PARAMETERS
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)
        
        for properties in self.model_part.Properties:            
            self.ModifyProperties(properties)
            
        for properties in self.inlet_model_part.Properties:            
            self.ModifyProperties(properties)
                    
        # CONTINUUM

        self.model_part.ProcessInfo.SetValue(SEARCH_TOLERANCE, self.search_tolerance)
        self.model_part.ProcessInfo.SetValue(AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION, self.amplified_continuum_search_radius_extension)
        self.model_part.ProcessInfo.SetValue(LOCAL_RESOLUTION_METHOD, self.local_resolution_method)

        self.model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, self.contact_mesh_option)
        self.model_part.ProcessInfo.SetValue(FAILURE_CRITERION_OPTION, self.failure_criterion_option)

        if (self.force_calculation_type_id == 3):
            self.model_part.ProcessInfo.SetValue(DONZE_G1, self.donze_g1)
            self.model_part.ProcessInfo.SetValue(DONZE_G2, self.donze_g2)
            self.model_part.ProcessInfo.SetValue(DONZE_G3, self.donze_g3)
            self.model_part.ProcessInfo.SetValue(DONZE_MAX_DEF, self.donze_max_def)

        if ( (self.test_type == "Triaxial") or (self.test_type == "Hydrostatic")):
            self.model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 1)
            
        self.model_part.ProcessInfo.SetValue(FIXED_VEL_TOP, self.fixed_vel_top)
        self.model_part.ProcessInfo.SetValue(FIXED_VEL_BOT, self.fixed_vel_bot)
            
        # OTHERS

        self.model_part.ProcessInfo.SetValue(DUMMY_SWITCH, self.dummy_switch)        

        # RESOLUTION METHODS AND PARAMETERS
        # Creating the solution strategy
        self.settings = ExplicitSolverSettings()
        self.settings.r_model_part = self.model_part
        self.settings.contact_model_part = self.contact_model_part
        self.settings.fem_model_part = self.fem_model_part
        self.settings.inlet_model_part = self.inlet_model_part   
        self.settings.cluster_model_part = self.cluster_model_part
                
        self.cplusplus_strategy = ContinuumExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                      self.MoveMeshFlag, self.delta_option, self.search_tolerance, self.coordination_number, self.creator_destructor, self.time_integration_scheme, self.search_strategy)
                                                    
                                                      
        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy Initialize function (initializes all elements and performs other necessary tasks before iterating)

        #Setting the constitutive LAWS

    def Initial_Critical_Time(self):
        (self.cplusplus_strategy).InitialTimeStepCalculation()

    #

    def Solve(self):
        (self.cplusplus_strategy).Solve()
        
    def PrepareContactElementsForPrinting(self):
        (self.cplusplus_strategy).PrepareContactElementsForPrinting()
        
    
    def AddAdditionalVariables(self, model_part, Param):

        model_part.AddNodalSolutionStepVariable(COHESIVE_GROUP)  # Continuum group        
        model_part.AddNodalSolutionStepVariable(SKIN_SPHERE)

        # ONLY VISUALIZATION

        if (Var_Translator(Param.PostExportSkinSphere) ):
            model_part.AddNodalSolutionStepVariable(EXPORT_SKIN_SPHERE)
        #    model_part.AddNodalSolutionStepVariable(PREDEFINED_SKIN)

    def AddDofs(self, model_part):

        for node in model_part.Nodes:
            node.AddDof(DISPLACEMENT_X, REACTION_X)
            node.AddDof(DISPLACEMENT_Y, REACTION_Y)
            node.AddDof(DISPLACEMENT_Z, REACTION_Z)
            node.AddDof(VELOCITY_X, REACTION_X)
            node.AddDof(VELOCITY_Y, REACTION_Y)
            node.AddDof(VELOCITY_Z, REACTION_Z)
            node.AddDof(ANGULAR_VELOCITY_X, REACTION_X)
            node.AddDof(ANGULAR_VELOCITY_Y, REACTION_Y)
            node.AddDof(ANGULAR_VELOCITY_Z, REACTION_Z)

        print("DOFs for the DEM solution added correctly")


    def AddClusterVariables(self, model_part, Param):
        pass
    
    
    def ModifyProperties(self, properties):
        BaseExplicitStrategy.ModifyProperties(self, properties)
        
        ContinuumConstitutiveLawString = properties[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]            
        ContinuumConstitutiveLaw = globals().get(ContinuumConstitutiveLawString)()        
        ContinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties)
        
        
        
