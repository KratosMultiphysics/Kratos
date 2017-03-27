from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import math
import time
import cluster_file_reader

class ExplicitStrategy:

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, scheme, Param, procedures):

        # Initialization of member variables        

        self.spheres_model_part = all_model_parts.spheres_model_part
        self.inlet_model_part = all_model_parts.DEM_inlet_model_part
        self.fem_model_part = all_model_parts.rigid_face_model_part
        self.cluster_model_part = all_model_parts.cluster_model_part
        self.contact_model_part = all_model_parts.contact_model_part
        
        self.Parameters = Param

        if not (hasattr(Param, "ComputeStressTensorOption")):
            self.compute_stress_tensor_option = 0
        else:
            self.compute_stress_tensor_option = self.Var_Translator(Param.ComputeStressTensorOption)

        if (hasattr(Param, "PostStressStrainOption") and self.Var_Translator(Param.PostStressStrainOption)):
            self.compute_stress_tensor_option = 1
            self.print_stress_tensor_option = 1
        else:
            self.print_stress_tensor_option = 0

        if not hasattr(Param, "AutomaticTimestep"):
            self.critical_time_option = 0
        else:
            self.critical_time_option = self.Var_Translator(Param.AutomaticTimestep)
                   
        self.trihedron_option        = self.Var_Translator(Param.PostEulerAngles)
        self.rotation_option         = self.Var_Translator(Param.RotationOption)
        self.bounding_box_option     = self.Var_Translator(Param.BoundingBoxOption)
        self.fix_velocities_flag     = 0
        self.Procedures              = procedures
        self.time_integration_scheme = scheme
        #self.time_integration_scheme.SetRotationOption(self.rotation_option)

        self.clean_init_indentation_option = self.Var_Translator(Param.CleanIndentationsOption)
        self.contact_mesh_option           = 0
        if (hasattr(Param, "ContactMeshOption")):
            self.contact_mesh_option      = self.Var_Translator(Param.ContactMeshOption)
        self.automatic_bounding_box_option = self.Var_Translator(Param.AutomaticBoundingBoxOption)

        self.delta_option = self.Var_Translator(Param.DeltaOption)

        self.search_tolerance = 0.0
        self.coordination_number = 10.0
        self.case_option = 3
        self.search_control = 1

        if (hasattr(Param, "LocalResolutionMethod")):
            if (Param.LocalResolutionMethod == "hierarchical"):
                self.local_resolution_method = 1
            elif (Param.LocalResolutionMethod == "area_distribution"):
                self.local_resolution_method = 2
            else:
                self.local_resolution_method = 1
        else:
            self.local_resolution_method = 1

        if (Param.DeltaOption == "None"):
            self.delta_option = 0

        elif (Param.DeltaOption == "Absolute"):
            self.delta_option = 1
            self.search_tolerance = Param.SearchTolerance

        elif (Param.DeltaOption == "Coordination_Number"):
            self.delta_option = 2
            self.coordination_number = Param.CoordinationNumber
            self.search_tolerance = 0.01 * 0.0001 #Param.MeanRadius

        # TIME RELATED PARAMETERS
        self.delta_time = Param.MaxTimeStep
        self.max_delta_time = Param.MaxTimeStep
        self.final_time = Param.FinalTime

        # BOUNDING_BOX
        self.enlargement_factor = Param.BoundingBoxEnlargementFactor
        self.top_corner = Array3()
        self.bottom_corner = Array3()
        self.top_corner[0] = Param.BoundingBoxMaxX
        self.top_corner[1] = Param.BoundingBoxMaxY
        self.top_corner[2] = Param.BoundingBoxMaxZ
        self.bottom_corner[0] = Param.BoundingBoxMinX
        self.bottom_corner[1] = Param.BoundingBoxMinY
        self.bottom_corner[2] = Param.BoundingBoxMinZ

        if not (hasattr(Param, "BoundingBoxStartTime")):
            self.bounding_box_start_time  = 0.0
        else:
            self.bounding_box_start_time  = Param.BoundingBoxStartTime

        if not (hasattr(Param, "BoundingBoxStopTime")):
            self.bounding_box_stop_time  = self.final_time
        else:
            self.bounding_box_stop_time  = Param.BoundingBoxStopTime

        # GLOBAL PHYSICAL ASPECTS
        self.gravity = Vector(3)
        self.gravity[0] = Param.GravityX
        self.gravity[1] = Param.GravityY
        self.gravity[2] = Param.GravityZ

        self.virtual_mass_option = 0
        self.nodal_mass_coeff = Param.VirtualMassCoefficient

        if (self.nodal_mass_coeff != 1.0):
            self.virtual_mass_option = 1

        self.rolling_friction_option = self.Var_Translator(Param.RollingFrictionOption)

        if not (hasattr(Param, "GlobalDamping")):
            self.global_damping = 0.0
            print("\nGlobal Damping parameter not found! No damping will be applied...\n")
        else:
            self.global_damping = Param.GlobalDamping

        # PRINTING VARIABLES
        self.print_export_id = self.Var_Translator(Param.PostExportId)
        self.print_export_skin_sphere = 0
        self.poisson_ratio_option = 0

        # RESOLUTION METHODS AND PARAMETERS
        self.n_step_search = int(Param.NeighbourSearchFrequency)
        self.safety_factor = Param.DeltaTimeSafetyFactor  # For critical time step @53214

        # CREATOR-DESTRUCTOR
        self.creator_destructor = creator_destructor
        self.dem_fem_search = dem_fem_search

        # STRATEGIES
        self.search_strategy = OMP_DEMSearch(-1.0, -1.0, -1.0)
        if (hasattr(Param, "PeriodicDomainOption")):
            if self.Var_Translator(Param.PeriodicDomainOption):
                self.search_strategy = OMP_DEMSearch(Param.BoundingBoxMaxX-Param.BoundingBoxMinX, Param.BoundingBoxMaxY-Param.BoundingBoxMinY, Param.BoundingBoxMaxZ-Param.BoundingBoxMinZ)

        self.SetContinuumType()

    def SetContinuumType(self):
        self.continuum_type = False

    def Var_Translator(self, variable):

        if (variable == "OFF" or variable == "0" or variable == 0 or variable == "No"):
            variable = 0
        else:
            variable = 1

        return variable

    def SetVariablesAndOptions(self):

        # Setting ProcessInfo variables

        # SIMULATION FLAGS
        self.spheres_model_part.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, self.virtual_mass_option)
        self.spheres_model_part.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, self.critical_time_option)
        self.spheres_model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_option)
        self.spheres_model_part.ProcessInfo.SetValue(TRIHEDRON_OPTION, self.trihedron_option)
        self.spheres_model_part.ProcessInfo.SetValue(ROTATION_OPTION, self.rotation_option)
        self.spheres_model_part.ProcessInfo.SetValue(BOUNDING_BOX_OPTION, self.bounding_box_option)
        self.spheres_model_part.ProcessInfo.SetValue(SEARCH_CONTROL, self.search_control)
        self.spheres_model_part.ProcessInfo.SetValue(FIX_VELOCITIES_FLAG, self.fix_velocities_flag)
        self.spheres_model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0)
        self.spheres_model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_option)
        self.spheres_model_part.ProcessInfo.SetValue(BOUNDING_BOX_START_TIME, self.bounding_box_start_time)
        self.spheres_model_part.ProcessInfo.SetValue(BOUNDING_BOX_STOP_TIME, self.bounding_box_stop_time)
        self.spheres_model_part.ProcessInfo.SetValue(COMPUTE_STRESS_TENSOR_OPTION, self.compute_stress_tensor_option)
        self.spheres_model_part.ProcessInfo.SetValue(PRINT_STRESS_TENSOR_OPTION, self.print_stress_tensor_option)
        self.spheres_model_part.ProcessInfo.SetValue(CONTINUUM_OPTION, self.continuum_type)

        # GLOBAL PHYSICAL ASPECTS
        self.spheres_model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)

        # GLOBAL MATERIAL PROPERTIES
        self.spheres_model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)
        self.spheres_model_part.ProcessInfo.SetValue(ROLLING_FRICTION_OPTION, self.rolling_friction_option)
        self.spheres_model_part.ProcessInfo.SetValue(GLOBAL_DAMPING, self.global_damping)

        # SEARCH-RELATED
        self.do_search_neighbours = True # Hard-coded until needed as an option
        self.spheres_model_part.ProcessInfo.SetValue(SEARCH_TOLERANCE, self.search_tolerance)
        self.spheres_model_part.ProcessInfo.SetValue(COORDINATION_NUMBER, self.coordination_number)
        self.spheres_model_part.ProcessInfo.SetValue(LOCAL_RESOLUTION_METHOD, self.local_resolution_method)

        # PRINTING VARIABLES

        self.spheres_model_part.ProcessInfo.SetValue(PRINT_EXPORT_ID, self.print_export_id)

        # TIME RELATED PARAMETERS
        self.spheres_model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)

        os.chdir("..")
        for properties in self.spheres_model_part.Properties:
            self.ModifyProperties(properties)

        for properties in self.inlet_model_part.Properties:
            self.ModifyProperties(properties)

        for properties in self.cluster_model_part.Properties:
            self.ModifyProperties(properties)
        
        # RESOLUTION METHODS AND PARAMETERS
        # Creating the solution strategy
        self.settings = ExplicitSolverSettings()
        self.settings.r_model_part = self.spheres_model_part
        self.settings.contact_model_part = self.contact_model_part
        self.settings.fem_model_part = self.fem_model_part
        self.settings.inlet_model_part = self.inlet_model_part
        self.settings.cluster_model_part = self.cluster_model_part

    def CheckMomentumConservation(self):

        previous_discontinuum_constitutive_law_string = ""
        current_discontinuum_constitutive_law_string = ""
        output_message = "\n*********************************************************************************************\n"+\
                           "*********************************************************************************************\n"+\
                           "WARNING:                                                                                     \n"+\
                           "A mix of constitutive laws is being used. The momentum conservation law will not be fulfilled\n"+\
                           "in the interaction of particles of different constitutive laws. Please, cancel the simulation\n"+\
                           "if that fact represents a problem. Otherwise, computations will resume in a few seconds.     \n"+\
                           "*********************************************************************************************\n"+\
                           "*********************************************************************************************\n"
        delay = 1  # Inserting the desired delay in seconds
        counter = 0
        for properties in self.spheres_model_part.Properties:
            current_discontinuum_constitutive_law_string = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]
            if ((counter > 0) and (previous_discontinuum_constitutive_law_string != current_discontinuum_constitutive_law_string)):
                self.Procedures.KRATOSprint(output_message)
                time.sleep(delay) # Inserting a delay so the user has ample time to read the message
                break
            previous_discontinuum_constitutive_law_string = current_discontinuum_constitutive_law_string
            counter += 1


    def CreateCPlusPlusStrategy(self):

        self.SetVariablesAndOptions()

        if (self.Parameters.IntegrationScheme == 'Verlet_Velocity'):
            self.cplusplus_strategy = IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.time_integration_scheme, self.search_strategy, self.do_search_neighbours) 
                                                              #TODO: remove time_integration_scheme. no longer necessary and maybe safety_factor
        else:
            self.cplusplus_strategy = ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                             self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                             self.time_integration_scheme, self.search_strategy, self.do_search_neighbours)
                                                             #TODO: remove time_integration_scheme. no longer necessary
                                
    def BeforeInitialize(self):
        self.CreateCPlusPlusStrategy()
        self.RebuildListOfDiscontinuumSphericParticles()
        self.SetNormalRadiiOnAllParticles()
        self.SetSearchRadiiOnAllParticles()
        
    def Initialize(self):                                                                     
        self.CheckMomentumConservation()
        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy (C++) Initialize function (initializes all elements and performs other necessary tasks before starting the time loop in Python)


    def Solve(self):
        time = self.spheres_model_part.ProcessInfo[TIME]
        self.FixDOFsManually(time)
        (self.cplusplus_strategy).ResetPrescribedMotionFlagsRespectingImposedDofs()
        self.FixExternalForcesManually(time)
        (self.cplusplus_strategy).Solve()
        
    def SetNormalRadiiOnAllParticles(self):
        (self.cplusplus_strategy).SetNormalRadiiOnAllParticles(self.spheres_model_part)
        
    def SetSearchRadiiOnAllParticles(self):
        (self.cplusplus_strategy).SetSearchRadiiOnAllParticles(self.spheres_model_part, self.search_tolerance, 1.0)
        
    def RebuildListOfDiscontinuumSphericParticles(self):
        (self.cplusplus_strategy).RebuildListOfDiscontinuumSphericParticles()

    def Compute_RigidFace_Movement(self):
        (self.cplusplus_strategy).Compute_RigidFace_Movement()


    def FixDOFsManually(self,time):
    #if time>1.0:
        #for node in self.spheres_model_part.Nodes:
            #node.Fix(VELOCITY_X)
            #node.SetSolutionStepValue(VELOCITY_X, 0.0)
        pass

    def FixExternalForcesManually(self,time):
        pass

    def AddAdditionalVariables(self, balls_model_part, DEM_parameters):
        pass

    def AddClusterVariables(self, spheres_model_part, Param):
        pass

    def AddDofs(self, spheres_model_part):

        for node in spheres_model_part.Nodes:
            node.AddDof(VELOCITY_X, REACTION_X)
            node.AddDof(VELOCITY_Y, REACTION_Y)
            node.AddDof(VELOCITY_Z, REACTION_Z)
            node.AddDof(ANGULAR_VELOCITY_X, REACTION_X)
            node.AddDof(ANGULAR_VELOCITY_Y, REACTION_Y)
            node.AddDof(ANGULAR_VELOCITY_Z, REACTION_Z)

        print("DOFs for the DEM solution added correctly")

    def PrepareElementsForPrinting(self):
        (self.cplusplus_strategy).PrepareElementsForPrinting()

    def coeff_of_rest_diff(self, gamma, desired_coefficient_of_restit):

        if gamma <= 1.0/math.sqrt(2.0) :
            return math.exp(-gamma/math.sqrt(1.0-gamma*gamma)*(math.pi-math.atan(2.0*gamma*math.sqrt(1.0-gamma*gamma)/(-2.0*gamma*gamma+1.0))))-desired_coefficient_of_restit
        elif gamma < 1.0 :
            return math.exp(-gamma/math.sqrt(1.0-gamma*gamma)*math.atan(2.0*gamma*math.sqrt(1.0-gamma*gamma)/(2.0*gamma*gamma-1.0)))-desired_coefficient_of_restit
        elif gamma == 1.0 :
            return 0.135335283 - desired_coefficient_of_restit
        else:
            return math.exp(-gamma/math.sqrt(gamma*gamma-1.0)*math.log((gamma/math.sqrt(gamma*gamma-1.0)+1.0)/(gamma/math.sqrt(gamma*gamma-1.0)-1.0)))-desired_coefficient_of_restit


    def RootByBisection(self, f, a, b, tol, maxiter, coefficient_of_restitution):

        if coefficient_of_restitution < 0.001 :
            coefficient_of_restitution = 0.001

        if coefficient_of_restitution > 0.999 :
            return 0.0
        k=0
        gamma = 0.5 * (a + b)

        while b - a > tol and k <= maxiter:
            coefficient_of_restitution_trial = self.coeff_of_rest_diff(gamma, coefficient_of_restitution)

            if self.coeff_of_rest_diff(a, coefficient_of_restitution) * coefficient_of_restitution_trial < 0:
                b = gamma

            elif coefficient_of_restitution_trial == 0:
                return gamma

            else:
                a = gamma

            gamma = 0.5 * (a + b)
            k += 1

        return gamma

    def GammaForHertzThornton(self, e):

        if e < 0.001:
            e = 0.001

        if e > 0.999:
            return 0.0

        h1  = -6.918798
        h2  = -16.41105
        h3  =  146.8049
        h4  = -796.4559
        h5  =  2928.711
        h6  = -7206.864
        h7  =  11494.29
        h8  = -11342.18
        h9  =  6276.757
        h10 = -1489.915

        alpha = e*(h1+e*(h2+e*(h3+e*(h4+e*(h5+e*(h6+e*(h7+e*(h8+e*(h9+e*h10)))))))))

        return math.sqrt(1.0/(1.0 - (1.0+e)*(1.0+e) * math.exp(alpha)) - 1.0)

    def IntegrationSchemeTranslator(self, name):
        class_name = None

        if name == 'Forward_Euler':
            class_name = 'ForwardEulerScheme'
        elif name == 'Symplectic_Euler':
            class_name = 'SymplecticEulerScheme'
        elif name == 'Taylor_Scheme':
            class_name = 'TaylorScheme'        
        elif name == 'Newmark_Beta_Method':
            class_name = 'NewmarkBetaScheme'
        elif name == 'Verlet_Velocity':
            class_name = 'VerletVelocityScheme'

        return class_name

    def GetSchemeInstance(self, class_name):
        return globals().get(class_name)()

    def GetScheme(self, name):
        class_name = self.IntegrationSchemeTranslator(name)
        scheme = None
        error_status = 0
        summary = ''

        if not class_name == None:
            try:
                scheme = self.GetSchemeInstance(class_name)
                return scheme, error_status, summary
            except:
                error_status = 1
                summary = 'The class corresponding to the scheme name (' + name + ') has not been added to python. Please, select a different name or add the required class.'
        else:
            error_status = 2
            summary = 'The scheme name (' + name + ') does not designate any available scheme. Please, select a different one'

        return scheme, error_status, summary

    def ModifyProperties(self, properties):
        DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]
        DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
        DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties)

        coefficient_of_restitution = properties[COEFFICIENT_OF_RESTITUTION]

        type_of_law = DiscontinuumConstitutiveLaw.GetTypeOfLaw()

        write_gamma = False
        write_AlphaFunction = False

        if (type_of_law == 'Linear'):
            gamma = self.RootByBisection(self.coeff_of_rest_diff, 0.0, 16.0, 0.0001, 300, coefficient_of_restitution)
            write_gamma = True

        elif (type_of_law == 'Hertz' or type_of_law == 'Dependent_friction'):
            gamma = self.GammaForHertzThornton(coefficient_of_restitution)
            write_gamma = True

        else:
            pass

        if write_gamma == True:
            properties[DAMPING_GAMMA] = gamma
            
        if properties.Has(CLUSTER_FILE_NAME):
            cluster_file_name = properties[CLUSTER_FILE_NAME]
            [name, list_of_coordinates, list_of_radii, size, volume, inertias] = cluster_file_reader.ReadClusterFile(cluster_file_name)
            pre_utils = PreUtilities(self.spheres_model_part)
            pre_utils.SetClusterInformationInProperties(name, list_of_coordinates, list_of_radii, size, volume, inertias, properties)
            self.Procedures.KRATOSprint(properties)
            
        if properties.Has(DEM_INTEGRATION_SCHEME_NAME):  
            scheme_name = properties[DEM_INTEGRATION_SCHEME_NAME]
        else:
            scheme_name = self.Parameters.IntegrationScheme
            
        scheme, error_status, summary_mssg = self.GetScheme(scheme_name)
        scheme.SetIntegrationSchemeInProperties(properties)
        
