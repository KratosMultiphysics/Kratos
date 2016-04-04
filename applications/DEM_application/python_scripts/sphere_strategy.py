from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import math
import time

class ExplicitStrategy:

    def __init__(self, model_part, fem_model_part, cluster_model_part, inlet_model_part, creator_destructor, dem_fem_search, scheme, Param, procedures):

        # Initialization of member variables
        # SIMULATION FLAGS
        self.Parameters              = Param        
        self.self_strain_option      = self.Var_Translator(Param.StressStrainOption)
        self.critical_time_option    = self.Var_Translator(Param.AutoReductionOfTimeStepOption)
        self.trihedron_option        = self.Var_Translator(Param.PostEulerAngles)
        self.rotation_option         = self.Var_Translator(Param.RotationOption)
        self.bounding_box_option     = self.Var_Translator(Param.BoundingBoxOption)
        self.fix_velocities_flag     = 0        
        self.Procedures              = procedures
        self.time_integration_scheme = scheme
        self.time_integration_scheme.SetRotationOption(self.rotation_option)

        self.clean_init_indentation_option = self.Var_Translator(Param.CleanIndentationsOption)
        self.contact_mesh_option           = self.Var_Translator(Param.ContactMeshOption)
        self.automatic_bounding_box_option = self.Var_Translator(Param.AutomaticBoundingBoxOption)

        self.delta_option = self.Var_Translator(Param.DeltaOption)

        self.search_tolerance = 0.0
        self.coordination_number = 10.0
        self.case_option = 3
        self.search_control                 = 1

        if (hasattr(Param, "LocalResolutionMethod")):
            if(Param.LocalResolutionMethod == "hierarchical"):
              self.local_resolution_method = 1
            elif(Param.LocalResolutionMethod == "area_distribution"):
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
            self.search_tolerance = 0.01 * Param.MeanRadius

        # TIME RELATED PARAMETERS
        self.delta_time = Param.MaxTimeStep
        self.max_delta_time = Param.MaxTimeStep
        self.final_time = Param.FinalTime

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

        if not (hasattr(Param, "BoundingBoxStartTime")):
            self.bounding_box_start_time  = 0.0
        else:
            self.bounding_box_start_time  = Param.BoundingBoxStartTime

        if not (hasattr(Param, "BoundingBoxStopTime")):
            self.bounding_box_stop_time  = self.final_time
        else:
            self.bounding_box_stop_time  = Param.BoundingBoxStopTime

        # MODEL
        self.model_part = model_part
        self.fem_model_part = fem_model_part
        self.cluster_model_part = cluster_model_part
        self.inlet_model_part = inlet_model_part
        self.contact_model_part = ModelPart("ContactModelPart")

        # GLOBAL PHYSICAL ASPECTS
        self.gravity = Vector(3)
        self.gravity[0] = Param.GravityX
        self.gravity[1] = Param.GravityY
        self.gravity[2] = Param.GravityZ

        self.virtual_mass_option = 0
        self.nodal_mass_coeff = Param.VirtualMassCoefficient

        if (self.nodal_mass_coeff != 1.00):
            self.virtual_mass_option = 1

        self.rolling_friction_option = self.Var_Translator(Param.RollingFrictionOption)

        # PRINTING VARIABLES
        self.print_export_id = self.Var_Translator(Param.PostExportId)
        self.print_export_skin_sphere = 0

        # RESOLUTION METHODS AND PARAMETERS
        self.n_step_search = int(Param.NeighbourSearchFrequency)
        self.safety_factor = Param.DeltaTimeSafetyFactor  # For critical time step

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
        self.model_part.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, self.virtual_mass_option)
        self.model_part.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, self.critical_time_option)
        self.model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_option)
        self.model_part.ProcessInfo.SetValue(TRIHEDRON_OPTION, self.trihedron_option)
        self.model_part.ProcessInfo.SetValue(ROTATION_OPTION, self.rotation_option)
        self.model_part.ProcessInfo.SetValue(BOUNDING_BOX_OPTION, self.bounding_box_option)
        self.model_part.ProcessInfo.SetValue(SEARCH_CONTROL, self.search_control)
        self.model_part.ProcessInfo.SetValue(FIX_VELOCITIES_FLAG, self.fix_velocities_flag)
        self.model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0)
        self.model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_option)
        self.model_part.ProcessInfo.SetValue(STRESS_STRAIN_OPTION, self.self_strain_option)
        self.model_part.ProcessInfo.SetValue(BOUNDING_BOX_START_TIME, self.bounding_box_start_time)
        self.model_part.ProcessInfo.SetValue(BOUNDING_BOX_STOP_TIME, self.bounding_box_stop_time)

        # GLOBAL PHYSICAL ASPECTS
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)

        # GLOBAL MATERIAL PROPERTIES
        self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)

        # SEARCH-RELATED
        self.model_part.ProcessInfo.SetValue(SEARCH_TOLERANCE, self.search_tolerance)
        self.model_part.ProcessInfo.SetValue(COORDINATION_NUMBER, self.coordination_number)
        self.model_part.ProcessInfo.SetValue(LOCAL_RESOLUTION_METHOD, self.local_resolution_method)

        # PRINTING VARIABLES

        self.model_part.ProcessInfo.SetValue(ROLLING_FRICTION_OPTION, self.rolling_friction_option)
        self.model_part.ProcessInfo.SetValue(PRINT_EXPORT_ID, self.print_export_id)

        # TIME RELATED PARAMETERS
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)

        for properties in self.model_part.Properties:
            self.ModifyProperties(properties)

        for properties in self.inlet_model_part.Properties:
            self.ModifyProperties(properties)

        for properties in self.cluster_model_part.Properties:
            self.ModifyProperties(properties)

        # RESOLUTION METHODS AND PARAMETERS
        # Creating the solution strategy
        self.settings = ExplicitSolverSettings()
        self.settings.r_model_part = self.model_part
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
        counter = 0
        for properties in self.model_part.Properties:
            current_discontinuum_constitutive_law_string = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]
            if ((counter > 0) and (previous_discontinuum_constitutive_law_string != current_discontinuum_constitutive_law_string)):
                self.Procedures.KRATOSprint(output_message)
                time.sleep(20) # Inserting a delay of 20 seconds so the user has ample time to read the message
                break
            previous_discontinuum_constitutive_law_string = current_discontinuum_constitutive_law_string
            counter += 1


    def Initialize(self):

        self.SetVariablesAndOptions()

        if (self.Parameters.IntegrationScheme == 'Verlet_Velocity'):
            self.cplusplus_strategy = IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.time_integration_scheme, self.search_strategy)
        else:
            self.cplusplus_strategy = ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                             self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                             self.time_integration_scheme, self.search_strategy)

        self.CheckMomentumConservation()

        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy (C++) Initialize function (initializes all elements and performs other necessary tasks before starting the time loop in Python)


    def Solve(self):
        (self.cplusplus_strategy).ResetPrescribedMotionFlagsRespectingImposedDofs()
        self.FixDOFsManually()
        self.FixExternalForcesManually()
        (self.cplusplus_strategy).Solve()

    def Compute_RigidFace_Movement(self):
        (self.cplusplus_strategy).Compute_RigidFace_Movement()


    def FixDOFsManually(self):
        pass
    
    def FixExternalForcesManually(self):
        pass
        
    def AddAdditionalVariables(self, balls_model_part, DEM_parameters):
        pass

    def AddClusterVariables(self, model_part, Param):
        pass

    def AddDofs(self, model_part):

        for node in model_part.Nodes:
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

    def ModifyProperties(self, properties):
        DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]
        DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
        DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties)

        coefficient_of_restitution = properties[COEFFICIENT_OF_RESTITUTION]

        type_of_law = DiscontinuumConstitutiveLaw.GetTypeOfLaw()

        write_gamma = False

        if (type_of_law == 'Linear'):
            gamma = self.RootByBisection(self.coeff_of_rest_diff, 0.0, 16.0, 0.0001, 300, coefficient_of_restitution)
            write_gamma = True

        elif (type_of_law == 'Hertz'):
            gamma = self.GammaForHertzThornton(coefficient_of_restitution)
            write_gamma = True

        else:
            pass

        if write_gamma == True:
            properties[DAMPING_GAMMA] = gamma
