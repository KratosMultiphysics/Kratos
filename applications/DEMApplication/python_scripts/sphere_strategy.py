from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import math
import time
import cluster_file_reader

class ExplicitStrategy(object):

    #def __init__(self, all_model_parts, creator_destructor, dem_fem_search, scheme, DEM_parameters, procedures):
    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):

        # Initialization of member variables

        self.spheres_model_part = all_model_parts.Get("SpheresPart")
        self.inlet_model_part = all_model_parts.Get("DEMInletPart")
        self.fem_model_part = all_model_parts.Get("RigidFacePart")
        self.cluster_model_part = all_model_parts.Get("ClusterPart")
        self.contact_model_part = all_model_parts.Get("ContactPart")

        self.DEM_parameters = DEM_parameters

        if not "ComputeStressTensorOption" in DEM_parameters.keys():
            self.compute_stress_tensor_option = 0
        else:
            self.compute_stress_tensor_option = DEM_parameters["ComputeStressTensorOption"].GetBool()

        if "PostStressStrainOption" in DEM_parameters.keys() and DEM_parameters["PostStressStrainOption"].GetBool():
            self.compute_stress_tensor_option = 1
            self.print_stress_tensor_option = 1
        else:
            self.print_stress_tensor_option = 0

        if not "AutomaticTimestep" in DEM_parameters.keys():
            self.critical_time_option = 0
        else:
            self.critical_time_option = DEM_parameters["AutomaticTimestep"].GetBool() #TODO: add suffix option

        self.trihedron_option        = DEM_parameters["PostEulerAngles"].GetBool()
        self.rotation_option         = DEM_parameters["RotationOption"].GetBool()
        self.bounding_box_option     = DEM_parameters["BoundingBoxOption"].GetBool()
        self.fix_velocities_flag     = 0
        self.Procedures              = procedures
        #self.time_integration_scheme = scheme
        #self.time_integration_scheme.SetRotationOption(self.rotation_option)

        self.clean_init_indentation_option = DEM_parameters["CleanIndentationsOption"].GetBool()
        self.contact_mesh_option           = 0
        if "ContactMeshOption" in DEM_parameters.keys():
            self.contact_mesh_option      = DEM_parameters["ContactMeshOption"].GetBool()
        self.automatic_bounding_box_option = DEM_parameters["AutomaticBoundingBoxOption"].GetBool()

        self.delta_option = DEM_parameters["DeltaOption"].GetString() #TODO: this is not an option (bool) let's change the name to something including 'type'

        self.search_increment = 0.0
        self.coordination_number = 10.0
        self.case_option = 3
        self.search_control = 1

        if "LocalResolutionMethod" in DEM_parameters.keys():
            if (DEM_parameters["LocalResolutionMethod"].GetString() == "hierarchical"):
                self.local_resolution_method = 1
            elif (DEM_parameters["LocalResolutionMethod"].GetString() == "area_distribution"):
                self.local_resolution_method = 2
            else:
                self.local_resolution_method = 1
        else:
            self.local_resolution_method = 1

        if DEM_parameters["DeltaOption"].GetString() == "None":
            self.delta_option = 0

        elif DEM_parameters["DeltaOption"].GetString() == "Absolute":
            self.delta_option = 1
            self.search_increment = DEM_parameters["SearchTolerance"].GetDouble()

        elif DEM_parameters["DeltaOption"].GetString() == "Coordination_Number":
            self.delta_option = 2
            self.coordination_number = DEM_parameters["CoordinationNumber"].GetDouble()
            self.search_increment = 0.01 * 0.0001 #DEM_parameters-MeanRadius

        # TIME RELATED PARAMETERS
        self.delta_time = DEM_parameters["MaxTimeStep"].GetDouble()
        self.max_delta_time = DEM_parameters["MaxTimeStep"].GetDouble()
        self.end_time = DEM_parameters["FinalTime"].GetDouble()

        # BOUNDING_BOX
        self.enlargement_factor = DEM_parameters["BoundingBoxEnlargementFactor"].GetDouble()
        self.top_corner = Array3()
        self.bottom_corner = Array3()
        self.bottom_corner[0] = DEM_parameters["BoundingBoxMinX"].GetDouble()
        self.bottom_corner[1] = DEM_parameters["BoundingBoxMinY"].GetDouble()
        self.bottom_corner[2] = DEM_parameters["BoundingBoxMinZ"].GetDouble()
        self.top_corner[0] = DEM_parameters["BoundingBoxMaxX"].GetDouble()
        self.top_corner[1] = DEM_parameters["BoundingBoxMaxY"].GetDouble()
        self.top_corner[2] = DEM_parameters["BoundingBoxMaxZ"].GetDouble()

        if not "BoundingBoxStartTime" in DEM_parameters.keys():
            self.bounding_box_start_time  = 0.0
        else:
            self.bounding_box_start_time  = DEM_parameters["BoundingBoxStartTime"].GetDouble()

        if not "BoundingBoxStopTime" in DEM_parameters.keys():
            self.bounding_box_stop_time  = self.end_time
        else:
            self.bounding_box_stop_time  = DEM_parameters["BoundingBoxStopTime"].GetDouble()

        # GLOBAL PHYSICAL ASPECTS
        self.gravity = Vector(3)
        self.gravity[0] = DEM_parameters["GravityX"].GetDouble()
        self.gravity[1] = DEM_parameters["GravityY"].GetDouble()
        self.gravity[2] = DEM_parameters["GravityZ"].GetDouble()

        self.virtual_mass_option = 0
        self.nodal_mass_coeff = DEM_parameters["VirtualMassCoefficient"].GetDouble()

        if (self.nodal_mass_coeff != 1.0):
            self.virtual_mass_option = 1

        self.rolling_friction_option = DEM_parameters["RollingFrictionOption"].GetBool()

        if not "GlobalDamping" in DEM_parameters.keys():
            self.global_damping = 0.0
            Logger.PrintWarning("DEM", "\nGlobal Damping parameter not found! No damping will be applied...\n")
        else:
            self.global_damping = DEM_parameters["GlobalDamping"].GetDouble()

        # PRINTING VARIABLES
        self.print_export_id = DEM_parameters["PostExportId"].GetBool()
        self.print_export_skin_sphere = 0
        self.poisson_ratio_option = 0

        # RESOLUTION METHODS AND PARAMETERS
        self.n_step_search = DEM_parameters["NeighbourSearchFrequency"].GetInt() #TODO: NeighbourSearchFrequency change name to something that includes number of steps
        self.safety_factor = DEM_parameters["DeltaTimeSafetyFactor"].GetDouble()  # For critical time step @53214

        # CREATOR-DESTRUCTOR
        self.creator_destructor = creator_destructor
        self.dem_fem_search = dem_fem_search

        # STRATEGIES
        self.search_strategy = OMP_DEMSearch()
        if "PeriodicDomainOption" in DEM_parameters.keys():
            if DEM_parameters["PeriodicDomainOption"].GetBool():
                self.search_strategy = OMP_DEMSearch(DEM_parameters["BoundingBoxMinX"].GetDouble(),
                                                     DEM_parameters["BoundingBoxMinY"].GetDouble(),
                                                     DEM_parameters["BoundingBoxMinZ"].GetDouble(),
                                                     DEM_parameters["BoundingBoxMaxX"].GetDouble(),
                                                     DEM_parameters["BoundingBoxMaxY"].GetDouble(),
                                                     DEM_parameters["BoundingBoxMaxZ"].GetDouble())


        self.SetContinuumType()
        self.do_search_neighbours = True # Hard-coded until needed as an option

    def SetContinuumType(self):
        self.continuum_type = False

    def Var_Translator(self, variable):

        if (variable == "OFF" or variable == "0" or variable == 0 or variable == "No"):
            variable = 0
        else:
            variable = 1

        return variable

    def SetOneOrZeroInProcessInfoAccordingToBoolValue(self, model_part, variable, bool_value): #TODO: to be removed, because the Kratos variables should be bools already
        if bool_value:
            model_part.ProcessInfo.SetValue(variable, 1)
        else:
            model_part.ProcessInfo.SetValue(variable, 0)

    def SetVariablesAndOptions(self):

        # Setting ProcessInfo variables

        # SIMULATION FLAGS
        self.spheres_model_part.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, self.virtual_mass_option)
        self.spheres_model_part.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, self.critical_time_option)
        self.spheres_model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_option)
        self.spheres_model_part.ProcessInfo.SetValue(TRIHEDRON_OPTION, self.trihedron_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, ROTATION_OPTION, self.rotation_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, BOUNDING_BOX_OPTION, self.bounding_box_option)
        self.spheres_model_part.ProcessInfo.SetValue(SEARCH_CONTROL, self.search_control)
        self.spheres_model_part.ProcessInfo.SetValue(FIX_VELOCITIES_FLAG, self.fix_velocities_flag)
        self.spheres_model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0)
        self.spheres_model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_option)
        self.spheres_model_part.ProcessInfo.SetValue(BOUNDING_BOX_START_TIME, self.bounding_box_start_time)
        self.spheres_model_part.ProcessInfo.SetValue(BOUNDING_BOX_STOP_TIME, self.bounding_box_stop_time)
        self.spheres_model_part.ProcessInfo.SetValue(COMPUTE_STRESS_TENSOR_OPTION, self.compute_stress_tensor_option)
        self.spheres_model_part.ProcessInfo.SetValue(PRINT_STRESS_TENSOR_OPTION, self.print_stress_tensor_option)
        self.spheres_model_part.ProcessInfo.SetValue(CONTINUUM_OPTION, self.continuum_type)

        self.spheres_model_part.ProcessInfo.SetValue(DOMAIN_IS_PERIODIC, 0) #TODO: DOMAIN_IS_PERIODIC should be a bool, and should have the suffix option
        if "PeriodicDomainOption" in self.DEM_parameters.keys():
            if self.DEM_parameters["PeriodicDomainOption"].GetBool():
                self.spheres_model_part.ProcessInfo.SetValue(DOMAIN_IS_PERIODIC, 1) #TODO: DOMAIN_IS_PERIODIC should be a bool, and should have the suffix option

        self.spheres_model_part.ProcessInfo.SetValue(DOMAIN_MIN_CORNER, self.bottom_corner)
        self.spheres_model_part.ProcessInfo.SetValue(DOMAIN_MAX_CORNER, self.top_corner)
        self.spheres_model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)

        # GLOBAL MATERIAL PROPERTIES
        self.spheres_model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, ROLLING_FRICTION_OPTION, self.rolling_friction_option)
        self.spheres_model_part.ProcessInfo.SetValue(GLOBAL_DAMPING, self.global_damping)

        # SEARCH-RELATED
        self.search_increment_for_walls = self.search_increment # for the moment, until all bugs have been removed
        self.spheres_model_part.ProcessInfo.SetValue(SEARCH_RADIUS_INCREMENT, self.search_increment)
        self.spheres_model_part.ProcessInfo.SetValue(SEARCH_RADIUS_INCREMENT_FOR_WALLS, self.search_increment_for_walls)
        self.spheres_model_part.ProcessInfo.SetValue(COORDINATION_NUMBER, self.coordination_number)
        self.spheres_model_part.ProcessInfo.SetValue(LOCAL_RESOLUTION_METHOD, self.local_resolution_method)
        if self.contact_mesh_option:
            self.spheres_model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, 1)
        else:
            self.spheres_model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, 0)

        # PRINTING VARIABLES

        self.spheres_model_part.ProcessInfo.SetValue(PRINT_EXPORT_ID, self.print_export_id)

        # TIME RELATED PARAMETERS
        self.spheres_model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)

        #-----os.chdir('..')   # check functionality

        for properties in self.spheres_model_part.Properties:
            self.ModifyProperties(properties)

        for properties in self.inlet_model_part.Properties:
            self.ModifyProperties(properties)

        for properties in self.cluster_model_part.Properties:
            self.ModifyProperties(properties)

        for properties in self.fem_model_part.Properties:
            self.ModifyProperties(properties, 1)

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

        strategy_parameters = self.DEM_parameters["strategy_parameters"]

        if (self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Velocity_Verlet'):
            self.cplusplus_strategy = IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.search_strategy, strategy_parameters, self.do_search_neighbours)
        else:
            self.cplusplus_strategy = ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                             self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                             self.search_strategy, strategy_parameters, self.do_search_neighbours)

    def AddVariables(self):
        pass

    def BeforeInitialize(self):
        self.CreateCPlusPlusStrategy()
        self.RebuildListOfDiscontinuumSphericParticles()
        self.SetNormalRadiiOnAllParticles()
        self.SetSearchRadiiOnAllParticles()

    def Initialize(self):
        self.CheckMomentumConservation()
        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy (C++) Initialize function (initializes all elements and performs other necessary tasks before starting the time loop in Python)

    def SetDt(self, dt):
        self.dt = dt

    def Predict(self):
        pass

    def Check(self):
        pass

    def Solve(self): # deprecated
        self.SolveSolutionStep()

    def SolveSolutionStep(self):
        time = self.spheres_model_part.ProcessInfo[TIME]
        self.FixDOFsManually(time)
        (self.cplusplus_strategy).ResetPrescribedMotionFlagsRespectingImposedDofs()
        self.FixExternalForcesManually(time)
        (self.cplusplus_strategy).Solve()

    def AdvanceInTime(self, step, time):
        step += 1
        time = time + self.dt

        return step, time

    def FinalizeSolutionStep(self):
        pass

    def SetNormalRadiiOnAllParticles(self):
        (self.cplusplus_strategy).SetNormalRadiiOnAllParticles(self.spheres_model_part)

    def SetSearchRadiiOnAllParticles(self):
        (self.cplusplus_strategy).SetSearchRadiiOnAllParticles(self.spheres_model_part, self.search_increment, 1.0)

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

    def AddClusterVariables(self, spheres_model_part, DEM_parameters):
        pass

    def AddDofs(self, spheres_model_part):

        for node in spheres_model_part.Nodes:
            node.AddDof(VELOCITY_X)
            node.AddDof(VELOCITY_Y)
            node.AddDof(VELOCITY_Z)
            node.AddDof(ANGULAR_VELOCITY_X)
            node.AddDof(ANGULAR_VELOCITY_Y)
            node.AddDof(ANGULAR_VELOCITY_Z)

        Logger.Print("DOFs for the DEM solution added correctly", label="DEM")

    def PrepareElementsForPrinting(self):
        (self.cplusplus_strategy).PrepareElementsForPrinting()

    def PrepareContactElementsForPrinting(self):
        (self.cplusplus_strategy).PrepareContactElementsForPrinting()

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

    def TranslationalIntegrationSchemeTranslator(self, name):
        class_name = None

        if name == 'Forward_Euler':
            class_name = 'ForwardEulerScheme'
        elif name == 'Symplectic_Euler':
            class_name = 'SymplecticEulerScheme'
        elif name == 'Taylor_Scheme':
            class_name = 'TaylorScheme'
        elif name == 'Velocity_Verlet':
            class_name = 'VelocityVerletScheme'

        return class_name

    def RotationalIntegrationSchemeTranslator(self, name_translational, name_rotational):
        class_name = None

        if name_rotational == 'Direct_Integration':
            if name_translational == 'Forward_Euler':
                class_name = 'ForwardEulerScheme'
            elif name_translational == 'Symplectic_Euler':
                class_name = 'SymplecticEulerScheme'
            elif name_translational == 'Taylor_Scheme':
                class_name = 'TaylorScheme'
            elif name_translational == 'Velocity_Verlet':
                class_name = 'VelocityVerletScheme'
        elif name_rotational == 'Runge_Kutta':
            class_name = 'RungeKuttaScheme'
        elif name_rotational == 'Quaternion_Integration':
            class_name = 'QuaternionIntegrationScheme'

        return class_name

    def GetTranslationalSchemeInstance(self, class_name):
             return globals().get(class_name)()

    def GetRotationalSchemeInstance(self, class_name):
             return globals().get(class_name)()

    def GetTranslationalScheme(self, name):
        class_name = self.TranslationalIntegrationSchemeTranslator(name)
        translational_scheme = None
        error_status = 0
        summary = ''

        if not class_name == None:
            try:
                translational_scheme = self.GetTranslationalSchemeInstance(class_name)
                return translational_scheme, error_status, summary
            except:
                error_status = 1
                summary = 'The class corresponding to the translational integration scheme named ' + name + ' has not been added to python. Please, select a different name or add the required class.'
        else:
            error_status = 2
            summary = 'The translational integration scheme name ' + name + ' does not designate any available scheme. Please, select a different one'

        return translational_scheme, error_status, summary

    def GetRotationalScheme(self, name_translational, name_rotational):
        class_name = self.RotationalIntegrationSchemeTranslator(name_translational, name_rotational)
        rotational_scheme = None
        error_status = 0
        summary = ''

        if not class_name == None:
            try:
                rotational_scheme = self.GetRotationalSchemeInstance(class_name)
                return rotational_scheme, error_status, summary
            except:
                error_status = 1
                summary = 'The class corresponding to the rotational integration scheme name ' + name + ' has not been added to python. Please, select a different name or add the required class.'
        else:
            error_status = 2
            summary = 'The rotational integration scheme name ' + name + ' does not designate any available scheme. Please, select a different one'

        return rotational_scheme, error_status, summary

    def ModifyProperties(self, properties, param = 0):

        if not param:
            DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]
            DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
            DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties, True)

            coefficient_of_restitution = properties[COEFFICIENT_OF_RESTITUTION]

            type_of_law = DiscontinuumConstitutiveLaw.GetTypeOfLaw()

            write_gamma = False

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
                if not properties.Has(BREAKABLE_CLUSTER):
                    properties.SetValue(BREAKABLE_CLUSTER, False)

        if properties.Has(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_NAME):
            translational_scheme_name = properties[DEM_TRANSLATIONAL_INTEGRATION_SCHEME_NAME]
        else:
            translational_scheme_name = self.DEM_parameters["TranslationalIntegrationScheme"].GetString()

        if properties.Has(PARTICLE_FRICTION):
            self.Procedures.KRATOSprint("---------------------------------------------------")
            self.Procedures.KRATOSprint("  WARNING: Property PARTICLE_FRICTION is deprecated ")
            self.Procedures.KRATOSprint("  since April 11th, 2018, replace with FRICTION")
            self.Procedures.KRATOSprint("  Automatic replacement is done now.")
            self.Procedures.KRATOSprint("---------------------------------------------------")
            properties[FRICTION] = properties[PARTICLE_FRICTION]
        if properties.Has(WALL_FRICTION):
            self.Procedures.KRATOSprint("-------------------------------------------------")
            self.Procedures.KRATOSprint("  WARNING: Property WALL_FRICTION is deprecated")
            self.Procedures.KRATOSprint("  since April 11th, 2018, replace with FRICTION")
            self.Procedures.KRATOSprint("  Automatic replacement is done now.")
            self.Procedures.KRATOSprint("-------------------------------------------------")
            properties[FRICTION] = properties[WALL_FRICTION]

        translational_scheme, error_status, summary_mssg = self.GetTranslationalScheme(translational_scheme_name)

        translational_scheme.SetTranslationalIntegrationSchemeInProperties(properties, True)

        if properties.Has(DEM_ROTATIONAL_INTEGRATION_SCHEME_NAME):
            rotational_scheme_name = properties[DEM_ROTATIONAL_INTEGRATION_SCHEME_NAME]
        else:
            rotational_scheme_name = self.DEM_parameters["RotationalIntegrationScheme"].GetString()

        rotational_scheme, error_status, summary_mssg = self.GetRotationalScheme(translational_scheme_name, rotational_scheme_name)
        rotational_scheme.SetRotationalIntegrationSchemeInProperties(properties, True)

        if not properties.Has(ROLLING_FRICTION_WITH_WALLS):
            properties[ROLLING_FRICTION_WITH_WALLS] = properties[ROLLING_FRICTION]

