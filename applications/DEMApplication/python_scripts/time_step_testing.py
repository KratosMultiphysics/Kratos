import KratosMultiphysics as Kratos
from KratosMultiphysics.DEMApplication import *
import main_script as DEM_main_script

class TimeStepTester(object):
    def __init__(self):
        #self.schemes_list = ["Forward_Euler", "Taylor_Scheme", "Symplectic_Euler", "Velocity_Verlet"]
        self.schemes_list = ["Symplectic_Euler", "Velocity_Verlet"]
        self.stable_time_steps_list = []

    def Run(self):

        for scheme in self.schemes_list:
            self.RunForACertainScheme(scheme)

        self.Finalize()

    def RunForACertainScheme(self, scheme):
        print("Computing stable time step for scheme: "+ scheme)
        tolerance = 1e-7
        dt = 1e-2
        previous_dt = 0.0
        while dt > previous_dt + tolerance:
            try:
                self.RunTestCaseWithCustomizedDtAndScheme(dt, scheme)
            except SystemExit:
                factor = min(0.5, 0.5*(dt-previous_dt))
                dt = factor * dt
                print("decreasing dt by " + str(factor))
                continue

            previous_dt = dt
            dt = dt * 1.5
            print("increasing dt by 1.5")

        self.stable_time_steps_list.append(previous_dt)

    @classmethod
    def RunTestCaseWithCustomizedDtAndScheme(self, dt, scheme):
        model = Kratos.Model()
        CustomizedSolutionForTimeStepTesting(model, dt, scheme).Run()

    def Finalize(self):

        print("\n")
        print("#############################")
        print("List of tested schemes:")
        print(self.schemes_list)
        print("List of stable time steps:")
        print(self.stable_time_steps_list)
        print("#############################")
        print("\n")


class CustomizedSolutionForTimeStepTesting(DEM_main_script.Solution):

    def __init__(self, model, dt, scheme):
        self.customized_time_step = dt
        self.customized_scheme = scheme
        super(CustomizedSolutionForTimeStepTesting, self).__init__(model)

    def LoadParametersFile(self):
        self.DEM_parameters = Kratos.Parameters(
            """
            {
                "Dimension"                        : 3,
                "BoundingBoxOption"                : true,
                "BoundingBoxEnlargementFactor"     : 1.1,
                "AutomaticBoundingBoxOption"       : false,
                "BoundingBoxEnlargementFactor"     : 1.0,
                "BoundingBoxMaxX"                  : 1e3,
                "BoundingBoxMaxY"                  : 1e3,
                "BoundingBoxMaxZ"                  : 1e3,
                "BoundingBoxMinX"                  : -1e3,
                "BoundingBoxMinY"                  : -1e3,
                "BoundingBoxMinZ"                  : -1e3,
                "dem_inlet_option"                 : false,
                "GravityX"                         : 0.0,
                "GravityY"                         : 0.0,
                "GravityZ"                         : 0.0,
                "VelocityTrapOption"               : false,
                "RotationOption"                   : true,
                "CleanIndentationsOption"          : true,
                "RemoveBallsInEmbeddedOption"      : true,
                "DeltaOption"                      : "Absolute",
                "SearchTolerance"                  : 0.0,
                "CoordinationNumber"               : 10,
                "AmplifiedSearchRadiusExtension"   : 1.10000e+00,
                "ModelDataInfo"                    : false,
                "VirtualMassCoefficient"           : 1.0,
                "RollingFrictionOption"            : false,
                "DontSearchUntilFailure"           : false,
                "ContactMeshOption"                : false,
                "OutputFileType"                   : "Binary",
                "Multifile"                        : "multiple_files",
                "IntegrationScheme"                : "Forward_Euler",
                "AutomaticTimestep"                : false,
                "DeltaTimeSafetyFactor"            : 1.0,
                "MaxTimeStep"                      : 1e-4,
                "FinalTime"                        : 4.0,
                "ControlTime"                      : 100,
                "NeighbourSearchFrequency"         : 1,
                "PeriodicDomainOption"             : false,
                "ElementType"                      : "SphericPartDEMElement3D",
                "GraphExportFreq"                  : 1e-5,
                "VelTrapGraphExportFreq"           : 1e-3,
                "OutputTimeStep"                   : 1e-3,
                "PostDisplacement"                 : false,
                "PostVelocity"                     : false,
                "PostElasticForces"                : false,
                "PostContactForces"                : false,
                "PostRigidElementForces"           : false,
                "PostTangentialElasticForces"      : false,
                "PostPressure"                     : false,
                "PostTotalForces"                  : false,
                "PostShearStress"                  : false,
                "PostNonDimensionalVolumeWear"     : false,
                "PostNodalArea"                    : false,
                "PostRHS"                          : false,
                "PostDampForces"                   : false,
                "PostAppliedForces"                : false,
                "PostRadius"                       : false,
                "PostGroupId"                      : false,
                "PostExportId"                     : false,
                "PostAngularVelocity"              : false,
                "PostParticleMoment"               : false,
                "PostEulerAngles"                  : false,
                "PostContactSigma"                 : false,
                "PostContactTau"                   : false,
                "PostLocalContactForce"            : false,
                "PostFailureCriterionState"        : false,
                "PostContactFailureId"             : false,
                "PostMeanContactArea"              : false,
                "PostStressStrainOption"           : false,
                "PostRollingResistanceMoment"      : false,
                "TranslationalIntegrationScheme"   : "Taylor_Scheme",
                "RotationalIntegrationScheme"      :  "Direct_Integration",

                "problem_name"                     : "TimeStepTests"
                }
            """
            )

        self.DEM_parameters["TranslationalIntegrationScheme"].SetString(self.customized_scheme)

        self.DEM_parameters["MaxTimeStep"].SetDouble(self.customized_time_step)

        default_input_parameters = self.GetDefaultInputParameters()
        self.DEM_parameters.ValidateAndAssignDefaults(default_input_parameters)

    def ReadModelParts(self, max_node_Id=0, max_elem_Id=0, max_cond_Id=0):

        properties = Kratos.Properties(0)
        properties_walls = Kratos.Properties(0)
        self.SetHardcodedProperties(properties, properties_walls)
        self.spheres_model_part.AddProperties(properties)
        self.rigid_face_model_part.AddProperties(properties_walls)


        DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]
        DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
        DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties, False)

        translational_scheme = ForwardEulerScheme()
        translational_scheme.SetTranslationalIntegrationSchemeInProperties(properties, True)
        rotational_scheme = ForwardEulerScheme()
        rotational_scheme.SetRotationalIntegrationSchemeInProperties(properties, True)

        element_name = "SphericParticle3D"
        PropertiesProxiesManager().CreatePropertiesProxies(self.spheres_model_part)

        coordinates = Kratos.Array3()
        coordinates[0] = 0.0
        coordinates[1] = -0.1
        coordinates[2] = 0.0
        radius = 0.1
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)

        coordinates[0] = 0.0
        coordinates[1] = 0.2
        coordinates[2] = 0.0
        radius = 0.1
        self.creator_destructor.CreateSphericParticle(self.spheres_model_part, coordinates, properties, radius, element_name)


        for node in self.spheres_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY_X, 10.0)
            node.SetSolutionStepValue(Kratos.VELOCITY_Y, 5.0)
            node.SetSolutionStepValue(Kratos.VELOCITY_Z, 0.0)

        self.initial_test_energy = self.ComputeEnergy()

        self.rigid_face_model_part.CreateNewNode(11, -0.5, -0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(12, -0.5, -0.5, 0.5)

        self.rigid_face_model_part.CreateNewNode(13, 0.5, -0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(14, 0.5, -0.5, 0.5)


        self.rigid_face_model_part.CreateNewNode(15, 0.5, 0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(16, 0.5, 0.5, 0.5)

        self.rigid_face_model_part.CreateNewNode(17, -0.5, 0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(18, -0.5, 0.5, 0.5)

        condition_name = "RigidFace3D3N"
        self.rigid_face_model_part.CreateNewCondition(condition_name, 1, [11, 12, 13], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 2, [12, 13, 14], self.rigid_face_model_part.GetProperties()[0])

        self.rigid_face_model_part.CreateNewCondition(condition_name, 3, [13, 14, 15], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 4, [14, 15, 16], self.rigid_face_model_part.GetProperties()[0])

        self.rigid_face_model_part.CreateNewCondition(condition_name, 5, [15, 16, 17], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 6, [16, 17, 18], self.rigid_face_model_part.GetProperties()[0])

        self.rigid_face_model_part.CreateNewCondition(condition_name, 7, [17, 18, 11], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 8, [18, 11, 12], self.rigid_face_model_part.GetProperties()[0])


    def ComputeEnergy(self):
        this_test_total_energy = 0.0

        for element in self.spheres_model_part.Elements:
            this_test_total_energy += element.Calculate(PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY, self.spheres_model_part.ProcessInfo)
            this_test_total_energy += element.Calculate(PARTICLE_ROTATIONAL_KINEMATIC_ENERGY, self.spheres_model_part.ProcessInfo)
            this_test_total_energy += element.Calculate(PARTICLE_ELASTIC_ENERGY, self.spheres_model_part.ProcessInfo)

        return this_test_total_energy

    @classmethod
    def SetHardcodedProperties(self, properties, properties_walls):

        properties[PARTICLE_DENSITY] = 2650.0
        properties[Kratos.YOUNG_MODULUS] = 7.0e6
        properties[Kratos.POISSON_RATIO] = 0.30
        properties[PARTICLE_FRICTION] = 0.0
        properties[PARTICLE_COHESION] = 0.0
        properties[COEFFICIENT_OF_RESTITUTION] = 1.0
        properties[Kratos.PARTICLE_MATERIAL] = 1
        properties[ROLLING_FRICTION] = 0.0
        properties[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEMContinuumConstitutiveLaw"
        properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = "DEM_D_Hertz_viscous_Coulomb"

        properties_walls[FRICTION] = 0.0
        properties_walls[WALL_COHESION] = 0.0
        properties_walls[COMPUTE_WEAR] = 0
        properties_walls[SEVERITY_OF_WEAR] = 0.001
        properties_walls[IMPACT_WEAR_SEVERITY] = 0.001
        properties_walls[BRINELL_HARDNESS] = 200.0
        properties_walls[Kratos.YOUNG_MODULUS] = 7.0e10
        properties_walls[Kratos.POISSON_RATIO] = 0.30


    def FinalizeTimeStep(self, time):
        super(CustomizedSolutionForTimeStepTesting, self).FinalizeTimeStep(time)

        current_test_energy = self.ComputeEnergy()
        #if not self.step%200:
        #    print("Energy: "+str(current_test_energy))

        if current_test_energy/self.initial_test_energy > 1.5:
            print("GAINING ENERGY!!")
            print("time step is:" + str(self.customized_time_step))
            import sys
            sys.exit()

    def PrintResultsForGid(self, time):
        pass


if __name__ == '__main__':
    TimeStepTester().Run()
