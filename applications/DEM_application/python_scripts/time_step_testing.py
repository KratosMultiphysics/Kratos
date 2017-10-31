import KratosMultiphysics as Kratos
from KratosMultiphysics.DEMApplication import *
import main_script as DEM_main_script

import platform
import os

class TimeStepTester(object):    
    def __init__(self):
        if platform.system()=="Windows":
            os.system("setenv OMP_NUM_THREADS 1") #NOT WORKING!
        else:
            os.environ['OMP_NUM_THREADS']='1'

    def Run(self):
        list_of_time_step_values = [1e-5, 1e-4, 5e-4, 1e-3, 2e-3]
        
        for dt in list_of_time_step_values:
            self.RunTestCaseWithCustomizedDt(dt)
            
        self.Finalize()
            
    def RunTestCaseWithCustomizedDt(self, dt):
        CustomizedSolutionForTimeStepTesting(dt).Run()
        
    def Finalize(self):
        if platform.system()=="Windows":
            os.system("setenv OMP_NUM_THREADS ") # Trying to set a 'default' value
        else:
            os.system("export OMP_NUM_THREADS=") # Trying to set a 'default' value
        
        
class CustomizedSolutionForTimeStepTesting(DEM_main_script.Solution):
    
    def __init__(self, dt):
        self.customized_time_step = dt
        super(CustomizedSolutionForTimeStepTesting,self).__init__()
    
    def LoadParametersFile(self):
        self.DEM_parameters = Kratos.Parameters(
            """
            {
                "Dimension"                        : 3,
                "drag_modifier_type"               : 3,
                "project_from_particles_option"    : false,
                "consider_lift_force_option"       : false,
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
                "Dempack"                          : false,
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
                "HorizontalFixVel"                 : true,
                "IntegrationScheme"                : "Forward_Euler",
                "AutomaticTimestep"                : false,
                "DeltaTimeSafetyFactor"            : 1.0,
                "MaxTimeStep"                      : 1e-4,
                "FinalTime"                        : 5e-1,
                "ControlTime"                      : 100,
                "NeighbourSearchFrequency"         : 1,
                "PeriodicDomainOption"             : false,
                "MaterialModel"                    : "Hertz",
                "G1"                               : 0.0,
                "G2"                               : 0.0,
                "G3"                               : 0.0,
                "MaxDef"                           : 0.0,
                "FailureCriterionType"             : "Uncoupled",
                "AreaFactor"                       : 1,
                "LocalContactDamping"              : "Normal",
                "LocalDampingFactor"               : 1.0,
                "GlobalForceReduction"             : 0.0,
                "TestType"                         : "None",
                "ConfinementPressure"              : 0.0,
                "LoadingVelocityTop"               : 0.0,
                "LoadingVelocityBot"               : 0.0,
                "FemPlates"                        : false,
                "StressStrainOption"               : false,
                "MeshType"                         : "Current",
                "MeshPath"                         : "0",
                "SpecimenLength"                   : 0.30,
                "SpecimenDiameter"                 : 0.15,
                "MeasuringSurface"                 : 0.01767145867644375,
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
                "PostExportSkinSphere"             : false,
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
                "PredefinedSkinOption"             : false,
                "PostRollingResistanceMoment"      : false,
                "MeanRadius"                       : 0.0001,

                "problem_name"                     : "TimeStepTests"
                }
            """
            )
            
        self.DEM_parameters["MaxTimeStep"].SetDouble(self.customized_time_step)
            
    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        
        properties = Kratos.Properties(0)
        properties_walls = Kratos.Properties(0)
        self.SetHardcodedProperties(properties, properties_walls)                
        self.spheres_model_part.AddProperties(properties)
        self.rigid_face_model_part.AddProperties(properties_walls)
        
        
        DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]
        DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
        DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties, False)
        scheme = SymplecticEulerScheme()        
        scheme.SetIntegrationSchemeInProperties(properties, False)
                
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
            
        self.ComputeEnergy()
        
        self.rigid_face_model_part.CreateNewNode(11, -0.5, -0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(12, -0.5, -0.5,  0.5)
        
        self.rigid_face_model_part.CreateNewNode(13,  0.5, -0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(14,  0.5, -0.5,  0.5)        
       
        
        self.rigid_face_model_part.CreateNewNode(15,  0.5,  0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(16,  0.5,  0.5,  0.5)
        
        self.rigid_face_model_part.CreateNewNode(17, -0.5,  0.5, -0.5)
        self.rigid_face_model_part.CreateNewNode(18, -0.5,  0.5,  0.5)
        
        condition_name = "RigidFace3D3N"
        self.rigid_face_model_part.CreateNewCondition(condition_name, 1, [11,12,13], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 2, [12,13,14], self.rigid_face_model_part.GetProperties()[0])

        self.rigid_face_model_part.CreateNewCondition(condition_name, 3, [13,14,15], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 4, [14,15,16], self.rigid_face_model_part.GetProperties()[0])
        
        self.rigid_face_model_part.CreateNewCondition(condition_name, 5, [15,16,17], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 6, [16,17,18], self.rigid_face_model_part.GetProperties()[0])
        
        self.rigid_face_model_part.CreateNewCondition(condition_name, 7, [17,18,11], self.rigid_face_model_part.GetProperties()[0])
        self.rigid_face_model_part.CreateNewCondition(condition_name, 8, [18,11,12], self.rigid_face_model_part.GetProperties()[0])
            
        
    def ComputeEnergy(self):
        self.this_test_total_energy = 0.0
        
        for element in self.spheres_model_part.Elements:            
            self.this_test_total_energy += element.Calculate(PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY, self.spheres_model_part.ProcessInfo)
            self.this_test_total_energy += element.Calculate(PARTICLE_ROTATIONAL_KINEMATIC_ENERGY, self.spheres_model_part.ProcessInfo)
            self.this_test_total_energy += element.Calculate(PARTICLE_ELASTIC_ENERGY, self.spheres_model_part.ProcessInfo)
            
                    
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
        
        properties_walls[WALL_FRICTION] = 0.0
        properties_walls[WALL_COHESION] = 0.0
        properties_walls[COMPUTE_WEAR] = 0
        properties_walls[SEVERITY_OF_WEAR] = 0.001
        properties_walls[IMPACT_WEAR_SEVERITY] = 0.001
        properties_walls[BRINELL_HARDNESS] = 200.0
        properties_walls[Kratos.YOUNG_MODULUS] = 7.0e10
        properties_walls[Kratos.POISSON_RATIO] = 0.30                                                                                                                                                                     
        
        
    def FinalizeTimeStep(self, time):
        super(CustomizedSolutionForTimeStepTesting,self).FinalizeTimeStep(time)
        old_energy = self.this_test_total_energy
        self.ComputeEnergy()
        print(self.this_test_total_energy)
        if self.this_test_total_energy/old_energy > 1.01 :
            print("GAINING ENERGY!!")
            print("time step is:" + str(self.customized_time_step))
            lele
    
    def PrintResultsForGid(self, time):
        pass

if __name__ == '__main__':
    TimeStepTester().Run()