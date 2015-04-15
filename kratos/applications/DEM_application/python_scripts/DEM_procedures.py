from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import math
import DEM_explicit_solver_var as DEM_parameters
import DEM_material_test_script
import os
import shutil

def Var_Translator(variable):

    if (variable == "OFF" or variable == "0" or variable == 0):
        variable = 0
    else:
        variable = 1

    return variable    


class MdpaCreator(object):

    def __init__(self, path, DEM_parameters):
        
        self.problem_parameters = param
        self.current_path = path

        # Creating necessary directories

        self.post_mdpas = str(self.current_path) + '/' + str(self.problem_parameters.problem_name) + '_post_mdpas'
        os.chdir(self.current_path)
        if not os.path.isdir(self.post_mdpas):
            os.makedirs(str(self.post_mdpas))

    def WriteMdpa(self, model_part):
        
        os.chdir(self.post_mdpas)
        time = model_part.ProcessInfo.GetValue(TIME)
        mdpa = open(str(self.problem_parameters.problem_name) + '_post_' + str(time) + '.mdpa', 'w')
        mdpa.write('\n')
        mdpa.write('Begin ModelPartData')
        mdpa.write('\n')
        mdpa.write('//  VARIABLE_NAME value')
        mdpa.write('\n')
        mdpa.write('End ModelPartData')
        mdpa.write('\n')
        mdpa.write('\n')
        mdpa.write('\n')
        mdpa.write('\n')
        mdpa.write('Begin Nodes')
        mdpa.write('\n')
        
        for node in model_part.Nodes:
            mdpa.write(str(node.Id).rjust(12) + ' ' + str(node.X).rjust(8) + ' ' + str(node.Y).rjust(8) + ' ' + str(node.Z).rjust(8))
            mdpa.write('\n')

        mdpa.write('End Nodes')
        mdpa.write('\n')
        mdpa.write('\n')
        mdpa.write('Begin NodalData RADIUS')
        mdpa.write('\n')

        for node in model_part.Nodes:
            mdpa.write(str(node.Id) + ' ' + str(0) + ' ' + str(node.GetSolutionStepValue(RADIUS)))
            mdpa.write('\n')

        mdpa.write('End NodalData')
        mdpa.write('\n')


class GranulometryUtils(object):

    def __init__(self, domain_volume, model_part):

        if (domain_volume <= 0.0):
            raise ValueError("Error: The input domain volume must be strictly positive")

        self.spheres_model_part = model_part
        self.UpdateData(domain_volume)

    def UpdateData(self, domain_volume):
        
        self.physics_calculator = SphericElementGlobalPhysicsCalculator(self.spheres_model_part)
        self.number_of_spheres    = self.spheres_model_part.NumberOfElements(0)
        self.solid_volume       = self.physics_calculator.CalculateTotalVolume(self.spheres_model_part)
        self.d_50               = self.physics_calculator.CalculateD50(self.spheres_model_part)

        if (self.number_of_spheres == 0):
            self.spheres_per_area = 0.0

        else:
            self.spheres_per_area              = domain_volume / self.number_of_spheres

        self.voids_volume       = domain_volume - self.solid_volume
        self.global_porosity    = self.voids_volume / domain_volume

    def PrintCurrentData(self):
        
        print("number_of_spheres: ", self.number_of_spheres)
        print("solid volume: ", self.solid_volume)
        print("voids volume: ", self.voids_volume)
        print("global porosity: ", self.global_porosity)
        print("D50: ", self.d_50)
        print("spheres per area unit: ", self.spheres_per_area)


class PostUtils(object):

    def __init__(self, DEM_parameters, spheres_model_part):
        
        self.DEM_parameters = DEM_parameters
        self.spheres_model_part = spheres_model_part
        self.post_utilities = PostUtilities()

        self.vel_trap_graph_counter = 0
        self.dem_io = DEMIo()
        self.vel_trap_graph_frequency = int(self.dem_io.VelTrapGraphExportFreq/spheres_model_part.ProcessInfo.GetValue(DELTA_TIME))
        if self.vel_trap_graph_frequency < 1:
            self.vel_trap_graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time
    def ComputeMeanVelocitiesinTrap(self, file_name, time_dem):

        if (self.DEM_parameters.VelocityTrapOption):
            if (self.vel_trap_graph_counter == self.vel_trap_graph_frequency):
                self.vel_trap_graph_counter = 0
            average_velocity = Array3()
            low_point = Array3()

            low_point[0] = self.DEM_parameters.VelocityTrapMinX
            low_point[1] = self.DEM_parameters.VelocityTrapMinY
            low_point[2] = self.DEM_parameters.VelocityTrapMinZ
            high_point = Array3()
            high_point[0] = self.DEM_parameters.VelocityTrapMaxX
            high_point[1] = self.DEM_parameters.VelocityTrapMaxY
            high_point[2] = self.DEM_parameters.VelocityTrapMaxZ

            average_velocity = self.post_utilities.VelocityTrap(self.spheres_model_part, low_point, high_point)
            f = open(file_name, 'a')
            tmp = str(time_dem) + "   " + str(average_velocity[0]) + "   " + str(average_velocity[1]) + "   " + str(average_velocity[2]) + "\n"
            f.write(tmp)
            f.flush()
            f.close()

            self.vel_trap_graph_counter += 1

class Procedures(object):

    def __init__(self, DEM_parameters):

        # GLOBAL VARIABLES OF THE SCRIPT
        # Defining list of skin particles (For a test tube of height 30 cm and diameter 15 cm)

        # Initialization of member variables
        # SIMULATION FLAGS
        self.rotation_OPTION               = Var_Translator(DEM_parameters.RotationOption)
        self.bounding_box_OPTION           = Var_Translator(DEM_parameters.BoundingBoxOption)
        self.automatic_bounding_box_OPTION = Var_Translator(DEM_parameters.AutomaticBoundingBoxOption)
        self.contact_mesh_OPTION           = Var_Translator(DEM_parameters.ContactMeshOption)
        self.arlequin                      = 0
        if (hasattr(DEM_parameters, "arlequin")):
          self.arlequin                    = Var_Translator(DEM_parameters.arlequin)
        #self.solver = solver
        
        # SIMULATION SETTINGS
        self.b_box_minX = DEM_parameters.BoundingBoxMinX
        self.b_box_minY = DEM_parameters.BoundingBoxMinY
        self.b_box_minZ = DEM_parameters.BoundingBoxMinZ
        self.b_box_maxX = DEM_parameters.BoundingBoxMaxX
        self.b_box_maxY = DEM_parameters.BoundingBoxMaxY
        self.b_box_maxZ = DEM_parameters.BoundingBoxMaxZ
        self.bounding_box_enlargement_factor = DEM_parameters.BoundingBoxEnlargementFactor
       
        # MODEL
        self.domain_size = DEM_parameters.Dimension
 

    def AddCommonVariables(self, model_part, DEM_parameters):
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
        if(DEM_parameters.PostGroupId):
            model_part.AddNodalSolutionStepVariable(GROUP_ID)   
            
    def AddSpheresVariables(self, model_part, Param):

        # KINEMATIC
        model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
        model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
        
        # FORCES
        model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
        model_part.AddNodalSolutionStepVariable(CONTACT_FORCES)
        model_part.AddNodalSolutionStepVariable(RIGID_ELEMENT_FORCE)
        model_part.AddNodalSolutionStepVariable(DAMP_FORCES)
        model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_FORCE)

        # BASIC PARTICLE PROPERTIES
        model_part.AddNodalSolutionStepVariable(RADIUS)
        model_part.AddNodalSolutionStepVariable(NODAL_MASS)

        # ROTATION RELATED PROPERTIES
        if (Var_Translator(Param.RotationOption)):
            model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA)
            model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_DAMP_RATIO)
            if( Var_Translator(Param.RollingFrictionOption)):
                model_part.AddNodalSolutionStepVariable(ROLLING_FRICTION)

        # OTHER PROPERTIES
        model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD
                
        # LOCAL AXIS
        if (Param.PostEulerAngles == "1" or Param.PostEulerAngles == 1):
            model_part.AddNodalSolutionStepVariable(EULER_ANGLES)

        # FLAGS

        if(Var_Translator(Param.PostGroupId)):
            model_part.AddNodalSolutionStepVariable(GROUP_ID)            # Differenced groups for plotting, etc..
        # ONLY VISUALIZATION
        if (Var_Translator(Param.PostExportId)):
            model_part.AddNodalSolutionStepVariable(EXPORT_ID)

        if (Var_Translator(Param.PostGroupId)):
            model_part.AddNodalSolutionStepVariable(EXPORT_GROUP_ID)
            
    def AddRigidFaceVariables(self, model_part, Param):

        model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
        model_part.AddNodalSolutionStepVariable(CONTACT_FORCES)
        model_part.AddNodalSolutionStepVariable(PRESSURE)
        model_part.AddNodalSolutionStepVariable(TANGENTIAL_ELASTIC_FORCES)
        model_part.AddNodalSolutionStepVariable(SHEAR_STRESS)
        model_part.AddNodalSolutionStepVariable(NODAL_AREA)
        if (Var_Translator(Param.PostNonDimensionalVolumeWear)):
            model_part.AddNodalSolutionStepVariable(NON_DIMENSIONAL_VOLUME_WEAR)
            model_part.AddNodalSolutionStepVariable(IMPACT_WEAR)
        
    def AddElasticFaceVariables(self, model_part, Param): #Only used in CSM coupling
        self.AddRigidFaceVariables(model_part,Param)
        model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
    
    def AddMappingVariables(self, model_part, Param): 
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(DUMMY_1)
        model_part.AddNodalSolutionStepVariable(DUMMY_1)
        model_part.AddNodalSolutionStepVariable(DUMMY_2)
        model_part.AddNodalSolutionStepVariable(DUMMY_3)
        model_part.AddNodalSolutionStepVariable(DUMMY_3D_1)
        model_part.AddNodalSolutionStepVariable(DUMMY_3D_2)
        model_part.AddNodalSolutionStepVariable(DUMMY_3D_3)
     
    def AddClusterVariables(self, model_part, DEM_parameters):
        # KINEMATIC
        model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
        model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
        model_part.AddNodalSolutionStepVariable(EULER_ANGLES)
        
        # FORCES
        model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
        model_part.AddNodalSolutionStepVariable(RIGID_ELEMENT_FORCE)
        model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
        
        model_part.AddNodalSolutionStepVariable(PRINCIPAL_MOMENTS_OF_INERTIA)
        model_part.AddNodalSolutionStepVariable(NODAL_MASS)     
        model_part.AddNodalSolutionStepVariable(CHARACTERISTIC_LENGTH)
        model_part.AddNodalSolutionStepVariable(RADIUS)
        model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)
        # Not really necessary but it is filled by inlet
        model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)

    def AddMpiVariables(self, model_part):
        pass

    def ModelData(self, spheres_model_part, contact_model_part, solver):
        
        # Previous Calculations.
        Model_Data = open('Model_Data.txt', 'w')

        # mean radius, and standard deviation:
        i = 0.0
        sum_radi = 0.0
        partial_sum_squared = 0.0
        total_sum_squared = 0.0
        volume = 0.0
        area = 0.0
        mean = 0.0
        var  = 0.0
        rel_std_dev = 0.0

        for node in spheres_model_part.Nodes:

            sum_radi += node.GetSolutionStepValue(RADIUS)
            partial_sum_squared = node.GetSolutionStepValue(RADIUS) ** 2.0
            total_sum_squared += partial_sum_squared
            volume += 4 * 3.141592 / 3 * node.GetSolutionStepValue(RADIUS) ** 3.0
            area += 3.141592 * partial_sum_squared
            i += 1.0

        if (i>0.0):
            mean = sum_radi / i
            var = total_sum_squared / i - mean ** 2.0
        std_dev = 0.0

        if(abs(var) > 1e-9):
            std_dev = var ** 0.5

        if (i>0.0):
            rel_std_dev = std_dev / mean

        Model_Data.write("Radius Mean: " + str(mean) + '\n')
        Model_Data.write("Std Deviation: " + str(std_dev) + '\n')
        Model_Data.write("Relative Std Deviation: " + str(rel_std_dev) + '\n')
        Model_Data.write("Total Particle Volume 3D: " + str(volume) + '\n')
        Model_Data.write("Total Particle Area 2D: " + str(area) + '\n')
        Model_Data.write('\n')

        Total_Particles = len(spheres_model_part.Nodes)

        Total_Contacts = 0

        Coordination_Number = 0.0

        if(self.contact_mesh_OPTION):

            for bar in contact_model_part.Elements:

                Total_Contacts += 1.0

            Coordination_Number = 1.0 * (Total_Contacts * 2.0) / (Total_Particles)

        Model_Data.write("Total Number of Particles: " + str(Total_Particles) + '\n')
        Model_Data.write("Total Number of Contacts: " + str(Total_Contacts) + '\n')
        Model_Data.write("Coordination Number NC: " + str(Coordination_Number) + '\n')
        Model_Data.write('\n')

        #Model_Data.write("Volume Elements: " + str(total_volume) + '\n')

        Model_Data.close()

        return Coordination_Number

    def MeasureBOT(self, solver):
        
        tol = 2.0
        y_mean = 0.0
        counter = 0.0

        for node in self.BOT:
            r = node.GetSolutionStepValue(RADIUS)
            y = node.Y
            y_mean += (y - r) * r
            counter += r

        return (y_mean, counter)

    def MeasureTOP(self, solver):
        
        tol = 2.0
        y_mean = 0.0
        counter = 0.0

        for node in self.TOP:
            r = node.GetSolutionStepValue(RADIUS)
            y = node.Y

            y_mean += (y + r) * r
            counter += r

        return (y_mean, counter)

    def MonitorPhysicalProperties(self, model_part, physics_calculator, properties_list):
        
        # This function returns a list of arrays (also lists)
        # Each array contains the values of the physical properties at the current time
        time = model_part.ProcessInfo.GetValue(TIME)
        present_prop = []

        if (len(properties_list) == 0):  # The first array in the list only contains the entries names
            names = []
            names.append("time")
            names.append("mass")
            names.append("gravitational_energy")
            names.append("kinetic_energy")
            names.append("elastic_energy")
            names.append("momentum")
            names.append("angular_momentum")
            names.append("total_energy")

            properties_list.append(names)

        # Calculating current values
        mass = physics_calculator.CalculateTotalMass(model_part)
        center = physics_calculator.CalculateCenterOfMass(model_part)
        initial_center = physics_calculator.GetInitialCenterOfMass()
        gravity_energy = physics_calculator.CalculateGravitationalPotentialEnergy(model_part, initial_center)
        kinetic_energy = physics_calculator.CalculateKineticEnergy(model_part)
        elastic_energy = physics_calculator.CalculateElasticEnergy(model_part)
        momentum = physics_calculator.CalculateTotalMomentum(model_part)
        angular_momentum = physics_calculator.CalulateTotalAngularMomentum(model_part)
        total_energy = gravity_energy + kinetic_energy + elastic_energy

        # Filling in the entries values corresponding to the entries names above
        present_prop.append(time)
        present_prop.append(mass)
        present_prop.append(gravity_energy)
        present_prop.append(kinetic_energy)
        present_prop.append(elastic_energy)
        present_prop.append(momentum)
        present_prop.append(angular_momentum)
        present_prop.append(total_energy)

        properties_list.append(present_prop)

        return properties_list

    # def PlotPhysicalProperties(self, properties_list, path):

    # This function creates one graph for each physical property.
    # properties_list[0][0] = 'time'
    # properties_list[0][j] = 'property_j'
    # properties_list[i][j] = value of property_j at time properties_list[i][0]

        # n_measures     = len(properties_list)
        # entries        = properties_list[0]
        # n_entries      = len(entries)
        # time_vect      = []
        # os.chdir(path)

        # for j in range(1, n_measures):
            # time_vect.append(properties_list[j][0])

        # for i in range(1, n_entries):
            # prop_vect_i = []

            # for j in range(1, n_measures):
                # prop_i_j = properties_list[j][i]

                # if (hasattr(prop_i_j, '__getitem__')): # Checking if it is an iterable object (a vector). If yes, take the modulus
                    # mod_prop_i_j = 0.0

                    # for k in range(len(prop_i_j)):
                        # mod_prop_i_j += prop_i_j[k] * prop_i_j[k]

                    # prop_i_j = sqrt(mod_prop_i_j) # Euclidean norm

                # prop_vect_i.append(prop_i_j)

            # plt.figure(i)
            # plot = plt.plot(time_vect, prop_vect_i)
            # plt.xlabel(entries[0])
            # plt.ylabel(entries[i])
            # plt.title('Evolution of ' + entries[i] + ' in time')
            # plt.savefig(entries[i] + '.pdf')

    def SetPredefinedSkin(self, spheres_model_part):

        for element in spheres_model_part.Elements:

            if (element.GetNode(0).GetSolutionStepValue(PREDEFINED_SKIN) > 0.0):  # PREDEFINED_SKIN is a double

                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE, 1)
                
    def SetCustomSkin(self,spheres_model_part):
    
        for element in spheres_model_part.Elements:

            x = element.GetNode(0).X
            y = element.GetNode(0).Y
            #z = element.GetNode(0).Z
          
            if(x>21.1):
                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE,1)
            if(x<1.25):
                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE,1)
            if(y>1.9):
                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE,1)
            if(y<0.1):
                element.GetNode(0).SetSolutionStepValue(SKIN_SPHERE,1)

    def CreateDirectories(self, main_path, problem_name):
        
        root             = main_path + '/' + problem_name
        post_path        = root + '_Post_Files'
        list_path        = root + '_Post_Lists'
        data_and_results = root + '_Results_and_Data'
        graphs_path      = root + '_Graphs'
        MPI_results      = root + '_MPI_results'

        shutil.rmtree(main_path + '/' + problem_name + '_Post_Files', ignore_errors = True)
        shutil.rmtree(main_path + '/' + problem_name + '_Graphs', ignore_errors = True)
        
        for directory in [post_path, list_path, data_and_results, graphs_path, MPI_results]:
            if not os.path.isdir(directory):
                os.makedirs(str(directory))

        return [post_path,list_path,data_and_results,graphs_path,MPI_results]
    
    def FindMaxNodeIdInModelPart(self, model_part):

        maxid = 0

        for node in model_part.Nodes:
            if (node.Id > maxid):
                maxid = node.Id

        return maxid

    
    def SetBoundingBox(self, spheres_model_part, clusters_model_part, rigid_faces_model_part, creator_destructor):
        
        b_box_low = Array3()
        b_box_high = Array3()
        b_box_low[0] = self.b_box_minX
        b_box_low[1] = self.b_box_minY
        b_box_low[2] = self.b_box_minZ
        b_box_high[0] = self.b_box_maxX
        b_box_high[1] = self.b_box_maxY
        b_box_high[2] = self.b_box_maxZ
               
        creator_destructor.SetLowNode(b_box_low)
        creator_destructor.SetHighNode(b_box_high)
        
        creator_destructor.CalculateSurroundingBoundingBox(spheres_model_part, clusters_model_part, rigid_faces_model_part, self.bounding_box_enlargement_factor, self.automatic_bounding_box_OPTION)
            
            
    def PreProcessModel(self, DEM_parameters):
        pass

    def KRATOSprint(self,message):
        print(message)
        sys.stdout.flush()
        
        
                
# #~CHARLIE~# Aixo no ho entenc 
class DEMFEMProcedures(object):

    def __init__(self, DEM_parameters, graphs_path, spheres_model_part, RigidFace_model_part):

        # GLOBAL VARIABLES OF THE SCRIPT
        self.TestType = DEM_parameters.TestType

        # Initialization of member variables
        # SIMULATION FLAGS
        self.rotation_OPTION     = Var_Translator(DEM_parameters.RotationOption)
        self.bounding_box_OPTION = Var_Translator(DEM_parameters.BoundingBoxOption)
        self.contact_mesh_OPTION = Var_Translator(DEM_parameters.ContactMeshOption)

        self.graphs_path = graphs_path
        self.spheres_model_part = spheres_model_part
        self.RigidFace_model_part = RigidFace_model_part
        #self.solver = solver

        self.fem_mesh_nodes = []

        self.graph_counter = 0
        self.balls_graph_counter = 0

        self.graph_frequency        = int(DEM_parameters.GraphExportFreq/spheres_model_part.ProcessInfo.GetValue(DELTA_TIME))
        if self.graph_frequency < 1:
            self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time
        os.chdir(self.graphs_path)
        #self.graph_forces = open(DEM_parameters.problem_name +"_force_graph.grf", 'w')                
        
        def open_graph_files(self,RigidFace_model_part):
            #os.chdir(self.graphs_path)
            for mesh_number in range(1, self.RigidFace_model_part.NumberOfMeshes()):
                if(self.RigidFace_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]): 
                    self.graph_forces[self.RigidFace_model_part.GetMesh((mesh_number))[IDENTIFIER]] = open(str(DEM_parameters.problem_name) + "_" + str( self.RigidFace_model_part.GetMesh((mesh_number))[IDENTIFIER]) + "_force_graph.grf", 'w')
        
        self.graph_forces = {}  
        
        
        
        
        def open_balls_graph_files(self,spheres_model_part):
            #os.chdir(self.graphs_path)
            for mesh_number in range(1, self.spheres_model_part.NumberOfMeshes()):
                if(self.spheres_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]): 
                 self.particle_graph_forces[self.spheres_model_part.GetMesh((mesh_number))[TOP]] = open(str(DEM_parameters.problem_name) + "_" + str( self.spheres_model_part.GetMesh((mesh_number))[TOP]) + "_particle_force_graph.grf", 'w');
        
        self.particle_graph_forces = {}      
        
        
        
        
        if self.TestType == "None":  
            print (self.graphs_path)
            open_graph_files(self,RigidFace_model_part)            
            open_balls_graph_files(self,spheres_model_part)   

        # SIMULATION SETTINGS

        self.bounding_box_enlargement_factor = DEM_parameters.BoundingBoxEnlargementFactor
        # MODEL
        self.domain_size = DEM_parameters.Dimension

        # PRINTING VARIABLES
        self.print_radius           = Var_Translator(DEM_parameters.PostRadius)
        self.print_velocity         = Var_Translator(DEM_parameters.PostVelocity)
        self.print_angular_velocity = Var_Translator(DEM_parameters.PostAngularVelocity)
        self.print_displacement     = Var_Translator(DEM_parameters.PostDisplacement)
        self.print_total_forces     = Var_Translator(DEM_parameters.PostTotalForces)
        self.print_damp_forces      = Var_Translator(DEM_parameters.PostDampForces)
        self.print_applied_forces   = Var_Translator(DEM_parameters.PostAppliedForces)
        self.print_particle_moment  = Var_Translator(DEM_parameters.PostParticleMoment)
        self.print_euler_angles     = Var_Translator(DEM_parameters.PostEulerAngles)
        self.print_group_id         = Var_Translator(DEM_parameters.PostGroupId)
        self.print_export_id        = Var_Translator(DEM_parameters.PostExportId)
        self.self_strain_option     = Var_Translator(DEM_parameters.StressStrainOption)

        if ((DEM_parameters.ElementType == "SphericContPartDEMElement3D") or(DEM_parameters.ElementType == "CylinderContPartDEMElement3D")):
            self.print_export_skin_sphere = Var_Translator(DEM_parameters.PostExportSkinSphere)
            self.predefined_skin_option = Var_Translator(DEM_parameters.PredefinedSkinOption)
            if (self.contact_mesh_OPTION):
                self.print_local_contact_force = Var_Translator(DEM_parameters.PostLocalContactForce)
                self.print_failure_criterion_state = Var_Translator(DEM_parameters.PostFailureCriterionState)
                #self.print_unidimensional_damage = Var_Translator(DEM_parameters.PostUnidimensionalDamage)
                self.print_contact_failure = Var_Translator(DEM_parameters.PostContactFailureId)
                self.print_contact_tau = Var_Translator(DEM_parameters.PostContactTau)
                self.print_contact_sigma = Var_Translator(DEM_parameters.PostContactSigma)
                self.print_mean_contact_area = Var_Translator(DEM_parameters.PostMeanContactArea)
        else:
            self.print_export_skin_sphere = Var_Translator(DEM_parameters.PostExportSkinSphere)
            
    
                    
    def close_graph_files(self,RigidFace_model_part):
        for mesh_number in range(1, self.RigidFace_model_part.NumberOfMeshes()):
            if(self.RigidFace_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]): 
                self.graph_forces[self.RigidFace_model_part.GetMesh((mesh_number))[IDENTIFIER]].close()
                
                
    def close_balls_graph_files(self,spheres_model_part):
        for mesh_number in range(1, self.spheres_model_part.NumberOfMeshes()):
            if(self.spheres_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]): 
                self.particle_graph_forces[self.spheres_model_part.GetMesh((mesh_number))[TOP]].close()            
                
                

    def MeasureForces(self):    # not used atm
        
        if self.TestType == "None":
            if self.RigidFace_model_part.NumberOfMeshes() > 1:
                for mesh_number in range(1, self.RigidFace_model_part.NumberOfMeshes()):
                    if(self.RigidFace_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]):
                        self.fem_mesh_nodes = self.RigidFace_model_part.GetMesh(mesh_number).Nodes
                        self.total_force_x = 0.0
                        self.total_force_y = 0.0
                        self.total_force_z = 0.0

                        for node in self.fem_mesh_nodes:
                            force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                            force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                            force_node_z = node.GetSolutionStepValue(ELASTIC_FORCES)[2]
                            self.total_force_x += force_node_x
                            self.total_force_y += force_node_y
                            self.total_force_z += force_node_z


    def PrintGraph(self, time):

        if DEM_parameters.TestType == "None":            
            if(self.graph_counter == self.graph_frequency):
                self.graph_counter = 0
                if self.RigidFace_model_part.NumberOfMeshes() > 1:
                    
                    for mesh_number in range(1, self.RigidFace_model_part.NumberOfMeshes()):
                        if(self.RigidFace_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]):
                            mesh_nodes = self.RigidFace_model_part.GetMesh(mesh_number).Nodes
                            
                            total_force = Array3()

                            total_force[0] = 0.0
                            total_force[1] = 0.0
                            total_force[2] = 0.0
                            
                            PostUtilities().IntegrationOfForces(mesh_nodes, total_force)
                                                                                      
                            self.graph_forces[self.RigidFace_model_part.GetMesh((mesh_number))[IDENTIFIER]].write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_force[0]).rjust(13)+" "+str("%.6g"%total_force[1]).rjust(13)+" "+str("%.6g"%total_force[2]).rjust(13)+"\n")
                            self.graph_forces[self.RigidFace_model_part.GetMesh((mesh_number))[IDENTIFIER]].flush()

            self.graph_counter += 1

    def FinalizeGraphs(self,RigidFace_model_part):

        if DEM_parameters.TestType == "None":
            self.close_graph_files(RigidFace_model_part)
            
            
    def PrintBallsGraph(self, time):

        if self.TestType == "None":
            if(self.balls_graph_counter == self.graph_frequency):
                self.balls_graph_counter = 0
                if self.spheres_model_part.NumberOfMeshes() > 1:
                    
                    for mesh_number in range(1, self.spheres_model_part.NumberOfMeshes()):
                        if(self.spheres_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]):
                            self.mesh_nodes = self.spheres_model_part.GetMesh(mesh_number).Nodes
                            self.total_force_x = 0.0
                            self.total_force_y = 0.0
                            self.total_force_z = 0.0

                            for node in self.mesh_nodes:
                                force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                                force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                                force_node_z = node.GetSolutionStepValue(ELASTIC_FORCES)[2]
                                self.total_force_x += force_node_x
                                self.total_force_y += force_node_y
                                self.total_force_z += force_node_z                    

                            self.particle_graph_forces[self.spheres_model_part.GetMesh((mesh_number))[TOP]].write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_x).rjust(13)+" "+str("%.6g"%self.total_force_y).rjust(13)+" "+str("%.6g"%self.total_force_z).rjust(13)+"\n")
                            self.particle_graph_forces[self.spheres_model_part.GetMesh((mesh_number))[TOP]].flush()                    


            self.balls_graph_counter += 1        
            
    def FinalizeBallsGraphs(self,spheres_model_part):

        if DEM_parameters.TestType == "None":
            self.close_balls_graph_files(spheres_model_part)
            
            
            
    def ApplyNodalRotation(self,time):

            #if (time < 0.5e-2 ) :
            if (time < 0.5e-1 ) :    
                
                #while ( time < DEM_parameters.FinalTime):
                #print("TIME STEP BEGINS.  STEP:"+str(time)+"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            
                vz = 10

                for mesh_number in range(1, self.spheres_model_part.NumberOfMeshes()):
                    if(self.spheres_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]):
                        self.mesh_nodes = self.spheres_model_part.GetMesh(mesh_number).Nodes

                        for node in self.mesh_nodes:

                              node.SetSolutionStepValue(VELOCITY_Z, vz)
                              node.Fix(VELOCITY_Z)

            else: 
                d0 = 1.694
                avance = 0.50100
                distance = (d0 - avance)
                rpm = 628

                vx = distance * rpm * math.cos(rpm * time)
                #vx = distance * rpm * math.cos(1.0)
                vz = distance * rpm * math.sin(rpm * time)
                #vz = - distance * rpm * math.sin(1.0)

                for mesh_number in range(1, self.spheres_model_part.NumberOfMeshes()):
                    if(self.spheres_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]):
                        self.mesh_nodes = self.spheres_model_part.GetMesh(mesh_number).Nodes

                        for node in self.mesh_nodes:

                            node.SetSolutionStepValue(VELOCITY_X, vx)
                            node.SetSolutionStepValue(VELOCITY_Z, vz)
                            node.Fix(VELOCITY_X)
                            node.Fix(VELOCITY_Z)

           
    def ApplyMovementbySteps(self,time):

                    
            if (time < 0.001e-0 ) :    
                
                #while ( time < DEM_parameters.FinalTime):
                #print("TIME STEP BEGINS.  STEP:"+str(time)+"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            
                vy = 0.1
                for mesh_number in range(1, self.RigidFace_model_part.NumberOfMeshes()):
                    if(self.RigidFace_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]):
                        self.mesh_nodes = self.RigidFace_model_part.GetMesh(mesh_number).Nodes

                        for node in self.mesh_nodes:

                              node.SetSolutionStepValue(VELOCITY_Y, vy)
                              node.Fix(VELOCITY_Y)

            if (time > 0.001e-0 ) and (time < 1e-1 ) :    
                
                vy = 10
                for mesh_number in range(1, self.RigidFace_model_part.NumberOfMeshes()):
                    if(self.RigidFace_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]):
                        self.mesh_nodes = self.RigidFace_model_part.GetMesh(mesh_number).Nodes

                        for node in self.mesh_nodes:

                              node.SetSolutionStepValue(VELOCITY_Y, vy)
                              node.Fix(VELOCITY_Y)

            else :    
                
                vy = 0.1
                for mesh_number in range(1, self.RigidFace_model_part.NumberOfMeshes()):
                    if(self.RigidFace_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]):
                        self.mesh_nodes = self.RigidFace_model_part.GetMesh(mesh_number).Nodes

                        for node in self.mesh_nodes:

                              node.SetSolutionStepValue(VELOCITY_Y, vy)
                              node.Fix(VELOCITY_Y)




class Report(object):

    def __init__(self):
        pass

    def Prepare(self,timer,control_time):
        self.initial_pr_time      = timer.clock()
        self.initial_re_time      = timer.time()
        self.prev_time            = 0.0
        self.total_steps_expected = 0
        self.control_time         = control_time
        self.first_print          = True

    def BeginReport(self, timer):

        report = "Main loop starting..." + "\n" + "Total number of TIME STEPs expected in the calculation: " + str(self.total_steps_expected) + "\n"

        return report

    def StepiReport(self, timer, time, step):
        
        incremental_time = (timer.time() - self.initial_re_time) - self.prev_time

        report = ""

        if (incremental_time > self.control_time):
        
            percentage = 100 * (float(step) / self.total_steps_expected)

            report = report + "Real time calculation: " + str(timer.time() - self.initial_re_time) + "\n"\
                + "Simulation time: " + str(time) + "\n"\
                + "%s %.5f %s"%("Percentage Completed: ", percentage,"%") + "\n"\
                + "Time Step: " + str(step) + "\n"
        
            self.prev_time  = (timer.time() - self.initial_re_time)
        
        if ((timer.time() - self.initial_re_time > 60) and self.first_print == True and step != 0):
            
            self.first_print = False
            estimated_sim_duration = 60.0 * (self.total_steps_expected / step) # seconds

            report = report + "The calculation total estimated time is " + str(estimated_sim_duration) + "seconds" + "\n"\
                + "in minutes:" + str(estimated_sim_duration / 60.0)    + "min." + "\n"\
                + "in hours:"   + str(estimated_sim_duration / 3600.0)  + "hrs." + "\n"\
                + "in days:"    + str(estimated_sim_duration / 86400.0) + "days" + "\n" 

            if ((estimated_sim_duration / 86400) > 2.0):
                report = report +"WARNING!!!:       VERY LASTING CALCULATION" + "\n"

        return report

    def FinalReport(self, timer):
        elapsed_pr_time = timer.clock() - self.initial_pr_time
        elapsed_re_time = timer.time()  - self.initial_re_time

        report = "Calculation ends at instant: "              + str(timer.time()) + "\n"\
            + "Calculation ends at processing time instant: " + str(timer.clock()) + "\n"\
            + "Elapsed processing time: "                     + str(elapsed_pr_time) + "\n"\
            + "Elapsed real time: "                           + str(elapsed_re_time) +"\n"

        report = report + "ANALYSIS COMPLETED" + "\n"

        return report


class MaterialTest(object):

    def __init__(self):
        pass

    def Initialize(self, DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part):
        self.TestType = DEM_parameters.TestType

        if (self.TestType != "None"):
            self.script = DEM_material_test_script.MaterialTest(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)
            self.script.Initialize()
 
    def PrepareDataForGraph(self):
        if (self.TestType != "None"):
            self.script.PrepareDataForGraph()

    def MeasureForcesAndPressure(self):
        if (self.TestType != "None"):
            self.script.MeasureForcesAndPressure()
            
    def PrintGraph(self, time):
        if (self.TestType != "None"):
            self.script.PrintGraph(time)

    def FinalizeGraphs(self):
        if (self.TestType != "None"):
            self.script.FinalizeGraphs()
            
            
    def PrintChart(self):        
        if (self.TestType != "None"):
            self.script.PrintChart()  
            
    def ApplyMovementbySteps(self, time):        
        if (self.TestType != "None"):
            self.script.ApplyMovementbySteps(time)  
    
class MultifileList(object):

    def __init__(self,name,step):
        self.index = 0
        self.step = step
        self.name = name
        self.file = open("_list_"+self.name+"_"+str(step)+".post.lst","w")


class DEMIo(object):

    def __init__(self):
        # Printing variables
        self.global_variables            = []
        self.spheres_variables           = []
        self.sphere_local_axis_variables = []
        self.fem_boundary_variables      = []
        self.mapping_variables           = []
        self.cluster_variables           = []
        self.contact_variables           = []
        self.multifilelists              = []
        
        # Reading Post options from DEM_parameters  
        self.PostDisplacement             = getattr(DEM_parameters, "PostDisplacement", 0 )
        self.PostVelocity                 = getattr(DEM_parameters, "PostVelocity", 0 )
        self.PostTotalForces              = getattr(DEM_parameters, "PostTotalForces", 0 )
        self.PostNonDimensionalVolumeWear = getattr(DEM_parameters, "PostNonDimensionalVolumeWear", 0 )
        self.PostImpactWear               = getattr(DEM_parameters, "PostImpactWear", 0 )
        self.PostAppliedForces            = getattr(DEM_parameters, "PostAppliedForces", 0 )
        self.PostDampForces               = getattr(DEM_parameters, "PostDampForces", 0 )
        self.PostRadius                   = getattr(DEM_parameters, "PostRadius", 0 )
        self.PostExportId                 = getattr(DEM_parameters, "PostExportId", 0 )
        self.PostExportSkinSphere         = getattr(DEM_parameters, "PostExportSkinSphere", 0 )
        self.PostAngularVelocity          = getattr(DEM_parameters, "PostAngularVelocity", 0 )
        self.PostParticleMoment           = getattr(DEM_parameters, "PostParticleMoment", 0 )
        self.PostEulerAngles              = getattr(DEM_parameters, "PostEulerAngles", 0 )
        self.PostLocalContactForce        = getattr(DEM_parameters, "PostLocalContactForce", 0 )
        self.PostFailureCriterionState    = getattr(DEM_parameters, "PostFailureCriterionState", 0 )
        self.PostContactFailureId         = getattr(DEM_parameters, "PostContactFailureId", 0 )
        self.PostContactTau               = getattr(DEM_parameters, "PostContactTau", 0 )
        self.PostContactSigma             = getattr(DEM_parameters, "PostContactSigma", 0 )
        self.PostMeanContactArea          = getattr(DEM_parameters, "PostMeanContactArea", 0 )
        self.PostElasticForces            = getattr(DEM_parameters, "PostElasticForces", 0 )
        self.PostContactForces            = getattr(DEM_parameters, "PostContactForces", 0 )
        self.PostRigidElementForces       = getattr(DEM_parameters, "PostRigidElementForces", 0 )
        self.PostPressure                 = getattr(DEM_parameters, "PostPressure", 0 )
        self.PostTangentialElasticForces  = getattr(DEM_parameters, "PostTangentialElasticForces", 0 )
        self.PostShearStress              = getattr(DEM_parameters, "PostShearStress", 0 )
        self.PostNodalArea                = getattr(DEM_parameters, "PostNodalArea", 0 )
        self.VelTrapGraphExportFreq       = getattr(DEM_parameters, "VelTrapGraphExportFreq", 0 )
        self.PostTemperature              = getattr(DEM_parameters, "PostTemperature", 0 )
        self.PostHeatFlux                 = getattr(DEM_parameters, "PostHeatFlux", 0 )
  
    
    def PushPrintVar(self,variable,name,print_list):
        if (Var_Translator(variable)):
            print_list.append(name)

    def AddGlobalVariables(self):
        # Global Variables
        self.PushPrintVar(self.PostDisplacement,             DISPLACEMENT,                self.global_variables)
        self.PushPrintVar(self.PostVelocity,                 VELOCITY,                    self.global_variables)
        self.PushPrintVar(self.PostTotalForces,              TOTAL_FORCES,                self.global_variables)
        
    def AddSpheresVariables(self):
        # Spheres Variables
        self.PushPrintVar(self.PostAppliedForces,    EXTERNAL_APPLIED_FORCE,       self.spheres_variables)
        self.PushPrintVar(self.PostDampForces,       DAMP_FORCES,                  self.spheres_variables)
        self.PushPrintVar(self.PostRadius,           RADIUS,                       self.spheres_variables)
        self.PushPrintVar(self.PostExportId,         EXPORT_ID,                    self.spheres_variables)
        self.PushPrintVar(self.PostExportSkinSphere, EXPORT_SKIN_SPHERE,           self.spheres_variables)
        self.PushPrintVar(self.PostRigidElementForces, RIGID_ELEMENT_FORCE,          self.spheres_variables)
        self.PushPrintVar(self.PostTemperature,      TEMPERATURE,                  self.spheres_variables)
        self.PushPrintVar(self.PostHeatFlux,         HEATFLUX,                     self.spheres_variables)
        #self.PushPrintVar( 1,                                               DUMMY_1, self.spheres_variables) # miquel mapping

        # Spheres Rotation
        if (Var_Translator(DEM_parameters.RotationOption)):  # xapuza
            self.PushPrintVar(self.PostAngularVelocity, ANGULAR_VELOCITY, self.spheres_variables)
            self.PushPrintVar(self.PostParticleMoment,  PARTICLE_MOMENT,  self.spheres_variables)
            self.PushPrintVar(self.PostEulerAngles,     EULER_ANGLES,     self.sphere_local_axis_variables)

        # Spheres Strain
        if ((DEM_parameters.ElementType == "SphericContPartDEMElement3D") or (DEM_parameters.ElementType == "CylinderContPartDEMElement3D")):
            self.PushPrintVar(DEM_parameters.StressStrainOption, REPRESENTATIVE_VOLUME, self.spheres_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_XX,         self.spheres_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_XY,         self.spheres_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_XZ,         self.spheres_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_YX,         self.spheres_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_YY,         self.spheres_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_YZ,         self.spheres_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_ZX,         self.spheres_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_ZY,         self.spheres_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_ZZ,         self.spheres_variables)
    
    
    def AddArlequinVariables(self):
         self.PushPrintVar( 1, DISTANCE, self.global_variables)
         self.PushPrintVar( 1, BORDER, self.global_variables)
         self.PushPrintVar( 1, SOLUTION, self.global_variables)
         self.PushPrintVar( 1, DUMMY_1, self.global_variables)
         self.PushPrintVar( 1, DUMMY_2, self.global_variables)
         self.PushPrintVar( 1, DUMMY_3, self.global_variables)
         self.PushPrintVar( 1, ALPHA_ARLEQUIN, self.global_variables)
         self.PushPrintVar( 1, LUMPED_PROJECTION_NODAL_MASS, self.global_variables)
     
    def AddFEMBoundaryVariables(self):
        
        self.PushPrintVar(self.PostElasticForces,            ELASTIC_FORCES, self.fem_boundary_variables)
        self.PushPrintVar(self.PostContactForces,            CONTACT_FORCES, self.fem_boundary_variables)
        self.PushPrintVar(self.PostPressure,                 PRESSURE, self.fem_boundary_variables)
        self.PushPrintVar(self.PostTangentialElasticForces,  TANGENTIAL_ELASTIC_FORCES, self.fem_boundary_variables)
        self.PushPrintVar(self.PostShearStress,              SHEAR_STRESS, self.fem_boundary_variables)
        self.PushPrintVar(self.PostNodalArea,                NODAL_AREA, self.fem_boundary_variables)
        if (Var_Translator(DEM_parameters.PostNonDimensionalVolumeWear)):
            self.PushPrintVar(1,                             NON_DIMENSIONAL_VOLUME_WEAR, self.fem_boundary_variables)
            self.PushPrintVar(1,                             IMPACT_WEAR,                 self.fem_boundary_variables)
        
    def AddMappingVariables(self):
        self.PushPrintVar( 1,                                               DUMMY_1, self.mapping_variables)

    def AddClusterVariables(self):
        pass
    
    def AddContactVariables(self):
        # Contact Elements Variables
        if ((DEM_parameters.ElementType == "SphericContPartDEMElement3D") or (DEM_parameters.ElementType == "CylinderContPartDEMElement3D")):
            if (Var_Translator(DEM_parameters.ContactMeshOption)):
                self.PushPrintVar(self.PostLocalContactForce,     LOCAL_CONTACT_FORCE,     self.contact_variables)
                self.PushPrintVar(self.PostFailureCriterionState, FAILURE_CRITERION_STATE, self.contact_variables)
                self.PushPrintVar(self.PostContactFailureId,      CONTACT_FAILURE,         self.contact_variables)
                self.PushPrintVar(self.PostContactTau,            CONTACT_TAU,             self.contact_variables)
                self.PushPrintVar(self.PostContactSigma,          CONTACT_SIGMA,           self.contact_variables)
                self.PushPrintVar(self.PostMeanContactArea,       MEAN_CONTACT_AREA,       self.contact_variables)

    def AddMpiVariables(self):
        pass

    def EnableMpiVariables(self):
        pass

    def Configure(self,problem_name,encoding,file_system,contact_mesh_option):
        self.problem_name = problem_name

        if (encoding == "Binary"):
            self.encoding = GiDPostMode.GiD_PostBinary
        else:
            self.encoding = GiDPostMode.GiD_PostAscii

        if (DEM_parameters.Multifile == "multiple_files"):
            self.filesystem = MultiFileFlag.MultipleFiles
        else:
            self.filesystem = MultiFileFlag.SingleFile

        self.deformed_mesh_flag  = WriteDeformedMeshFlag.WriteDeformed
        self.write_conditions    = WriteConditionsFlag.WriteConditions
        self.contact_mesh_option = contact_mesh_option

        self.gid_io = GidIO(self.problem_name, 
                            self.encoding, 
                            self.filesystem, 
                            self.deformed_mesh_flag,
                            self.write_conditions)

        self.post_utility = PostUtilities()

    def SetOutputName(self,name):
        self.gid_io.ChangeOutputName(name)

    def SetMultifileLists(self,multifile_list):
        for mfilelist in multifile_list:
            self.multifilelists.append(mfilelist)

        for mfilelist in self.multifilelists:
            mfilelist.file.write("Multiple\n")
            mfilelist.index = 1

    def PrintMultifileLists(self,time, post_path):
        for mfilelist in self.multifilelists:
            
            if mfilelist.index == mfilelist.step:
                
                if (self.encoding == GiDPostMode.GiD_PostBinary):
                    mfilelist.file.write(post_path+"/"+mfilelist.name+"_"+"%.12g"%time+".post.bin\n")
                else:
                    mfilelist.file.write(post_path+"/"+mfilelist.name+"_"+"%.12g"%time+".post.msh\n")
                    mfilelist.file.write(post_path+"/"+mfilelist.name+"_"+"%.12g"%time+".post.res\n")
                mfilelist.file.flush()
                mfilelist.index = 0
            mfilelist.index += 1

    def CloseMultifiles(self):
        for mfilelist in self.multifilelists:
            mfilelist.file.close()

    def InitializeMesh(self, mixed_model_part, spheres_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, mapping_model_part): #MIQUEL MAPPING
        if (self.filesystem == MultiFileFlag.SingleFile):

            self.post_utility.AddModelPartToModelPart(mixed_model_part, spheres_model_part)

            if (self.contact_mesh_option == "ON"):
                self.post_utility.AddModelPartToModelPart(mixed_model_part, contact_model_part)

            self.post_utility.AddModelPartToModelPart(mixed_model_part, rigid_face_model_part)
            self.post_utility.AddModelPartToModelPart(mixed_model_part, mapping_model_part) #MIQUEL MAPPING

            self.gid_io.InitializeMesh(0.0) 
            self.gid_io.WriteMesh(rigid_face_model_part.GetCommunicator().LocalMesh())
            self.gid_io.WriteMesh(mapping_model_part.GetCommunicator().LocalMesh())         #MIQUEL MAPPING
            self.gid_io.WriteSphereMesh(spheres_model_part.GetCommunicator().LocalMesh())

            if (self.contact_mesh_option == "ON"):
                self.gid_io.WriteMesh(contact_model_part.GetCommunicator().LocalMesh())

            self.gid_io.FinalizeMesh()
            self.gid_io.InitializeResults(0.0, mixed_model_part.GetCommunicator().LocalMesh())

    def InitializeResults(self, mixed_model_part, spheres_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, mapping_model_part, time): #MIQUEL MAPPING
        if (self.filesystem == MultiFileFlag.MultipleFiles):
            mixed_model_part.Elements.clear()
            mixed_model_part.Nodes.clear()

            self.post_utility.AddModelPartToModelPart(mixed_model_part, spheres_model_part)
            if (self.contact_mesh_option == "ON"):                
                self.post_utility.AddModelPartToModelPart(mixed_model_part, contact_model_part)
            self.post_utility.AddModelPartToModelPart(mixed_model_part, rigid_face_model_part)
            self.post_utility.AddModelPartToModelPart(mixed_model_part, mapping_model_part) #MIQUEL MAPPING

            self.gid_io.InitializeMesh(time) 
            self.gid_io.WriteSphereMesh(spheres_model_part.GetCommunicator().LocalMesh())
            if (self.contact_mesh_option == "ON"):
                self.gid_io.WriteMesh(contact_model_part.GetCommunicator().LocalMesh())
            self.gid_io.WriteMesh(rigid_face_model_part.GetCommunicator().LocalMesh())
            self.gid_io.WriteMesh(mapping_model_part.GetCommunicator().LocalMesh())#MIQUEL MAPPING
            self.gid_io.FinalizeMesh()            
            self.gid_io.InitializeResults(time, mixed_model_part.GetCommunicator().LocalMesh())

    def FinalizeMesh(self):
        if (self.filesystem == MultiFileFlag.SingleFile):
            self.gid_io.FinalizeResults()

    def FinalizeResults(self):
        if (self.filesystem == MultiFileFlag.MultipleFiles):
            self.gid_io.FinalizeResults() 

    def PrintingGlobalVariables(self, export_model_part, time):
        for variable in self.global_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)

    def PrintingSpheresVariables(self, export_model_part, time):
        for variable in self.spheres_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)
        for variable in self.sphere_local_axis_variables:
            self.gid_io.WriteLocalAxesOnNodes(variable, export_model_part.Nodes, time, 0)            

    def PrintingFEMBoundaryVariables(self, export_model_part, time):
        for variable in self.fem_boundary_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)
    
    def PrintingClusterVariables(self, export_model_part, time):
        for variable in self.cluster_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)

    def PrintingContactElementsVariables(self, export_model_part, time):
        if (self.contact_mesh_option == "ON"):
            for variable in self.contact_variables:
                self.gid_io.PrintOnGaussPoints(variable, export_model_part, time)

    def PrintingMappingVariables(self, export_model_part, time):
        for variable in self.mapping_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)
    
    #def PrintingArlequinVariables(self, export_model_part, time):
    #    self.gid_io.PrintOnGaussPoints(IN_ARLEQUIN, export_model_part, time)                

    def PrintResults(self, mixed_model_part, spheres_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, mapping_model_part, time):
        if (self.filesystem == MultiFileFlag.MultipleFiles):
            self.InitializeResults(mixed_model_part,
                                   spheres_model_part,
                                   rigid_face_model_part,
                                   cluster_model_part,
                                   contact_model_part,
                                   mapping_model_part,
                                   time)

        self.PrintingGlobalVariables(mixed_model_part, time)
        self.PrintingSpheresVariables(spheres_model_part, time)
        self.PrintingFEMBoundaryVariables(rigid_face_model_part, time)
        self.PrintingClusterVariables(cluster_model_part, time)
        self.PrintingContactElementsVariables(contact_model_part, time)
        self.PrintingMappingVariables(mapping_model_part, time)
        #self.PrintingArlequinVariables(rigid_face_model_part, time)

        if (self.filesystem == MultiFileFlag.MultipleFiles):
            self.FinalizeResults()


class ParallelUtils(object):

    def __init__(self):
        pass

    def Repart(self, spheres_model_part):
        pass

    def CalculateModelNewIds(self, spheres_model_part):
        pass

    def PerformInitialPartition(self, model_part, model_part_io_solid, input_file_name):
        return [model_part_io_solid, model_part, '']

    def GetSearchStrategy(self, solver, model_part):
        return solver.search_strategy
