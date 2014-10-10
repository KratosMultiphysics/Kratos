from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import DEM_explicit_solver_var as DEM_parameters
import DEM_material_test_script
import os

def Var_Translator(variable):

    if (variable == "OFF" or variable == "0" or variable == 0):
        variable = 0
    else:
        variable = 1

    return variable
    
def FindMaxNodeIdInModelPart(model_part):

    maxid = 0

    for node in model_part.Nodes:
        if (node.Id > maxid):
            maxid = node.Id

    return maxid


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
            mdpa.write(str(node.Id) + '   ' + str(node.X) + '  ' + str(node.Y) + '  ' + str(node.Z))
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

        self.balls_model_part = model_part
        self.UpdateData(domain_volume)

    def UpdateData(self, domain_volume):
        
        self.physics_calculator = SphericElementGlobalPhysicsCalculator(self.balls_model_part)
        self.number_of_balls    = self.balls_model_part.NumberOfElements(0)
        self.solid_volume       = self.physics_calculator.CalculateTotalVolume(self.balls_model_part)
        self.d_50               = self.physics_calculator.CalculateD50(self.balls_model_part)
        self.balls_per_area     = domain_volume / self.number_of_balls
        self.voids_volume       = domain_volume - self.solid_volume
        self.global_porosity    = self.voids_volume / domain_volume

    def PrintCurrentData(self):
        
        print("number_of_balls: ", self.number_of_balls)
        print("solid volume: ", self.solid_volume)
        print("voids volume: ", self.voids_volume)
        print("global porosity: ", self.global_porosity)
        print("D50: ", self.d_50)
        print("balls per area unit: ", self.balls_per_area)


class PostUtils(object):

    def __init__(self, DEM_parameters, balls_model_part):
        
        self.DEM_parameters = DEM_parameters
        self.balls_model_part = balls_model_part
        self.post_utilities = PostUtilities()

    def ComputeMeanVelocitiesinTrap(self, file_name, time_dem):

        if (self.DEM_parameters.VelocityTrapOption):
            average_velocity = Array3()
            low_point = Array3()

            low_point[0] = self.DEM_parameters.VelocityTrapMinX
            low_point[1] = self.DEM_parameters.VelocityTrapMinY
            low_point[2] = self.DEM_parameters.VelocityTrapMinZ
            high_point = Array3()
            high_point[0] = self.DEM_parameters.VelocityTrapMaxX
            high_point[1] = self.DEM_parameters.VelocityTrapMaxY
            high_point[2] = self.DEM_parameters.VelocityTrapMaxZ

            average_velocity = self.post_utilities.VelocityTrap(self.balls_model_part, low_point, high_point)
            f = open(file_name, 'a')
            tmp = str(time_dem) + "   " + str(average_velocity[0]) + "   " + str(average_velocity[1]) + "   " + str(average_velocity[2]) + "\n"
            f.write(tmp)
            f.flush()
            f.close()


class Procedures(object):

    def __init__(self, DEM_parameters):

        # GLOBAL VARIABLES OF THE SCRIPT
        # Defining list of skin particles (For a test tube of height 30 cm and diameter 15 cm)

        # Initialization of member variables
        # SIMULATION FLAGS
        self.rotation_OPTION        = Var_Translator(DEM_parameters.RotationOption)
        self.bounding_box_OPTION    = Var_Translator(DEM_parameters.BoundingBoxOption)
        self.contact_mesh_OPTION    = Var_Translator(DEM_parameters.ContactMeshOption)
        #self.solver = solver

        # SIMULATION SETTINGS
        self.bounding_box_enlargement_factor = DEM_parameters.BoundingBoxEnlargementFactor
       
        # MODEL
        self.domain_size = DEM_parameters.Dimension

    def AddCommonVariables(self, model_part, DEM_parameters):
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
        model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)

        if(DEM_parameters.PostGroupId):
            model_part.AddNodalSolutionStepVariable(GROUP_ID)

    def AddMpiVariables(self, model_part):
        pass

    def ModelData(self, balls_model_part, contact_model_part, solver):
        
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

        for node in balls_model_part.Nodes:

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

        Total_Particles = len(balls_model_part.Nodes)

        Total_Contacts = 0

        Coordination_Number = 0.0

        if(self.contact_mesh_OPTION):

            for bar in contact_model_part.Elements:

                Total_Contacts += 1.0

            Coordination_Number = 1.0 * (double(Total_Contacts) * 2.0) / double(Total_Particles)

        Model_Data.write("Total Number of Particles: " + str(Total_Particles) + '\n')
        Model_Data.write("Total Number of Contacts: " + str(Total_Contacts) + '\n')
        Model_Data.write("Coordination Number NC: " + str(Coordination_Number) + '\n')
        Model_Data.write('\n')

        # Model_Data.write("Volume Elements: " + str(total_volume) + '\n')

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

    def SetPredefinedSkin(self, balls_model_part):

        for element in balls_model_part.Elements:

            if (element.GetNode(0).GetSolutionStepValue(PREDEFINED_SKIN) > 0.0):  # PREDEFINED_SKIN is a double

                element.SetValue(SKIN_SPHERE, 1)
                
    def SetCustomSkin(self,balls_model_part):
    
        for element in balls_model_part.Elements:

            x = element.GetNode(0).X
            y = element.GetNode(0).Y
            #z = element.GetNode(0).Z
          
            if(x>21.1):
              element.SetValue(SKIN_SPHERE,1)
            if(x<1.25):
              element.SetValue(SKIN_SPHERE,1)
            if(y>1.9):
              element.SetValue(SKIN_SPHERE,1)
            if(y<0.1):
              element.SetValue(SKIN_SPHERE,1)

    def CreateDirectories(self, main_path,problem_name):

        root             = main_path + '/' + problem_name

        post_path        = root + '_Post_Files'
        list_path        = root + '_Post_Lists'
        data_and_results = root + '_Results_and_Data'
        graphs_path      = root + '_Graphs'
        MPI_results      = root + '_MPI_results'

        for directory in [post_path, list_path, data_and_results, graphs_path, MPI_results]:
            if not os.path.isdir(directory):
                os.makedirs(str(directory))

        return [post_path,list_path,data_and_results,graphs_path,MPI_results]

    def PreProcessModel(self, DEM_parameters):
        pass

    def KRATOSprint(self,message):
        print(message)
        sys.stdout.flush()

# #~CHARLIE~# Aixo no ho entenc 
class DEMFEMProcedures(object):

    def __init__(self, DEM_parameters, graphs_path, balls_model_part, RigidFace_model_part):

        # GLOBAL VARIABLES OF THE SCRIPT
        self.TestType = DEM_parameters.TestType

        # Initialization of member variables
        # SIMULATION FLAGS
        self.rotation_OPTION     = Var_Translator(DEM_parameters.RotationOption)
        self.bounding_box_OPTION = Var_Translator(DEM_parameters.BoundingBoxOption)
        self.contact_mesh_OPTION = Var_Translator(DEM_parameters.ContactMeshOption)

        self.graphs_path = graphs_path
        self.balls_model_part = balls_model_part
        self.RigidFace_model_part = RigidFace_model_part
        #self.solver = solver

        self.top_mesh_nodes = []

        self.graph_counter = 0;
        self.graph_frequency        = int(DEM_parameters.GraphExportFreq/balls_model_part.ProcessInfo.GetValue(DELTA_TIME))
        os.chdir(self.graphs_path)
        self.graph_forces = open(DEM_parameters.problem_name +"_force_graph.grf", 'w')
        self.total_force_top = 0.0
        self.total_force_x = 0.0
        self.total_force_y = 0.0
        self.total_force_z = 0.0

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

    def MeasureForces(self):
        
        if self.TestType == "None":
            if self.RigidFace_model_part.NumberOfMeshes() > 1:
                for mesh_number in range(1, self.RigidFace_model_part.NumberOfMeshes()):
                    if(self.RigidFace_model_part.GetMesh(mesh_number)[FORCE_INTEGRATION_GROUP]):
                        self.top_mesh_nodes = self.RigidFace_model_part.GetMesh(mesh_number).Nodes
                        self.total_force_x = 0.0
                        self.total_force_y = 0.0
                        self.total_force_z = 0.0

                        for node in self.top_mesh_nodes:
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
                self.graph_forces.write(str(time)+"    "+str(self.total_force_x)+"    "+str(self.total_force_y)+"    "+str(self.total_force_z)+'\n')
                self.graph_forces.flush()                

            self.graph_counter += 1

    def FinalizeGraphs(self):

        if DEM_parameters.TestType == "None":
            os.chdir(self.graphs_path)
            self.graph_forces.close()


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

    def Initialize(self, DEM_parameters, procedures, solver, graphs_path, post_path, balls_model_part, rigid_face_model_part):
        self.type = DEM_parameters.TestType

        if (self.type != "None"):
            self.script = DEM_material_test_script.MaterialTest(DEM_parameters, procedures, solver, graphs_path, post_path, balls_model_part, rigid_face_model_part)
            self.script.Initialize()
 
    def PrepareDataForGraph(self):
        if (self.type != "None"):
            self.script.PrepareDataForGraph()

    def MeasureForcesAndPressure(self):
        if (self.type != "None"):
            self.script.MeasureForcesAndPressure()
            
    def PrintGraph(self, step):
        if (self.type != "None"):
            self.script.PrintGraph(step)

    def FinalizeGraphs(self):
        if (self.type != "None"):
            self.script.FinalizeGraphs()


class MultifileList(object):

    def __init__(self,name,step):
        self.index = 0
        self.step = step
        self.name = name
        self.file = open(self.name+"_"+str(step)+".post.lst","w")


class DEMIo(object):

    def __init__(self):
        # Printing variables
        self.global_variables          = []
        self.ball_variables            = []
        self.ball_local_axis_variables = []
        self.contact_variables         = []
        self.multifilelists            = []

    def PushPrintVar(self,variable,name,print_list):
        if (Var_Translator(variable)):
            print_list.append(name)

    def AddGlobalVariables(self):
        # Global Variables
        self.PushPrintVar(DEM_parameters.PostDisplacement, DISPLACEMENT,    self.global_variables)
        self.PushPrintVar(DEM_parameters.PostVelocity,     VELOCITY,        self.global_variables)
        self.PushPrintVar(DEM_parameters.PostTotalForces,  TOTAL_FORCES,    self.global_variables)

    def AddBallVariables(self):
        # Balls Variables
        self.PushPrintVar(DEM_parameters.PostAppliedForces,    EXTERNAL_APPLIED_FORCE, self.ball_variables)
        self.PushPrintVar(DEM_parameters.PostDampForces,       DAMP_FORCES,            self.ball_variables)
        self.PushPrintVar(DEM_parameters.PostRadius,           RADIUS,                 self.ball_variables)
        self.PushPrintVar(DEM_parameters.PostExportId,         EXPORT_ID,              self.ball_variables)
        self.PushPrintVar(DEM_parameters.PostExportSkinSphere, EXPORT_SKIN_SPHERE,     self.ball_variables)

        # Balls Rotation
        if (Var_Translator(DEM_parameters.RotationOption)):  # xapuza
            self.PushPrintVar(DEM_parameters.PostAngularVelocity, ANGULAR_VELOCITY, self.ball_variables)
            self.PushPrintVar(DEM_parameters.PostParticleMoment,  PARTICLE_MOMENT,  self.ball_variables)
            self.PushPrintVar(DEM_parameters.PostEulerAngles,     EULER_ANGLES,     self.ball_local_axis_variables)

        # Balls Strain
        if ((DEM_parameters.ElementType == "SphericContPartDEMElement3D") or(DEM_parameters.ElementType == "CylinderContPartDEMElement3D")):
            self.PushPrintVar(DEM_parameters.StressStrainOption, REPRESENTATIVE_VOLUME, self.ball_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_XX,         self.ball_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_XY,         self.ball_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_XZ,         self.ball_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_YX,         self.ball_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_YY,         self.ball_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_YZ,         self.ball_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_ZX,         self.ball_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_ZY,         self.ball_variables)
            self.PushPrintVar(DEM_parameters.StressStrainOption, DEM_STRESS_ZZ,         self.ball_variables)

    def AddContactVariables(self):
        # Contact Elements Variables
        if ((DEM_parameters.ElementType == "SphericContPartDEMElement3D") or(DEM_parameters.ElementType == "CylinderContPartDEMElement3D")):
            if (Var_Translator(DEM_parameters.ContactMeshOption)):
                self.PushPrintVar(DEM_parameters.PostLocalContactForce,     LOCAL_CONTACT_FORCE,     self.contact_variables)
                self.PushPrintVar(DEM_parameters.PostFailureCriterionState, FAILURE_CRITERION_STATE, self.contact_variables)
                self.PushPrintVar(DEM_parameters.PostContactFailureId,      CONTACT_FAILURE,         self.contact_variables)
                self.PushPrintVar(DEM_parameters.PostContactTau,            CONTACT_TAU,             self.contact_variables)
                self.PushPrintVar(DEM_parameters.PostContactSigma,          CONTACT_SIGMA,           self.contact_variables)
                self.PushPrintVar(DEM_parameters.PostMeanContactArea,       MEAN_CONTACT_AREA,       self.contact_variables)

    def AddMpiVariables(self):
        pass

    def EnableMpiVariables(self):
        pass

    def Configure(self,problem_name,encoding,file_system,contact_mesh_option):
        self.problem_name = problem_name;

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

    def PrintMultifileLists(self,time):
        for mfilelist in self.multifilelists:
            if mfilelist.index == mfilelist.step:
                mfilelist.file.write(mfilelist.name+"_"+str(time)+".post.bin\n")
                mfilelist.index = 0
            mfilelist.index += 1

    def CloseMultifiles(self):
        for mfilelist in self.multifilelists:
            mfilelist.file.close()

    def InitializeMesh(self,mixed_model_part,balls_model_part,rigid_face_model_part,contact_model_part):
        if (self.filesystem == MultiFileFlag.SingleFile):

            self.post_utility.AddModelPartToModelPart(mixed_model_part, balls_model_part)

            if (self.contact_mesh_option == "ON"):
                self.post_utility.AddModelPartToModelPart(mixed_model_part, contact_model_part)

            self.post_utility.AddModelPartToModelPart(mixed_model_part, rigid_face_model_part)
            self.gid_io.InitializeMesh(0.0) 
            self.gid_io.WriteMesh(rigid_face_model_part.GetMesh())
            self.gid_io.WriteSphereMesh(balls_model_part.GetMesh())

            if (self.contact_mesh_option == "ON"):
                self.gid_io.WriteMesh(contact_model_part.GetMesh())

            self.gid_io.FinalizeMesh()
            self.gid_io.InitializeResults(0.0, mixed_model_part.GetMesh())

    def InitializeResults(self,mixed_model_part,balls_model_part,rigid_face_model_part,contact_model_part,time):
        if (self.filesystem == MultiFileFlag.MultipleFiles):
            mixed_model_part.Elements.clear()
            mixed_model_part.Nodes.clear()

            self.post_utility.AddModelPartToModelPart(mixed_model_part, balls_model_part)
            if (self.contact_mesh_option == "ON"):
                self.post_utility.AddModelPartToModelPart(mixed_model_part, contact_model_part)
            self.post_utility.AddModelPartToModelPart(mixed_model_part, rigid_face_model_part)

            self.gid_io.InitializeMesh(time) 
            self.gid_io.WriteSphereMesh(balls_model_part.GetMesh())
            if (self.contact_mesh_option == "ON"):
                self.gid_io.WriteMesh(contact_model_part.GetMesh())
            self.gid_io.WriteMesh(rigid_face_model_part.GetMesh())
            self.gid_io.FinalizeMesh()
            
            self.gid_io.InitializeResults(time, mixed_model_part.GetMesh())

    def FinalizeMesh(self):
        if (self.filesystem == MultiFileFlag.SingleFile):
            self.gid_io.FinalizeResults()

    def FinalizeResults(self):
        if (self.filesystem == MultiFileFlag.MultipleFiles):
            self.gid_io.FinalizeResults() 

    def PrintingGlobalVariables(self, export_model_part, time):
        for variable in self.global_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)

    def PrintingBallsVariables(self, export_model_part, time):
        for variable in self.ball_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)
        for variable in self.ball_local_axis_variables:
            self.gid_io.WriteLocalAxesOnNodes(variable, export_model_part.Nodes, time, 0)            

    def PrintingContactElementsVariables(self, export_model_part, time):
        if (self.contact_mesh_option == "ON"):
            for variable in self.contact_variables:
                self.gid_io.PrintOnGaussPoints(variable, export_model_part, time)
                #self.gid_io.PrintOnGaussPoints(CONTACT_ORIENTATION, export_model_part, time)

    def PrintResults(self,mixed_model_part,balls_model_part,rigid_face_model_part,contact_model_part,time):
        if (self.filesystem == MultiFileFlag.MultipleFiles):
            self.InitializeResults(mixed_model_part,
                                   balls_model_part,
                                   rigid_face_model_part,
                                   contact_model_part,
                                   time)

        self.PrintingGlobalVariables(mixed_model_part, time)
        self.PrintingBallsVariables(balls_model_part, time)
        self.PrintingContactElementsVariables(rigid_face_model_part, time)

        if (self.filesystem == MultiFileFlag.MultipleFiles):
            self.FinalizeResults()


class ParallelUtils(object):

    def __init__(self):
        pass

    def Repart(self, balls_model_part):
        pass

    def CalculateModelNewIds(self, balls_model_part):
        pass

    def PerformInitialPartition(self, model_part, model_part_io_solid, input_file_name):
        return [model_part_io_solid, model_part, '']

    def GetSearchStrategy(self, solver, model_part):
        return solver.search_strategy
