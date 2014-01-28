from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import os
# from numpy import *

# from KratosMultiphysics.mpi import * #CARLOS


def Var_Translator(variable):

    if (variable == "OFF" or variable == "0" or variable == 0):
        variable = 0
    else:
        variable = 1

    return variable


class MdpaCreator:

    def __init__(self, path, param):
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


class GranulometryUtils:

    def __init__(self, domain_volume, model_part):

        self.balls_model_part = model_part
        self.physics_calculator = SphericElementGlobalPhysicsCalculator(self.balls_model_part)
        self.UpdateData(domain_volume)

    def UpdateData(self, domain_volume):

        self.number_of_balls = self.balls_model_part.NumberOfElements(0)
        self.solid_volume = self.physics_calculator.CalculateTotalVolume(self.balls_model_part)
        self.d_50 = self.physics_calculator.CalculateD50(self.balls_model_part)
        self.balls_per_area = domain_volume / self.number_of_balls
        self.voids_volume = domain_volume - self.solid_volume
        self.global_porosity = self.voids_volume / domain_volume

    def PrintCurrentData(self):

        print("solid volume: ", self.solid_volume)
        print("voids volume: ", self.voids_volume)
        print("D50: ", self.d_50)
        print("global porosity: ", self.global_porosity)
        print("number_of_balls: ", self.number_of_balls)
        print("balls per area unit: ", self.balls_per_area)


class PostUtils:

    def __init__(self, param, balls_model_part):
        self.param = param
        self.balls_model_part = balls_model_part
        self.post_utilities = PostUtilities()

    def ComputeMeanVelocitiesinTrap(self, file_name, time_dem):

        if (self.param.VelocityTrapOption):
            average_velocity = Array3()
            low_point = Array3()

            low_point[0] = self.param.VelocityTrapMinX
            low_point[1] = self.param.VelocityTrapMinY
            low_point[2] = self.param.VelocityTrapMinZ
            high_point = Array3()
            high_point[0] = self.param.VelocityTrapMaxX
            high_point[1] = self.param.VelocityTrapMaxY
            high_point[2] = self.param.VelocityTrapMaxZ

            average_velocity = self.post_utilities.VelocityTrap(self.balls_model_part, low_point, high_point)
            f = open(file_name, 'a')
            tmp = str(time_dem) + "   " + str(average_velocity[0]) + "   " + str(average_velocity[1]) + "   " + str(average_velocity[2]) + "\n"
            f.write(tmp)
            f.flush()
            f.close()


class Procedures:

    def __init__(self, param):

        # GLOBAL VARIABLES OF THE SCRIPT
        # Defining list of skin particles (For a test tube of height 30 cm and diameter 15 cm)

        self.sup_layer_fm = list()
        self.inf_layer_fm = list()
        self.sup_plate_fm = list()
        self.inf_plate_fm = list()
        self.special_selection = list()
        self.others = list()
        self.SKIN = list()
        self.LAT = list()
        self.BOT = list()
        self.TOP = list()
        self.XLAT = list()  # only lat, not the corner ones
        self.XTOP = list()  # only top, not corner ones...
        self.XBOT = list()
        self.XTOPCORNER = list()
        self.XBOTCORNER = list()

        # Initialization of member variables
        # SIMULATION FLAGS
        self.rotation_OPTION = Var_Translator(param.RotationOption)
        self.bounding_box_OPTION = Var_Translator(param.BoundingBoxOption)
        self.continuum_OPTION = Var_Translator(param.ContinuumOption)
        self.contact_mesh_OPTION = Var_Translator(Var_Translator(param.ContactMeshOption) & Var_Translator(param.ContinuumOption))

        # SIMULATION SETTINGS

        self.bounding_box_enlargement_factor = param.BoundingBoxEnlargementFactor
       # MODEL
        self.domain_size = param.Dimension

        # PRINTING VARIABLES

        self.print_radius = Var_Translator(param.PostRadius)
        self.print_velocity = Var_Translator(param.PostVelocity)
        self.print_angular_velocity = Var_Translator(param.PostAngularVelocity)
        self.print_displacement = Var_Translator(param.PostDisplacement)
        self.print_radial_displacement = Var_Translator(param.PostRadialDisplacement)
        self.print_total_forces = Var_Translator(param.PostTotalForces)
        self.print_damp_forces = Var_Translator(param.PostDampForces)
        self.print_applied_forces = Var_Translator(param.PostAppliedForces)
        self.print_particle_moment = Var_Translator(param.PostParticleMoment)
        self.print_particle_cohesion = Var_Translator(param.PostParticleCohesion)
        self.print_particle_tension = Var_Translator(param.PostParticleTension)
        self.print_euler_angles = Var_Translator(param.PostEulerAngles)
        self.print_group_id = Var_Translator(param.PostGroupId)
        self.print_export_id = Var_Translator(param.PostExportId)

        self.total_volume = param.TotalElementsVolume

        if (Var_Translator(param.ContinuumOption)):
            self.print_export_skin_sphere = Var_Translator(param.PostExportSkinSphere)
            self.predefined_skin_option = Var_Translator(param.PredefinedSkinOption)
            if (self.contact_mesh_OPTION):
                self.print_local_contact_force = Var_Translator(param.PostLocalContactForce)
                self.print_failure_criterion_state = Var_Translator(param.PostFailureCriterionState)
                self.print_unidimensional_damage = Var_Translator(param.PostUnidimensionalDamage)
                self.print_contact_failure = Var_Translator(param.PostContactFailureId)
                self.print_contact_tau = Var_Translator(param.PostContactTau)
                self.print_contact_sigma = Var_Translator(param.PostContactSigma)
                self.print_mean_contact_area = Var_Translator(param.PostMeanContactArea)

    def AddMpiVariables(self, model_part):

        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        model_part.AddNodalSolutionStepVariable(INTERNAL_ENERGY)
        model_part.AddNodalSolutionStepVariable(OSS_SWITCH)

    def PerformInitialPartition(self, model_part, model_part_io_solid, input_file_name):

        self.domain_size = 3

        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "before performing the division")
        number_of_partitions = mpi.size  # we set it equal to the number of processors

        if mpi.rank == 0:
            print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "start partition process")
            partitioner = MortonDivideInputToPartitionsProcess(model_part_io_solid, number_of_partitions, data.domain_size)
            partitioner.Execute()

        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "division performed")
        mpi.world.barrier()

        MPICommSetup = SetMPICommunicatorProcess(model_part)
        MPICommSetup.Execute()

        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "Comunicator Set")

        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "Reading: " + input_file_name + "_" + str(mpi.rank))

        my_input_filename = input_file_name + "_" + str(mpi.rank)
        model_part_io_solid = ModelPartIO(my_input_filename)

        return model_part_io_solid

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

        for node in balls_model_part.Nodes:

            sum_radi += node.GetSolutionStepValue(RADIUS)
            partial_sum_squared = node.GetSolutionStepValue(RADIUS) ** 2.0
            total_sum_squared += partial_sum_squared
            volume += 4 * 3.141592 / 3 * node.GetSolutionStepValue(RADIUS) ** 3.0
            area += 3.141592 * partial_sum_squared
            i += 1.0

        mean = sum_radi / i
        var = total_sum_squared / i - mean ** 2.0
        std_dev = 0.0

        if(abs(var) > 1e-9):
            std_dev = var ** 0.5

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

    def CylinderSkinDetermination(self, model_part, solver, param):

        # SKIN DETERMINATION
        total_cross_section = 0.0

        # Cylinder dimensions

        h = param.SpecimenHeight
        d = param.SpecimenWidth

        eps = 2.0

        surface = 2 * (3.141592 * d * d * 0.25) + (3.141592 * d * h)

        xlat_area = 0.0
        xbot_area = 0.0
        xtop_area = 0.0
        xbotcorner_area = 0.0
        xtopcorner_area = 0.0

        for element in model_part.Elements:

            element.SetValue(SKIN_SPHERE, 0)

            node = element.GetNode(0)
            r = node.GetSolutionStepValue(RADIUS)
            x = node.X
            y = node.Y
            z = node.Z
            node_group = node.GetSolutionStepValue(GROUP_ID)
            cross_section = 3.141592 * r * r

            if ((x * x + z * z) >= ((d / 2 - eps * r) * (d / 2 - eps * r))):

                element.SetValue(SKIN_SPHERE, 1)
                self.LAT.append(node)

                if ((y > eps * r) and (y < (h - eps * r))):

                    self.SKIN.append(element)
                    self.XLAT.append(node)

                    xlat_area = xlat_area + cross_section

            if ((y <= eps * r) or (y >= (h - eps * r))):

                element.SetValue(SKIN_SPHERE, 1)
                self.SKIN.append(element)

                if (y <= eps * r):

                    self.BOT.append(node)

                elif (y >= (h - eps * r)):

                    self.TOP.append(node)

                if ((x * x + z * z) >= ((d / 2 - eps * r) * (d / 2 - eps * r))):

                    if (y > h / 2):

                        self.XTOPCORNER.append(node)
                        xtopcorner_area = xtopcorner_area + cross_section

                    else:

                        self.XBOTCORNER.append(node)
                        xbotcorner_area = xbotcorner_area + cross_section
                else:

                    if (y <= eps * r):

                        self.XBOT.append(node)
                        xbot_area = xbot_area + cross_section

                    elif (y >= (h - eps * r)):

                        self.XTOP.append(node)
                        xtop_area = xtop_area + cross_section

        print("End ", h, "x", d, "Cylinder Skin Determination", "\n")

        return (xtop_area, xbot_area, xlat_area, xtopcorner_area, xbotcorner_area)

    def BtsSkinDetermination(self, model_part, solver, param):

        # SKIN DETERMINATION

        # Cylinder dimensions

        h = 0.086
        d = 0.15
        eps = 2.0

        for element in model_part.Elements:

            element.SetValue(SKIN_SPHERE, 0)

            node = element.GetNode(0)
            r = node.GetSolutionStepValue(RADIUS)
            x = node.X
            y = node.Y
            z = node.Z

            if ((x * x + y * y) >= ((d / 2 - eps * r) * (d / 2 - eps * r))):

                element.SetValue(SKIN_SPHERE, 1)

            if ((z <= eps * r) or (z >= (h - eps * r))):

                element.SetValue(SKIN_SPHERE, 1)

        print("End 30x15 Bts Skin Determination", "\n")

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

    def PrintingGlobalVariables(self, gid_io, export_model_part, time):

        if (self.print_displacement):
            gid_io.WriteNodalResults(DISPLACEMENT, export_model_part.Nodes, time, 0)
        if (self.print_velocity):
            gid_io.WriteNodalResults(VELOCITY, export_model_part.Nodes, time, 0)
        if (self.print_total_forces):
            gid_io.WriteNodalResults(TOTAL_FORCES, export_model_part.Nodes, time, 0)
        if (self.print_group_id):
            gid_io.WriteNodalResults(EXPORT_GROUP_ID, export_model_part.Nodes, time, 0)
        gid_io.Flush()
        sys.stdout.flush()

    def PrintingBallsVariables(self, gid_io, export_model_part, time):

        if (self.print_radial_displacement):
            gid_io.WriteNodalResults(RADIAL_DISPLACEMENT, export_model_part.Nodes, time, 0)
        if (self.print_applied_forces):
            gid_io.WriteNodalResults(EXTERNAL_APPLIED_FORCE, export_model_part.Nodes, time, 0)
        if (self.print_damp_forces):
            gid_io.WriteNodalResults(DAMP_FORCES, export_model_part.Nodes, time, 0)
        if (self.print_radius):
            gid_io.WriteNodalResults(RADIUS, export_model_part.Nodes, time, 0)
        if (self.print_particle_cohesion):
            gid_io.WriteNodalResults(PARTICLE_COHESION, export_model_part.Nodes, time, 0)
        if (self.print_particle_tension):
            gid_io.WriteNodalResults(PARTICLE_TENSION, export_model_part.Nodes, time, 0)
        if (self.print_export_id):
            gid_io.WriteNodalResults(EXPORT_ID, export_model_part.Nodes, time, 0)

        if (self.continuum_OPTION):

            #if (self.print_export_particle_failure_id):
                #gid_io.WriteNodalResults(EXPORT_PARTICLE_FAILURE_ID, export_model_part.Nodes, time, 0)
            if (self.print_export_skin_sphere):
                gid_io.WriteNodalResults(EXPORT_SKIN_SPHERE, export_model_part.Nodes, time, 0)

        # Aixo sempre per que si no hi ha manera de debugar
        # gid_io.WriteNodalResults(PARTITION_INDEX, export_model_part.Nodes, time, 0)
        # gid_io.WriteNodalResults(INTERNAL_ENERGY, export_model_part.Nodes, time, 0)

        if (self.rotation_OPTION):  # xapuza
            if (self.print_angular_velocity):
                gid_io.WriteNodalResults(ANGULAR_VELOCITY, export_model_part.Nodes, time, 0)
            if (self.print_particle_moment):
                gid_io.WriteNodalResults(PARTICLE_MOMENT, export_model_part.Nodes, time, 0)
            if (self.print_euler_angles):
                gid_io.WriteLocalAxesOnNodes(EULER_ANGLES, export_model_part.Nodes, time, 0)

        gid_io.Flush()
        sys.stdout.flush()

    def PrintingContactElementsVariables(self, gid_io, export_model_part, time):
        if (self.contact_mesh_OPTION):  # xapuza
            if (self.print_local_contact_force):
                gid_io.PrintOnGaussPoints(LOCAL_CONTACT_FORCE, export_model_part, time)
            if (self.print_mean_contact_area):
                gid_io.PrintOnGaussPoints(MEAN_CONTACT_AREA, export_model_part, time)
            if (self.print_contact_failure):
                gid_io.PrintOnGaussPoints(CONTACT_FAILURE, export_model_part, time)
            if (self.print_failure_criterion_state):
                gid_io.PrintOnGaussPoints(FAILURE_CRITERION_STATE, export_model_part, time)
            if (self.print_unidimensional_damage):
                gid_io.PrintOnGaussPoints(UNIDIMENSIONAL_DAMAGE, export_model_part, time)
            if (self.print_contact_tau):
                gid_io.PrintOnGaussPoints(CONTACT_TAU, export_model_part, time)
            if (self.print_contact_sigma):
                gid_io.PrintOnGaussPoints(CONTACT_SIGMA, export_model_part, time)
           # gid_io.PrintOnGaussPoints(NON_ELASTIC_STAGE,export_model_part,time)
        gid_io.Flush()
        sys.stdout.flush()

# DEM CONTINUUM # # #

    def SetPredefinedSkin(self, balls_model_part):

        for element in balls_model_part.Elements:

            if (element.GetNode(0).GetSolutionStepValue(PREDEFINED_SKIN) > 0.0):  # PREDEFINED_SKIN is a double

                element.SetValue(SKIN_SPHERE, 1)

    def ListDefinition(self, model_part, solver):

        for node in model_part.Nodes:
            if (node.GetSolutionStepValue(GROUP_ID) == 1):  # reserved for speciment particles with imposed displacement and strain-stress measurement (superior). Doesn't recive pressure
                self.sup_layer_fm.append(node)
            elif (node.GetSolutionStepValue(GROUP_ID) == 2):  # reserved for speciment particles with imposed displacement and strain-stress measurement (superior). Doesn't recive pressure
                self.inf_layer_fm.append(node)
            elif (node.GetSolutionStepValue(GROUP_ID) == 3):  # reserved for auxiliar strain-stress measurement plate (superior)
                self.sup_plate_fm.append(node)
            elif (node.GetSolutionStepValue(GROUP_ID) == 4):  # reserved for auxiliar strain-stress measurement plate (inferior)
                self.inf_plate_fm.append(node)
            elif (node.GetSolutionStepValue(GROUP_ID) == 5):
                self.special_selection.append(node)
            else:
                self.others.append(node)

        return (self.sup_layer_fm, self.inf_layer_fm, self.sup_plate_fm, self.inf_plate_fm)
