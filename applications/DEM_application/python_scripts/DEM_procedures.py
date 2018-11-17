# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import math
import os
import shutil
import sys
import weakref
from glob import glob

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import DEM_material_test_script

def Flush(a):
    a.flush()


def KratosPrint(*args):
    Logger.Print(*args, label="DEM")
    Logger.Flush()


def Var_Translator(variable):

    if variable == "OFF" or variable == "0" or variable == 0:
        variable = 0
    else:
        variable = 1

    return variable


def GetBoolParameterIfItExists(set_of_parameters, parameter_key):
    if parameter_key in set_of_parameters.keys():
        return set_of_parameters[parameter_key].GetBool()
    return False


class MdpaCreator(object):

    def __init__(self, path, DEM_parameters):

        self.DEM_parameters = DEM_parameters
        self.current_path = path

        # Creating necessary directories
        self.post_mdpas = os.path.join(str(self.current_path), str(self.DEM_parameters["problem_name"].GetString()) + '_post_mdpas')
        if not os.path.isdir(self.post_mdpas):
            os.makedirs(str(self.post_mdpas))

    def WriteMdpa(self, model_part):
        time = model_part.ProcessInfo.GetValue(TIME)
        absolute_path_to_file = os.path.join(self.post_mdpas, str(self.DEM_parameters["problem_name"].GetString()) + '_post_' + str(time) + '.mdpa')
        mdpa = open(absolute_path_to_file, 'w')
        mdpa.write('Begin ModelPartData' + '\n')
        mdpa.write('//  VARIABLE_NAME value')
        mdpa.write('End ModelPartData' + '\n' + '\n' + '\n' + '\n')
        mdpa.write('Begin Nodes' + '\n')

        for node in model_part.Nodes:
            mdpa.write(str(node.Id) + ' ' + str(node.X) + ' ' + str(node.Y) + ' ' + str(node.Z) + '\n')
        mdpa.write('End Nodes' + '\n' + '\n')

        mdpa.write('Begin Elements SphericParticle3D' + '\n')
        for element in model_part.Elements:
            mdpa.write(str(element.Id) + ' ' + '1' + ' ' + str(element.GetNode(0).Id) + '\n')
        mdpa.write('End Elements' + '\n' + '\n')

        self.WriteVariableData(RADIUS, mdpa, model_part)
        #self.WriteVariableData(VELOCITY_X, mdpa, model_part)
        #self.WriteVariableData(VELOCITY_Y, mdpa, model_part)
        #self.WriteVariableData(VELOCITY_Z, mdpa, model_part)

    @classmethod
    def WriteVariableData(self, variable_name, mdpa, model_part):

        mdpa.write('Begin NodalData ' + str(variable_name) + '\n')
        for node in model_part.Nodes:
            mdpa.write(str(node.Id) + ' ' + str(0) + ' ' + str(node.GetSolutionStepValue(variable_name)) + '\n')
        mdpa.write('End NodalData' + '\n' + '\n')


class SetOfModelParts(object):
    def __init__(self, model_parts_list):
        self.MaxNodeId = 0
        self.MaxElemId = 0
        self.MaxCondId = 0

        self.model_parts = dict()
        self.mp_list = []
        for mp in model_parts_list:
            self.model_parts[mp.Name] = mp
            self.mp_list.append(mp)

        self.spheres_model_part = self.Get("SpheresPart")
        self.rigid_face_model_part = self.Get("RigidFacePart")
        self.cluster_model_part = self.Get("ClusterPart")
        self.DEM_inlet_model_part = self.Get("DEMInletPart")
        self.mapping_model_part = self.Get("MappingPart")
        self.contact_model_part = self.Get("ContactPart")

    def ComputeMaxIds(self):

        for mp in self.mp_list:
            self.GetMaxIds(mp)

    def GetMaxIds(self, model_part):

        for node in model_part.Nodes:
            self.MaxNodeId = max(self.MaxNodeId, node.Id)

        for elem in model_part.Elements:
            self.MaxElemId = max(self.MaxElemId, elem.Id)

        for cond in model_part.Conditions:
            self.MaxCondId = max(self.MaxCondId, cond.Id)

    def Get(self, name):
        return self.model_parts[name]

    def Add(self, model_part, name=None):
        if name != None:
            self.model_parts[name] = model_part
        else:
            self.model_parts[model_part.Name] = model_part

        self.mp_list.append(model_part)


class GranulometryUtils(object):

    def __init__(self, domain_volume, model_part):

        if domain_volume <= 0.0:
            raise ValueError(
                "Error: The input domain volume must be strictly positive!")

        self.spheres_model_part = model_part
        self.UpdateData(domain_volume)

    def UpdateData(self, domain_volume):

        self.physics_calculator = SphericElementGlobalPhysicsCalculator(self.spheres_model_part)
        self.number_of_spheres = self.spheres_model_part.NumberOfElements(0)
        self.solid_volume = self.physics_calculator.CalculateTotalVolume(self.spheres_model_part)
        self.d_50 = self.physics_calculator.CalculateD50(self.spheres_model_part)

        if self.number_of_spheres == 0:
            self.spheres_per_area = 0.0
        else:
            self.spheres_per_area = domain_volume / self.number_of_spheres

        self.voids_volume = domain_volume - self.solid_volume
        self.global_porosity = self.voids_volume / domain_volume

    def PrintCurrentData(self):

        Logger.Print("number_of_spheres: ", self.number_of_spheres, label="")
        Logger.Print("solid volume: ", self.solid_volume, label="")
        Logger.Print("voids volume: ", self.voids_volume, label="")
        Logger.Print("global porosity: ", self.global_porosity, label="")
        Logger.Print("D50: ", self.d_50, label="")
        Logger.Print("spheres per area unit: ", self.spheres_per_area, label="")


class PostUtils(object):

    def __init__(self, DEM_parameters, spheres_model_part):

        self.DEM_parameters = DEM_parameters
        self.spheres_model_part = spheres_model_part
        self.post_utilities = PostUtilities()

        self.vel_trap_graph_counter = 0
        # TODO: change the name of VelTrapGraphExportFreq to VelTrapGraphExportTimeInterval
        self.vel_trap_graph_frequency = int(self.DEM_parameters["VelTrapGraphExportFreq"].GetDouble() / spheres_model_part.ProcessInfo.GetValue(DELTA_TIME))
        if self.vel_trap_graph_frequency < 1:
            self.vel_trap_graph_frequency = 1 # that means it is not possible to print results with a higher frequency than the computations delta time


        self.previous_vector_of_inner_nodes = []
        self.previous_time = 0.0

    @classmethod
    def Flush(self, a):
        a.flush()

    def ComputeMeanVelocitiesInTrap(self, file_name, time_dem, graphs_path):

        if self.DEM_parameters["VelocityTrapOption"].GetBool():
            compute_flow = "ComputeFlow" in self.DEM_parameters.keys() and self.DEM_parameters["ComputeFlow"].GetBool()

            self.vel_trap_graph_counter += 1

            if self.vel_trap_graph_counter == self.vel_trap_graph_frequency:
                self.vel_trap_graph_counter = 0
                average_velocity = Array3()
                low_point = Array3()

                low_point[0] = self.DEM_parameters["VelocityTrapMinX"].GetDouble()
                low_point[1] = self.DEM_parameters["VelocityTrapMinY"].GetDouble()
                low_point[2] = self.DEM_parameters["VelocityTrapMinZ"].GetDouble()
                high_point = Array3()
                high_point[0] = self.DEM_parameters["VelocityTrapMaxX"].GetDouble()
                high_point[1] = self.DEM_parameters["VelocityTrapMaxY"].GetDouble()
                high_point[2] = self.DEM_parameters["VelocityTrapMaxZ"].GetDouble()

                average_velocity = self.post_utilities.VelocityTrap(self.spheres_model_part, low_point, high_point)

                if compute_flow:
                    vector_of_inner_nodes = []
                    for node in self.spheres_model_part.Nodes:
                        if (node.X > low_point[0]) & (node.Y > low_point[1]) & (node.Z > low_point[2]) & (node.X < high_point[0]) & (node.Y < high_point[1]) & (node.Z < high_point[2]):
                            vector_of_inner_nodes.append(node)

                    crossing_spheres = 0
                    crossing_volume = 0.0

                    for node in vector_of_inner_nodes:
                        id_found = False
                        for previous_node in self.previous_vector_of_inner_nodes:
                            if node.Id == previous_node.Id:
                                id_found = True
                                break
                        # This only happens if None of the previous nodes were capable of setting id_found = True.
                        if id_found is False:
                            crossing_spheres = crossing_spheres + 1
                            radius = node.GetSolutionStepValue(RADIUS)
                            crossing_volume = crossing_volume + 4.0 / 3.0 * math.pi * radius * radius * radius

                    time_between_measures = self.spheres_model_part.ProcessInfo.GetValue(TIME) - self.previous_time
                    number_of_spheres_flow = float(crossing_spheres) / time_between_measures
                    net_volume_flow = crossing_volume / time_between_measures

                    self.previous_time = self.spheres_model_part.ProcessInfo.GetValue(TIME)
                    self.previous_vector_of_inner_nodes = vector_of_inner_nodes

                absolute_path_to_file = os.path.join(graphs_path, filename)
                f = open(absolute_path_to_file, 'a')
                tmp = str(time_dem) + "   " + str(average_velocity[0]) + "   " + str(average_velocity[1]) + "   " + str(average_velocity[2])
                if compute_flow:
                    tmp = tmp + "   " + str(net_volume_flow) + "   " + str(number_of_spheres_flow)
                tmp = tmp + "\n"

                f.write(tmp)
                self.Flush(f)

    @classmethod
    def PrintEulerAngles(self, spheres_model_part, cluster_model_part):
        PostUtilities().ComputeEulerAngles(spheres_model_part, cluster_model_part)


class DEMEnergyCalculator(object):

    def __init__(self, DEM_parameters, spheres_model_part, cluster_model_part, graphs_path, energy_plot):

        self.calculate_option = False

        if "EnergyCalculationOption" in DEM_parameters.keys():
            if DEM_parameters["EnergyCalculationOption"].GetBool():
                self.calculate_option = True
                self.DEM_parameters = DEM_parameters
                self.SpheresModelPart = spheres_model_part
                self.ClusterModelPart = cluster_model_part
                absolute_path_to_file = os.path.join(graphs_path, energy_plot)
                self.energy_plot = open(absolute_path_to_file, 'w')
                self.SpheresEnergyUtil = SphericElementGlobalPhysicsCalculator(spheres_model_part)
                self.ClusterEnergyUtil = SphericElementGlobalPhysicsCalculator(cluster_model_part)
                self.PotentialEnergyReferencePoint = Array3()
                self.PotentialEnergyReferencePoint[0] = self.DEM_parameters["PotentialEnergyReferencePointX"].GetDouble()
                self.PotentialEnergyReferencePoint[1] = self.DEM_parameters["PotentialEnergyReferencePointY"].GetDouble()
                self.PotentialEnergyReferencePoint[2] = self.DEM_parameters["PotentialEnergyReferencePointZ"].GetDouble()
                self.translational_kinematic_energy = 0.0
                self.rotational_kinematic_energy = 0.0
                self.kinematic_energy = 0.0
                self.gravitational_energy = 0.0
                self.elastic_energy = 0.0
                self.inelastic_frictional_energy = 0.0
                self.inelastic_viscodamping_energy = 0.0
                self.external_energy = 0.0
                self.total_energy = 0.0
                self.graph_frequency = int(self.DEM_parameters["GraphExportFreq"].GetDouble() / spheres_model_part.ProcessInfo.GetValue(DELTA_TIME))  # TODO: change the name GraphExportFreq to GraphExportTimeInterval
                self.energy_graph_counter = 0
                self.energy_plot.write(str("Time").rjust(9) + "   " + str("Trans kinematic energy").rjust(22) + "   " + str("Rot kinematic energy").rjust(20) + "   " + str("Kinematic energy").rjust(16) + "   " + str("Gravitational energy").rjust(20) + "   " + str("Elastic energy").rjust(14) + "   " + str("Frictional energy").rjust(16) + "   " + str("Viscodamping energy").rjust(19) + "   " + str("Total energy").rjust(12) + "\n")

    def CalculateEnergyAndPlot(self, time):
        if self.calculate_option:
            if self.energy_graph_counter == self.graph_frequency:
                self.energy_graph_counter = 0
                self.CalculateEnergy()
                self.PlotEnergyGraph(time)

            self.energy_graph_counter += 1

    def CalculateEnergy(self):

        self.translational_kinematic_energy = self.SpheresEnergyUtil.CalculateTranslationalKinematicEnergy(self.SpheresModelPart) + self.ClusterEnergyUtil.CalculateTranslationalKinematicEnergy(self.ClusterModelPart)
        self.rotational_kinematic_energy = self.SpheresEnergyUtil.CalculateRotationalKinematicEnergy(self.SpheresModelPart) + self.ClusterEnergyUtil.CalculateRotationalKinematicEnergy(self.ClusterModelPart)
        self.kinematic_energy = self.translational_kinematic_energy + self.rotational_kinematic_energy
        self.gravitational_energy = self.SpheresEnergyUtil.CalculateGravitationalPotentialEnergy(self.SpheresModelPart, self.PotentialEnergyReferencePoint) + self.ClusterEnergyUtil.CalculateGravitationalPotentialEnergy(self.ClusterModelPart, self.PotentialEnergyReferencePoint)
        self.elastic_energy = self.SpheresEnergyUtil.CalculateElasticEnergy(self.SpheresModelPart) + self.ClusterEnergyUtil.CalculateElasticEnergy(self.ClusterModelPart)
        self.inelastic_frictional_energy = self.SpheresEnergyUtil.CalculateInelasticFrictionalEnergy(self.SpheresModelPart) + self.ClusterEnergyUtil.CalculateInelasticFrictionalEnergy(self.ClusterModelPart)
        self.inelastic_viscodamping_energy = self.SpheresEnergyUtil.CalculateInelasticViscodampingEnergy(self.SpheresModelPart) + self.ClusterEnergyUtil.CalculateInelasticViscodampingEnergy(self.ClusterModelPart)
        self.total_energy = self.kinematic_energy + self.gravitational_energy + self.elastic_energy + self.inelastic_frictional_energy + self.inelastic_viscodamping_energy

    def PlotEnergyGraph(self, time):

        plot_kinematic = self.kinematic_energy
        plot_translational_kinematic = self.translational_kinematic_energy
        plot_rotational_kinematic = self.rotational_kinematic_energy
        plot_gravitational = self.gravitational_energy
        plot_elastic = self.elastic_energy
        plot_inelastic_frictional = self.inelastic_frictional_energy
        plot_inelastic_viscodamping = self.inelastic_viscodamping_energy
        plot_total = self.total_energy
        self.energy_plot.write(str("%.8g" % time).rjust(9) + "   " + str("%.6g" % plot_translational_kinematic).rjust(22) + "   " + str("%.6g" % plot_rotational_kinematic).rjust(20) + "   " + str("%.6g" % plot_kinematic).rjust(16) + "   " + str("%.6g"%plot_gravitational).rjust(20) + "   " + str("%.6g" % plot_elastic).rjust(14) + "   " + str("%.6g" % plot_inelastic_frictional).rjust(16) + "   " + str("%.6g" % plot_inelastic_viscodamping).rjust(19) + "   " + str("%.6g" % plot_total).rjust(12) + '\n')
        self.energy_plot.flush()

    def FinalizeEnergyPlot(self):
        if self.calculate_option:
            self.energy_plot.close


class Procedures(object):

    def __init__(self, DEM_parameters):

        # GLOBAL VARIABLES OF THE SCRIPT
        # Defining list of skin particles (For a test tube of height 30 cm and diameter 15 cm)

        # Initialization of member variables
        self.DEM_parameters = DEM_parameters

        # SIMULATION FLAGS
        self.rotation_OPTION = self.DEM_parameters["RotationOption"].GetBool()
        self.bounding_box_OPTION = self.DEM_parameters["BoundingBoxOption"].GetBool()
        self.automatic_bounding_box_OPTION = self.DEM_parameters["AutomaticBoundingBoxOption"].GetBool()

        self.contact_mesh_OPTION = False
        if "ContactMeshOption" in self.DEM_parameters.keys():
            self.contact_mesh_OPTION = self.DEM_parameters["ContactMeshOption"].GetBool()

        # SIMULATION SETTINGS
        self.b_box_minX = self.DEM_parameters["BoundingBoxMinX"].GetDouble()
        self.b_box_minY = self.DEM_parameters["BoundingBoxMinY"].GetDouble()
        self.b_box_minZ = self.DEM_parameters["BoundingBoxMinZ"].GetDouble()
        self.b_box_maxX = self.DEM_parameters["BoundingBoxMaxX"].GetDouble()
        self.b_box_maxY = self.DEM_parameters["BoundingBoxMaxY"].GetDouble()
        self.b_box_maxZ = self.DEM_parameters["BoundingBoxMaxZ"].GetDouble()
        self.bounding_box_enlargement_factor = self.DEM_parameters["BoundingBoxEnlargementFactor"].GetDouble()

        # MODEL
        self.domain_size = self.DEM_parameters["Dimension"].GetInt()
        self.aux = AuxiliaryUtilities()

    def Barrier(self):
        pass

    def SetTranslationalScheme(self):
        if self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Forward_Euler':
            translational_scheme = ForwardEulerScheme()
        elif self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Symplectic_Euler':
            translational_scheme = SymplecticEulerScheme()
        elif self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Taylor_Scheme':
            translational_scheme = TaylorScheme()
        elif (self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Velocity_Verlet'):
            translational_scheme = VelocityVerletScheme()
        else:
            self.KRATOSprint('Error: selected translational integration scheme not defined. Please select a different scheme')
            sys.exit("\nExecution was aborted.\n")
        return translational_scheme

    def SetRotationalScheme(self):
        if self.DEM_parameters["RotationalIntegrationScheme"].GetString() == 'Direct_Integration':
            if self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Forward_Euler':
                rotational_scheme = ForwardEulerScheme()
            elif self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Symplectic_Euler':
                rotational_scheme = SymplecticEulerScheme()
            elif self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Taylor_Scheme':
                rotational_scheme = TaylorScheme()
            elif (self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Velocity_Verlet'):
                rotational_scheme = VelocityVerletScheme()
        elif self.DEM_parameters["RotationalIntegrationScheme"].GetString() == 'Runge_Kutta':
            rotational_scheme = RungeKuttaScheme()
        elif self.DEM_parameters["RotationalIntegrationScheme"].GetString() == 'Quaternion_Integration':
            rotational_scheme = QuaternionIntegrationScheme()
        else:
            self.KRATOSprint('Error: selected rotational integration scheme not defined. Please select a different scheme')
            sys.exit("\nExecution was aborted.\n")
        return rotational_scheme

    def AddAllVariablesInAllModelParts(self, solver, translational_scheme, rotational_scheme, all_model_parts, DEM_parameters):

        spheres_model_part = all_model_parts.Get('SpheresPart')
        cluster_model_part = all_model_parts.Get('ClusterPart')
        DEM_inlet_model_part = all_model_parts.Get('DEMInletPart')
        rigid_face_model_part = all_model_parts.Get('RigidFacePart')

        self.solver = weakref.proxy(solver)
        self.translational_scheme = weakref.proxy(translational_scheme)
        self.rotational_scheme = weakref.proxy(rotational_scheme)
        self.AddCommonVariables(spheres_model_part, DEM_parameters)
        self.AddSpheresVariables(spheres_model_part, DEM_parameters)
        self.AddMpiVariables(spheres_model_part)
        self.solver.AddAdditionalVariables(spheres_model_part, DEM_parameters)
        self.AddCommonVariables(cluster_model_part, DEM_parameters)
        self.AddClusterVariables(cluster_model_part, DEM_parameters)
        self.AddMpiVariables(cluster_model_part)
        self.AddCommonVariables(DEM_inlet_model_part, DEM_parameters)
        self.AddSpheresVariables(DEM_inlet_model_part, DEM_parameters)
        self.solver.AddAdditionalVariables(DEM_inlet_model_part, DEM_parameters)
        self.AddCommonVariables(rigid_face_model_part, DEM_parameters)
        self.AddRigidFaceVariables(rigid_face_model_part, DEM_parameters)
        self.AddMpiVariables(rigid_face_model_part)

    @classmethod
    def AddCommonVariables(self, model_part, DEM_parameters):
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)

    def AddSpheresVariables(self, model_part, DEM_parameters):

        # KINEMATIC
        # TODO: only if self.DEM_parameters-RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(DELTA_ROTATION)
        # TODO: only if self.DEM_parameters-RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
        # TODO: only if self.DEM_parameters-RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
        model_part.AddNodalSolutionStepVariable(NORMAL_IMPACT_VELOCITY)
        model_part.AddNodalSolutionStepVariable(TANGENTIAL_IMPACT_VELOCITY)
        # TODO: only if self.DEM_parameters-RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(ORIENTATION)
        # JIG: SHOULD BE REMOVED IN THE FUTURE
        model_part.AddNodalSolutionStepVariable(ORIENTATION_REAL)
        # JIG: SHOULD BE REMOVED IN THE FUTURE
        model_part.AddNodalSolutionStepVariable(ORIENTATION_IMAG)
        # TODO: only if self.DEM_parameters-RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(ANGULAR_MOMENTUM)
        model_part.AddNodalSolutionStepVariable(FACE_NORMAL_IMPACT_VELOCITY)
        model_part.AddNodalSolutionStepVariable(FACE_TANGENTIAL_IMPACT_VELOCITY)
        model_part.AddNodalSolutionStepVariable(LINEAR_IMPULSE)
        # ****************** Quaternion Integration BEGIN ******************
        # TODO: only if self.DEM_parameters-RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(LOCAL_AUX_ANGULAR_VELOCITY)
        # TODO: only if self.DEM_parameters-RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(AUX_ORIENTATION)
        # ******************* Quaternion Integration END *******************

        # FORCES
        model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
        model_part.AddNodalSolutionStepVariable(LOCAL_CONTACT_FORCE)
        model_part.AddNodalSolutionStepVariable(CONTACT_FORCES)
        model_part.AddNodalSolutionStepVariable(RIGID_ELEMENT_FORCE)
        model_part.AddNodalSolutionStepVariable(DAMP_FORCES)
        # TODO: only if self.DEM_parameters-RotationOption! Check that no one accesses them in c++ without checking the rotation option
        model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_FORCE)
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_MOMENT)
        model_part.AddNodalSolutionStepVariable(FORCE_REACTION)
        model_part.AddNodalSolutionStepVariable(MOMENT_REACTION)

        # BASIC PARTICLE PROPERTIES
        model_part.AddNodalSolutionStepVariable(RADIUS)
        model_part.AddNodalSolutionStepVariable(NODAL_MASS)
        model_part.AddNodalSolutionStepVariable(REPRESENTATIVE_VOLUME)
        model_part.AddNodalSolutionStepVariable(NEIGHBOUR_SIZE)
        model_part.AddNodalSolutionStepVariable(NEIGHBOUR_RATIO)

        # ROTATION RELATED PROPERTIES
        if self.DEM_parameters["RotationOption"].GetBool():
            # TODO: only if self.DEM_parameters-RotationOption! Check that no one accesses them in c++ without checking the rotation option
            model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA)
            # TODO: only if self.DEM_parameters-RotationOption! Check that no one accesses them in c++ without checking the rotation option
            model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_DAMP_RATIO)
            if self.DEM_parameters["RollingFrictionOption"].GetBool():
                model_part.AddNodalSolutionStepVariable(ROLLING_RESISTANCE_MOMENT)

        # OTHER PROPERTIES
        model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD

        if "PostSkinSphere" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostSkinSphere"].GetBool():
                model_part.AddNodalSolutionStepVariable(SKIN_SPHERE)

        # LOCAL AXIS
        if DEM_parameters["PostEulerAngles"].GetBool():
            model_part.AddNodalSolutionStepVariable(EULER_ANGLES)

        if "PostStressStrainOption" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostStressStrainOption"].GetBool():
                model_part.AddNodalSolutionStepVariable(DEM_STRESS_TENSOR)

        if self.solver.poisson_ratio_option:
            model_part.AddNodalSolutionStepVariable(POISSON_VALUE)

        # Nano Particle
        if self.DEM_parameters["ElementType"].GetString() == "SwimmingNanoParticle":
            model_part.AddNodalSolutionStepVariable(CATION_CONCENTRATION)
            model_part.AddNodalSolutionStepVariable(DRAG_COEFFICIENT)

        # ONLY VISUALIZATION
        if self.DEM_parameters["PostExportId"].GetBool():  # TODO: add suffix Option
            model_part.AddNodalSolutionStepVariable(EXPORT_ID)

        #model_part.AddNodalSolutionStepVariable(SPRAYED_MATERIAL)

    @classmethod
    def AddRigidFaceVariables(self, model_part, DEM_parameters):

        model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
        model_part.AddNodalSolutionStepVariable(CONTACT_FORCES)
        model_part.AddNodalSolutionStepVariable(DEM_PRESSURE)
        model_part.AddNodalSolutionStepVariable(TANGENTIAL_ELASTIC_FORCES)
        model_part.AddNodalSolutionStepVariable(SHEAR_STRESS)
        model_part.AddNodalSolutionStepVariable(DEM_NODAL_AREA)
        model_part.AddNodalSolutionStepVariable(NON_DIMENSIONAL_VOLUME_WEAR)
        model_part.AddNodalSolutionStepVariable(IMPACT_WEAR)
        model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
        model_part.AddNodalSolutionStepVariable(DELTA_ROTATION)
        model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
        model_part.AddNodalSolutionStepVariable(LOCAL_ANGULAR_VELOCITY)
        model_part.AddNodalSolutionStepVariable(LOCAL_AUX_ANGULAR_VELOCITY)
        model_part.AddNodalSolutionStepVariable(ORIENTATION_REAL) # JIG: SHOULD BE REMOVED IN THE FUTURE
        model_part.AddNodalSolutionStepVariable(ORIENTATION_IMAG) # JIG: SHOULD BE REMOVED IN THE FUTURE
        model_part.AddNodalSolutionStepVariable(ORIENTATION)
        model_part.AddNodalSolutionStepVariable(AUX_ORIENTATION)
        model_part.AddNodalSolutionStepVariable(ANGULAR_MOMENTUM)

        # FORCES
        model_part.AddNodalSolutionStepVariable(RIGID_ELEMENT_FORCE)
        model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_FORCE)
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_MOMENT)

        # PHYSICAL PROPERTIES
        model_part.AddNodalSolutionStepVariable(PRINCIPAL_MOMENTS_OF_INERTIA)
        model_part.AddNodalSolutionStepVariable(CLUSTER_VOLUME)
        model_part.AddNodalSolutionStepVariable(NODAL_MASS)
        model_part.AddNodalSolutionStepVariable(CHARACTERISTIC_LENGTH)
        model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)

    def AddElasticFaceVariables(self, model_part, DEM_parameters): #Only used in CSM coupling
        self.AddRigidFaceVariables(model_part,self.DEM_parameters)
        model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)

    @classmethod
    def AddClusterVariables(self, model_part, DEM_parameters):
        # KINEMATIC
        model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)
        model_part.AddNodalSolutionStepVariable(DELTA_ROTATION)
        model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
        model_part.AddNodalSolutionStepVariable(LOCAL_ANGULAR_VELOCITY)
        model_part.AddNodalSolutionStepVariable(ORIENTATION)
        # JIG: SHOULD BE REMOVED IN THE FUTURE
        model_part.AddNodalSolutionStepVariable(ORIENTATION_REAL)
        # JIG: SHOULD BE REMOVED IN THE FUTURE
        model_part.AddNodalSolutionStepVariable(ORIENTATION_IMAG)
        model_part.AddNodalSolutionStepVariable(ANGULAR_MOMENTUM)
        # ****************** Quaternion Integration BEGIN ******************
        model_part.AddNodalSolutionStepVariable(LOCAL_AUX_ANGULAR_VELOCITY)
        model_part.AddNodalSolutionStepVariable(AUX_ORIENTATION)
        # ******************* Quaternion Integration END *******************

        # FORCES
        model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
        model_part.AddNodalSolutionStepVariable(RIGID_ELEMENT_FORCE)
        model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_FORCE)
        model_part.AddNodalSolutionStepVariable(EXTERNAL_APPLIED_MOMENT)

        # PHYSICAL PROPERTIES
        model_part.AddNodalSolutionStepVariable(PRINCIPAL_MOMENTS_OF_INERTIA)
        model_part.AddNodalSolutionStepVariable(CLUSTER_VOLUME)
        model_part.AddNodalSolutionStepVariable(NODAL_MASS)
        model_part.AddNodalSolutionStepVariable(CHARACTERISTIC_LENGTH)
        model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)

        # OTHER PROPERTIES
        model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)   # Colour defined in GiD

        # LOCAL AXIS
        if DEM_parameters["PostEulerAngles"].GetBool():
            model_part.AddNodalSolutionStepVariable(EULER_ANGLES)

    def AddMpiVariables(self, model_part):
        pass

    def SetInitialNodalValues(self, spheres_model_part, cluster_model_part, DEM_inlet_model_part, rigid_face_model_part):
        pass

    @classmethod
    def SetUpBufferSizeInAllModelParts(self, spheres_model_part, spheres_b_size, cluster_model_part, clusters_b_size, DEM_inlet_model_part, inlet_b_size, rigid_face_model_part, rigid_b_size):
        spheres_model_part.SetBufferSize(spheres_b_size)
        cluster_model_part.SetBufferSize(clusters_b_size)
        DEM_inlet_model_part.SetBufferSize(inlet_b_size)
        rigid_face_model_part.SetBufferSize(rigid_b_size)

    @classmethod
    def FindMaxNodeIdAccrossModelParts(self, creator_destructor, all_model_parts):

        max_candidates = []
        max_candidates.append(creator_destructor.FindMaxNodeIdInModelPart(all_model_parts.Get("SpheresPart")))
        max_candidates.append(creator_destructor.FindMaxElementIdInModelPart(all_model_parts.Get("SpheresPart")))
        max_candidates.append(creator_destructor.FindMaxNodeIdInModelPart(all_model_parts.Get("RigidFacePart")))
        max_candidates.append(creator_destructor.FindMaxNodeIdInModelPart(all_model_parts.Get("ClusterPart")))

        return max(max_candidates)

    def ModelData(self, spheres_model_part, solver):

        contact_model_part = solver.contact_model_part
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
        var = 0.0
        rel_std_dev = 0.0

        for node in spheres_model_part.Nodes:

            sum_radi += node.GetSolutionStepValue(RADIUS)
            partial_sum_squared = node.GetSolutionStepValue(RADIUS) ** 2.0
            total_sum_squared += partial_sum_squared
            volume += 4 * 3.141592 / 3 * node.GetSolutionStepValue(RADIUS) ** 3.0
            area += 3.141592 * partial_sum_squared
            i += 1.0

        if i > 0.0:
            mean = sum_radi / i
            var = total_sum_squared / i - mean ** 2.0
        std_dev = 0.0

        if abs(var) > 1e-9:
            std_dev = var ** 0.5

        if i > 0.0:
            rel_std_dev = std_dev / mean

        Model_Data.write("Radius Mean: " + str(mean) + '\n')
        Model_Data.write("Std Deviation: " + str(std_dev) + '\n')
        Model_Data.write("Relative Std Deviation: " + str(rel_std_dev) + '\n')
        Model_Data.write("Total Particle Volume 3D: " + str(volume) + '\n')
        Model_Data.write("Total Particle Area 2D: " + str(area) + '\n')
        Model_Data.write('\n')

        Total_Particles = len(spheres_model_part.Nodes)

        if solver.continuum_type:
            Coordination_Number = 0.0

            if self.contact_mesh_OPTION:
                Total_Contacts = contact_model_part.NumberOfElements(0)

                if Total_Particles:
                    Coordination_Number = 2.0 * Total_Contacts / Total_Particles

            Model_Data.write("Total Number of Particles: " + str(Total_Particles) + '\n')
            Model_Data.write("Total Number of Bonds: " + str(Total_Contacts) + '\n')
            Model_Data.write("Bonded Coordination Number NC: " + str(Coordination_Number) + '\n')
            Model_Data.write('\n')
            #Model_Data.write("Volume Elements: " + str(total_volume) + '\n')
            self.KRATOSprint("Coordination Number: " + str(Coordination_Number) + "\n")

        Model_Data.close()

    @classmethod
    def MonitorPhysicalProperties(self, model_part, physics_calculator, properties_list):

        # This function returns a list of arrays (also lists)
        # Each array contains the values of the physical properties at the current time
        time = model_part.ProcessInfo.GetValue(TIME)
        present_prop = []

        if not properties_list:  # The first array in the list only contains the entries names
            names = []
            names.append("time")
            names.append("mass")
            names.append("gravitational_energy")
            names.append("kinematic_energy")
            #names.append("elastic_energy")
            names.append("momentum")
            names.append("angular_momentum")
            names.append("total_energy")

            properties_list.append(names)

        # Calculating current values
        mass = physics_calculator.CalculateTotalMass(model_part)
        center = physics_calculator.CalculateCenterOfMass(model_part)
        initial_center = physics_calculator.GetInitialCenterOfMass()
        gravity_energy = physics_calculator.CalculateGravitationalPotentialEnergy(model_part, initial_center)
        kinematic_energy = physics_calculator.CalculateKinematicEnergy(model_part)
        #elastic_energy = physics_calculator.CalculateElasticEnergy(model_part)
        momentum = physics_calculator.CalculateTotalMomentum(model_part)
        angular_momentum = physics_calculator.CalulateTotalAngularMomentum(model_part)
        total_energy = gravity_energy + kinematic_energy  # + elastic_energy

        # Filling in the entries values corresponding to the entries names above
        present_prop.append(time)
        present_prop.append(mass)
        present_prop.append(gravity_energy)
        present_prop.append(kinematic_energy)
        #present_prop.append(elastic_energy)
        present_prop.append(momentum)
        present_prop.append(angular_momentum)
        present_prop.append(total_energy)

        properties_list.append(present_prop)

        return properties_list

    @classmethod
    def RemoveFoldersWithResults(self, main_path, problem_name, run_code=''):
        shutil.rmtree(os.path.join(main_path, problem_name + '_Post_Files' + run_code), ignore_errors=True)
        shutil.rmtree(os.path.join(main_path, problem_name + '_Graphs'), ignore_errors=True)
        shutil.rmtree(os.path.join(main_path, problem_name + '_Results_and_Data'), ignore_errors=True)
        shutil.rmtree(os.path.join(main_path, problem_name + '_MPI_results'), ignore_errors=True)

        try:
            file_to_remove = os.path.join(main_path, problem_name)+"DEM.time"
            os.remove(file_to_remove)
        except OSError:
            pass
        try:
            file_to_remove = os.path.join(main_path, problem_name)+"DEM_Inlet.time"
            os.remove(file_to_remove)
        except OSError:
            pass

        try:
            file_to_remove = os.path.join(main_path, problem_name)+"DEM_FEM_boundary.time"
            os.remove(file_to_remove)
        except OSError:
            pass

        try:
            file_to_remove = os.path.join(main_path, problem_name)+"DEM_Clusters.time"
            os.remove(file_to_remove)
        except OSError:
            pass

        try:
            file_to_remove = os.path.join(main_path, "TimesPartialRelease")
            os.remove(file_to_remove)
        except OSError:
            pass

        try:
            file_to_remove = os.path.join(main_path, problem_name)+".post.lst"
            os.remove(file_to_remove)
        except OSError:
            pass


    @classmethod
    def CreateDirectories(self, main_path, problem_name, run_code=''):

        root = os.path.join(main_path, problem_name)
        post_path = root + '_Post_Files' + run_code
        data_and_results = root + '_Results_and_Data'
        graphs_path = root + '_Graphs'
        MPI_results = root + '_MPI_results'

        self.RemoveFoldersWithResults(main_path, problem_name, run_code)

        for directory in [post_path, data_and_results, graphs_path, MPI_results]:
            if not os.path.isdir(directory):
                os.makedirs(str(directory))

        return [post_path, data_and_results, graphs_path, MPI_results]

    @classmethod
    def FindMaxNodeIdInModelPart(self, model_part):

        maxid = 0

        for node in model_part.Nodes:
            if node.Id > maxid:
                maxid = node.Id

        return maxid

    def SetBoundingBoxLimits(self, all_model_parts, creator_destructor):

        bounding_box_time_limits = []
        if self.DEM_parameters["BoundingBoxOption"].GetBool():
            self.SetBoundingBox(all_model_parts.Get("SpheresPart"), all_model_parts.Get("ClusterPart"), all_model_parts.Get("RigidFacePart"), all_model_parts.Get("DEMInletPart"), creator_destructor)
            bounding_box_time_limits = [self.solver.bounding_box_start_time, self.solver.bounding_box_stop_time]
            return bounding_box_time_limits

    def SetBoundingBox(self, spheres_model_part, clusters_model_part, rigid_faces_model_part, DEM_inlet_model_part, creator_destructor):

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
        creator_destructor.CalculateSurroundingBoundingBox(spheres_model_part, clusters_model_part, rigid_faces_model_part, DEM_inlet_model_part, self.bounding_box_enlargement_factor, self.automatic_bounding_box_OPTION)

    @classmethod
    def DeleteFiles(self):
        files_to_delete_list = glob('*.time')
        for to_erase_file in files_to_delete_list:
            try:
                os.remove(to_erase_file)
            except OSError:
                pass

    def PreProcessModel(self, DEM_parameters):
        pass

    def CheckVariableType(self, var, expected_type, msg):
        actual_type = type(var)
        if actual_type is int and expected_type is float:
            return
        if actual_type is not expected_type:
            self.KRATOSprint(
                "**************************************************************************")
            self.KRATOSprint(
                "ERROR: Input parameter of wrong type in file 'DEM_explicit_solver_var.py'.")
            a = str(expected_type)
            b = str(var)
            self.KRATOSprint("The type expected was " +
                             a + " but " + b + " was read.")
            self.KRATOSprint(
                "**************************************************************************")
            sys.exit()

    @classmethod
    def Flush(self, a):
        a.flush()

    def KRATOSprint(self, message):
        Logger.Print(message, label="DEM")
        Logger.Flush()


class DEMFEMProcedures(object):

    def __init__(self, DEM_parameters, graphs_path, spheres_model_part, RigidFace_model_part):

        # GLOBAL VARIABLES OF THE SCRIPT
        self.DEM_parameters = DEM_parameters

        if not "TestType" in DEM_parameters.keys():
            self.TestType = "None"
        else:
            self.TestType = self.DEM_parameters["TestType"].GetString()

        # Initialization of member variables
        # SIMULATION FLAGS
        # TODO: Why is this in DEM FEM Procs also?
        self.rotation_OPTION = self.DEM_parameters["RotationOption"].GetBool()
        self.bounding_box_OPTION = self.DEM_parameters["BoundingBoxOption"].GetBool()

        # TODO: This is already in the Procedures object. why to repeat it?
        self.contact_mesh_OPTION = False
        if "ContactMeshOption" in self.DEM_parameters.keys():
            self.contact_mesh_OPTION = self.DEM_parameters["ContactMeshOption"].GetBool()

        self.graphs_path = graphs_path
        self.spheres_model_part = spheres_model_part
        self.RigidFace_model_part = RigidFace_model_part
        #self.solver = solver
        self.aux = AuxiliaryUtilities()

        self.fem_mesh_nodes = []

        self.graph_counter = 1
        self.balls_graph_counter = 0

        self.graph_frequency = int((self.DEM_parameters["GraphExportFreq"].GetDouble() / spheres_model_part.ProcessInfo.GetValue(DELTA_TIME))+1.0)
        if self.graph_frequency < 1:
            # that means it is not possible to print results with a higher frequency than the computations delta time
            self.graph_frequency = 1

        self.mesh_motion = DEMFEMUtilities()

        def Flush(self, a):
            a.flush()

        def open_graph_files(self, RigidFace_model_part):
            for smp in self.RigidFace_model_part.SubModelParts:
                if smp[FORCE_INTEGRATION_GROUP]:
                    identifier = smp[IDENTIFIER]
                    absolute_path_to_file = os.path.join(self.graphs_path, str(self.DEM_parameters["problem_name"].GetString()) + "_" + str(identifier) + "_force_graph.grf")
                    self.graph_forces[identifier] = open(absolute_path_to_file, 'w')
                    self.graph_forces[identifier].write(str("#time").rjust(12) + " " + str("total_force[0]").rjust(13) + " " + str("total_force[1]").rjust(13) + " " + str("total_force[2]").rjust(13) + " " + str("total_moment[0]").rjust(13) + " " + str("total_moment[1]").rjust(13) + " " + str("total_moment[2]").rjust(13) + "\n")

        self.graph_forces = {}

        def open_balls_graph_files(self, spheres_model_part):
            for smp in self.spheres_model_part.SubModelParts:
                if smp[FORCE_INTEGRATION_GROUP]:
                    identifier = smp[IDENTIFIER]
                    absolute_path_to_file = os.path.join(self.graphs_path, str(self.DEM_parameters["problem_name"].GetString()) + "_" + str(identifier) + "_particle_force_graph.grf")
                    self.particle_graph_forces[identifier] = open(absolute_path_to_file, 'w')
                    self.particle_graph_forces[identifier].write(str("#time").rjust(12) + " " + str("total_force_x").rjust(13) + " " + str("total_force_y").rjust(13) + " " + str("total_force_z").rjust(13) + "\n")

        def evaluate_computation_of_fem_results():

            self.spheres_model_part.ProcessInfo.SetValue(COMPUTE_FEM_RESULTS_OPTION, 0)
            elastic_forces = self.DEM_parameters["PostElasticForces"].GetBool()
            tangential_elastic_forces = self.DEM_parameters["PostTangentialElasticForces"].GetBool()
            dem_pressure = self.DEM_parameters["PostPressure"].GetBool()

            if not "PostContactForces" in self.DEM_parameters.keys():
                contact_forces = 0
            else:
                contact_forces = self.DEM_parameters["PostContactForces"].GetBool()

            if not "PostShearStress" in self.DEM_parameters.keys():
                shear_stress = 0
            else:
                shear_stress = self.DEM_parameters["PostShearStress"].GetBool()

            if not "PostNodalArea" in self.DEM_parameters.keys():
                dem_nodal_area = 0
            else:
                dem_nodal_area = self.DEM_parameters["PostNodalArea"].GetBool()

            integration_groups = False

            if self.RigidFace_model_part.NumberOfSubModelParts() > 0:
                for smp in self.RigidFace_model_part.SubModelParts:
                    if smp[FORCE_INTEGRATION_GROUP]:
                        integration_groups = True
                        break
            if elastic_forces or contact_forces or dem_pressure or tangential_elastic_forces or shear_stress or dem_nodal_area or integration_groups:
                self.spheres_model_part.ProcessInfo.SetValue(COMPUTE_FEM_RESULTS_OPTION, 1)

        self.particle_graph_forces = {}

        if self.TestType == "None":
            open_graph_files(self, RigidFace_model_part)
            open_balls_graph_files(self, spheres_model_part)

        # SIMULATION SETTINGS
        self.bounding_box_enlargement_factor = self.DEM_parameters["BoundingBoxEnlargementFactor"].GetDouble()

        # MODEL
        self.domain_size = self.DEM_parameters["Dimension"].GetInt()
        evaluate_computation_of_fem_results()

    def MoveAllMeshes(self, all_model_parts, time, dt):

        spheres_model_part = all_model_parts.Get("SpheresPart")
        DEM_inlet_model_part = all_model_parts.Get("DEMInletPart")
        rigid_face_model_part = all_model_parts.Get("RigidFacePart")
        cluster_model_part = all_model_parts.Get("ClusterPart")

        self.mesh_motion.MoveAllMeshes(rigid_face_model_part, time, dt)
        self.mesh_motion.MoveAllMeshes(spheres_model_part, time, dt)
        self.mesh_motion.MoveAllMeshes(DEM_inlet_model_part, time, dt)
        self.mesh_motion.MoveAllMeshes(cluster_model_part, time, dt)

    # def MoveAllMeshesUsingATable(self, model_part, time, dt):

    #     self.mesh_motion.MoveAllMeshesUsingATable(model_part, time, dt)

        # for smp in model_part.SubModelParts:

        #     if not smp[TABLE_NUMBER]:
        #         continue

        #     Logger.Print("Info:", label="")
        #     Logger.Print(smp[IDENTIFIER], label="")
        #     Logger.Print(smp[TABLE_NUMBER], label="")

        #     for node in smp.Nodes:

        #         old_coords = Vector(3)
        #         old_coords[0] = node.X
        #         old_coords[1] = node.Y
        #         old_coords[2] = node.Z

        #         velocity = Vector(3)
        #         velocity[0] = model_part.GetTable(smp[TABLE_NUMBER]).GetValue(time)
        #         velocity[1] = 0.0
        #         velocity[2] = 0.0
        #         node.SetSolutionStepValue(VELOCITY, velocity)

        #         node.X = old_coords[0] + velocity[0] * dt
        #         node.Y = old_coords[1] + velocity[1] * dt
        #         node.Z = old_coords[2] + velocity[2] * dt

        #         displacement = Vector(3)
        #         displacement[0] = node.X - node.X0
        #         displacement[1] = node.Y - node.Y0
        #         displacement[2] = node.Z - node.Z0
        #         node.SetSolutionStepValue(DISPLACEMENT, displacement)

    @classmethod
    def UpdateTimeInModelParts(self, all_model_parts, time, dt, step, is_time_to_print):

        spheres_model_part = all_model_parts.Get("SpheresPart")
        cluster_model_part = all_model_parts.Get("ClusterPart")
        DEM_inlet_model_part = all_model_parts.Get("DEMInletPart")
        rigid_face_model_part = all_model_parts.Get("RigidFacePart")

        self.UpdateTimeInOneModelPart(spheres_model_part, time, dt, step, is_time_to_print)
        self.UpdateTimeInOneModelPart(cluster_model_part, time, dt, step, is_time_to_print)
        self.UpdateTimeInOneModelPart(DEM_inlet_model_part, time, dt, step, is_time_to_print)
        self.UpdateTimeInOneModelPart(rigid_face_model_part, time, dt, step, is_time_to_print)

    @classmethod
    def UpdateTimeInOneModelPart(self, model_part, time, dt, step, is_time_to_print):
        model_part.ProcessInfo[TIME] = time
        model_part.ProcessInfo[DELTA_TIME] = dt
        model_part.ProcessInfo[TIME_STEPS] = step
        model_part.ProcessInfo[IS_TIME_TO_PRINT] = is_time_to_print

    def close_graph_files(self, RigidFace_model_part):

        for smp in self.RigidFace_model_part.SubModelParts:
            if smp[FORCE_INTEGRATION_GROUP]:
                identifier = smp[IDENTIFIER]
                self.graph_forces[identifier].close()

    def close_balls_graph_files(self, spheres_model_part):

        for mesh_number in range(0, self.spheres_model_part.NumberOfSubModelParts()):
            if smp[FORCE_INTEGRATION_GROUP]:
                identifier = smp[IDENTIFIER]
                self.particle_graph_forces[identifier].close()

    @classmethod
    def PrintPoisson(self, model_part, DEM_parameters, filename, time):

        if DEM_parameters["Dimension"].GetInt() == 3:
            poisson, dummy, _ = PostUtilities().ComputePoisson(model_part)
        else:
            poisson, dummy, _ = PostUtilities().ComputePoisson2D(model_part)

        file_open = open(filename, 'a')
        data = str(time) + "  " + str(poisson) + "\n"
        file_open.write(data)

    def PrintGraph(self, time):

        if self.TestType == "None":

            if self.graph_counter == self.graph_frequency:
                self.graph_counter = 0

                for smp in self.RigidFace_model_part.SubModelParts:
                    if smp[FORCE_INTEGRATION_GROUP]:
                        mesh_nodes = smp.Nodes

                        total_force = Array3()
                        total_force[0] = 0.0
                        total_force[1] = 0.0
                        total_force[2] = 0.0

                        total_moment = Array3()
                        total_moment[0] = 0.0
                        total_moment[1] = 0.0
                        total_moment[2] = 0.0

                        rotation_center = smp[ROTATION_CENTER]

                        PostUtilities().IntegrationOfForces(mesh_nodes, total_force, rotation_center, total_moment)

                        identifier = smp[IDENTIFIER]

                        self.graph_forces[identifier].write(str("%.8g" % time).rjust(12) +
                                                            " " + str("%.6g" % total_force[0]).rjust(13) + " " + str("%.6g" % total_force[1]).rjust(13) +
                                                            " " + str("%.6g" % total_force[2]).rjust(13) + " " + str("%.6g" % total_moment[0]).rjust(13) +
                                                            " " + str("%.6g" % total_moment[1]).rjust(13) + " " + str("%.6g" % total_moment[2]).rjust(13) + "\n")
                        self.graph_forces[identifier].flush()

            self.graph_counter += 1

    def FinalizeGraphs(self, RigidFace_model_part):

        if not "TestType" in self.DEM_parameters.keys():
            self.close_graph_files(RigidFace_model_part)

    def PrintBallsGraph(self, time):

        if not "TestType" in self.DEM_parameters.keys():

            if self.balls_graph_counter == self.graph_frequency:
                self.balls_graph_counter = 0

                for smp in self.spheres_model_part.SubModelParts:
                    if smp[FORCE_INTEGRATION_GROUP]:
                        mesh_nodes = smp.Nodes

                        total_force = Array3()
                        total_force[0] = 0.0
                        total_force[1] = 0.0
                        total_force[2] = 0.0

                        PostUtilities().IntegrationOfElasticForces(mesh_nodes, total_force)

                        identifier = smp[IDENTIFIER]
                        self.particle_graph_forces[identifier].write(str("%.8g" % time).rjust(12) + " " + str("%.6g" % total_force[0]).rjust(13) + " " + str("%.6g" % total_force[1]).rjust(13) + " " + str("%.6g" % total_force[2]).rjust(13) + "\n")
                        self.particle_graph_forces[identifier].flush()

            self.balls_graph_counter += 1

    def FinalizeBallsGraphs(self, spheres_model_part):

        if not "TestType" in self.DEM_parameters.keys():
            self.close_balls_graph_files(spheres_model_part)

    def ApplyNodalRotation(self, time):

            #if (time < 0.5e-2 ) :
        if time < 3.8e-5:

            #while (time < self.DEM_parameters["FinalTime"].GetDouble()):
            #print("TIME STEP BEGINS.  STEP:"+str(time)+"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

            d0 = 1.694
            avance = 0.50100
            distance = (d0 - avance)
            w = 62.8
            distance = 2

            vx = -distance * w * math.sin(w * time)
            #vx = distance * rpm * math.cos(1.0)
            vy = distance * w * math.cos(w * time)
            #vz = - distance * rpm * math.sin(1.0)

            for smp in self.spheres_model_part.SubModelParts:
                if smp[FORCE_INTEGRATION_GROUP]:
                    self.mesh_nodes = smp.Nodes

                    for node in self.mesh_nodes:
                        node.SetSolutionStepValue(VELOCITY_X, vx)
                        node.SetSolutionStepValue(VELOCITY_Y, vy)
                        node.Fix(VELOCITY_X)
                        node.Fix(VELOCITY_Y)
        else:

            d0 = 1.694
            avance = 0.50100
            distance = (d0 - avance)
            w = 62.8
            distance = 2

            vx = -distance * w * math.sin(w * time)
            #vx = distance * rpm * math.cos(1.0)
            vy = distance * w * math.cos(w * time)
            #vz = - distance * rpm * math.sin(1.0)
            radius = 1.0001

            for smp in self.spheres_model_part.SubModelParts:
                if smp[FORCE_INTEGRATION_GROUP]:
                    self.mesh_nodes = smp.Nodes

                for node in self.mesh_nodes:
                    node.SetSolutionStepValue(RADIUS, radius)

            for smp in self.spheres_model_part.SubModelParts:
                if smp[FORCE_INTEGRATION_GROUP]:
                    self.mesh_nodes = smp.Nodes

                    for node in self.mesh_nodes:
                        node.SetSolutionStepValue(VELOCITY_X, vx)
                        node.SetSolutionStepValue(VELOCITY_Y, vy)
                        node.Fix(VELOCITY_X)
                        node.Fix(VELOCITY_Y)


class Report(object):

    def __init__(self):
        pass

    def Prepare(self, timer, control_time):
        self.initial_pr_time = timer.clock()
        self.initial_re_time = timer.time()
        self.prev_time = 0.0
        self.total_steps_expected = 0
        self.control_time = control_time
        self.first_print = True

    def BeginReport(self, timer):
        label = "DEM: "
        report = "Main loop starting..." + "\n" + \
            label + "Total number of TIME STEPs expected in the calculation: " + \
            str(self.total_steps_expected) + "\n" + label

        return report

    def StepiReport(self, timer, time, step):

        incremental_time = (timer.time() - self.initial_re_time) - self.prev_time
        report = ""
        label = "DEM: "

        if incremental_time > self.control_time:
            percentage = 100 * (float(step) / self.total_steps_expected)
            elapsed_time = timer.time() - self.initial_re_time

            report = report + "Real time calculation: " + str(elapsed_time) + " seconds" + "\n"\
                            + label + "In minutes: " + str(elapsed_time / 60.0) + " minutes" + "\n"\
                            + label + "In hours: " + str(elapsed_time / 3600.0) + " hours" + "\n"\
                            + label + "Simulation time: " + str(time) + " seconds" + "\n"\
                            + label + "%s %.5f %s" % ("Percentage Completed: ", percentage, "%") + "\n"\
                            + label + "Computed time steps: " + str(step) + " out of " + str(self.total_steps_expected) + "\n" + label

            self.prev_time = (timer.time() - self.initial_re_time)

        if (timer.time() - self.initial_re_time > 60) and self.first_print and step != 0:
            self.first_print = False
            estimated_sim_duration = 60.0 * (self.total_steps_expected / step)  # seconds

            report = report + "\n" + label + "The total estimated computation time is " + str(estimated_sim_duration) + " seconds" + "\n"\
                + label + "In minutes: " + str(estimated_sim_duration / 60.0) + " minutes" + "\n"\
                + label + "In hours:   " + str(estimated_sim_duration / 3600.0) + " hours" + "\n"\
                + label + "In days:    " + str(estimated_sim_duration / 86400.0) + " days" + "\n" + label

        return report

    def FinalReport(self, timer):
        elapsed_pr_time = timer.clock() - self.initial_pr_time
        elapsed_re_time = timer.time() - self.initial_re_time
        label = "DEM: "

        report = "Calculation ends at instant: " + str(timer.time()) + "\n"\
            + label + "Calculation ends at processing time instant: " + str(timer.clock()) + "\n"\
            + label + "Elapsed processing time: " + str(elapsed_pr_time) + "\n"\
            + label + "Elapsed real time: " + str(elapsed_re_time) + "\n" + label

        report = report + "\n" + label + "ANALYSIS COMPLETED"

        return report


class PreUtils(object):

    def __init__(self):
        pass

class MaterialTest(object):

    def __init__(self):
        pass

    def Initialize(self, DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part):

        if not "TestType" in DEM_parameters.keys():
            self.TestType = "None"
        else:
            self.TestType = DEM_parameters["TestType"].GetString()

        if self.TestType != "None":
            self.script = DEM_material_test_script.MaterialTest(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)
            self.script.Initialize()

            #self.PreUtils = DEM_material_test_script.PreUtils(spheres_model_part)
            #self.PreUtils.BreakBondUtility(spheres_model_part)

    def PrepareDataForGraph(self):
        if self.TestType != "None":
            self.script.PrepareDataForGraph()

    def MeasureForcesAndPressure(self):
        if self.TestType != "None":
            self.script.MeasureForcesAndPressure()

    def PrintGraph(self, time):
        if self.TestType != "None":
            self.script.PrintGraph(time)

    def FinalizeGraphs(self):
        if self.TestType != "None":
            self.script.FinalizeGraphs()

    def PrintChart(self):
        if self.TestType != "None":
            self.script.PrintChart()

    def GenerateGraphics(self):
        if self.TestType != "None":
            self.script.GenerateGraphics()


class MultifileList(object):

    def __init__(self, post_path, name, step, which_folder):
        self.index = 0
        self.step = step
        self.name = name
        self.which_folder = which_folder
        if which_folder == "inner":
            absolute_path_to_file = os.path.join(post_path, "_list_" + self.name + "_" + str(step) + ".post.lst")
        else:
            absolute_path_to_file = os.path.join(post_path, self.name + ".post.lst")

        self.file = open(absolute_path_to_file, "w")


class DEMIo(object):

    def __init__(self, model, DEM_parameters, post_path, all_model_parts):

        self.model = model

        self.post_path = post_path
        self.mixed_model_part = model.CreateModelPart("Mixed_Part")
        self.mixed_spheres_and_clusters_model_part = model.CreateModelPart("MixedSpheresAndClustersPart")
        self.mixed_spheres_not_in_cluster_and_clusters_model_part = model.CreateModelPart("MixedSpheresNotInClusterAndClustersPart")

        self.spheres_model_part = all_model_parts.Get("SpheresPart")
        self.cluster_model_part = all_model_parts.Get("ClusterPart")
        self.rigid_face_model_part = all_model_parts.Get("RigidFacePart")
        self.contact_model_part = all_model_parts.Get("ContactPart")
        self.mapping_model_part = all_model_parts.Get("MappingPart")

        # Printing variables
        self.DEM_parameters = DEM_parameters
        self.global_variables = []
        self.spheres_and_clusters_variables = []
        self.spheres_and_clusters_local_axis_variables = []
        self.spheres_not_in_cluster_and_clusters_variables = []
        self.spheres_not_in_cluster_and_clusters_local_axis_variables = []
        self.spheres_variables = []
        self.spheres_local_axis_variables = []
        self.fem_boundary_variables = []
        self.clusters_variables = []
        self.rigid_body_variables = []
        self.contact_variables = []
        self.multifilelists = []

        # Reading Post options from DEM_parameters
        self.PostDisplacement = self.DEM_parameters["PostDisplacement"].GetBool()
        self.PostVelocity = self.DEM_parameters["PostVelocity"].GetBool()
        self.PostTotalForces = self.DEM_parameters["PostTotalForces"].GetBool()
        self.PostNonDimensionalVolumeWear = self.DEM_parameters["PostNonDimensionalVolumeWear"].GetBool()
        self.PostAppliedForces = self.DEM_parameters["PostAppliedForces"].GetBool()
        self.PostDampForces = self.DEM_parameters["PostDampForces"].GetBool()
        self.PostRadius = self.DEM_parameters["PostRadius"].GetBool()
        self.PostExportId = self.DEM_parameters["PostExportId"].GetBool()
        self.PostSkinSphere = GetBoolParameterIfItExists(self.DEM_parameters, "PostSkinSphere")
        self.PostAngularVelocity = self.DEM_parameters["PostAngularVelocity"].GetBool()
        self.PostParticleMoment = self.DEM_parameters["PostParticleMoment"].GetBool()
        self.PostEulerAngles = self.DEM_parameters["PostEulerAngles"].GetBool()
        self.PostRollingResistanceMoment = self.DEM_parameters["PostRollingResistanceMoment"].GetBool()
        self.PostLocalContactForce = GetBoolParameterIfItExists(self.DEM_parameters, "PostLocalContactForce")
        self.PostFailureCriterionState = GetBoolParameterIfItExists(self.DEM_parameters, "PostFailureCriterionState")
        self.PostContactFailureId = GetBoolParameterIfItExists(self.DEM_parameters, "PostContactFailureId")
        self.PostContactTau = GetBoolParameterIfItExists(self.DEM_parameters, "PostContactTau")
        self.PostContactSigma = GetBoolParameterIfItExists(self.DEM_parameters, "PostContactSigma")
        self.PostMeanContactArea = GetBoolParameterIfItExists(self.DEM_parameters, "PostMeanContactArea")
        self.PostElasticForces = self.DEM_parameters["PostElasticForces"].GetBool()
        self.PostContactForces = self.DEM_parameters["PostContactForces"].GetBool()
        self.PostRigidElementForces = self.DEM_parameters["PostRigidElementForces"].GetBool()
        self.PostPressure = self.DEM_parameters["PostPressure"].GetBool()
        self.PostTangentialElasticForces = self.DEM_parameters["PostTangentialElasticForces"].GetBool()
        self.PostShearStress = self.DEM_parameters["PostShearStress"].GetBool()
        self.PostNodalArea = self.DEM_parameters["PostNodalArea"].GetBool()
        self.PostTemperature = GetBoolParameterIfItExists(self.DEM_parameters, "PostTemperature")
        self.PostHeatFlux = GetBoolParameterIfItExists(self.DEM_parameters, "PostHeatFlux")
        self.PostNeighbourSize = GetBoolParameterIfItExists(self.DEM_parameters, "PostNeighbourSize")
        self.PostBrokenRatio = GetBoolParameterIfItExists(self.DEM_parameters, "PostBrokenRatio")
        self.PostNormalImpactVelocity = GetBoolParameterIfItExists(self.DEM_parameters, "PostNormalImpactVelocity")
        self.PostTangentialImpactVelocity = GetBoolParameterIfItExists(self.DEM_parameters, "PostTangentialImpactVelocity")
        self.VelTrapGraphExportFreq = self.DEM_parameters["VelTrapGraphExportFreq"].GetDouble()
        if not "PostCharacteristicLength" in self.DEM_parameters.keys():
            self.PostCharacteristicLength = 0
        else:
            self.PostCharacteristicLength = self.DEM_parameters["PostCharacteristicLength"].GetBool()

        #self.PostFaceNormalImpactVelocity = GetBoolParameterIfItExists(self.DEM_parameters, "PostFaceNormalImpactVelocity", 0)
        #self.PostFaceTangentialImpactVelocity = GetBoolParameterIfItExists(self.DEM_parameters, "PostFaceTangentialImpactVelocity", 0)

        if not "PostBoundingBox" in self.DEM_parameters.keys():
            self.PostBoundingBox = 0
        else:
            self.PostBoundingBox = self.DEM_parameters["PostBoundingBox"].GetBool()

        #self.automatic_bounding_box_option = Var_Translator(self.DEM_parameters["AutomaticBoundingBoxOption"].GetBool())
        #self.b_box_minX = self.DEM_parameters["BoundingBoxMinX"].GetDouble()
        #self.b_box_minY = self.DEM_parameters["BoundingBoxMinY"].GetDouble()
        #self.b_box_minZ = self.DEM_parameters["BoundingBoxMinZ"].GetDouble()
        #self.b_box_maxX = self.DEM_parameters["BoundingBoxMaxX"].GetDouble()
        #self.b_box_maxY = self.DEM_parameters["BoundingBoxMaxY"].GetDouble()
        #self.b_box_maxZ = self.DEM_parameters["BoundingBoxMaxZ"].GetDouble()

        self.continuum_element_types = ["SphericContPartDEMElement3D", "CylinderContPartDEMElement2D", "IceContPartDEMElement3D"]

        self.OpenMultiFileLists()

        #Analytic
        if not "PostNormalImpactVelocity" in self.DEM_parameters.keys():
            self.PostNormalImpactVelocity = 0
        else:
            self.PostNormalImpactVelocity = self.DEM_parameters["PostNormalImpactVelocity"].GetBool()

        if not "PostTangentialImpactVelocity" in self.DEM_parameters.keys():
            self.PostTangentialImpactVelocity = 0
        else:
            self.PostTangentialImpactVelocity = self.DEM_parameters["PostTangentialImpactVelocity"].GetBool()

        if not "PostFaceNormalImpactVelocity" in self.DEM_parameters.keys():
            self.PostFaceNormalImpactVelocity = 0
        else:
            self.PostFaceNormalImpactVelocity = self.DEM_parameters["PostFaceNormalImpactVelocity"].GetBool()

        if not "PostFaceTangentialImpactVelocity" in self.DEM_parameters.keys():
            self.PostFaceTangentialImpactVelocity = 0
        else:
            self.PostFaceTangentialImpactVelocity = self.DEM_parameters["PostFaceTangentialImpactVelocity"].GetBool()

        # Ice
        self.sea_settings = self.DEM_parameters["virtual_sea_surface_settings"]

        if self.sea_settings["print_sea_surface"].GetBool():
            self.SeaSurfaceX1 = self.sea_settings["PostVirtualSeaSurfaceX1"].GetDouble()
            self.SeaSurfaceY1 = self.sea_settings["PostVirtualSeaSurfaceY1"].GetDouble()
            self.SeaSurfaceX2 = self.sea_settings["PostVirtualSeaSurfaceX2"].GetDouble()
            self.SeaSurfaceY2 = self.sea_settings["PostVirtualSeaSurfaceY2"].GetDouble()
            self.SeaSurfaceX3 = self.sea_settings["PostVirtualSeaSurfaceX3"].GetDouble()
            self.SeaSurfaceY3 = self.sea_settings["PostVirtualSeaSurfaceY3"].GetDouble()
            self.SeaSurfaceX4 = self.sea_settings["PostVirtualSeaSurfaceX4"].GetDouble()
            self.SeaSurfaceY4 = self.sea_settings["PostVirtualSeaSurfaceY4"].GetDouble()

    def OpenMultiFileLists(self):
        one_level_up_path = os.path.join(self.post_path, "..")
        self.multifiles = (
            #MultifileList(one_level_up_path, self.DEM_parameters["problem_name"].GetString(), 1, "outer"),
            MultifileList(self.post_path, self.DEM_parameters["problem_name"].GetString(), 1, "inner"),
            MultifileList(self.post_path, self.DEM_parameters["problem_name"].GetString(), 2, "inner"),
            MultifileList(self.post_path, self.DEM_parameters["problem_name"].GetString(), 5, "inner"),
            MultifileList(self.post_path, self.DEM_parameters["problem_name"].GetString(), 10, "inner"),
            MultifileList(self.post_path, self.DEM_parameters["problem_name"].GetString(), 20, "inner"),
            MultifileList(self.post_path, self.DEM_parameters["problem_name"].GetString(), 50, "inner"),
        )
        self.SetMultifileLists(self.multifiles)


    def KRATOSprint(self, message):
        Logger.Print(message,label="DEM")
        Logger.Flush()

    @classmethod
    def Flush(self, a):
        a.flush()

    def ShowPrintingResultsOnScreen(self, all_model_parts):
        self.KRATOSprint("*******************  PRINTING RESULTS FOR GID  ***************************")
        self.KRATOSprint("                        (" + str(all_model_parts.Get("SpheresPart").NumberOfElements(0)) + " elements)")
        self.KRATOSprint("                        (" + str(all_model_parts.Get("SpheresPart").NumberOfNodes(0)) + " nodes)")
        self.KRATOSprint("")

    def Initialize(self, DEM_parameters):
        self.AddGlobalVariables()
        self.AddSpheresVariables()
        self.AddSpheresAndClustersVariables()
        self.AddSpheresNotInClusterAndClustersVariables()
        self.AddFEMBoundaryVariables()
        self.AddClusterVariables()
        self.AddRigidBodyVariables()
        self.AddContactVariables()
        self.AddMpiVariables()
        self.Configure(DEM_parameters["problem_name"].GetString(), DEM_parameters["OutputFileType"].GetString(), DEM_parameters["Multifile"].GetString(), DEM_parameters["ContactMeshOption"].GetBool())
        self.SetOutputName(DEM_parameters["problem_name"].GetString())

    @classmethod
    def PushPrintVar(self, variable, name, print_list):
        if Var_Translator(variable):
            print_list.append(name)

    def AddGlobalVariables(self):
        self.PushPrintVar(self.PostDisplacement, DISPLACEMENT, self.global_variables)
        self.PushPrintVar(self.PostVelocity, VELOCITY, self.global_variables)
        self.PushPrintVar(self.PostTotalForces, TOTAL_FORCES, self.global_variables)
        self.PushPrintVar(self.PostAppliedForces, EXTERNAL_APPLIED_FORCE,  self.global_variables)
        self.PushPrintVar(self.PostAppliedForces, EXTERNAL_APPLIED_MOMENT, self.global_variables)
        if self.DEM_parameters["PostAngularVelocity"].GetBool():
            self.PushPrintVar(self.PostAngularVelocity, ANGULAR_VELOCITY, self.global_variables)
        if self.DEM_parameters["PostParticleMoment"].GetBool():
            self.PushPrintVar(self.PostParticleMoment, PARTICLE_MOMENT, self.global_variables)

    def AddSpheresAndClustersVariables(self):  # variables common to spheres and clusters
        self.PushPrintVar(self.PostRigidElementForces,  RIGID_ELEMENT_FORCE,     self.spheres_and_clusters_variables)

    # variables common to spheres and clusters
    def AddSpheresNotInClusterAndClustersVariables(self):
        if self.DEM_parameters["PostEulerAngles"].GetBool():
            self.PushPrintVar(self.PostEulerAngles, EULER_ANGLES, self.spheres_not_in_cluster_and_clusters_local_axis_variables)

    def AddSpheresVariables(self):
        self.PushPrintVar(self.PostDampForces, DAMP_FORCES, self.spheres_variables)
        self.PushPrintVar(self.PostRadius, RADIUS, self.spheres_variables)
        self.PushPrintVar(self.PostExportId, EXPORT_ID, self.spheres_variables)
        self.PushPrintVar(self.PostTemperature, TEMPERATURE, self.spheres_variables)
        self.PushPrintVar(self.PostHeatFlux, HEATFLUX, self.spheres_variables)
        self.PushPrintVar(self.PostNormalImpactVelocity, NORMAL_IMPACT_VELOCITY, self.spheres_variables)
        self.PushPrintVar(self.PostTangentialImpactVelocity, TANGENTIAL_IMPACT_VELOCITY, self.spheres_variables)
        self.PushPrintVar(self.PostFaceNormalImpactVelocity, FACE_NORMAL_IMPACT_VELOCITY, self.spheres_variables)
        self.PushPrintVar(self.PostFaceTangentialImpactVelocity, FACE_TANGENTIAL_IMPACT_VELOCITY, self.spheres_variables)
        #self.PushPrintVar(self.PostLinearImpulse, LINEAR_IMPULSE, self.spheres_variables)
        #self.PushPrintVar(1, DELTA_DISPLACEMENT, self.spheres_variables)  # Debugging
        #self.PushPrintVar(1, PARTICLE_ROTATION_ANGLE, self.spheres_variables)  # Debugging

        if "PostRollingResistanceMoment" in self.DEM_parameters.keys():
            if self.DEM_parameters["RotationOption"].GetBool():
                if self.DEM_parameters["RollingFrictionOption"].GetBool():
                    self.PushPrintVar(self.PostRollingResistanceMoment, ROLLING_RESISTANCE_MOMENT, self.spheres_variables)

        if "PostSkinSphere" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostSkinSphere"].GetBool():
                self.PushPrintVar(self.PostSkinSphere, SKIN_SPHERE, self.spheres_variables)

        if "PostNeighbourSize" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostNeighbourSize"].GetBool():
                self.PushPrintVar(self.PostNeighbourSize, NEIGHBOUR_SIZE, self.spheres_variables)

        if "PostBrokenRatio" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostBrokenRatio"].GetBool():
                self.PushPrintVar(self.PostBrokenRatio, NEIGHBOUR_RATIO, self.spheres_variables)

        # NANO (TODO: must be removed from here.)
        if self.DEM_parameters["ElementType"].GetString() == "SwimmingNanoParticle":
            self.PushPrintVar(self.PostHeatFlux, CATION_CONCENTRATION, self.spheres_variables)

        if "PostStressStrainOption" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostStressStrainOption"].GetBool():
                self.PushPrintVar(1, REPRESENTATIVE_VOLUME, self.spheres_variables)
                self.PushPrintVar(1, DEM_STRESS_TENSOR, self.spheres_variables)

        if "PostReactions" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostReactions"].GetBool():
                self.PushPrintVar(1, FORCE_REACTION, self.spheres_variables)
                self.PushPrintVar(1, MOMENT_REACTION, self.spheres_variables)

        if "PostPoissonRatio" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostPoissonRatio"].GetBool():
                self.PushPrintVar(1, POISSON_VALUE, self.spheres_variables)

    def AddFEMBoundaryVariables(self):
        self.PushPrintVar(self.PostElasticForces, ELASTIC_FORCES, self.fem_boundary_variables)
        self.PushPrintVar(self.PostContactForces, CONTACT_FORCES, self.fem_boundary_variables)
        self.PushPrintVar(self.PostPressure, DEM_PRESSURE, self.fem_boundary_variables)
        self.PushPrintVar(self.PostTangentialElasticForces, TANGENTIAL_ELASTIC_FORCES, self.fem_boundary_variables)
        self.PushPrintVar(self.PostShearStress, SHEAR_STRESS, self.fem_boundary_variables)
        self.PushPrintVar(self.PostNodalArea, DEM_NODAL_AREA, self.fem_boundary_variables)
        if Var_Translator(self.PostNonDimensionalVolumeWear):
            self.PushPrintVar(1, NON_DIMENSIONAL_VOLUME_WEAR, self.fem_boundary_variables)
            self.PushPrintVar(1, IMPACT_WEAR, self.fem_boundary_variables)

    def AddClusterVariables(self):

        if self.PostCharacteristicLength:
            self.PushPrintVar(self.PostCharacteristicLength, CHARACTERISTIC_LENGTH, self.clusters_variables)

        if self.DEM_parameters["PostEulerAngles"].GetBool():
            # JIG: SHOULD BE REMOVED IN THE FUTURE
            self.PushPrintVar(self.PostEulerAngles, ORIENTATION_REAL, self.clusters_variables)
            # JIG: SHOULD BE REMOVED IN THE FUTURE
            self.PushPrintVar(self.PostEulerAngles, ORIENTATION_IMAG, self.clusters_variables)
            #self.PushPrintVar(self.PostEulerAngles, ORIENTATION, self.clusters_variables)

    def AddRigidBodyVariables(self):
        pass
        #self.PushPrintVar(1,                         PARTICLE_MOMENT,              self.rigid_body_variables)
        #self.PushPrintVar(1,                         DELTA_DISPLACEMENT,           self.rigid_body_variables)
        #self.PushPrintVar(1,                         DELTA_ROTATION,               self.rigid_body_variables)
        #self.PushPrintVar(1,                         EXTERNAL_APPLIED_FORCE,       self.rigid_body_variables)
        #self.PushPrintVar(1,                         EXTERNAL_APPLIED_MOMENT,      self.rigid_body_variables)

    def AddContactVariables(self):
        # Contact Elements Variables
        if self.DEM_parameters["ContactMeshOption"].GetBool():
            self.PushPrintVar(self.PostLocalContactForce, LOCAL_CONTACT_FORCE, self.contact_variables)
            if self.DEM_parameters["ElementType"].GetString() in self.continuum_element_types:
                self.PushPrintVar(self.PostFailureCriterionState, FAILURE_CRITERION_STATE, self.contact_variables)
                self.PushPrintVar(self.PostContactFailureId, CONTACT_FAILURE, self.contact_variables)
                self.PushPrintVar(self.PostContactTau, CONTACT_TAU, self.contact_variables)
                self.PushPrintVar(self.PostContactSigma, CONTACT_SIGMA, self.contact_variables)
                self.PushPrintVar(self.PostMeanContactArea, MEAN_CONTACT_AREA, self.contact_variables)

    def AddMpiVariables(self):
        pass

    def Configure(self, problem_name, encoding, file_system, contact_mesh_option):
        self.problem_name = problem_name

        if encoding == "Binary":
            self.encoding = GiDPostMode.GiD_PostBinary
        else:
            self.encoding = GiDPostMode.GiD_PostAscii

        if self.DEM_parameters["Multifile"].GetString() == "multiple_files":
            self.filesystem = MultiFileFlag.MultipleFiles
        else:
            self.filesystem = MultiFileFlag.SingleFile

        self.deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
        self.write_conditions = WriteConditionsFlag.WriteConditions
        self.contact_mesh_option = contact_mesh_option

        problem_name = os.path.join(self.post_path, self.problem_name)
        self.gid_io = GidIO(problem_name,
                            self.encoding,
                            self.filesystem,
                            self.deformed_mesh_flag,
                            self.write_conditions)

        self.post_utility = PostUtilities()

    def SetOutputName(self, name):
        problem_name = os.path.join(self.post_path, self.problem_name)
        self.gid_io.ChangeOutputName(problem_name)

    def SetMultifileLists(self, multifile_list):
        for mfilelist in multifile_list:
            self.multifilelists.append(mfilelist)

        for mfilelist in self.multifilelists:
            mfilelist.file.write("Multiple\n")
            mfilelist.index = 1

    def PrintMultifileLists(self, time, post_path):
        for mfilelist in self.multifilelists:

            if mfilelist.index == mfilelist.step:

                if self.encoding == GiDPostMode.GiD_PostBinary:
                    text_to_print = self.GetMultiFileListName(mfilelist.name) + "_" + "%.12g" % time + ".post.bin\n"
                    if mfilelist.which_folder == "outer":
                        path_of_file = os.path.dirname(mfilelist.file.name)
                        text_to_print = os.path.join(os.path.relpath(
                            post_path, path_of_file), text_to_print)
                    mfilelist.file.write(text_to_print)
                else:
                    text_to_print1 = self.GetMultiFileListName(mfilelist.name) + "_" + "%.12g" % time + ".post.msh\n"
                    text_to_print2 = self.GetMultiFileListName(mfilelist.name) + "_" + "%.12g" % time + ".post.res\n"
                    if mfilelist.which_folder == "outer":
                        path_of_file = os.path.dirname(mfilelist.file.name)
                        text_to_print1 = os.path.join(os.path.relpath(post_path, path_of_file), text_to_print1)
                        text_to_print2 = os.path.join(os.path.relpath(post_path, path_of_file), text_to_print2)
                    mfilelist.file.write(text_to_print1)
                    mfilelist.file.write(text_to_print2)
                self.Flush(mfilelist.file)
                mfilelist.index = 0

            mfilelist.index += 1

    @classmethod
    def GetMultiFileListName(self, name):
        return name

    def CloseMultifiles(self):
        for mfilelist in self.multifilelists:
            mfilelist.file.close()

    def AddModelPartsToMixedModelPart(self):
        self.post_utility.AddModelPartToModelPart(self.mixed_model_part, self.spheres_model_part)
        if self.contact_mesh_option:
            self.post_utility.AddModelPartToModelPart(self.mixed_model_part, self.contact_model_part)
        self.post_utility.AddModelPartToModelPart(self.mixed_model_part, self.rigid_face_model_part)
        self.post_utility.AddModelPartToModelPart(self.mixed_model_part, self.cluster_model_part)
        self.post_utility.AddModelPartToModelPart(self.mixed_spheres_and_clusters_model_part, self.spheres_model_part)
        self.post_utility.AddModelPartToModelPart(self.mixed_spheres_and_clusters_model_part, self.cluster_model_part)

        self.post_utility.AddSpheresNotBelongingToClustersToMixModelPart(self.mixed_spheres_not_in_cluster_and_clusters_model_part, self.spheres_model_part)
        self.post_utility.AddModelPartToModelPart(self.mixed_spheres_not_in_cluster_and_clusters_model_part, self.cluster_model_part)


    def InitializeMesh(self, all_model_parts):
        if self.filesystem == MultiFileFlag.SingleFile:
            self.AddModelPartsToMixedModelPart()
            self.gid_io.InitializeMesh(0.0)
            self.gid_io.WriteMesh(rigid_face_model_part.GetCommunicator().LocalMesh())
            self.gid_io.WriteClusterMesh(cluster_model_part.GetCommunicator().LocalMesh())
            if self.DEM_parameters["ElementType"].GetString() == "CylinderContPartDEMElement2D":
                self.gid_io.WriteCircleMesh(spheres_model_part.GetCommunicator().LocalMesh())
            else:
                self.gid_io.WriteSphereMesh(spheres_model_part.GetCommunicator().LocalMesh())

            if self.contact_mesh_option:
                self.gid_io.WriteMesh(contact_model_part.GetCommunicator().LocalMesh())

            self.gid_io.FinalizeMesh()
            self.gid_io.InitializeResults(0.0, self.mixed_model_part.GetCommunicator().LocalMesh())
            #self.gid_io.InitializeResults(0.0, mixed_spheres_and_clusters_model_part.GetCommunicator().LocalMesh())

    def InitializeResults(self, spheres_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, mapping_model_part, creator_destructor, dem_fem_search, time, bounding_box_time_limits):  # MIQUEL MAPPING

        if self.filesystem == MultiFileFlag.MultipleFiles:
            self.RemoveElementsAndNodes()
            self.AddModelPartsToMixedModelPart()
            self.gid_io.InitializeMesh(time)
            if self.DEM_parameters["ElementType"].GetString() == "CylinderContPartDEMElement2D":
                self.gid_io.WriteCircleMesh(spheres_model_part.GetCommunicator().LocalMesh())
            else:
                self.gid_io.WriteSphereMesh(spheres_model_part.GetCommunicator().LocalMesh())
            if self.contact_mesh_option:
                #We overwrite the Id of the properties 0 not to overlap with other entities that use layer 0 for PRINTING
                contact_model_part.GetProperties(0)[0].Id = 9184
                self.gid_io.WriteMesh(contact_model_part.GetCommunicator().LocalMesh())

            self.gid_io.WriteMesh(rigid_face_model_part.GetCommunicator().LocalMesh())
            self.gid_io.WriteClusterMesh(cluster_model_part.GetCommunicator().LocalMesh())

            #Compute and print the graphical bounding box if active in time
            if self.DEM_parameters["BoundingBoxOption"].GetBool() and (time >= bounding_box_time_limits[0] and time <= bounding_box_time_limits[1]):
                self.ComputeAndPrintBoundingBox(spheres_model_part, rigid_face_model_part, contact_model_part, creator_destructor)

            # Ice. Printing a virtual sea surface
            if self.sea_settings["print_sea_surface"].GetBool():
                self.ComputeAndPrintSeaSurface(spheres_model_part, rigid_face_model_part)

            #self.ComputeAndPrintDEMFEMSearchBinBoundingBox(spheres_model_part, rigid_face_model_part, dem_fem_search)#MSIMSI

            self.gid_io.FinalizeMesh()
            self.gid_io.InitializeResults(time, self.mixed_model_part.GetCommunicator().LocalMesh())
            #self.gid_io.InitializeResults(time, mixed_spheres_and_clusters_model_part.GetCommunicator().LocalMesh())

    def FinalizeMesh(self):
        if self.filesystem == MultiFileFlag.SingleFile:
            self.gid_io.FinalizeResults()

    def FinalizeResults(self):
        if self.filesystem == MultiFileFlag.MultipleFiles:
            self.gid_io.FinalizeResults()

    def PrintingGlobalVariables(self, export_model_part, time):
        for variable in self.global_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)

    def PrintingSpheresAndClustersVariables(self, export_model_part, time):
        for variable in self.spheres_and_clusters_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)
        for variable in self.spheres_and_clusters_local_axis_variables:
            self.gid_io.WriteLocalAxesOnNodes(variable, export_model_part.Nodes, time, 0)

    def PrintingSpheresNotInClusterAndClustersVariables(self, export_model_part, time):
        for variable in self.spheres_not_in_cluster_and_clusters_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)
        for variable in self.spheres_not_in_cluster_and_clusters_local_axis_variables:
            self.gid_io.WriteLocalAxesOnNodes(variable, export_model_part.Nodes, time, 0)

    def PrintingSpheresVariables(self, export_model_part, time):
        for variable in self.spheres_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)
        for variable in self.spheres_local_axis_variables:
            self.gid_io.WriteLocalAxesOnNodes(variable, export_model_part.Nodes, time, 0)

    def PrintingFEMBoundaryVariables(self, export_model_part, time):
        for variable in self.fem_boundary_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)

    def PrintingClusterVariables(self, export_model_part, time):
        for variable in self.clusters_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)

    def PrintingRigidBodyVariables(self, export_model_part, time):
        for variable in self.rigid_body_variables:
            self.gid_io.WriteNodalResults(variable, export_model_part.Nodes, time, 0)

    def PrintingContactElementsVariables(self, export_model_part, time):
        if self.contact_mesh_option:
            for variable in self.contact_variables:
                self.gid_io.PrintOnGaussPoints(variable, export_model_part, time)

    def PrintResults(self, all_model_parts, creator_destructor, dem_fem_search, time, bounding_box_time_limits):

        #TODO: move these definitions to the constructor! (__init__) moved!
        # self.spheres_model_part = spheres_model_part = all_model_parts.Get("SpheresPart")
        # self.cluster_model_part = cluster_model_part = all_model_parts.Get("ClusterPart")
        # self.rigid_face_model_part = rigid_face_model_part = all_model_parts.Get("RigidFacePart")
        # self.contact_model_part = contact_model_part = all_model_parts.Get("ContactPart")
        # self.mapping_model_part = mapping_model_part = all_model_parts.Get("MappingPart")

        if self.filesystem == MultiFileFlag.MultipleFiles:
            self.InitializeResults(self.spheres_model_part,
                                   self.rigid_face_model_part,
                                   self.cluster_model_part,
                                   self.contact_model_part,
                                   self.mapping_model_part,
                                   creator_destructor,
                                   dem_fem_search,
                                   time,
                                   bounding_box_time_limits)

        self.PrintingGlobalVariables(self.mixed_model_part, time)
        self.PrintingSpheresAndClustersVariables(self.mixed_spheres_and_clusters_model_part, time)
        self.PrintingSpheresNotInClusterAndClustersVariables(self.mixed_spheres_not_in_cluster_and_clusters_model_part, time)
        self.PrintingSpheresVariables(self.spheres_model_part, time)
        self.PrintingFEMBoundaryVariables(self.rigid_face_model_part, time)
        self.PrintingRigidBodyVariables(self.rigid_face_model_part, time)
        self.PrintingClusterVariables(self.cluster_model_part, time)
        self.PrintingContactElementsVariables(self.contact_model_part, time)

        self.RemoveElementsAndNodes()

        if self.filesystem == MultiFileFlag.MultipleFiles:
            self.FinalizeResults()

    def RemoveElementsAndNodes(self):
        self.mixed_model_part.Elements.clear()
        self.mixed_model_part.Nodes.clear()
        self.mixed_spheres_and_clusters_model_part.Elements.clear()
        self.mixed_spheres_and_clusters_model_part.Nodes.clear()
        self.mixed_spheres_not_in_cluster_and_clusters_model_part.Elements.clear()
        self.mixed_spheres_not_in_cluster_and_clusters_model_part.Nodes.clear()

    def ComputeAndPrintBoundingBox(self, spheres_model_part, rigid_face_model_part, contact_model_part, creator_destructor):

        if self.PostBoundingBox:
            # Creation of bounding box's model part
            bounding_box_model_part = self.model.CreateModelPart("BoundingBoxPart")

            max_node_Id = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(spheres_model_part)
            max_FEM_node_Id = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(rigid_face_model_part)
            max_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(spheres_model_part)
            max_FEM_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(rigid_face_model_part)
            max_contact_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(contact_model_part)

            if max_FEM_node_Id > max_node_Id:
                max_node_Id = max_FEM_node_Id

            if max_FEM_element_Id > max_element_Id:
                max_element_Id = max_FEM_element_Id

            if max_contact_element_Id > max_element_Id:
                max_element_Id = max_contact_element_Id

            BBMaxX = creator_destructor.GetHighNode()[0]
            BBMaxY = creator_destructor.GetHighNode()[1]
            BBMaxZ = creator_destructor.GetHighNode()[2]
            BBMinX = creator_destructor.GetLowNode()[0]
            BBMinY = creator_destructor.GetLowNode()[1]
            BBMinZ = creator_destructor.GetLowNode()[2]

            self.BuildGraphicalBoundingBox(bounding_box_model_part, max_node_Id, max_element_Id, BBMinX, BBMinY, BBMinZ, BBMaxX, BBMaxY, BBMaxZ)

            self.gid_io.WriteMesh(bounding_box_model_part.GetCommunicator().LocalMesh())

            self.model.DeleteModelPart("BoundingBoxPart")

    def ComputeAndPrintSeaSurface(self, spheres_model_part, rigid_face_model_part):

        # Creation of sea surface model part
        sea_surface_model_part = ModelPart("SeaSurfacePart")

        max_node_Id = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(spheres_model_part)
        max_FEM_node_Id = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(rigid_face_model_part)
        max_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(spheres_model_part)
        max_FEM_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(rigid_face_model_part)

        if max_FEM_node_Id > max_node_Id:
            max_node_Id = max_FEM_node_Id

        if max_FEM_element_Id > max_element_Id:
            max_element_Id = max_FEM_element_Id

        # Z = 0.0 as sea level. We will always assume this value
        node1 = sea_surface_model_part.CreateNewNode(max_node_Id + 9, self.SeaSurfaceX1, self.SeaSurfaceY1, 0.0)
        node2 = sea_surface_model_part.CreateNewNode(max_node_Id + 10, self.SeaSurfaceX2, self.SeaSurfaceY2, 0.0)
        node3 = sea_surface_model_part.CreateNewNode(max_node_Id + 11, self.SeaSurfaceX3, self.SeaSurfaceY3, 0.0)
        node4 = sea_surface_model_part.CreateNewNode(max_node_Id + 12, self.SeaSurfaceX4, self.SeaSurfaceY4, 0.0)

        ''' Properties colours: 0 -> grey,        1 -> dark blue, 2 -> pink,       3 -> light blue,       4 -> dark red,    5 -> light green
                                6 -> light brown, 7 -> red-brown, 8 -> dark brown, 9 -> dark green/blue, 10 -> dark purple'''

        # Sea Surface Element, consisting in a quadrilateral. Property 3 corresponds to a light blue for water
        sea_surface_model_part.CreateNewCondition("RigidFace3D4N", max_element_Id + 1, [node1.Id, node2.Id, node3.Id, node4.Id], Properties(3))

        self.gid_io.WriteMesh(
            sea_surface_model_part.GetCommunicator().LocalMesh())

    def ComputeAndPrintDEMFEMSearchBinBoundingBox(self, spheres_model_part, rigid_face_model_part, dem_fem_search):

        bounding_box_model_part = self.model.CreateModelPart("BoundingBoxPart")

        max_node_Id = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(spheres_model_part)
        max_FEM_node_Id = ParticleCreatorDestructor().FindMaxNodeIdInModelPart(rigid_face_model_part)
        max_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(spheres_model_part)
        max_FEM_element_Id = ParticleCreatorDestructor().FindMaxElementIdInModelPart(rigid_face_model_part)

        if max_FEM_node_Id > max_node_Id:
            max_node_Id = max_FEM_node_Id

        if max_FEM_element_Id > max_element_Id:
            max_element_Id = max_FEM_element_Id

        BBMaxX = dem_fem_search.GetBBHighPoint()[0]
        BBMaxY = dem_fem_search.GetBBHighPoint()[1]
        BBMaxZ = dem_fem_search.GetBBHighPoint()[2]
        BBMinX = dem_fem_search.GetBBLowPoint()[0]
        BBMinY = dem_fem_search.GetBBLowPoint()[1]
        BBMinZ = dem_fem_search.GetBBLowPoint()[2]

        DX = (BBMaxX - BBMinX)
        DY = (BBMaxY - BBMinY)
        DZ = (BBMaxZ - BBMinZ)

        #The cases with 0 thickness in one direction, a 10% of the shortest other two is given to the 0-thickness direction.
        if DX == 0:
            height = min(DY, DZ)
            BBMinX = BBMinX - 0.05 * height
            BBMaxX = BBMaxX + 0.05 * height
        if DY == 0:
            height = min(DX, DZ)
            BBMinY = BBMinY - 0.05 * height
            BBMaxY = BBMaxY + 0.05 * height
        if DZ == 0:
            height = min(DX, DY)
            BBMinZ = BBMinZ - 0.05 * height
            BBMaxZ = BBMaxZ + 0.05 * height

        volume = DX * DY * DZ

        if abs(volume) > 1e21:
            BBMaxX = 0.0
            BBMaxY = 0.0
            BBMaxZ = 0.0
            BBMinX = 0.0
            BBMinY = 0.0
            BBMinZ = 0.0

        self.BuildGraphicalBoundingBox(bounding_box_model_part, max_node_Id, max_element_Id, BBMinX, BBMinY, BBMinZ, BBMaxX, BBMaxY, BBMaxZ)

        self.model.DeleteModelPart("BoundingBoxPart")

        #self.gid_io.WriteMesh(bounding_box_model_part.GetCommunicator().LocalMesh()) #BOUNDING BOX IMPLEMENTATION

    @classmethod
    def BuildGraphicalBoundingBox(self, bounding_box_model_part, max_node_Id, max_element_Id, BBMinX, BBMinY, BBMinZ, BBMaxX, BBMaxY, BBMaxZ):
        # BB Nodes:
        node1 = bounding_box_model_part.CreateNewNode(max_node_Id + 1, BBMinX, BBMinY, BBMinZ)
        node2 = bounding_box_model_part.CreateNewNode(max_node_Id + 2, BBMaxX, BBMinY, BBMinZ)
        node3 = bounding_box_model_part.CreateNewNode(max_node_Id + 3, BBMaxX, BBMaxY, BBMinZ)
        node4 = bounding_box_model_part.CreateNewNode(max_node_Id + 4, BBMinX, BBMaxY, BBMinZ)
        node5 = bounding_box_model_part.CreateNewNode(max_node_Id + 5, BBMinX, BBMinY, BBMaxZ)
        node6 = bounding_box_model_part.CreateNewNode(max_node_Id + 6, BBMaxX, BBMinY, BBMaxZ)
        node7 = bounding_box_model_part.CreateNewNode(max_node_Id + 7, BBMaxX, BBMaxY, BBMaxZ)
        node8 = bounding_box_model_part.CreateNewNode(max_node_Id + 8, BBMinX, BBMaxY, BBMaxZ)

        # BB Elements:
        props = Properties(10000)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 1, [node1.Id, node4.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 2, [node4.Id, node8.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 3, [node8.Id, node5.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 4, [node5.Id, node1.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 5, [node1.Id, node2.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 6, [node3.Id, node4.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 7, [node7.Id, node8.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 8, [node5.Id, node6.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 9, [node6.Id, node2.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 10, [node2.Id, node3.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 11, [node3.Id, node7.Id], props)
        bounding_box_model_part.CreateNewCondition("RigidEdge3D", max_element_Id + 12, [node7.Id, node6.Id], props)



class ParallelUtils(object):

    def __init__(self):
        pass

    def Repart(self, spheres_model_part):
        pass

    def CalculateModelNewIds(self, spheres_model_part):
        pass

    def PerformInitialPartition(self, model_part):
        pass

    @classmethod
    def SetCommunicator(self, spheres_model_part, model_part_io_spheres, spheres_mp_filename):
        MPICommSetup = 0
        return [model_part_io_spheres, spheres_model_part, MPICommSetup]

    @classmethod
    def GetSearchStrategy(self, solver, model_part):
        return solver.search_strategy
