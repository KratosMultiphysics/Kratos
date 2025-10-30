import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
import beam_mapper_test_case
import mapper_test_case
from math import cos
import math
import os
import numpy as np
import json

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)

class BladeMappingTests(beam_mapper_test_case.BeamMapperTestCase):
    '''This class contains basic tests for mapping on real geometries
    In this case it is a remodeled NREL Phase VI wind turbine blade
    It also serves as a showcase on how to use the Mapper in FSI
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters):
        structure_mdpa_file_name = "blade_line"
        fluid_mdpa_file_name     = "blade_quad"
        super(BladeMappingTests, cls).setUpModelParts(structure_mdpa_file_name, fluid_mdpa_file_name)

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

        mapper_parameters_p = mapper_parameters.Clone()
        mapper_parameters_s = mapper_parameters.Clone()

        mapper_parameters_p.AddEmptyValue("interface_submodel_part_destination").SetString("pressure_side_quad")
        mapper_parameters_s.AddEmptyValue("interface_submodel_part_destination").SetString("suction_side_quad")

        cls.model_part_structure = cls.model_part_origin
        cls.model_part_fluid = cls.model_part_destination

        cls.mapper_pressure_side = KratosMapping.MapperFactory.CreateMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_p)
        cls.mapper_suction_side = KratosMapping.MapperFactory.CreateMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_s)

        cls.print_output = False # this can be overridden in derived classes to print the output

    def test_map_displacements(self):
        SetDisplacements(self.model_part_structure)
        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.DISPLACEMENT, "Blade_" + self.mapper_type + "_Structure_prescr_disp")
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.ROTATION, "Blade_" + self.mapper_type + "_Structure_prescr_rot")

        self.mapper_pressure_side.Map(KM.DISPLACEMENT, KM.ROTATION, KM.MESH_DISPLACEMENT)
        self.mapper_suction_side.Map(KM.DISPLACEMENT, KM.ROTATION, KM.MESH_DISPLACEMENT)
        
        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.MESH_DISPLACEMENT, "Blade_" + self.mapper_type + "_Fluid_mapped_disp")

        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_fluid, KM.MESH_DISPLACEMENT, GetFilePath(self.__GetFileName("blade_map_disp")))

    def test_map_forces(self):
        SetReactions(self.model_part_fluid)
        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.REACTION, "Blade_" + self.mapper_type + "_Fluid_prescr_force")

        self.mapper_pressure_side.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION)
        self.mapper_suction_side.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION)

        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.FORCE, "Blade_" + self.mapper_type + "_Structure_mapped_force")
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.MOMENT, "Blade_" + self.mapper_type + "_Structure_mapped_moment")

        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_structure, KM.FORCE, GetFilePath(self.__GetFileName("blade_map_force")))
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_structure, KM.MOMENT, GetFilePath(self.__GetFileName("blade_map_moment")))

    def __GetFileName(self, file_appendix):
        return os.path.join("result_files", self.mapper_type, file_appendix)

def SetDisplacements(model_part_structure):
    for node in model_part_structure.Nodes:
        lenght_beam = 4.521
        alfa = 0.000020929 
        beta = 0.0
        r = lenght_beam / alfa

        theta_X = - node.Z / r
        theta_Y = 0.0
        theta_Z = (beta * node.Z) / lenght_beam

        e_x = np.array([1.0, 0.0, 0.0])
        e_y = np.array([0.0, 1.0, 0.0])
        e_z = np.array([0.0, 0.0, 1.0])

        Rx = CalculateRotationMatrixWithAngle(e_x, theta_X)
        Ry = CalculateRotationMatrixWithAngle(e_y, theta_Y)
        Rz = CalculateRotationMatrixWithAngle(e_z, theta_Z)

        R_temp = np.dot(Ry, Rz)
        R = np.dot(Rx, R_temp)

        ROTATION = GetRotationVector(R)

        node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, r - r*math.cos(-theta_X))
        node.SetSolutionStepValue(KM.DISPLACEMENT_Z, r * math.sin(-theta_X) - node.Z )
        node.SetSolutionStepValue(KM.ROTATION_X, ROTATION[0] )
        node.SetSolutionStepValue(KM.ROTATION_Y, ROTATION[1] )
        node.SetSolutionStepValue(KM.ROTATION_Z, ROTATION[2] )

def SetReactions(model_part_fluid):
    for node in model_part_fluid.Nodes:
        react_x = cos(node.X*5)
        react_y = 0.0
        react_z = cos(node.Z*2)
        node.SetSolutionStepValue(KM.REACTION, KM.Vector([react_x, react_y, react_z]))

def CalculateRotationMatrixWithAngle(Axis, Angle):
    rotation_matrix = np.zeros((3, 3))

    rotation_matrix[0][0] = math.cos( Angle ) + Axis[0]**2 * (1 - math.cos( Angle ))
    rotation_matrix[0][1] = Axis[0] * Axis[1] * (1 - math.cos( Angle )) - Axis[2] * math.sin( Angle )
    rotation_matrix[0][2] = Axis[0] * Axis[2] * (1 - math.cos( Angle )) + Axis[1] * math.sin( Angle )

    rotation_matrix[1][0] = Axis[0] * Axis[1] * (1 - math.cos( Angle )) + Axis[2] * math.sin( Angle )
    rotation_matrix[1][1] = math.cos( Angle ) + Axis[1]**2 * (1 - math.cos( Angle ))
    rotation_matrix[1][2] = Axis[1] * Axis[2] * (1 - math.cos( Angle )) - Axis[0] * math.sin( Angle )

    rotation_matrix[2][0] = Axis[0] * Axis[2] * (1 - math.cos( Angle )) - Axis[1] * math.sin( Angle )
    rotation_matrix[2][1] = Axis[1] * Axis[2] * (1 - math.cos( Angle )) + Axis[0] * math.sin( Angle )
    rotation_matrix[2][2] = math.cos( Angle ) + Axis[2]**2 * (1 - math.cos( Angle ))

    return rotation_matrix

def GetRotationVector(rotation_matrix):
    # see Non-linear Modeling and Analysis of Solids and Structures (Steen Krenk 2009) P52
    rotation_vector = np.array([0.0, 0.0, 0.0])
    angle = rotation_matrix[0][0] + rotation_matrix[1][1] + rotation_matrix[2][2] - 1.0

    angle = angle/2.0
    if (angle > 1.0):
        angle = 1.0
    elif (angle < -1.0):
        angle = -1.0

    angle = math.acos(angle) # between 0 and pi

    EPS = 1E-6
    M_PI = math.pi
    if (angle < EPS):
        rotation_vector[0] = 0.0
        rotation_vector[1] = 0.0
        rotation_vector[2] = 0.0
        return rotation_vector
    elif ((M_PI - angle) < EPS):
        product11 = (rotation_matrix[0][0] + 1.0) / 2.0
        product22 = (rotation_matrix[1][1] + 1.0) / 2.0
        product33 = (rotation_matrix[2][2] + 1.0) / 2.0
        product12 = (rotation_matrix[0][1] + 1.0) / 2.0
        product23 = (rotation_matrix[1][2] + 1.0) / 2.0
        product13 = (rotation_matrix[0][2] + 1.0) / 2.0
        tmp1 = math.sqrt(product11)
        tmp2 = math.sqrt(product22)
        tmp3 = math.sqrt(product33)

        rotation_vector[0] = tmp1
        rotation_vector[1] = tmp2
        rotation_vector[2] = tmp3
        tmp12 = rotation_vector[0] * rotation_vector[1]
        tmp13 = rotation_vector[0] * rotation_vector[2]
        tmp23 = rotation_vector[1] * rotation_vector[2]
        if (math.fabs(tmp12) < EPS or math.fabs(tmp12 - product12) < math.fabs(tmp12 + product12)):
            if (math.fabs(tmp13) < EPS or math.fabs(tmp13 - product13) < math.fabs(tmp13 + product13)):
                if (math.fabs(tmp23) < EPS or math.fabs(tmp23 - product23) < math.fabs(tmp23 + product23)):
                    rotation_vector[0] *= M_PI
                    rotation_vector[1] *= M_PI
                    rotation_vector[2] *= M_PI
                    return rotation_vector

        rotation_vector[0] =  tmp1
        rotation_vector[1] = -tmp2
        rotation_vector[2] = -tmp3
        tmp12 = rotation_vector[0] * rotation_vector[1]
        tmp13 = rotation_vector[0] * rotation_vector[2]
        tmp23 = rotation_vector[1] * rotation_vector[2]
        if (math.fabs(tmp12) < EPS or math.fabs(tmp12 - product12) < math.fabs(tmp12 + product12)):
            if (math.fabs(tmp13) < EPS or math.fabs(tmp13 - product13) < math.fabs(tmp13 + product13)):
                if (math.fabs(tmp23) < EPS or math.fabs(tmp23 - product23) < math.fabs(tmp23 + product23)):
                    rotation_vector[0] *= M_PI
                    rotation_vector[1] *= M_PI
                    rotation_vector[2] *= M_PI
                    return rotation_vector

        rotation_vector[0] = -tmp1
        rotation_vector[1] =  tmp2
        rotation_vector[2] = -tmp3
        tmp12 = rotation_vector[0] * rotation_vector[1]
        tmp13 = rotation_vector[0] * rotation_vector[2]
        tmp23 = rotation_vector[1] * rotation_vector[2]
        if (math.fabs(tmp12) < EPS or math.fabs(tmp12 - product12) < math.fabs(tmp12 + product12)):
            if (math.fabs(tmp13) < EPS or math.fabs(tmp13 - product13) < math.fabs(tmp13 + product13)):
                if (math.fabs(tmp23) < EPS or math.fabs(tmp23 - product23) < math.fabs(tmp23 + product23)):
                    rotation_vector[0] *= M_PI
                    rotation_vector[1] *= M_PI
                    rotation_vector[2] *= M_PI
                    return rotation_vector

        rotation_vector[0] = -tmp1
        rotation_vector[1] = -tmp2
        rotation_vector[2] =  tmp3
        tmp12 = rotation_vector[0] * rotation_vector[1]
        tmp13 = rotation_vector[0] * rotation_vector[2]
        tmp23 = rotation_vector[1] * rotation_vector[2]
        if (math.fabs(tmp12) < EPS or math.fabs(tmp12 - product12) < math.fabs(tmp12 + product12)):
            if (math.fabs(tmp13) < EPS or math.fabs(tmp13 - product13) < math.fabs(tmp13 + product13)):
                if (math.fabs(tmp23) < EPS or math.fabs(tmp23 - product23) < math.fabs(tmp23 + product23)):
                    rotation_vector[0] *= M_PI
                    rotation_vector[1] *= M_PI
                    rotation_vector[2] *= M_PI
                    return rotation_vector
        assert(False)
    
    tmp = angle / 2.0 / math.sin(angle)
    rotation_vector[0] = -(rotation_matrix[1][2] - rotation_matrix[2][1]) * tmp
    rotation_vector[1] =  (rotation_matrix[0][2] - rotation_matrix[2][0]) * tmp
    rotation_vector[2] = -(rotation_matrix[0][1] - rotation_matrix[1][0]) * tmp

    return rotation_vector