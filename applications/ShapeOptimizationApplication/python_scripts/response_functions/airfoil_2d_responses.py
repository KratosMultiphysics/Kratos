import KratosMultiphysics
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import math

def DotProduct(A,B):
    result = 0
    for i,j in zip(A,B):
        result += i*j
    return result

class AngleOfAttackResponseFunction(ResponseFunctionInterface):
    def __init__(self, response_id, response_settings, model):
        self.model = model
        if not response_settings.Has("trailing_edge_model_part")  or not response_settings.Has("leading_edge_model_part"):
            raise(Exception("Please define a model part with the trailing edge node and a model part with the leading edge node"))

        self.trailing_edge_model_part_name = response_settings["trailing_edge_model_part"].GetString()
        self.leading_edge_sub_model_part_name = response_settings["leading_edge_model_part"].GetString()

    def Initialize(self):
        self.te_model_part = self.model[self.trailing_edge_model_part_name]
        self.le_model_part = self.model[self.leading_edge_sub_model_part_name]
        self.main_model_part = self.te_model_part.GetRootModelPart()
        self.reference_direction = self.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.FREE_STREAM_VELOCITY_DIRECTION)

        for node in self.te_model_part.Nodes:
            self.te_node = node
            break
        for node in self.le_model_part.Nodes:
            self.le_node = node
            break

    def CalculateValue(self):
        self.aoa = self._CalculateAOA(self.te_node.X,self.te_node.Y, self.le_node.X, self.le_node.Y)

    def CalculateGradient(self):
        x_diff = self.le_node.X - self.te_node.X
        y_diff = self.le_node.Y - self.te_node.Y

        cshape_sensitivity_0 = x_diff
        cshape_sensitivity_1 = y_diff
        cshape_sensitivity_2 = cshape_sensitivity_0**2 + cshape_sensitivity_1**2
        cshape_sensitivity_3 = 1/cshape_sensitivity_2
        cshape_sensitivity_4 = cshape_sensitivity_0*self.reference_direction[0] + cshape_sensitivity_1*self.reference_direction[1]
        cshape_sensitivity_5 = cshape_sensitivity_3*cshape_sensitivity_4
        cshape_sensitivity_6 = self.reference_direction[0]**2 + self.reference_direction[1]**2
        cshape_sensitivity_7 = 1/(math.sqrt(cshape_sensitivity_2)*math.sqrt(cshape_sensitivity_6)*math.sqrt(-cshape_sensitivity_3*cshape_sensitivity_4**2/cshape_sensitivity_6 + 1)+1e-9)

        self.te_x_gradient = -cshape_sensitivity_7*(-cshape_sensitivity_0*cshape_sensitivity_5 + self.reference_direction[0])
        self.te_y_gradient = -cshape_sensitivity_7*(-cshape_sensitivity_1*cshape_sensitivity_5 + self.reference_direction[1])
        self.le_x_gradient = cshape_sensitivity_7*(-cshape_sensitivity_0*cshape_sensitivity_5 + self.reference_direction[0])
        self.le_y_gradient = cshape_sensitivity_7*(-cshape_sensitivity_1*cshape_sensitivity_5 + self.reference_direction[1])

    def GetValue(self):
        return self.aoa

    def GetNodalGradient(self, variable):
        zero_vector = KratosMultiphysics.Vector(3, 0.0)
        gradient ={node.Id : zero_vector for node in self.main_model_part.Nodes}

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.te_x_gradient
        shape_gradient[1] = self.te_y_gradient
        gradient[self.te_node.Id] = shape_gradient

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.le_x_gradient
        shape_gradient[1] = self.le_y_gradient
        gradient[self.le_node.Id] = shape_gradient

        return gradient

    def _CalculateAOA(self, te_x, te_y, le_x, le_y):
        ref_dir_norm = KratosMultiphysics.Vector(self.reference_direction).norm_2()
        chord_vector = KratosMultiphysics.Vector(3, 0.0)
        chord_vector[0] = te_x-le_x
        chord_vector[1] = te_y-le_y
        chord_norm = chord_vector.norm_2()
        aoa = math.acos(DotProduct(self.reference_direction, chord_vector)/(ref_dir_norm*chord_norm))
        aoa_deg = math.degrees(aoa) # Angle in degrees!!!
        return aoa_deg

class ChordLengthResponseFunction(ResponseFunctionInterface):


    def __init__(self, response_id, response_settings, model):
        self.model = model

        if not response_settings.Has("trailing_edge_model_part")  or not response_settings.Has("leading_edge_model_part"):
            raise(Exception("Please define a model part with the trailing edge node and a model part with the leading edge node"))

        self.trailing_edge_model_part_name = response_settings["trailing_edge_model_part"].GetString()
        self.leading_edge_sub_model_part_name = response_settings["leading_edge_model_part"].GetString()

    def Initialize(self):
        self.te_model_part = self.model[self.trailing_edge_model_part_name]
        self.le_model_part = self.model[self.leading_edge_sub_model_part_name]
        self.main_model_part = self.te_model_part.GetRootModelPart()

        for node in self.te_model_part.Nodes:
            self.te_node = node
            break
        for node in self.le_model_part.Nodes:
            self.le_node = node
            break

    def _ComputeChord(self, te_x, te_y, le_x, le_y):
        chord = math.sqrt((te_x-le_x)**2+(te_y-le_y)**2)

        return chord

    def CalculateValue(self):
        self.chord = self._ComputeChord(self.te_node.X,self.te_node.Y, self.le_node.X,self.le_node.Y)

    def CalculateGradient(self):
        x_diff = self.le_node.X - self.te_node.X
        y_diff = self.le_node.Y - self.te_node.Y

        self.te_x_gradient = -x_diff/math.sqrt(x_diff**2 + y_diff**2)
        self.te_y_gradient = -y_diff/math.sqrt(x_diff**2 + y_diff**2)
        self.le_x_gradient = x_diff/math.sqrt(x_diff**2 + y_diff**2)
        self.le_y_gradient = y_diff/math.sqrt(x_diff**2 + y_diff**2)

    def GetValue(self):
        return self.chord

    def GetNodalGradient(self, variable):
        zero_vector = KratosMultiphysics.Vector(3, 0.0)
        gradient ={node.Id : zero_vector for node in self.main_model_part.Nodes}

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.te_x_gradient
        shape_gradient[1] = self.te_y_gradient
        gradient[self.te_node.Id] = shape_gradient

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.le_x_gradient
        shape_gradient[1] = self.le_y_gradient
        gradient[self.le_node.Id] = shape_gradient

        return gradient



class PerimeterResponseFunction(ResponseFunctionInterface):

    def __init__(self, response_id, response_settings, model):
        self.model = model
        self.model_part_name = response_settings["model_part_name"].GetString()
        self.step_size = response_settings["step_size"].GetDouble() if response_settings.Has("step_size") else 1e-6

    def Initialize(self):
        self.main_model_part = self.model[self.model_part_name]

    def _ComputePerimeter(self,  model_part):
        return KSO.GeometryUtilities(model_part).CalculateLength(model_part.Conditions)

    def CalculateValue(self):
        pass

    def CalculateGradient(self):
        pass

    def GetValue(self):

        return self._ComputePerimeter(self.main_model_part)

    def GetNodalGradient(self, variable):
        gradient = {}
        initial_perimeter =  self._ComputePerimeter(self.main_model_part)

        for node in self.main_model_part.Nodes:
            shape_gradient = KratosMultiphysics.Vector(3, 0.0)

            node.X += self.step_size
            current_perimeter =  self._ComputePerimeter(self.main_model_part)
            shape_gradient[0] = (current_perimeter-initial_perimeter)/self.step_size
            node.X -= self.step_size

            node.Y += self.step_size
            current_perimeter =  self._ComputePerimeter(self.main_model_part)
            shape_gradient[1] = (current_perimeter-initial_perimeter)/self.step_size
            node.Y -= self.step_size

            gradient[node.Id] = shape_gradient

        return gradient


#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

class WingAngleOfAttackResponseFunction(ResponseFunctionInterface):

    def __init__(self, response_id, response_settings, model):
        self.model = model

        if not response_settings.Has("trailing_edge_model_part") or not response_settings.Has("leading_edge_model_part"):
            raise Exception("Please define 'trailing_edge_model_part' and 'leading_edge_model_part' in the settings.")

        self.trailing_edge_model_part_name = response_settings["trailing_edge_model_part"].GetString()
        self.leading_edge_model_part_name = response_settings["leading_edge_model_part"].GetString()

    def Initialize(self):
        self.te_model_part = self.model[self.trailing_edge_model_part_name]
        self.le_model_part = self.model[self.leading_edge_model_part_name]
        self.main_model_part = self.te_model_part.GetRootModelPart()

        ref_dir = self.main_model_part.ProcessInfo.GetValue(
            KratosMultiphysics.CompressiblePotentialFlowApplication.FREE_STREAM_VELOCITY_DIRECTION)
        ref_norm = (ref_dir[0]**2 + ref_dir[1]**2 + ref_dir[2]**2)**0.5
        self.reference_direction = KratosMultiphysics.Vector([ref_dir[0]/(ref_norm + 1e-12),
                                           ref_dir[1]/(ref_norm + 1e-12),
                                           ref_dir[2]/(ref_norm + 1e-12)])

        self.le_centroid = self._ComputeCentroid(self.le_model_part)
        self.te_centroid = self._ComputeCentroid(self.te_model_part)

    def _ComputeCentroid(self, model_part):
        n = len(model_part.Nodes)
        if n == 0:
            raise Exception(f"Model part {model_part.Name} has no nodes.")
        cx = sum(node.X for node in model_part.Nodes) / n
        cy = sum(node.Y for node in model_part.Nodes) / n
        cz = sum(node.Z for node in model_part.Nodes) / n
        return KratosMultiphysics.Vector([cx, cy, cz])

    def CalculateValue(self):
        self.aoa = self._ComputeAngleOfAttackProjectedXZ()
        return self.aoa

    def _ComputeAngleOfAttackProjectedXZ(self):
        chord_vec = KratosMultiphysics.Vector([self.te_centroid[0] - self.le_centroid[0],
                            0.0,
                            self.te_centroid[2] - self.le_centroid[2]])

        chord_norm = (chord_vec[0]**2 + chord_vec[2]**2)**0.5 + 1e-12
        chord_unit = KratosMultiphysics.Vector([chord_vec[0]/chord_norm, 0.0, chord_vec[2]/chord_norm])

        ref_vec = KratosMultiphysics.Vector([self.reference_direction[0], 0.0, self.reference_direction[2]])
        ref_norm = (ref_vec[0]**2 + ref_vec[2]**2)**0.5 + 1e-12
        ref_unit = KratosMultiphysics.Vector([ref_vec[0]/ref_norm, 0.0, ref_vec[2]/ref_norm])

        dot_prod = chord_unit[0]*ref_unit[0] + chord_unit[2]*ref_unit[2]
        dot_prod = max(min(dot_prod, 1.0), -1.0)  
        angle = math.acos(dot_prod)

        cross_z = chord_unit[0]*ref_unit[2] - chord_unit[2]*ref_unit[0]
        if cross_z > 0.0:
            angle *= -1.0

        aoa_deg = math.degrees(angle) # Angle in degrees!!!
        return aoa_deg

    def CalculateGradient(self):
        chord_vec = KratosMultiphysics.Vector([self.te_centroid[0] - self.le_centroid[0],
                            0.0,
                            self.te_centroid[2] - self.le_centroid[2]])
        chord_norm = (chord_vec[0]**2 + chord_vec[2]**2)**0.5 + 1e-12
        chord_unit = KratosMultiphysics.Vector([chord_vec[0]/chord_norm, 0.0, chord_vec[2]/chord_norm])

        ref_vec = KratosMultiphysics.Vector([self.reference_direction[0], 0.0, self.reference_direction[2]])
        ref_norm = (ref_vec[0]**2 + ref_vec[2]**2)**0.5 + 1e-12
        ref_unit = KratosMultiphysics.Vector([ref_vec[0]/ref_norm, 0.0, ref_vec[2]/ref_norm])

        dtheta_dchord = KratosMultiphysics.Vector([
            -(ref_unit[2]) / (chord_norm),
            0.0,
            ref_unit[0] / (chord_norm)
        ])

        self.dAngle_dx_te = dtheta_dchord
        self.dAngle_dx_le = KratosMultiphysics.Vector([-dtheta_dchord[0], 0.0, -dtheta_dchord[2]])

    def GetValue(self):
        return self.aoa

    def GetNodalGradient(self, variable):
        zero_vector = KratosMultiphysics.Vector(3, 0.0)
        gradient = {node.Id : zero_vector for node in self.main_model_part.Nodes}

        le_nodes = list(self.le_model_part.Nodes)
        te_nodes = list(self.te_model_part.Nodes)

        n_le = len(le_nodes)
        n_te = len(te_nodes)

        for node in te_nodes:
            gradient[node.Id] = self.dAngle_dx_te / n_te
        for node in le_nodes:
            gradient[node.Id] = self.dAngle_dx_le / n_le
        return gradient
    
class WingChordLengthResponseFunction(ResponseFunctionInterface):

    def __init__(self, response_id, response_settings, model):
        self.model = model

        if not response_settings.Has("trailing_edge_model_part") or not response_settings.Has("leading_edge_model_part"):
            raise Exception("Please define 'trailing_edge_model_part' and 'leading_edge_model_part' in the settings.")

        self.trailing_edge_model_part_name = response_settings["trailing_edge_model_part"].GetString()
        self.leading_edge_model_part_name = response_settings["leading_edge_model_part"].GetString()

    def Initialize(self):
        self.te_model_part = self.model[self.trailing_edge_model_part_name]
        self.le_model_part = self.model[self.leading_edge_model_part_name]
        self.main_model_part = self.te_model_part.GetRootModelPart()

        self.le_centroid = self._ComputeCentroid(self.le_model_part)
        self.te_centroid = self._ComputeCentroid(self.te_model_part)

    def _ComputeCentroid(self, model_part):
        n = len(model_part.Nodes)
        if n == 0:
            raise Exception(f"Model part {model_part.Name} has no nodes.")
        cx = sum(node.X for node in model_part.Nodes) / n
        cy = sum(node.Y for node in model_part.Nodes) / n
        cz = sum(node.Z for node in model_part.Nodes) / n
        return KratosMultiphysics.Vector([cx, cy, cz])

    def CalculateValue(self):
        self.chord = self._ComputeChordProjectedXZ()
        return self.chord

    def _ComputeChordProjectedXZ(self):
        chord_vec = KratosMultiphysics.Vector([self.te_centroid[0] - self.le_centroid[0],
                            0.0,
                            self.te_centroid[2] - self.le_centroid[2]])
        chord_length = (chord_vec[0]**2 + chord_vec[2]**2)**0.5
        return chord_length

    def CalculateGradient(self):
        chord_proj = KratosMultiphysics.Vector([self.te_centroid[0] - self.le_centroid[0],
                             0.0,
                             self.te_centroid[2] - self.le_centroid[2]])
        chord_norm = (chord_proj[0]**2 + chord_proj[2]**2)**0.5 + 1e-12

        self.dC_dx_te = KratosMultiphysics.Vector([chord_proj[0] / chord_norm, 0.0, chord_proj[2] / chord_norm])
        self.dC_dx_le = KratosMultiphysics.Vector([-chord_proj[0] / chord_norm, 0.0, -chord_proj[2] / chord_norm])

    def GetValue(self):
        return self.chord

    def GetNodalGradient(self, variable):
        zero_vector = KratosMultiphysics.Vector(3, 0.0)
        gradient = {node.Id : zero_vector for node in self.main_model_part.Nodes}

        le_nodes = list(self.le_model_part.Nodes)
        te_nodes = list(self.te_model_part.Nodes)

        n_le = len(le_nodes)
        n_te = len(te_nodes)

        for node in te_nodes:
            gradient[node.Id] = self.dC_dx_te / n_te
        for node in le_nodes:
            gradient[node.Id] = self.dC_dx_le / n_le
        return gradient
    
class WingPerimeterResponseFunction(ResponseFunctionInterface):

    def __init__(self, response_id, response_settings, model):
        self.model = model
        self.model_part_name = response_settings["model_part_name"].GetString()
        self.step_size = response_settings["step_size"].GetDouble() if response_settings.Has("step_size") else 1e-6

    def Initialize(self):
        self.main_model_part = self.model[self.model_part_name]

    def _ComputePerimeter(self,  model_part):
        return KSO.GeometryUtilities(model_part).CalculateArea(model_part.Conditions)/2

    def CalculateValue(self):
        pass

    def CalculateGradient(self):
        pass

    def GetValue(self):

        return self._ComputePerimeter(self.main_model_part)

    def GetNodalGradient(self, variable):
        gradient = {}
        initial_perimeter =  self._ComputePerimeter(self.main_model_part)

        for node in self.main_model_part.Nodes:
            shape_gradient = KratosMultiphysics.Vector(3, 0.0)

            node.X += self.step_size
            current_perimeter =  self._ComputePerimeter(self.main_model_part)
            shape_gradient[0] = (current_perimeter-initial_perimeter)/self.step_size
            node.X -= self.step_size

            node.Y += self.step_size
            current_perimeter =  self._ComputePerimeter(self.main_model_part)
            shape_gradient[1] = (current_perimeter-initial_perimeter)/self.step_size
            node.Y -= self.step_size

            node.Z += self.step_size
            current_perimeter =  self._ComputePerimeter(self.main_model_part)
            shape_gradient[2] = (current_perimeter-initial_perimeter)/self.step_size
            node.Z -= self.step_size

            gradient[node.Id] = shape_gradient

        return gradient