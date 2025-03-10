from typing import Optional
import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator


def Factory(model: KM.Model, parameters: KM.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"SubdivisionSurfaceControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SubdivisionSurfaceControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return SubdivisionSurfaceControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class SubdivisionSurfaceControl(Control):
    """
    CAD-based shape control using subdivision surfaces. 
    The kind of subdivision scheme defines the shape of the limit surface.

    This allows the optimization of open and closed shell geometries.
    """
    def __init__(self, name: str, model: KM.Model, parameters: KM.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        self.parameters = parameters
        self.model = model
        self.optimization_problem = optimization_problem
        self.use_controlled_input_mesh = False      # name suggestion: self.use_discrete_limit = False

        default_settings = KM.Parameters("""
        {
            "control_polygon_model_part_name"   : "control_polygon",
            "controlled_model_part_names"       : [],
            "subdivision_type"                  : "catmull_clark",
            "fix_boundary"                      : true
        }
        """)

        self.parameters.ValidateAndAssignDefaults(default_settings)

        self.subdivision_type = self.parameters["subdivision_type"].GetString()

        self.fix_boundary = self.parameters["fix_boundary"].GetBool()

        self.supported_subdivision_types = ["catmull_clark"]    # only catmull_clark for now, but possibility to add others such as Loop scheme
        if self.subdivision_type not in self.supported_subdivision_types:
            raise RuntimeError(f"The specified subdivision type \"{self.subdivision_type}\" is not supported. Supported variables are:\n\t" + "\n\t".join([k for k in self.supported_subdivision_types]))

        control_polygon_model_part_name = parameters["control_polygon_model_part_name"].GetString()
        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            KM.Logger.PrintWarning(f"No controlled model parts are provided for SubdivisionSurfaceControl. Limit surface of the control polygon will be used. [ control name = \"{self.GetName()}\"]")
        else:
            self.use_controlled_input_mesh = True
            self.controlled_model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"controlled_{self.GetName()}", controlled_model_names_parts, False)
            # self.model_part: 'Optional[KM.ModelPart]' = None    # needed?

        self.control_polygon_model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_polygon_{self.GetName()}", [control_polygon_model_part_name], False)

        self.sds_mapper = 0
        self.is_initialized = False

    @time_decorator(methodName="GetName")
    def Initialize(self) -> None:

        self.control_polygon_model_part = self.control_polygon_model_part_operation.GetModelPart()
        if self.use_controlled_input_mesh:
            self.controlled_model_part = self.controlled_model_part_operation.GetModelPart()
        # else create controlled model part by subdivision

        self.control_field = KM.Expression.NodalExpression(self.control_polygon_model_part)

        control_polygon = np.zeros((self.control_polygon_model_part.NumberOfNodes(), 3))
        for idx, node in enumerate(self.control_polygon_model_part.Nodes):
            control_polygon[idx,:] = np.array([node.X, node.Y, node.Z])
        
        KM.Expression.CArrayExpressionIO.Read(
            self.control_field, control_polygon
        )
        
        # only need to do it once, updating the control variables should suffice
        # create the mapping relation between fem mesh and control polygon
        # mapping_relation = KM.Expression.NodalExpression(self.controlled_model_part)
        mapping_relation_array = np.array(
            KOA.ControlUtils.SubdivisionSurfaceProjectionUtils.CalculateMappingRelation(
                # mapping_relation,
                self.control_polygon_model_part,
                self.controlled_model_part,
                self.subdivision_type,
                self.fix_boundary
            )
        )

        num_fe_nodes = self.controlled_model_part.NumberOfNodes()
        num_control_points = self.control_polygon_model_part.NumberOfNodes()
        mapping_relation_data = np.zeros((num_fe_nodes, num_control_points))
        for i in range(num_control_points):
            mapping_relation_data[:,i] = mapping_relation_array[i::num_control_points]

        self.mapping_relation = KM.Expression.NodalExpression(self.control_polygon_model_part)
        KM.Expression.LiteralExpressionIO.SetData(self.mapping_relation, mapping_relation_data)

        self.forward_matrix = KM.Matrix(mapping_relation_data)
        self.backward_matrix = KM.Matrix()

        # inverse mapping in case of matrix: M_inv_map = M.transpose()
        # backward mapping with pseudoinverse is not consistent
        # self.inverse_mapping_relation = np.transpose(self.mapping_relation)
        KOA.ExpressionUtils.Transpose(self.backward_matrix, self.forward_matrix)
        # self.inverse_mapping_relation = KM.Expression.NodalExpression(self.control_polygon_model_part)
        # KM.Expression.LiteralExpressionIO.SetData(self.inverse_mapping_relation, backward_matrix)
        # import pdb
        # pdb.set_trace()

        self.is_initialized = True
    
    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KOA.SHAPE]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = KM.Expression.NodalExpression(self.controlled_model_part)
        KM.Expression.LiteralExpressionIO.SetData(field, [0,0,0])
        return field

    def GetEmptyControlField(self):
        field = KM.Expression.NodalExpression(self.control_polygon_model_part)
        KM.Expression.LiteralExpressionIO.SetData(field, [0,0,0])
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_field   # characteristic values on the control polygon
    
    def GetPhysicalField(self) -> ContainerExpressionTypes:
        physical_shape_field = KM.Expression.NodalExpression(self.controlled_model_part)
        KM.Expression.NodalPositionExpressionIO.Read(physical_shape_field, KM.Configuration.Initial)
        return physical_shape_field

    @time_decorator(methodName="GetName")
    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 1:
            raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if KOA.SHAPE not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {KOA.SHAPE.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        physical_gradient = physical_gradient_variable_container_expression_map[KOA.SHAPE]
        if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
            raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.controlled_model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

        # print("MapGradient :: physical_gradient (KOA.SHAPE):\n", physical_gradient.Evaluate().shape)
        # import pdb
        # pdb.set_trace()
        control_gradient = self.ProjectBackward(physical_gradient)
        # KOA.ExpressionUtils.ComputeNodalVariableProductWithEntityMatrix(self.GetControlField(), physical_gradient, self.inverse_mapping_relation_matrix, self.control_polygon_model_part.Nodes)
        # KratosOA.ExpressionUtils.ComputeNodalVariableProductWithEntityMatrix(output_values, nodal_values, KratosOA.HELMHOLTZ_MASS_MATRIX, self.model_part.Elements)
        # print("MapGradient :: physical_gradient:\n", physical_gradient.Evaluate())
        # print("MapGradient :: control_gradient:\n", control_gradient.Evaluate())

        # filtered_gradient = self.filter.BackwardFilterIntegratedField(KOA.ExpressionUtils.ExtractData(physical_gradient, self.model_part))

        return KOA.ExpressionUtils.ExtractData(control_gradient, self.control_polygon_model_part.GetRootModelPart())

    def ProjectBackward(self, gradient):
        # gradient_unit = np.zeros(self.inverse_mapping_relation.shape[1])
        # gradient_unit[0] = 1.0
        mapped_gradient = KM.Expression.NodalExpression(self.control_polygon_model_part)
        KOA.ExpressionUtils.ProductWithEntityMatrix(mapped_gradient, self.backward_matrix, gradient)
        return mapped_gradient

    @time_decorator(methodName="GetName")
    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        # control_polygon_update = np.zeros((self.control_polygon_model_part.NumberOfNodes(), 3))
        # control_polygon_update[0,:] = np.array([0,0,0.25])    # artificial update for test

        # control_field_array = self.control_field.Evaluate()

        # new_control_field_array = control_field_array + control_polygon_update

        # print("Update :: new_control_field:\n", new_control_field.Evaluate())
        # print("Update :: self.control_field.PrintData():\n", self.control_field.Evaluate())
        # print("Update :: NormL2(self.control_field - new_control_field): ", KM.Expression.Utils.NormL2(self.control_field - new_control_field))
        
        if not IsSameContainerExpression(new_control_field, self.GetEmptyControlField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.control_polygon_model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")
        # import pdb
        # pdb.set_trace()
        
        if KM.Expression.Utils.NormL2(self.control_field - new_control_field) > 1e-15:
            # update the control SHAPE field
            control_update = new_control_field - self.control_field
            self.control_field = new_control_field
            # now update the physical field
            # new_physical_field = KM.Expression.NodalExpression(self.controlled_model_part)
            new_physical_field = self.GetPhysicalField()
            physical_field_update = KM.Expression.NodalExpression(self.controlled_model_part)
            KOA.ExpressionUtils.ProductWithEntityMatrix(physical_field_update, self.forward_matrix, control_update)
            new_physical_field += physical_field_update
            self._UpdateControlPolygon(self.control_field)
            self._UpdateMesh(new_physical_field)

            return True
        return False  
    
    def _UpdateControlPolygon(self, new_control_field) -> None:
        KM.Expression.NodalPositionExpressionIO.Write(new_control_field, KM.Configuration.Initial)
        KM.Expression.NodalPositionExpressionIO.Write(new_control_field, KM.Configuration.Current)

    
    def _UpdateMesh(self, new_physical_field: ContainerExpressionTypes) -> None:
        # physical_field = self.GetPhysicalField()
        KM.Expression.NodalPositionExpressionIO.Write(new_physical_field, KM.Configuration.Initial)
        KM.Expression.NodalPositionExpressionIO.Write(new_physical_field, KM.Configuration.Current)

    def GetControllingObjects(self):
        return self.controlling_objects
            


