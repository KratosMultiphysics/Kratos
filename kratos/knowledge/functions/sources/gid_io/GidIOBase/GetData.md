---
type: C++ Method
title: GetData
resource: sources/gid_io.cpp#L35-L37
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/conditions/mesh_condition/MeshCondition/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/constraints/linear_master_slave_constraint/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/mesh_element/MeshElement/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/hexahedra_3d_20/Hexahedra3D20/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/hexahedra_3d_27/Hexahedra3D27/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/hexahedra_3d_8/Hexahedra3D8/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_2/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_3/Line2D3/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_4/Line2D4/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_5/Line2D5/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_3d_2/Line3D2/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_3d_3/Line3D3/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/point_2d/Point2D/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/point_3d/Point3D/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/prism_3d_15/Prism3D15/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/prism_3d_6/Prism3D6/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/pyramid_3d_13/Pyramid3D13/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/pyramid_3d_5/Pyramid3D5/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrature_point_geometry/QuadraturePointGeometry/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_2d_4/Quadrilateral2D4/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_2d_8/Quadrilateral2D8/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_2d_9/Quadrilateral2D9/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_3d_4/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_3d_8/Quadrilateral3D8/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_3d_9/Quadrilateral3D9/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/sphere_3d_1/Sphere3D1/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/tetrahedra_3d_10/Tetrahedra3D10/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_10/Triangle2D10/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_15/Triangle2D15/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_3/Triangle2D3/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_6/Triangle2D6/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_6/Triangle3D6/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/condition/Condition/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/element/Element/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/geometrical_object/Has
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/geometrical_object/SetValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/GidIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/master_slave_constraint/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/ExecuteFinalize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_model_part_to_python/GetNonHistoricalVariablesNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/ProcessMultifreedomConstraint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/MakeRelationTopology
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/AssembleRelationMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/linear_multifreedom_constraint/LinearMultifreedomConstraint/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/DeepCopyEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/properties_utilities/CopyPropertiesValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/EraseNonHistoricalVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/ClearNonHistoricalData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int GidIOBase::GetData()`

# Called by

- [Clone](../../../../functions/conditions/mesh_condition/MeshCondition/Clone.md)
- [Clone](../../../../functions/constraints/linear_master_slave_constraint/Clone.md)
- [Clone](../../../../functions/elements/mesh_element/MeshElement/Clone.md)
- [Create](../../../../functions/geometries/geometry/Create.md)
- [Create](../../../../functions/geometries/hexahedra_3d_20/Hexahedra3D20/Create.md)
- [Create](../../../../functions/geometries/hexahedra_3d_27/Hexahedra3D27/Create.md)
- [Create](../../../../functions/geometries/hexahedra_3d_8/Hexahedra3D8/Create.md)
- [Create](../../../../functions/geometries/line_2d_2/Create.md)
- [Create](../../../../functions/geometries/line_2d_3/Line2D3/Create.md)
- [Create](../../../../functions/geometries/line_2d_4/Line2D4/Create.md)
- [Create](../../../../functions/geometries/line_2d_5/Line2D5/Create.md)
- [Create](../../../../functions/geometries/line_3d_2/Line3D2/Create.md)
- [Create](../../../../functions/geometries/line_3d_3/Line3D3/Create.md)
- [Create](../../../../functions/geometries/point_2d/Point2D/Create.md)
- [Create](../../../../functions/geometries/point_3d/Point3D/Create.md)
- [Create](../../../../functions/geometries/prism_3d_15/Prism3D15/Create.md)
- [Create](../../../../functions/geometries/prism_3d_6/Prism3D6/Create.md)
- [Create](../../../../functions/geometries/pyramid_3d_13/Pyramid3D13/Create.md)
- [Create](../../../../functions/geometries/pyramid_3d_5/Pyramid3D5/Create.md)
- [Create](../../../../functions/geometries/quadrature_point_geometry/QuadraturePointGeometry/Create.md)
- [Create](../../../../functions/geometries/quadrilateral_2d_4/Quadrilateral2D4/Create.md)
- [Create](../../../../functions/geometries/quadrilateral_2d_8/Quadrilateral2D8/Create.md)
- [Create](../../../../functions/geometries/quadrilateral_2d_9/Quadrilateral2D9/Create.md)
- [Create](../../../../functions/geometries/quadrilateral_3d_4/Create.md)
- [Create](../../../../functions/geometries/quadrilateral_3d_8/Quadrilateral3D8/Create.md)
- [Create](../../../../functions/geometries/quadrilateral_3d_9/Quadrilateral3D9/Create.md)
- [Create](../../../../functions/geometries/sphere_3d_1/Sphere3D1/Create.md)
- [Create](../../../../functions/geometries/tetrahedra_3d_10/Tetrahedra3D10/Create.md)
- [Create](../../../../functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/Create.md)
- [Create](../../../../functions/geometries/triangle_2d_10/Triangle2D10/Create.md)
- [Create](../../../../functions/geometries/triangle_2d_15/Triangle2D15/Create.md)
- [Create](../../../../functions/geometries/triangle_2d_3/Triangle2D3/Create.md)
- [Create](../../../../functions/geometries/triangle_2d_6/Triangle2D6/Create.md)
- [Create](../../../../functions/geometries/triangle_3d_3/Create.md)
- [Create](../../../../functions/geometries/triangle_3d_6/Triangle3D6/Create.md)
- [Clone](../../../../functions/includes/condition/Condition/Clone.md)
- [Clone](../../../../functions/includes/element/Element/Clone.md)
- [Has](../../../../functions/includes/geometrical_object/Has.md)
- [SetValue](../../../../functions/includes/geometrical_object/SetValue.md)
- [GidIO](../../../../functions/includes/gid_io/GidIO.md)
- [Clone](../../../../functions/includes/master_slave_constraint/Clone.md)
- [ExecuteFinalize](../../../../functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/ExecuteFinalize.md)
- [Execute](../../../../functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/Execute.md)
- [GetNonHistoricalVariablesNames](../../../../functions/python/add_model_part_to_python/GetNonHistoricalVariablesNames.md)
- [ProcessMultifreedomConstraint](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/ProcessMultifreedomConstraint.md)
- [MakeRelationTopology](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/MakeRelationTopology.md)
- [AssembleRelationMatrix](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/constraint_utilities/AssembleRelationMatrix.md)
- [CalculateLocalSystem](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/linear_multifreedom_constraint/LinearMultifreedomConstraint/CalculateLocalSystem.md)
- [WriteProperties](../../../../functions/sources/model_part_io/ModelPartIO/WriteProperties.md)
- [WriteDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/WriteDataBlock.md)
- [DeepCopyEntities](../../../../functions/utilities/auxiliar_model_part_utilities/DeepCopyEntities.md)
- [CopyPropertiesValues](../../../../functions/utilities/properties_utilities/CopyPropertiesValues.md)
- [EraseNonHistoricalVariable](../../../../functions/utilities/variable_utils/EraseNonHistoricalVariable.md)
- [ClearNonHistoricalData](../../../../functions/utilities/variable_utils/ClearNonHistoricalData.md)