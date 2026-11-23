---
type: C++ Method
title: GetClosestContactColorFromNeighbours
resource: modeler/operation/find_contacts_in_skin_model_part.cpp#L162-L197
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetCenterOfElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/FindConditionContact
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindContactsInSkinModelPart::GetClosestContactColorFromNeighbours( const std::vector<int>& rContactColors, const std::vector<int>& rNeighbourColors, const std::vector<array_1d<std::size_t, 3>>& rNeighbourIndexes, const array_1d<double,3>& rConditionCenter, const array_1d<double,3>& rConditionNormal, const CartesianMeshColors& rColors, int& rPairingStatus, double& rMinDistance, int& rTemporaryColor) const`

# Calls

- [GetCenterOfElement](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetCenterOfElement.md)
- [norm_2](../../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [FindConditionContact](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/FindConditionContact.md)