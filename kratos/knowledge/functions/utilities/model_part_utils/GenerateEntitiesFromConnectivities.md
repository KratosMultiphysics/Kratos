---
type: C++ Function
title: GenerateEntitiesFromConnectivities
resource: utilities/model_part_utils.h#L144-L184
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void GenerateEntitiesFromConnectivities( const TEntity& rPrototypeEntity, const std::vector<std::vector<std::size_t>>& rEntitiesConnectivities, ModelPart::NodesContainerType& rThisNodes, PointerVectorSet<TEntity, IndexedObject, std::less<typename IndexedObject::result_type>, std::equal_to<typename IndexedObject::result_type>, typename TEntity::Pointer, std::vector<typename TEntity::Pointer>>& rThisEntities, const Properties::Pointer pProperties = nullptr )`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)