---
type: C++ Function
title: GenerateEntitiesFromConnectivities
resource: utilities/model_part_utils.h#L288-L301
generated:
  by: okf-rs/0.7.0
---

# Signature

`static void GenerateEntitiesFromConnectivities( const std::string& rEntityName, const std::vector<std::size_t>& rEntitiesIds, const std::vector<std::size_t>& rPropertiesIds, const std::vector<std::vector<std::size_t>>& rEntitiesConnectivities, ModelPart::NodesContainerType& rThisNodes, ModelPart::PropertiesContainerType& rThisProperties, PointerVectorSet<TEntity, IndexedObject, std::less<typename IndexedObject::result_type>, std::equal_to<typename IndexedObject::result_type>, typename TEntity::Pointer, std::vector<typename TEntity::Pointer>>& rThisEntities )`