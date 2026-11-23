---
type: C++ Method
title: NotFoundError
resource: sources/registry.cpp#L158-L171
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/registry_item/RegistryItem/GetSubItemAvailableList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Registry::NotFoundError( const std::string& rFullName, const std::string& rItemName, RegistryItem* pCurrentItem )`

# Calls

- [GetSubItemAvailableList](../../../../functions/sources/registry_item/RegistryItem/GetSubItemAvailableList.md)