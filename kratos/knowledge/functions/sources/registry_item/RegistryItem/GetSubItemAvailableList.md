---
type: C++ Method
title: GetSubItemAvailableList
resource: sources/registry_item.cpp#L94-L104
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/registry/Registry/NotFoundError
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/registry_item/RegistryItem/NotFoundError
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> RegistryItem::GetSubItemAvailableList() const`

# Called by

- [NotFoundError](../../../../functions/sources/registry/Registry/NotFoundError.md)
- [NotFoundError](../../../../functions/sources/registry_item/RegistryItem/NotFoundError.md)