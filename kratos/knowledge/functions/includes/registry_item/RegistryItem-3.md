---
type: C++ Function
title: RegistryItem
resource: includes/registry_item.h#L187-L192
generated:
  by: okf-rs/0.7.0
---

# Signature

`RegistryItem( const std::string &rName, const std::function<std::shared_ptr<TItemType>(TArgs...)> &rValue) : mName(rName), mpValue(rValue()), mGetValueStringMethod(&RegistryItem::GetItemString<TItemType>)`