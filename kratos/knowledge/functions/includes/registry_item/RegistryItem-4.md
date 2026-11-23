---
type: C++ Function
title: RegistryItem
resource: includes/registry_item.h#L201-L206
generated:
  by: okf-rs/0.7.0
---

# Signature

`RegistryItem( const std::string& rName, const TItemType& rValue) : mName(rName), mpValue(Kratos::make_shared<TItemType>(rValue)), mGetValueStringMethod(&RegistryItem::GetItemString<TItemType>)`