---
type: C++ Function
title: RegistryItem
resource: includes/registry_item.h#L215-L220
generated:
  by: okf-rs/0.7.0
---

# Signature

`RegistryItem( const std::string& rName, const shared_ptr<TItemType>& pValue) : mName(rName), mpValue(pValue), mGetValueStringMethod(&RegistryItem::GetItemString<TItemType>)`