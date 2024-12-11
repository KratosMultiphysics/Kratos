---
title: `Registry` and `RegistryItem` in Kratos Multiphysics
keywords: Registry
tags: [Registry]  
sidebar: kratos_for_developers
summary: Introduction to registry framework
---

## Overview

The `Registry` and `RegistryItem` classes in *Kratos Multiphysics* provide a flexible, hierarchical, and type-safe registry system. This system enables efficient storage and retrieval of shared objects, prototypes, and configurations across the application. The design is optimized for extensibility, thread-safety, and integration with Kratos's infrastructure.

## Key Features

- **Hierarchical Structure**: A tree-like organization where each node is a `RegistryItem` that may store a value or sub-items.
- **Type-Safe Access**: Supports type-erased storage using `std::any` and allows type-safe retrieval.
- **Thread-Safe Operations**: Thread safety is ensured using locks for global operations.
- **Prototype Registration**: Facilitates the registration of class prototypes, enabling factory-like usage for dynamic object creation.
- **Integration with Python**: Python bindings expose key features for scripting and automation.
- **Custom Serialization**: Items can be serialized into JSON for debugging and reporting.

## Usage

### **Registering Prototypes**

The `KRATOS_REGISTRY_ADD_PROTOTYPE` macro simplifies the registration of prototypes. For example:

```cpp
KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.KratosMultiphysics", Modeler, ConnectivityPreserveModeler)
```

This registers the `ConnectivityPreserveModeler` as a prototype under the specified path. Once registered, objects can be dynamically created using the prototype.

### **Adding and Retrieving Items**

#### **Add a Value**

```cpp
Registry::AddItem<Variable<double>>("variables.all.DISPLACEMENT", DISPLACEMENT);
```

#### **Retrieve a Value**

```cpp
const auto& r_displacement_variable = Registry::GetValue<Variable<double>>("variables.all.DISPLACEMENT");
```

#### **Retrieve a Prototype**

```cpp
auto& prototype = Registry::GetValue<Modeler>("Modelers.All.ConnectivityPreserveModeler.Prototype");
auto p_instance = prototype.Create(model, parameters);
```

### **Python Integration**

The `Registry` is accessible in Python, allowing dynamic access to the registry's contents.

#### Accessing Registry Items in Python:

```python
from KratosMultiphysics import CppRegistry

# Check for an item
if CppRegistry.HasItem("variables.all.DISPLACEMENT"):
    displacement = CppRegistry.GetValue("variables.all.DISPLACEMENT")
    print("Displacement variable:", displacement)

# Iterate over registry keys
keys = CppRegistry.keys()
for key in keys:
    print("Registry Key:", key)
```

## Class API

### **Class: Registry**

#### Methods:

1. **AddItem**
   Registers a new item or prototype.
   ```cpp
   template<typename TItemType, class... TArgumentsList>
   static RegistryItem& AddItem(const std::string& rItemFullName, TArgumentsList&&... Arguments);
   ```

2. **GetValue**
   Retrieves a value stored in a registry item.
   ```cpp
   template<typename TDataType>
   static TDataType const& GetValue(const std::string& rItemFullName);
   ```

3. **RemoveItem**
   Removes an item from the registry.
   ```cpp
   static void RemoveItem(const std::string& rItemFullName);
   ```

4. **HasItem**
   Checks if an item exists.
   ```cpp
   static bool HasItem(const std::string& rItemFullName);
   ```

### **Class: RegistryItem**

#### Methods:

1. **AddItem**
   Adds a sub-item to the current registry item.
   ```cpp
   template<typename TItemType, class... TArgumentsList>
   RegistryItem& AddItem(const std::string& ItemName, TArgumentsList&&... Arguments);
   ```

2. **GetItem**
   Retrieves a sub-item by name.
   ```cpp
   RegistryItem& GetItem(const std::string& rItemName);
   ```

3. **GetValue**
   Retrieves the value stored in the item.
   ```cpp
   template<typename TDataType>
   TDataType const& GetValue() const;
   ```

4. **ToJson**
   Serializes the item and its sub-items to JSON.
   ```cpp
   std::string ToJson(const std::string& rTabSpacing = "", const std::size_t Level = 0) const;
   ```

## Examples

### **Prototype Registration**

Define and register a prototype:
```cpp
class CustomModeler : public Modeler {
    KRATOS_CLASS_POINTER_DEFINITION(CustomModeler);
    KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.KratosMultiphysics", Modeler, CustomModeler)
};
```

Retrieve and use the prototype:
```cpp
auto& prototype = Registry::GetValue<Modeler>("Modelers.KratosMultiphysics.CustomModeler.Prototype");
auto p_instance = prototype.Create(model, settings);
```

### **Testing**

The following tests validate the functionality of the registry and registry items:

1. **Adding and Retrieving Items**
   ```cpp
   KRATOS_TEST_CASE_IN_SUITE(RegistryItemGetValue, KratosCoreFastSuite) {
       Registry::AddItem<double>("Physics.Gravity", 9.81);
       double gravity = Registry::GetValue<double>("Physics.Gravity");
       KRATOS_EXPECT_DOUBLE_EQ(gravity, 9.81);
   }
   ```

2. **Prototype Access**
   ```cpp
   KRATOS_TEST_CASE_IN_SUITE(RegistryItemGetValueDerivedAsBase, KratosCoreFastSuite) {
       auto pProcess = Registry::GetValue<Process>("Processes.KratosMultiphysics.OutputProcess.Prototype");
       KRATOS_EXPECT_TRUE((std::is_base_of<Process, decltype(pProcess)>::value));
   }
   ```

3. **JSON Serialization**
   ```cpp
   KRATOS_TEST_CASE_IN_SUITE(RegistryJsonValue, KratosCoreFastSuite) {
       RegistryItem value_item("value_item", 3.14);
       std::string expected_json = R"({"value_item": "3.14"})";
       KRATOS_EXPECT_EQ(value_item.ToJson(), expected_json);
   }
   ```