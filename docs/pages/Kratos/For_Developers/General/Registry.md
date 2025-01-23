---
title: Registry and RegistryItem in Kratos Multiphysics
keywords: Registry
tags: [Registry]
sidebar: kratos_for_developers
summary: Introduction to registry framework
---

## Overview

The `Registry` and `RegistryItem` classes in *Kratos Multiphysics* provide a flexible, hierarchical, and type-safe registry system. This system enables efficient storage and retrieval of information, prototypes, and configurations across Kratos.

The Kratos `Registry` is internally structured in two layers. The deeper one is the C++ layer, which contains all the items that are registered at the C++ level. Complementary to this, there is the Python layer that adds all the Python items (i.e. modules) to the `Registry`. Note that the Python layer is deliberately design as an extension of the C++ one, meaning that all the C++ items can be accessed from Python as well. It is important to note that when accessing from Python, the users of the `Registry` do not require to know which is the layer containing the item that is to be accessed as the `Registry` handles this internally. Also note that, as expected, only C++ items can be retrieved when accessing from C++.

## Key Features

- **Hierarchical Structure**: A tree-like organization where each node is a `RegistryItem` that may store a value or sub-items.
- **Type-Safe Access**: Supports type-erased storage using `std::any` and allows type-safe retrieval.
- **Thread-Safe Operations**: Thread safety is ensured using locks for global operations.
- **Prototype Registration**: Facilitates the registration of class prototypes, enabling factory-like usage for dynamic object creation.
- **Integration with both C++ and Python**: Items from both C++ and Python layers can be registered.

## Usage

### **Registering Prototypes**

The `KRATOS_REGISTRY_ADD_PROTOTYPE` macro simplifies the registration of prototypes. Note that besides other purposes, it handles the key creation based on the class name. For example:

```cpp
KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.KratosMultiphysics", Modeler, ConnectivityPreserveModeler)
```

It creates two items for the `ConnectivityPreserveModeler`, one accessible with the key `Modelers.KratosMultiphysics.ConnectivityPreserveModeler.Prototype` and so for the other one with `Modelers.All.ConnectivityPreserveModeler.Prototype`. Once registered, objects can be dynamically created using the prototype. The macro can (and should) be called in custom C++ modules as well. This can be done as shown below for `CustomModeler` class:

```cpp
class CustomModeler : public Modeler {
    KRATOS_CLASS_POINTER_DEFINITION(CustomModeler);
    KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.KratosMultiphysics", Modeler, CustomModeler)
};
```

To then retrieve the corresponding prototype and use it as a factory as shown below

```cpp
auto& r_prototype = Registry::GetValue<Modeler>("Modelers.KratosMultiphysics.CustomModeler.Prototype");
auto p_instance = r_prototype.Create(model, settings);
```

### **Adding and Retrieving Items at C++ level**

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

### **Adding and Retrieving Items at Python level**

The `Registry` is accessible in Python, allowing dynamic registry and access to the registry's contents.

#### Registering Python items

Python items are registered when loading the corresponding module. The items to be registered are specified in the corresponding Python registry lists. Find [here](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_interface/python_registry_lists.py) example of the `KratosCore` python registry list or [here](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/FluidDynamicsApplication/python_registry_lists.py) that of the `FluidDynamicsApplication`.

Custom items can also be registered dynamically during the execution. Find an example below.

```python
class MyCustomStage(AnalysisStage):
    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

# Registering the above custom stage class
KratosMultiphysics.Registry.AddItem("Stages.KratosMultiphysics.MyCustomStage.ClassName", "MyCustomStage")
```

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