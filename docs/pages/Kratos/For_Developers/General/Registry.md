# Kratos Registry

The Kratos registry is a mechanism that enables the access to any class (process, utility, etc...) from anywhere in kratos (similar to how KratosComponents work)

In a nutshell, registry is composed of a **name** (string), which identifies the class, and a **prototype** (your class), which will be given to you once you acces it.

## Basics

### Register

In order to make your class avaialbe to use through the registry, you only have to add a macro to your class definition with the key.

To keep the registry clean and well organized and iterable, we have a couple rules on how you need to register classes

- Classes are registered **at the private section** of your class.
- Classes are registered **in the category they belong**: `ApplyInletProcess` is a `Process`, so it needs to be registered into `"Process"`.
- Classes are registered **in the application they are defined**, and in the "All" category. Following the last example. `ApplyInletProcess` should be in `"All"` and `"KratosMultiphyscis"`, because it belongs to the core.

You are also free to register things in custom places / categories.

Example:

```C++
class ApplyInletProcess {
    public: 
        ...
    private:

        KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.KratosMultiphysics", ApplyInletProcess)
        KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.All", ApplyInletProcess)

        ...
}
```

If the process were to beling to an application:

```C++
class FluidDynamicsProcess {
    public: 
        ...
    private:

        KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.FluidDynamicsApplication", FluidDynamicsProcess)
        KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.All", FluidDynamicsProcess)

        ...
}
```

### Access

In order to retrieve a proptotype or a value from the registry one simply hasmto use the `GetItem(std::string item_name)` function and get its value with `GetValue<TValueType>`.

For example, if I registered a custom Process called `MyCustomProcess` and I want to create an instance of it:

Example:
```C++
// Register
class MyCustomProcess {
    public: 
        ...
    private:

        KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.MyApplication", MyCustomProcess)
        KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.All", MyCustomProcess)

        ...
}
```

```C++
// Usage (Note that I don't have to be in MyApplication anymore)
auto my_model = Model();
auto my_params = KratosParameters();

auto my_custom_process = Registry::GetItem("Processes.All.MyCustomProcess").GetValue<Process>().Create(my_model, my_params);

my_custom_process.Execute();
```

## Advanced

### Dynamic Register

If More complex operations are desired, you can registry items dinamycally by calling the `Registry::AddItem<T>()`.

`AddItem` matches a `RegistryItem` with a `string` parameter if trying to add an intermediate node, or a `TDataType`, a `string`


### Value Register

Aside from classes, it is possible to directly literals in the registry

### Iterating

Registry allows iterating over sets of items 