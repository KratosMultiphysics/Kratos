---
title: Porting to PyBind11   common steps
keywords: 
tags: [Porting-to-PyBind11---common-steps.md]
sidebar: kratos_sidebar
summary: 
---

# Modifications in the CMakeLists.txt
*Pybind11* comes with a module to help compilation. Such module should be included at the beginning of the CMakeLists.txt by adding:

    include(pybind11Tools)

other than that one needs to change

    add_library(KratosFluidDynamicsApplication SHARED  ......)

to

    pybind11_add_module(KratosFluidDynamicsApplication MODULE ............. )

# Moving to module interface in pybind11
The definition of new python modules is slightly different in *Pybind* than in *Boost*. Python
the first change is that the python interfaces should be modified as in

![Change 1](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Porting-to-PyBind11---common-steps/change1.png)

[Link](https://github.com/KratosMultiphysics/Kratos/blob/eb4f48418e9b28e4a920fd1c524842ac42a84b44/applications/FluidDynamicsApplication/custom_python/add_custom_constitutive_laws_to_python.h)

and 

![Change 2](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Porting-to-PyBind11---common-steps/change2.png)

[Link](https://github.com/KratosMultiphysics/Kratos/blob/eb4f48418e9b28e4a920fd1c524842ac42a84b44/applications/FluidDynamicsApplication/custom_python/add_custom_constitutive_laws_to_python.cpp)

The changes include:
* Removing includes of *Boost* python
* `using namespace boost::python` ----> `namespace py = pybind11;`
* Adding the argument `pybind11::module& m`

Aside of this the class interface is slightly changed:
* A new argument `m` (the module) has to be passed as first argument.
* `bases` is not any longer required. Only its template argument must remain, in the same position of what was there before.
* The `noncopyable` template argument should not be provided (everything is `noncopyable` unless specified) - if something is to me made copyable, a copy constructor should be provided to python
* `init` has to be made as a .def (even the `first` one)
* a Pointer definition should be explicitly provided **if the base class provides it**

for example:

    class_< Euler2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "Euler2DLaw",  init<>() );

becomes

    py::class_< Euler2DLaw, Euler2DLaw::Pointer, ConstitutiveLaw >(m,"Euler2DLaw")
        .def(  py::init<>() );

One last important remark is to update the CMakeLists.txt

![Change 3](https://github.com/KratosMultiphysics/Documentation/blob/master/Wiki_files/Porting-to-PyBind11---common-steps/change_3.png?raw=true)

## Internal references in python
When references are to be used

     return_internal_reference<>() ---> py::return_value_policy::reference_internal

## Using python lists
*Pybind* list interface is similar to the one of `boost::python`, however the extraction of data from a list is slightly different

For example:

     boost::python::extract<ShellCrossSection::Pointer>(seclist[i])

Becomes:

     pybind11::cast<ShellCrossSection::Pointer >( seclist[i] );

## `staticmethod`

Definition of static methods is changed

     .def("AddToProperties", &TClassName::AddToProperties< Variable< TClassName >>).staticmethod("AddToProperties")

Becomes

     .def_static("AddToProperties", &TClassName::AddToProperties< Variable< TClassName > >)

## self_ns::str

This construct no longer exists in pybind11. Its previous use was to print a stringified version of your class. In pybind you must provide it yourself by defining the __str__ function 

```c++
.def(self_ns::str(self))
```

becomes 

```c++
.def("__str__", &Class::ToStrFunction)
```
Note that a default implementation of a `__str__` called `PrintObject`, which should be useful for most Kratos classes, is provided in `define_python.h`. It can be used as:
```c++
py::class_< SomeClass, SomeClass::Pointer, SomeBaseClass>
// other methods
.def("__str__", PrintObject<SomeClass>)
```

## Modifications in the "configure" file
modify the following line of the file:

-DPYTHON_LIBRARY="C:\Python34\libs\python34.lib"

with

-DPYTHON_EXECUTABLE="c:/python/python34.exe"

## Windows compatibility issues

Originally posted in [this PR](https://github.com/KratosMultiphysics/Kratos/pull/1830#issuecomment-380099368).


_**The problem**_
The problem appears on Windows when someone tries to link against an application. With boost, both the application code and application interface where compiled inside the same object, and their symbols exported when needed using the `KRATOS_API` macro. This ended up generating a dynamic object that was the python module and contained the application itself, altogether in the same package.

Pybind changed how things work a bit while creating modules, and one of the changes made impossible to to have a python module as a dependency of one application, drawing derived apps impossible to compile with the old system (as the symbols were inside the module).

In order to solve this, we agreed to divide the applications that are dependencies into `FooCore` and `FooApplication` objects, much as the `KratosCore` & `Kratos` objects work.

This has increased the number of symbols that need to be visible from the "interface", all files that with `add_xxxxx_to_python`, as they are no longer in the same place, and this is what is causing the linker errors.

_**Solution**_
Just compile your application on Windows and see if there are linking errors. I identified myself the errors under 4 categories:

**Missing exported classes**:

```c++
class foo {}
```
to
```c++
class KRATOS_API(APPLICATION) foo {}
```

**Variables**:

```c++
KRATOS_DEFINE_VARIABLE(type, FOO)
```
to
```c++
KRATOS_DEFINE_APPLICATION_VARIABLE(APPLICATION, type, FOO)
```

**Local Flags**:

```c++
KRATOS_DEFINE_LOCAL_FLAG(FOO)
```
to
```c++
KRATOS_DEFINE_LOCAL_APPLICATION_FLAG(APPLICATION, FOO)
```
**Template instantiation in the interface of non header classes**

This is more tricky. If you happen to expose a templated class and such class in non header-only, you must explicitly instantiate that in your "core" application, for example in the foo_application.cpp file.

```c++
#foo.h
template<typename Ta, typename Tb> class KRATOS_API(FOO_APPLICATION) MyTemplateClass{
}

#foo.cpp
MyTemplateClass::MyTemplateClass() {
    std::cout << "I am implemented in a cpp file" << std::endl;
}
```
and
```c++
# add_utilities_to_python.cpp
class_<MyTemplateClass<int, int>>(m, "MyTemplateClassInts").def(init<>);
class_<MyTemplateClass<double, double>>(m, "MyTemplateClassDoubles").def(init<>);
```

This will crash, and you need to add:

```c++
# foo_application.cpp
template class MyTemplateClass<int, int>;
template class MyTemplateClass<double, double>;
```
