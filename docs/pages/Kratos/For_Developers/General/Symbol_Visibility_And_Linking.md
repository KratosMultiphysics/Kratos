---
title: KRATOS_API for Windows Compilation
keywords: 
tags: [KRATOS_API Windows Compilation]
sidebar: kratos_for_developers
summary: 
---

Kratos Windows compilation is done via Microsoft Visual C++ compiler. Therefore when developement is done in C++ level, if some one has declarations in header file and definitions in `cpp` files, then followings need to considered.

# Classes
If a class is having seperated header (`.h, .hpp`) and source (`.c, .cpp`) files, add `KRATOS_API(${APPLICATION_NAME})` before class name in the header file as explained in following code snippet.
```cpp
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) TestClass
{
    // Rest of the code is unchanged
};
```
{: data-lang="C++"}

# Methods
If someone implemnets methods in a namespace and it has declaration in header and definition in cpp, then add `KRATOS_API(${APPLICATION_NAME})` to the front of the method declaration in the header file as explained in following code snippet.
```cpp
void KRATOS_API(FLUID_DYNAMICS_APPLICATION) TestMethod(int Value);
```
{: data-lang="C++"}

If the method is templated, and you have explicit template instantiations or definitions, then you need to again add `KRATOS_API(${APPLICATION_NAME})` to all of the template instantiations as explained in following code snippets.

In the header file:
```cpp
template<class T>
void KRATOS_API(FLUID_DYNAMICS_APPLICATION) TestTemplatedMethod(const T& Value);
```
{: data-lang="C++"}

In the cpp file:
```cpp
// Template definition
template<class T>
void TestTemplatedMethod<T>(const T& value)
{
    // Function body unchanged
}

// Explicit template instantiation and definition
template<>
void KRATOS_API(FLUID_DYNAMICS_APPLICATION) TestTemplatedMethod<int>(const int& Value)
{
    // Function body unchanged
}

// Explicit template instantiation
template void KRATOS_API(FLUID_DYNAMICS_APPLICATION) TestTemplatedMethod<double>(const double&);
```
{: data-lang="C++"}