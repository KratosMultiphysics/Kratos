---
title: KRATOS_API for Windows Compilation
keywords: 
tags: [KRATOS_API-for-Windows-Compilation.md]
sidebar: kratos_sidebar
summary: 
---

Kratos windows compilation is done via Microsoft Visual C++ compiler. Therefore when developement is done in C++ level, if some one has declarations in header file and definitions in `cpp` files, then followings need to considered.

# Classes
If a class is having seperated header file and cpp, add **KRATOS_API(<APPLICATION_NAME>)** before class name in the header file as explained in following code snippet.
```
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) TestClass
{
// rest of the code is unchanged
};
```

# Methods
If someone implemnets methods in a namespace and it has declaration in header and definition in cpp, then add **KRATOS_API(<APPLICATION_NAME>)** to the front of the method declaration in the header file as explained in following code snippet.
```
void KRATOS_API(FLUID_DYNAMICS_APPLICATION) TestMethod(int Value);
```

If the method is templated, and you have explicit template instantiations or definitions, then you need to again add **KRATOS_API(<APPLICATION_NAME>)** to all of the template instantiations as explained in following code snippets.

In the header file:
```
template<class T>
void KRATOS_API(FLUID_DYNAMICS_APPLICATION) TestTemplatedMethod(const T& Value);
```

In the cpp file:
```
// template definition
template<class T>
void TestTemplatedMethod<T>(const T& value)
{
// function body unchanged
}

// explicit template instantiation and definition
template<>
void KRATOS_API(FLUID_DYNAMICS_APPLICATION) TestTemplatedMethod<int>(const int& Value)
{
// function body unchanged
}

// explicit template instantiation
template void KRATOS_API(FLUID_DYNAMICS_APPLICATION) TestTemplatedMethod<double>(const double&);
```