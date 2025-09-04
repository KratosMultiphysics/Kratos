---
title: Smart Pointers within Kratos
keywords: 
tags: [How-to-use-Smart-Pointers-within-Kratos.md]
sidebar: kratos_for_developers
summary: 
---

## Overview

All the memory management within Kratos is done through "Shared Pointers". Essentially a shared pointer _is an entity which holds a counter with the number of existing instances of an object_. **EVERY TIME** a new Shared Pointer is created such counter is incremented or decremented to represent the number of living instances of the object pointed to.

A good description of the design and behaviour of `shared_ptr` can be found at the links:

* [Wikipedia](http://en.wikipedia.org/wiki/Smart_pointer) 
* [Boost](http://www.boost.org/doc/libs/1_63_0/libs/smart_ptr/smart_ptr.htm) 

In the practice when using Shared Pointers within Kratos one should be aware of their performance pitfalls so to be able to avoid them.

The use of an existing pointer does not imply any performance penalty with respect to a traditional pointer, **HOWEVER** while creating and destroying a traditional "_c-style_" pointer is a cheap operation, the creation or destruction of a shared point is relatively time consuming. This is so, as at the moment of creating/destructing a Shared Pointer the number of references to the object should be incremented or decremented, which implies an atomic operation if **OpenMP** parallelism is employed.

Within Kratos the great majority of shared_ptrs is stored in vectors, typically in the classes:

```cpp
PointerVector (used as a basis for the "geometry" class)
PointerVectorSet (NodesContainerType, ElementsContainerType, ConditionsContainerType)
PointerVectorMap
``` 

All of such objects provide two distinct operators, [] and (). 

### OPERATOR [] 

It returns a _reference_ to the object pointed by the underlying pointer. This operator should be used in most cases.

### OPERATOR () 

It returns a _pointer_ to the object, which will often require creating a copy of the pointer. This is often **slow**, only use it if you **really** need a pointer.

Similarly, **avoid** using GetGeometry().pGetNode(...) unless you really need a pointer, since this will imply allocating a new pointer and is hence an expensive operation.  
