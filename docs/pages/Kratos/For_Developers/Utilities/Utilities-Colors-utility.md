---
title: Utilities Colors utility
keywords: 
tags: [Utilities-Colors-utility.md]
sidebar: kratos_for_developers
summary: 
---

# Content
* [What is the colors utility][what]
* [Functions provided by this utility][functions]

[what]: https://github.com/KratosMultiphysics/Kratos/wiki/Colors-utility#what-is-the-colors-utility
[functions]: https://github.com/KratosMultiphysics/Kratos/wiki/Colors-utility#functions-provided-by-this-utility

# What is the colors utility

A model part contains entities like nodes, elements and conditions, and even, another model part. The structure of the main model part and it's sub model parts is very helpfull to define the problem we are about to solve.

While a sub model part (a model part) gives us the entities which are stored in, those entities doesn't know about the sub model parts they are stored in. Sometimes we need to know the sub model parts to which an entity belong to, e.g., when we create new entities. This problem is solved across the colors: *a **color** is a key and its corresponding value gives us the sub model part names*. The color map is defined as:

```cpp
typedef std::unordered_map<int,std::vector<std::string>> IntStringMapType;
```

When an entity belongs to more than one sub model part, there will be a combination of sub model parts and a new color will be defined. Using this concept of ***color**, every entity will have a unique color*. 

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Colors-Utility/model-part.png)

:warning: A color is an integer, not an RGB color

The previous picture show the structure of a model part:

```
MainModelPart
	BodySubModelPart
	SkinSubModelPart
		LeftSubModelPart
```

And the possible colors of this model part are:

| Color | Names |
| ------- | --------- |
| 0 | MainModelPart |
| 1 | BodySubModelPart |
| 2 | SkinSubModelPart|
| 3 | LeftSubModelPart|
| 4 | BodySubModelPart, SkinSubModelPart |
| 5 | BodySubModelPArt, SkinSubModelPart, LeftSubModelPart |

Finally, the colors of the entities are maps with the Id of an entity and its corresponding color:

```cpp
typedef std::unordered_map<IndexType,int> IndexIntMapType;
```


# Functions provided by this utility

