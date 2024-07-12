---
title: Overview
keywords: 
tags: [overview.md]
sidebar: documentation_guide
summary: 
---

## Overview

This documentation is the central location to provide all the necessary information to Kratos **users** and **developers**. This has two parts, first one is this documentation which is more focused on showcasing abilities and applications with examples tailored mostly for **users**. While this is tailored for **users**, documentors are allowed to have in here code snippets, hints, guides for developers. A seperate documentation is there which created using Doxygen for the **C++** and **python** layer of the code. Developers are advised to follow proper documentation guide lines in their **C++** or **python** source code so it can be linked/reused/viewed from this documentation.

This documentation has three guides.

* How to create documentation content for this documentation guide.
* How to properly document in **C++** source code.
* How to properly document in **python** source code.

## Latex support

MarkDown language is used to create documentation guide. This is enabled with latex support to illustrate equations if necessary. Following is an example latex code.
```latex
<p align="center">$$ \underline{r} = \sum_{i=1}^N{\underline{x}_i}$$</p>
```

This is showed as following in the respective content page.
<p align="center">$$ \underline{r} = \sum_{i=1}^N{\underline{x}_i}$$</p>

## Automation

This documentation uses the folder structure and special ```menu_info.json``` files in each folder structure to automate creation of navigation bars. And the MarkDown files starts with ```remote_*.md``` are automatically generated MarkDown files therefore, it is advised not to change them since they will be overwritten when the folder structure is processed.