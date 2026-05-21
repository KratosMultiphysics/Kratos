---
title: copy single to distributed vectorial
keywords: 
tags: [copy_single_to_distributed_vectorial.md]
sidebar: cosimulation_application
summary: 
---
## Overview
DataTransferOperator to take one single value and set it to all values on another interface. It is a version of CopySingleToDistributed which does the same but for vectorial variables. Used e.g. for FSI with RigidBody, where the RigidBody has one value and the fluid interface has many.
## Available Transfer Options
- **swap_sign**:
- **add_values**:
- **distribute_values**: