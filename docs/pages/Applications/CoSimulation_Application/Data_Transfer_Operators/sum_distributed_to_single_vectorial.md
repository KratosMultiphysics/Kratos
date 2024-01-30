---
title: sum distributed to single vectorial
keywords: 
tags: [sum_distributed_to_single_vectorial.md]
sidebar: cosimulation_application
summary: 
---
## Overview
DataTransferOperator to sum values on one interface and put it to another interface. It is a version of SumDistributedToSingle which does the same but for vectorial variables. Used e.g. for FSI with RigidBody, where the loads on the fluid interface are summed up and set to the RigidBody interface
## Available Transfer Options
- **swap_sign**:
- **add_values**: