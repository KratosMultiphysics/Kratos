---
title: copy with empty ranks
keywords: 
tags: [copy_with_empty_ranks.md]
sidebar: cosimulation_application
summary: 
---
## Overview
DataTransferOperator that copies values from one interface to another, without any checks. It is an extended version of CopyDataTransferOperation which considers solvers that only run in rank 0 (e.g. the RigidBodySolver).
## Available Transfer Options
- **swap_sign**:
- **add_values**: