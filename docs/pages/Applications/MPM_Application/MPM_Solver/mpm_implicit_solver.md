---
title: MPM Implicit Dynamic Solver
keywords: mpm constitutive laws
tags: [mpm constitutive laws]
sidebar: mpm_application
summary: 
---

The `MPMImplicitDynamicSolver` implements both the **Newmark** and **Bossak-Newmark** time schemes.

The parameters to be included in the `solver_settings` section of the `ProjectParameters.json` input file are the following.

## Newmark Time Scheme

```json
{
    "solver_type"             : "Dynamic",
    "time_integration_scheme" : "implicit",
    "scheme_type"             : "newmark",
    "newmark_beta"            : 0.25
}
```

## Bossak Time Scheme

```json
{
    "solver_type"             : "Dynamic",
    "time_integration_scheme" : "implicit",
    "scheme_type"             : "bossak",
    "damp_factor_m"           : -0.3,
    "newmark_beta"            : 0.25
}
```
