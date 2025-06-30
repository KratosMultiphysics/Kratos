---
title: Getting Kratos Binaries for Linux
keywords:
tags: [Getting-Kratos-Binaries-for-Linux.md]
sidebar: kratos_for_users
summary:
---

# Pip Install

Kratos binaries are available for Linux and Windows and can be optained with `pip`:

```console
pip install KratosMultiphysics-all
```

# Pip Packages

The following Kratos Pakcages ara avaialable

## Core

- KratosMultiphysics-all: Kratos core and all its dependencies
- KratosMultiphysics: Kratos Core

## Applications

- KratosStructuralMechanicsApplication
- KratosFluidDynamicsApplication
- KratosDEMApplication
- KratosContactStructuralMechanicsApplication
- KratosParticleMechanicsApplication;
- KratosConvectionDiffusionApplication;
- KratosDamApplication;
- KratosPoromechanicsApplication;
- KratosFSIApplication;
- KratosSwimmingDEMApplication;
- KratosLinearSolversApplication;
- KratosConstitutiveLawsApplication;
- KratosDelaunayMeshingApplication;
- KratosMeshingApplication;
- KratosDemStructuresCouplingApplication;
- KratosMeshMovingApplication;
- KratosCSharpWrapperApplication;
- KratosShapeOptimizationApplication;
- KratosCoSimulationApplication;
- KratosCableNetApplication;
- KratosRANSApplication;
- KratosMappingApplication;
- KratosCompressiblePotentialFlowApplication;

add_app ${KRATOS_APP_DIR}/IgaApplication;
add_app ${KRATOS_APP_DIR}/ChimeraApplication;
add_app ${KRATOS_APP_DIR}/StatisticsApplication;
add_app ${KRATOS_APP_DIR}/RomApplication;
add_app ${KRATOS_APP_DIR}/ShallowWaterApplication;
add_app ${KRATOS_APP_DIR}/OptimizationApplication;
add_app ${KRATOS_APP_DIR}/GeoMechanicsApplication;

## MPI

Please note that this packages are only available for linux:

- KratosTrilinosApplication

Please note that Kratos binaries are not available for MacOS.

