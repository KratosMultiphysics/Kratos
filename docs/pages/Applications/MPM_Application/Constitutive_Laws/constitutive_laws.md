---
title: Constitutive Laws
keywords: mpm constitutive laws
tags: [mpm constitutive laws]
sidebar: mpm_application
summary: 
---

## Linear Elasticity

The **linear elastic** constitutive law is identified by the following labels:

- `LinearElasticIsotropicPlaneStrain2DLaw`: two-dimensional problem, plane strain formulation;
- `LinearElasticIsotropicPlaneStress2DLaw`: two-dimensional problem, plane stress formulation;
- `LinearElasticIsotropic3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `YOUNG_MODULUS`
- `POISSON_RATIO`

## Hyperelastic NeoHookean

The **hyperelastic NeoHookean** constitutive law is identified by the following labels:

- `HyperElasticNeoHookeanPlaneStrain2DLaw`: two-dimensional problem, plane strain and irreducible formulation;
- `HyperElasticNeoHookeanPlaneStrainUP2DLaw`: two-dimensional problem, plane strain and mixed formulation;
- `HyperElasticNeoHookean3DLaw`: three-dimensional problem, irreducible formulation;
- `HyperElasticNeoHookeanUP3DLaw`: three-dimensional problem, mixed formulation.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `YOUNG_MODULUS`
- `POISSON_RATIO`

## Mohr Coulomb

The **plastic Mohr Coulomb** constitutive law is identified by the following labels:

- `HenckyMCPlasticPlaneStrain2DLaw`: two-dimensional problem, plane-strain formulation;
- `HenckyMCPlastic3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `YOUNG_MODULUS`
- `POISSON_RATIO`
- `COHESION`
- `INTERNAL_FRICTION_ANGLE`
- `INTERNAL_DILATANCY_ANGLE`

## Mohr Coulomb Strain Softening

The **Mohr Coulomb** with **Strain Softening** constitutive law is identified by the following labels:

- `HenckyMCStrainSofteningPlasticPlaneStrain2DLaw`: two-dimensional problem, plane-strain formulation;
- `HenckyMCStrainSofteningPlastic3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `YOUNG_MODULUS`
- `POISSON_RATIO`
- `COHESION`: cohesion (peak)
- `COHESION_RESIDUAL`: cohesion (residual)
- `INTERNAL_FRICTION_ANGLE`: internal friction angle (peak)
- `INTERNAL_FRICTION_ANGLE_RESIDUAL`: internal friction angle (residual)
- `INTERNAL_DILATANCY_ANGLE`: internal dilatancy angle (peak)
- `INTERNAL_DILATANCY_ANGLE_RESIDUAL`: internal dilatancy angle (residual)
- `SHAPE_FUNCTION_BETA`: exponential softening beta coefficient

## Modified Cam Clay

The **modified cam clay** constitutive law is identified by the following labels:

- `HenckyBorjaCamClayPlasticPlaneStrain2DLaw`: two-dimensional problem, plane-strain formulation;
- `HenckyBorjaCamClayPlastic3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `PRE_CONSOLIDATION_STRESS`: preconsolidation pressure
- `OVER_CONSOLIDATION_RATIO`: over Consolidation Ratio (OCR)
- `SWELLING_SLOPE`: slope of swelling line
- `NORMAL_COMPRESSION_SLOPE`: slope of Normal Consolidation Line (NCL)
- `CRITICAL_STATE_LINE`: slope of Critical State Line (CSL)
- `INITIAL_SHEAR_MODULUS`: initial Shear Modulus
- `ALPHA_SHEAR`: volumetric-deviatoric coupling constant

## Newtonian Fluid

The **displacement-based Newtonian fluid** constitutive law is identified by the following labels:

- `DispNewtonianFluidPlaneStrain2DLaw`: two-dimensional problem, plane-strain formulation;
- `DispNewtonianFluid3DLaw`: three-dimensional problem.

These are the admissible (string) values that can be assigned to the field `"name"` of the `"constitutive_law"` section in the file `ParticleMaterials.json`.
More details about the input file `ParticleMaterials.json` can be found [here](../Input_Files/json#particlematerialsjson).

The variables that must be included in the `"Variables"` section of the input file `ParticleMaterials.json` are:
- `DENSITY`
- `BULK_MODULUS`
- `DYNAMIC_VISCOSITY`
