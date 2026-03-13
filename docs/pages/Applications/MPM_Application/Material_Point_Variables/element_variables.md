---
title: List of MP Element Variables
keywords: mpm element variables
tags: [mpm element variables]
sidebar: mpm_application
summary: 
---

## MPMUpdatedLagrangian

The following variables are defined for each material point **element** of type [`MPMUpdatedLagrangian`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/custom_elements/mpm_updated_lagrangian.hpp).

* `MP_MATERIAL_ID` (`int`): id identifying the material (see the `properties_id` field in the `ParticleMaterials.json` [input file](../Input_Files/json#particlematerialsjson))
* `MP_MASS` (`double`): material point element **mass**
* `MP_DENSITY` (`double`): material point element **density**
* `MP_VOLUME` (`double`): material point element **volume**
* `MP_COORD` (`array_1d<double,3>`): material point element **coordinates**
* `MP_DISPLACEMENT` (`array_1d<double,3>`): material point element **displacement**
* `MP_VELOCITY` (`array_1d<double,3>`): material point element **velocity**
* `MP_ACCELERATION` (`array_1d<double,3>`): material point element **acceleration**
* `MP_CAUCHY_STRESS_VECTOR` (`Vector`): material point element **Cauchy stress vector**
* `MP_ALMANSI_STRAIN_VECTOR` (`Vector`): material point element **Almansi strain vector**

## MPMUpdatedLagrangianUP

For the material point elements of type [`MPMUpdatedLagrangianUP`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/custom_elements/mpm_updated_lagrangian_UP.hpp), implementing the mixed UP (displacement/pressure) formulation, the following additional variable is defined.

* `MP_PRESSURE` (`double`): material point element **pressure**
