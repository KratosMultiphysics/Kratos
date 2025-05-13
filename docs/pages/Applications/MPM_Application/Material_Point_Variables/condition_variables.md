---
title: List of MP Condition Variables
keywords: mpm element variables
tags: [mpm element variables]
sidebar: mpm_application
summary: 
---

## MPMParticleBaseCondition

The following variables are defined for each material point **condition** of type `MPMParticleBaseCondition`

* `MPC_AREA` (`double`): material point condition **integration area**
* `MPC_COORD` (`array_1d<double,3>`): material point condition **coordinates**
* `MPC_DISPLACEMENT` (`array_1d<double,3>`): material point condition **displacement**
* `MPC_VELOCITY` (`array_1d<double,3>`): material point condition **velocity**
* `MPC_ACCELERATION` (`array_1d<double,3>`): material point condition **acceleration**
* `MPC_NORMAL` (`array_1d<double,3>`): material point condition **normal vector**
* `MPC_DENSITY` (`double`): material point element **density**
* `MPC_VOLUME` (`double`): material point element **volume**
* `MPC_CAUCHY_STRESS_VECTOR` (`Vector`): material point element **Cauchy stress vector**
* `MPC_ALMANSI_STRAIN_VECTOR` (`Vector`): material point element **Almansi strain vector**

## MPMParticleBaseDirichletCondition

In addition to the variables defined for material point conditions of type `MPMParticleBaseCondition`,
the conditions of type `MPMParticleBaseDirichletCondition`, which are used to impose non-conforming
Dirichlet boundary conditions, define the following variables.

* `MPC_IMPOSED_DISPLACEMENT` (`array_1d<double,3>`): material point condition **imposed displacement**
* `MPC_IMPOSED_VELOCITY` (`array_1d<double,3>`): material point condition **imposed velocity**
* `MPC_IMPOSED_ACCELERATION` (`array_1d<double,3>`): material point condition **imposed acceleration**
* `MPC_CONTACT_FORCE` (`array_1d<double,3>`): material point condition **contact force**

## MPMParticlePenaltyDirichletCondition

In addition to the variables defined for material point conditions of type `MPMParticleBaseDirichletCondition`,
the conditions of type `MPMParticlePenaltyDirichletCondition`, which impose Dirichlet boundary conditions
by means of the penalty method, define the following variable

* `MPC_PENALTY_FACTOR` (`double`): material point condition **penalty factor**

## MPMParticlePointLoadCondition

In addition to the variables defined for material point conditions of type `MPMParticleBaseCondition`,
the conditions of type `MPMParticlePointLoadCondition`, which implement the action of a point load,
define the following variables

* `POINT_LOAD` (`array_1d<double,3>`): material point condition **load**
* `MPC_DELTA_DISPLACEMENT` (`array_1d<double,3>`): material point condition **displacement increment**
