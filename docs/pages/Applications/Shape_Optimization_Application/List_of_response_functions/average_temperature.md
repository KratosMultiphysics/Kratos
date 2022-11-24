---
title: Average temperature
keywords: 
tags: [average_temperature.md]
sidebar: shape_optimization_application
summary: 
---
## Introduction

This reponse function computes average temperature of a given model part
## Formulation

Following formulation is used.
<p align="center">$$ J   = \frac{1}{\sum_{\forall \Omega_i \in \Omega} 1}\sum_{\forall \Omega_i \in \Omega} T_i $$</p>

## Source

Location: ["applications/ConvectionDiffusionApplication/custom_response_functions/local_temperature_average_response_function.h"](https://github.com/KratosMultiphysics/Kratos/blob/shapeopt/kreisselmeier_aggregation/applications/ConvectionDiffusionApplication/custom_response_functions/local_temperature_average_response_function.h)