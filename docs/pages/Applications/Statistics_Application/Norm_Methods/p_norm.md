---
title: P Norm
keywords: norm, p_norm, statistics
tags: [p_norm.md]
sidebar: statistics_application
summary: 
---
P norm is applicable for `Array 3D`, `Vector` and `Matrix` variables. The p-norm is used as `pnorm_p` where, `p` represents the value to be used as `p` in the p-norm equation given below. `p` should be greater than or equal to 1.0.

For `Array 3D` and `Vector`:

<p align="center">$$ ||\underline{V}||_{p} = \left(\sum_{i=0}^{n-1}|v_{i}|^p \right )^{1/p}$$</p>

For `Matrix`:

<p align="center">$$ ||\underline{A}||_{p} = \left(\sum_{i=0}^{n-1}\sum_{j=0}^{n-1}|a_{ij}|^p \right )^{1/p}$$</p>
