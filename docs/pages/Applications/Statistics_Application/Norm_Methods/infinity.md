---
title: Infinity Norm
keywords: norm, infinity, statistics
tags: [infinity.md]
sidebar: statistics_application
summary: 
---
This is infinity (i.e. max norm) which is available for `Array 3D`, `Vector` and `Matrix` type.

For an `Array 3D` variable following equation is used.

<p align="center">$$ ||\underline{V}||_\infty = \max\left \lbrace {|V_X|, |V_Y|, |V_Z|}\right \rbrace$$</p>

For a `Vector` variable following equation is used.

<p align="center">$$ ||\underline{V}||_\infty = \max\left \lbrace {|V_0|, |V_1|, |V_2|,...,{V_{n-1}}}\right \rbrace$$</p>

For a `Matrix` variable following equation is used.

<p align="center">$$ ||\underline{A}||_\infty = \max_{0 \leq i < n}\left( \sum_{j=0}^{n-1} |a_{ij}|\right )$$</p>

