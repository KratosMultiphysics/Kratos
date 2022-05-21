---
title: Face angle
keywords: 
tags: [Face_angle.md]
sidebar: shape_optimization_application
summary: 
---
    Face angle response function.
    It aggregates the deviation of the face angles of all surface conditions using sqrt(sum(g_i)),
    where g_i are the condition wise violations - feasible conditions do not contribute

    It requires surface conditions in the modelpart, since they are used to compute the face orientation.
    Ideally the design surface model part is used.