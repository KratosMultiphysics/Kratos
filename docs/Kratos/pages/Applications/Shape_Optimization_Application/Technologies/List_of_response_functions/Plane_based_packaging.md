---
title: Plane based packaging
keywords: 
tags: [Plane_based_packaging.md]
sidebar: shape_optimization_application
summary: 
---
    A base class for packaging response functions that agglomerate the nodal violations
    into a single response function.
    The agglomeration happens by summing up the square of each nodal violation.
    Nodes that are feasible do NOT contribute to the response value/gradient.
    This is why a prediction of the violation using the gradients is not possible,
    only correction of violations (e.g. from the last step) will happen.

    Derived classes need to implement the calculation of the nodal violations

        A class that defines the response function for plane-based packaging. The plane is defined by an origin point and a normal vector.
    By default the normal of the plane indicates the feasible side of the plane (see setting 'feasible_in_normal_direction')