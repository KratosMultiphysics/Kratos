---
title: Mesh based packaging
keywords: 
tags: [Mesh_based_packaging.md]
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

    A class for mesh packaging response function. The mesh should contain surface conditions only.
    By default the normals of the conditions indicate the feasible side of the mesh (see setting 'feasible_in_normal_direction')