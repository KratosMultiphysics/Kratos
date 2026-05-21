---
title: Local Axes Orientation
keywords: 
tags: [Local-axes-orientation.md]
sidebar: kratos_for_users
summary: 
---

The purpose of the current page is to describe the convention employed in the definition of local axes within _Kratos_.
This is designed in particular for the "structural" case, consisting in Beams, Shells and Solid elements, however the same conventions shall be employed in other fields when applicable.

The following variables, of type `array_1d<double,3>` are defined within the _KratosCore_, and shall be used in naming the axis of choice.

* `LOCAL_AXIS_1`
* `LOCAL_AXIS_2`
* `LOCAL_AXIS_3`

## BEAM case

For the case of a Beam, the axis `LOCAL_AXIS_1` is chosen as the unit vector tangent to the beam axis, oriented following the beam natural numbering (for a 2 node beam, from node 0 towards node 1).

`LOCAL_AXIS_2` is expected to be orthogonal to the direction identified by `LOCAL_AXIS_1`. This axis can either be provided by the user or computed automatically. `LOCAL_AXIS_3` is then computed to form an orthonormal basis with the first 2.

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Local_axes_orientation/local_axis_convention.png)

### CASE 1 - Axis is user prescribed

In this case, `LOCAL_AXIS_2` is assumed to be _approximately_ orthogonal to the beam axis. It will be made orthogonal and normalized as a very first step.
It shall be assigned to the element employing the `GetValue`/`SetValue` method, so that it is possible to query whether it was prescribed or not by employing the function `Has()`.

### CASE 2 - Axis not prescribed 

By default, `LOCAL_AXIS_2` is initialized to the global y-axis when `abs(LOCAL_AXIS_1[2]) > cos(10)`, otherwise it is initialized to the global z-axis. It will be made orthogonal and normalized as a very first step.
 
### Orthonormalization of `LOCAL_AXIS_2`

After `LOCAL_AXIS_2` has been initialized by either case 1 or case 2, it is projected onto the plane defined by `LOCAL_AXIS_1` and then normalized.

