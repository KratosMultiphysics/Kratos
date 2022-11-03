---
title: Kreisselmeier Aggregated Stress
keywords: 
tags: [kreisselmeier_aggregated_stress.md]
sidebar: shape_optimization_application
summary: 
---

## Settings

<p align="center">
    <img src="../images/objectives_kmas.png" alt="Kreisselmeier aggregated stress advanced settings"/>
</p>
<p align="center">Figure 1: Kreisselmeier aggregated stress advanced settings</p>

Figure 1 illustrates `Response Settings` and `Advanced Settings` related to `Kreisselmeier Aggregated Stress` response.

The response function is computed as follows where $$\sigma$$ is the Von Mises stress:
<p align="center">$$ J   = \frac{1}{\rho}\log\left({\sum{e^{\rho \frac{\sigma}{\alpha}}}}\right)$$</p>

The coefficients are listed in the following table:

|Coefficient | Setting | Description |
|------------|---------|-------------|
|$$\rho$$    | `Aggregation penalty` | This works as a smoothning factor |
|$$\alpha$$ | `Scaling factor` | This works as a scaling factor. Too law values will result in infinity and too large values will result in less weighting in the sensitivity computation. It also accepts `max` as a keyword. Then the solver will compute the max Von Mises stress in the initial geometry and then use it as a scaling factor. |

`Echo level` is used to determine the verbosity of the output created by this response. Higher the value, higher the verbosity.

### Caution:

User needs to be cautious about using this response as an objective along with `Mass` as another constraint or vice versa. Because these two objectives are contradicting each other, therfore if the constraint is violated, then there will be weighting of sensitivities computed from `Mass` and `Kreisselmeier Aggregated Stress`. This weighting can be influenced by `Scaling Factor`. If a higher value is used than the actual Von Mises stresses, then it will be having a less wieght, and the optimizer will be removing mass in the case where `Mass` minimization is the objective. The usual range for this would be 20% - 50% of the max stress. Depending on the `Scaling factor`, optimizer will compute different shape updates.

