# Settlement on a dry column with uniform load
This test consists of a rectangular soil domain, 1.0m wide and 50.0m deep. The mesh consists of SmallStrainUPwDiffOrderElement2D6N elements. The water pressure is kept at 0.0 [Pa] in the entire domain, resulting in a dry column. A schematic can be found in the figure below:
![Schematic](Schematic.svg)

## Setup

The test is performed in 5 stages:
1. A K0 stage with a linear elastic model, with a Young's modulus of 1e9 [Pa] and a Poisson ratio of 0.2.
2. In the second stage, the material is switched to the abc model, and a uniform line-load of 20 [kN] is applied in the negative Y direction to the top of the model. The stage defines a time duration of 1 second, such that the time step is smaller than 0.001 days to make sure the abc model distributes stresses horizontally. The load stepping scheme is used, to apply external forces gradually (starting with steps of 10% of the initial external forces and automatically adjusting the step size according to the number of iterations and convergence behavior). Note: reset_displacement is set to true here, such that the total displacements start counting from the start of this stage.
3. In this settlement stage (100 days or 8640000 seconds), the same uniform load is applied, but the material is switched to the abc model. The time-step is > 0.001 days, such that there is no horizontal stress distribution by the abc model. 
4. A second load stage, where the uniform load is increased to 25 [kN] in the negative Y direction. The stage defines a time duration of 1 second, such that the time step is smaller than 0.001 days to make sure the abc model distributes stresses horizontally.
5. A second settlement stage, such that the total time reaches 10000 days (or 864000000 seconds). The same uniform load of 25 [kN] is applied in the negative Y direction, but the time-step is > 0.001 days, such that there is no horizontal stress distribution by the abc model.

The following common conditions hold for all stages:
  - Displacements on the bottom are fixed in all directions.
  - Displacements on the sides are fixed in the X direction.
  - The uniform loads follow the values described above, the LineLoadDiffOrderCondition2D3N condition is used to apply the load on the top of the model.
  - Water pressure is kept at 0.0 [Pa] in the entire domain.
  - Gravity is applied to the entire domain (-9.81 [m/s^2] in the negative Y direction).

## Assertions
The following assertions are made in this test:
1. The total vertical displacement at the top of the column is expected to be within 1% of the analytical solution, -3.22m, after 100 days (stage 3).
2. The total vertical displacement at the top of the column is expected to be within 1% of the analytical solution, -8.02m, after 10000 days (stage 5).

## Results

The following plot shows the evolution of settlement obtained with Kratos (coarse mesh), and compares it with the result from D-Settlement 23.2.1.  Nodes 2, 3, and 104 are located at the top edge of the soil column:

![Settlement](coarse_mesh/test_case_1_settlement_plot.svg)

The following plots show the vertical effective stress distribution ($`\sigma_{\mathrm{eff, yy}}`$) and the vertical total stress distribution ($`\sigma_{\mathrm{tot, yy}}`$) along the left edge of the soil column at certain times.

At the start of the analysis (after 0.1 days):

![Stress plot after 0.1 days](coarse_mesh/test_case_1_stress_plot_after_0.1_days.svg)

After 100.1 days:

![Stress plot after 100.1 days](coarse_mesh/test_case_1_stress_plot_after_100.1_days.svg)

After 10,000 days of settlement:

![Stress plot after 10000 days](coarse_mesh/test_case_1_stress_plot_after_10000_days.svg)
