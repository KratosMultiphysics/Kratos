# 1D wave propagation in a drained soil column, with constant mass and damping

This tests verifies that the ResidualBasedBlockBuilderAnderSolverWithMassAndDamping is used and gives results.

## Setup

This test is a copy of the test 1D wave propagation in a drained soil column. The difference is that mass and damping matrices only are setup once. For that block_builder and prebuild_dynamics are both set to true in ProjectParameters.json.
[See 1D wave propagation in a drained soil column](../test_1d_wave_prop_drained_soil/README.md)

## Assertions

The test asserts that the vertical speed of node 41 ( at half the height of the column ) matches the predicted value at specified times and asserts that the intended block builder and solver is used.
