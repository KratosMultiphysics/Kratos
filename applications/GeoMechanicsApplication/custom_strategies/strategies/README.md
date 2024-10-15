# GeoMechanicsNewtonRaphsonErosionProcessStrategy class.

## Setup

This class implements a Newton-Raphson process for the erosion piping process. Details can be found in

## Algotrithm

The algorithm is implemented in FinalizeSolutionStep function that loops over all piping elements until the erosion
groth is finished or all piping elements are checked.

### FinalizeSolutionStep

![FinalizeSolutionStep.svg](FinalizeSolutionStep.svg)

This function calls other functions that are described below.

### GetPipingElements function

This function returns a list of all piping elements. To create the list, all elements are checked for  
PIPE_START_ELEMENT property that holds Id of the piping starting element.
Then elements on the list are sorted based on X-coordinate of the element begin.
The current limitation is the following. The starting piping element has to have minimum or maximum value of the X
coordinate. This is needed to find a direction for the piping growth.

### InitialiseNumActivePipeElements function

It checks the list of the elements and returns a number of elements with PIPE_ACTIVE property At the beginning the
number is zero.

### CalculateMaxPipeHeight function

It calculates a maximum diameter using CalculateParticleDiameter function then multiplies it by a magic factor of 100 to
get the height value.

### CalculateParticleDiameter

This is SteadyStatePwPipingElement class function. In case of PIPE_MODIFIED_D
$$ diameter = 2.08\times 10^{-4} * \Bigg(\frac{PIPE\_D\_70}{2.08\times 10^{-4}}\Bigg)^{0.4}$$

otherwise, $$diameter = PIPE\_D\_70$$

### check_pipe_equilibrium function

This function performs iterative process with a prescribed maximum number of iterations. For each iteration the
modelling is
performed. Then a pipe height is calculated for each open piping element.

![check_pipe_equilibrium.svg](check_pipe_equilibrium.svg)

### Recalculate function

This function performs the modeling using prescribed GeoMechanicsNewtonRaphsonStrategy class.

![Recalculate.svg](Recalculate.svg)

### CalculateEquilibriumPipeHeight function

This is a function from SteadyStatePwPipingElement. It calculates the equilibrium pipe height as

$$PIPE\_MODEL\_FACTOR * \pi / 3.0 * D_p * (SolidDensity / FluidDensity - 1) * \eta *
\tan(\theta + pipeSlope) / dhdx$$

where $D_p$ is the particle diameter calculated with CalculateParticleDiameter function,\
$dhdx$ is the head gradient obtained with CalculateHeadGradient function,\
Currently, the pipe is assumed horizontal $pipeSlope = 0$.

### check_status_tip_element function

This function checks the open piping elements and if their height is larger than a maximum value of the pipe height or
too small then the function changes the pipe properties: PIPE_EROSION and PIPE_ACTIVE to false.

![check_status_tip_element.svg](check_status_tip_element.svg)

### save_or_reset_pipe_heights function

Based on the grow value, this function sets PREV_PIPE_HEIGHT equal to PIPE_HEIGHT or reset PIPE_HEIGHT to
PREV_PIPE_HEIGHT value.

### Interaction with Interface classes.

SteadyStatePwPipingElement class sets JointWidth equal to PIPE_HEIGHT to call functions from Interface classes. 

