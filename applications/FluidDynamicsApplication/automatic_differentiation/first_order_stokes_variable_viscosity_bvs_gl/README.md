# Stokes Formulation to be used as first-order Elements, via BVS stabilizaton and GL momentum equation form 

## ELEMENT DESCRIPTION:
Implements the element formulation found in paper:

    Felipe Galarce, Douglas R.Q. Pacheco,
    Fully consistent lowest-order finite element methods for generalised Stokes flows with variable viscosity,
    Computers & Mathematics with Applications,
    Volume 188,
    2025,
    Pages 40-49,
    ISSN 0898-1221,
    https://doi.org/10.1016/j.camwa.2025.03.013.

Specifically, it uses the form in Equation 11 of the paper. That means that the momentum equation in its Generalized Laplacian (GL) form and the stabilization used is the Boundary Vorticity Stabilisation (BVS).

## FILES DESCRIPTION:
This element only needs the element files, no new conditions. Therefore, the element's .h file is already present in the custom_elements directory and besides this file you will find the template to use for the .cpp file.
This template needs the substitution of the RHS and LHS definitions in  the specified places. In order to obtain the full file, run the template through the KratosFECompiler software, using the fe_definition.json file as input.
In the python_scripts directory, you will also find the file first_order_stokes_variable_viscosity_solver.py, which implements the solver to be used in order to run this element. One can just copy the default settings to the Projectparameters.json file to simulate.