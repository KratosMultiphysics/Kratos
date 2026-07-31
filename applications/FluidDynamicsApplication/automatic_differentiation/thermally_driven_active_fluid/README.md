# Divergence Constrained Mixed Element for NonLinear Fourth-Order Thermally Driven Active Fluid Equations

## ELEMENT DESCRIPTION:
Implements the element formulation found in paper:

    Nan Zheng, Qingguang Guan, Wenlong Pei, Wenju Zhao,
    A divergence constrained mixed finite element method for thermally driven active fluid model,
    Applied Numerical Mathematics,
    Volume 227,
    2026,
    Pages 81-107,
    ISSN 0168-9274,
    https://doi.org/10.1016/j.apnum.2026.04.008.


Specifically, it uses the form in equations 3.11-3.15, using Taylor-Hood mixed elements and the Dahlquist-Liniger-Nevanlinna (DLN) time integration scheme.

## FILES DESCRIPTION:
This element only needs the element files, no new conditions. Therefore, the element's .h file is already present in the custom_elements directory and besides this file you will find the template to use for the .cpp file.
This template needs the substitution of the RHS and LHS definitions in  the specified places. In order to obtain the full file, run the template through the KratosFECompiler software, using the fe_definition.json file as input.
In the python_scripts directory, you will also find the file thermally_driven_active_fluid_solver.py, which implements the solver to be used in order to run this element. One can just copy the default settings to the Projectparameters.json file to simulate.