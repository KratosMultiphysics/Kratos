---
title: the FEAST Solver in Kratos?
keywords: 
tags: [How-to-use-the-FEAST-Solver-in-Kratos.md]
sidebar: kratos_for_developers
summary: 
---

The [FEAST Solver](http://www.ecs.umass.edu/~polizzi/feast/) is used for solving an eigenvalue problem. It is based on a new algorithm which deviates fundamentally from the _Krylov_ subspace based techniques (_Arnoldi_ and _Lanczos_ algorithms), Davidson-Jacobi techniques or other traditional subspace iteration techniques

   > E. Polizzi, Density-Matrix-Based Algorithms for Solving Eigenvalue Problems, Phys. Rev. B. Vol. 79, 115112 (2009)

A documentation about how FEAST 4.0 works internally is available [here](https://arxiv.org/abs/2002.04807).

### Compilation of the FEAST Library

In order to compile the FEAST Library, the LinearSolversApplication with Intel MKL is required. The following lines must be added to the configure.sh file:

```console
-DLINEAR_SOLVERS_APPLICATION=ON                                 \
-DUSE_EIGEN_MKL=ON                                              \
-DUSE_EIGEN_FEAST=ON                                            \
```

Once these lines are included, it is just needed the execution of the `configure.sh` file. However, only compilation on Linux is supported at the moment. For more information, please consult the [LinearSolversApplication readme](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/LinearSolversApplication/README.md).

### How to use FEAST 4.0
FEAST 4.0 can be used for solving an eigenvalue problem in a structural mechanics context by specifying `eigen_value` as `solver_type` in the problem parameters `json`. The default input parameters for solving a real-valued problem using FEAST are:

```json
"eigensolver_settings":{
	    "solver_type" : "feast",
            "symmetric" : true,
            "number_of_eigenvalues" : 0,
            "e_min": 0.0,
            "e_max": 0.2,
            "search_lowest_eigenvalues" : false,
            "search_highest_eigenvalues" : false,
            "sort_eigenvalues" : false,
            "sort_order" : "sr",
            "subspace_size" : 0,
            "max_iteration" : 20,
            "tolerance" : 1e-12,
            "echo_level" : 0
	},
```

The settings controlled with these parameters are:
- `solver_type`: If FEAST is used. Another possibility is `eigen_eigensystem`. The following settings are only applicable, if `feast` is chosen.
- `symmetric`: True, if the input matrices are symmetric.
- `number_of_eigenvalues`: Specifies the number of eigenvalues to be found.
- `e_min` and `e_max`: Defines the search interval, where eigenvalues are found as `e_min` < &lambda; < `e_max`. Note, that for strucural mechanics, the relation between eigenvalue &lambda; and angular eigenfrequency &omega; is &omega;²=&lambda;
- `search_lowest_eigenvalues` and `search_highest_eigenvalues`: Either search for the highest or lowest eigenvalues in the specified region. Only available for real, symmetric problems and in combination with a defined `number_of_eigenvalues`.
- `sort_eigenvalues`: True, if the resulting eigenvalues should be sorted-
- `sort_order`: Defines the sort order. Possible choices are `sr` (smallest real part), `sm` (smallest magnitude), `lr` (largest real part), `lm` (largest magnitude).
- `subspace_size`: Manually defines the size of the subspace used to find eigenvalues. The FEAST documentation recommends a size of 1.5 * numer of eigenvalues to be found. This value is overwritten, if `numer_of_eigenvalues` is specified.
- `max_iteration`: Maximum number of FEAST iterations.
- `tolerance`: Defines the FEAST solver tolerance.

### Printing Results

In the StructuralMechanicsApplication there is a process for of postprocessing eigenvalues. It animates the results and writtes all EigenValues into one file for easier postprocessing.

**PostprocessEigenvaluesProcess**

This process is written in C++ but can be called by adding the corresponding python process to the `list_other_processes` in the `ProjectParameters.json`

```json
"list_other_processes" : [{
    "python_module"   : "postprocess_eigenvalues_process",
    "kratos_module"   : "KratosMultiphysics.StructuralMechanicsApplication",
    "help"                  : "This process postprocces the eigen values for GiD",
    "process_name"          : "PostProcessEigenvaluesProcess",
    "Parameters"            : { 
        "result_file_name" : "Structure",
        "result_file_format_use_ascii" : false,
        "computing_model_part_name"   : "computing_domain",
        "animation_steps"   :  20,
        "list_of_result_variables" : ["DISPLACEMENT"],
        "label_type" : "frequency"
    }
}]
```
Certain Parameters can be chosen (but don't have to, the settings above are the defaults. If the defaults are to be used, simply omit the corresponding parameter):
- `result_file_name` : Name of the result file, will be followed by "_EigenResults"
- `result_file_format_use_ascii` : Usually the results are printed in binary format, this setting is only used for testing and debugging
- `computing_model_part_name`   : Name of the Computing Model Part that is used as input for the solver
- `animation_steps`   :  Number of steps used to animate the eigenvectors
- `list_of_result_variables` : Variables that are used for the computation of the eigen-solution. These are used for the animation. In Structural Analysis this is the `DISPLACEMENT`
- `label_type` : Output the results in either `angular_frequency` [rad/s] or `frequency` [Hz]
