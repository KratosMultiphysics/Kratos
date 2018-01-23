## External Solvers Application

External Solvers application is the place where solvers that do not directly to the Kratos project are localed. 
This application provides different interfaces to these solvers which make use of third party libraries or codes not directly integrated in Kratos.

In the following list you can find which solvers are currently in kratos and links to the library projects or algortihms in which are bases:

### Direct Solvers
* __SuperLUSolver__: Based on [SuperLu](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) project.
* __PastixSolver__: Based on [Pastix](http://pastix.gforge.inria.fr/files/README-txt.html) project.

### Iterative Solvers:
* __SuperLUIterativeSolver__: Based on [SuperLu](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) project.
* __GMRESSolver__: Based in the [GMRES](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method) method.

### Eigen Solvers:
* __FEASTSolver__: Based on [FEAST](http://www.feast-solver.org/) project.
