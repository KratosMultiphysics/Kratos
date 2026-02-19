@page p_multigrid P-Multigrid

P-Multigrid refers to a multigrid method whose coarsening strategy is based on shape functions of high order elements (quadratic elements already qualify in this case). This implementation also supports models with multifreedom constraints.

@section overview Overview

The system's entry point is the @ref Kratos::PMultigridBuilderAndSolver "PMultigridBuilderAndSolver". As the name suggests, the main responsibilities of this class are
- Allocating the left hand side matrix, as well as the right hand side vector and the solution vector.
- Assembling the left hand side matrix, right hand side vector.
- Solving the resulting linear system of equations.

Due to the current design of @ref Kratos::Scheme "Scheme" and @ref Kratos::BuilderAndSolver "BuilderAndSolver", there is a number of other tasks as well that either tie in to the primary purpose, or help in pre- and postprocessing.
- Allocation, assembly and imposition of @ref Kratos::MasterSlaveConstraint "multifreedom constraints".
- Applying Dirichlet conditions.
- Partial reassembly of the linear system.
- Computing reactions.

The other exposed family of classes from this system is @ref Kratos::MultifreedomConstraint "MultifreedomConstraint". Its main difference compared to @ref Kratos::MasterSlaveConstraint "MasterSlaveConstraint" (which it derives from) is that it only models the constraint equations and does not impose a partition of slave- and master @ref Kratos::Dof "DoFs".
@see @ref Kratos::LinearMultifreedomConstraint "LinearMultifreedomConstraint"
@see @ref Kratos::LinkConstraint "LinkConstraint"

@note The rest of the classes in this system are meant for internal use and must not be exposed to the rest of Kratos.

Although the primary purpose of this system is to exploit the structure arising from a model using high order elements,
the multigrid feature can be completely disabled to use it as a standard @ref Kratos::BuilderAndSolver "BuilderAndSolver". One reason for
doing this would be taking advantage of different constraint imposition methods, such as
@ref Kratos::AugmentedLagrangeConstraintAssembler "augmented Lagrange".

@section coarse_hierarchy Coarse Hierarchy

@note The current implementation only supports a two-grid method, since the selection of high order elements in Kratos is limited, and does not yet justify an arbitrary depth. That said, adding support for it should be possible without major interface changes, but would involve minor changes in @ref Kratos::PMultigridBuilderAndSolver "PMultigridBuilderAndSolver" and @ref Kratos::PGrid "PGrid", as well as major changes in @ref Kratos::MakePRestrictionOperator "MakePRestrictionOperator".

The root grid (i.e.: the finest level) is stored in and represented by @ref Kratos::PMultigridBuilderAndSolver "PMultigridBuilderAndSolver", while
coarse grids are represented by @ref Kratos::PGrid "PGrid" in a linked list. The reason for this difference is the additional
set of responsibilities of the root grid, namely the allocation and assembly of the finest level. Coarse grids
do not perform assembly, but construct restriction operators that they then apply on the parent grid to compute
their own system.

Another important distinction is that the coarse grids can have floating point types different than the root grid.
This can be useful when the user has access to accelerator hardware (i.e.: GPUs). Coarse grids need not solve their
own problems with high precision so they might as well use single precision floating point numbers to save VRAM,
for example.

@section constraints Constraints

Constraint assembly and imposition is extracted through @ref Kratos::ConstraintAssembler "a dedicated interface" that
currently supports @ref Kratos::MasterSlaveConstraintAssembler "master-slave elimination",
@ref Kratos::AugmentedLagrangeConstraintAssembler "augmented Lagrange", and @ref Kratos::NoOpConstraintAssembler "a dummy"
for debugging.

@subsection constraints_in_a_multigrid_setting Constraints in a Multigrid Setting

Applying the multigrid preconditioner on a problem including multifreedom constraints involves transforming state variables and their residuals between 4 spaces.

@code
                + ---------------- + ------------------ +
                |                  |                    |
                | Dependent Space  |  Independent Space |
                |       (d)        |         (i)        |
                + ---------------- + ------------------ +
+ ----------- + + ---------------- + ------------------ +
|             | |                  |                    |
| Fine Grid   | |        B       <--->        A         |
|    (q)      | |        |         |                    |
+ ----------- + + -----R | ------- + ------------------ +
|             | |        v         |                    |
| Coarse Grid | |        C       <--->        D         |
|    (p)      | |                  |                    |
+ ----------- + + ---------------- + ------------------ +
@endcode

In the graph above, the original problem @p A (finest grid with polynomial order @f$q@f$, and independent degrees-of-freedom @f$u_i^q@f$) we must solve is in the <b>top right corner</b>. This system was constructed by constraining the one in the <b>top left corner</b>, which is the dependent system on the fine grid @p B (polynomial order @f$q@f$, dependent DoFs @f$u_d^q@f$). Furthermore, we know how to construct a coarse equivalent of this problem (using the restriction operator @f$R@f$), which is in the <b> bottom left corner</b> @p C (polynomial order @f$p<q@f$, dependent DoFs @f$u_d^p@f$). However, the problem we need to solve to obtain a coarse correction to our original problem is in the <b>bottom right corner</b> @p D (polynomial order @f$p@f$, independent DoFs @f$u_i^p@f$).

While performing a multigrid iteration, residuals must be transformed from system @p A (@f$r_i^q@f$) to @p D (@f$r_i^p@f$), and state variables from @p D (@f$u^p_i@f$) to @p A (@f$u^q_i@f$). This requires transforming residuals and state variables between dependent and independent spaces. To support such transformations, @ref Kratos::ConstraintAssembler "ConstraintAssembler" exposes the following methods.
- @ref Kratos::ConstraintAssembler::ComputeDependentResidual "ComputeDependentResidual"
- @ref Kratos::ConstraintAssembler::ComputeIndependentResidual "ComputeIndependentResidual"
- @ref Kratos::ConstraintAssembler::ComputeDependentSolution "ComputeDependentSolution"
- @ref Kratos::ConstraintAssembler::ComputeIndependentSolution "ComputeIndependentSolution"

@note If the multigrid If the multigrid feature is enabled, only the top-level constraints should be imposed via exact methods (i.e. master-slave elimination). Trying to impose constraints on coarse grids exactly may render the system ill-posed due to overconstraining. This is due to the fact that coarsening only changes the degrees-of-freedom the constraints are applied on, but does not filter linearly dependent constraints (i.e. the number of constraint equation remains constant through the coarsening process).

@section linear_solvers Linear Solvers

Unlike other @ref Kratos::BuilderAndSolver "BuilderAndSolvers", @ref Kratos::PMultigridBuilderAndSolver "PMultigridBuilderAndSolver" does not use the linear
solver provided to it from the python layer. The reason is that it must construct a solver for each grid level
separately. These fall into two categories
- smoothers for the finer grids (including the root grid)
- linear solver for the coarsest grid (usually an AMG solver).

Instead of passing linear solver instances, the user must provide two sets of parameters for the two different solver categories, after which the grids take care of constructing their own instances.
