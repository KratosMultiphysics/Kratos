@defgroup adjoints Adjoint System

@brief System responsible for computing sensitivities.

@details

@section adjoints-purpose Purpose

The adjoint system is responsible for computing the sensitivities of
@ref AdjointResponseFunction "objective functions" with respect to
design variables by solving an adjoint problem on top of the primal
(forward/physical) one.

@section adjoints-definition Definition

This section provides a basic overview of the mathematical terminology,
definitions and derivation of the adjoint problem.

Given a minimization problem
@f[
    \begin{align}
        j(s, u^n .. u^b, \dot u^n .. \dot u^b, \ddot u^n .. \ddot u^b) \rightarrow min!&    \\
        \text{subject to} \quad&                                                            \\
        r^k(s, u^k, \dot u^k, \ddot u) = 0 \quad &\forall \quad k \in \{ b .. n \}          \\
        g^k(u^k .. u^{k-b}, \dot u^k .. \dot u^{k-b}, \ddot u^k .. \ddot u^{k-b}) = 0 \quad &\forall \quad k \in \{ b .. n \} \\
        h^k(\dot u^k .. \dot u^{k-b}, \ddot u^k .. \ddot u^{k-b}) = 0 \quad &\forall \quad k \in \{ b .. n \} \\
        i^k(u^b .. u^0, \dot u^b .. \dot u^0, \ddot u^b .. \ddot u^0) = 0 \quad &\forall \quad k \in \{ 0 .. b-1 \}
    \end{align}
@f]
where
- @f$j@f$ is the objective function to be minimized,
- @f$s@f$ is the set of design variables,
- @f$u@f$ is the set of state variables defined by
- @f$r@f$ the physical residual,
- @f$b@f$ is the temporal finite-differencing scheme's order,
- @f$g^k@f$ is the velocity discretization at time step @f$k@f$,
- @f$h^k@f$ is the acceleration discretization at time step @f$k@f$, and
- @f$i@f$ is the set of initial conditions,

the task is to compute @f$\frac{dj}{ds}@f$.

Unfortunately, this involves the Lagrangian of @f$j@f$.
@f[
    L = j + \sum_{k=b}^n \left[ \kappa^k r^k \right] + \sum_{k=b}^n \left[ \lambda^k g^k \right] + \sum_{k=b}^n \left[ \mu^k h^k \right] + \sum_{k=0}^{b-1} \left[ \nu^k i^k \right]
@f]
where @f$\kappa^k@f$, @f$\lambda^k@f$, @f$\mu^k@f$ and @f$\nu^k@f$ are all Lagrange multipliers for their respective
constraint equations they enforce.

@section adjoints-structure Structure

Assembling, solving and postprocessing the adjoint problem is slightly
more involved than the primal problem, but on the hand is at least linear
since we target discrete adjoints.

@subsection adjoints-structure-assembly Assembly

@subsection adjoints-structure-solution Solution

@subsection adjoints-structure-post-processing Postprocessing
