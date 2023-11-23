//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos {


/** @brief A solver similar to @ref AMGCLSolver but with a direct interface to AMGCL.
 *  @details This class has 2 main differences compared to @ref AMGCLSolver:
 *           - the AMGCL solver instance is stored and does not get reconstructed
 *             at every call to @ref Solve, leading to better performance in cases
 *             when the solver is called repeatedly for the same system.
 *           - if compiled with GPU support, a single precision GPU backend is
 *             supported in addition to the standard double precision backend.
 *
 *           Default Parameters:
 *           @code
 *           {
 *              "solver_type" : "amgcl_raw",
 *              "verbosity" : 0,
 *              "tolerance" : 1e-6,
 *              "gpgpu_backend" : "",
 *              "amgcl_settings" : {
 *                  "precond" : {
 *                      "class" : "amg",
 *                      "relax" : {
 *                          "type" : "ilu0"
 *                      },
 *                      "coarsening" : {
 *                          "type" : "aggregation",
 *                          "aggr" : {
 *                              "eps_strong" : 0.08,
 *                              "block_size" : 1
 *                          }
 *                      },
 *                      "coarse_enough" : 333,
 *                      "npre" : 1,
 *                      "npost" : 1
 *                  },
 *                  "solver" : {
 *                      "type" : "cg",
 *                      "maxiter" : 555,
 *                      "tol" : 1e-6
 *                  }
 *              }
 *          }
 *          @endcode
 *
 *          Parameters:
 *          - "solver_type": the name of this class referenced from the JSON interface
 *          - "verbosity": level of information printed from the solver. A higher value
 *                         will result in more verbose output. Level 4 and above will
 *                         print the system matrices and throw and terminate the program.
 *          - "tolerance": relative tolerance to check convergence after solving. Note that
 *                         this setting does not get passed on to AMGCL. Set the tolerance
 *                         directly in @a "amgcl_settings" to control AMGCL's tolerance.
 *          - "gpgpu_backend": [@a "" (default), @a "double", or "float"] control what
 *                             backend AMGCL should use. Available options are:
 *                             - @a "": use the built-in double precision backend that runs
 *                                      on the CPU. This is the default setting (no GPU).
 *                             - @a "double": use a double precision GPU backend. The exact
 *                                            type of backend depends on what Kratos was
 *                                            compiled with (controlled by @a AMGCL_GPGPU_BACKEND).
 *                             - @a "float": use a single precision GPU backend. The exact
 *                                           type of backend depends on what Kratos was
 *                                           compiled with (controlled by @a AMGCL_GPGPU_BACKEND).
 *          - "amgcl_settings": subparameter tree passed on directly to AMGCL. See AMGCL's
 *                              documentation for available options.
 */
template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer = Reorderer<TSparseSpace, TDenseSpace> >
class AMGCLRawSolver final : public LinearSolver<TSparseSpace,
                                                 TDenseSpace,
                                                 TReorderer>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(AMGCLRawSolver);

    using Base =  LinearSolver<TSparseSpace, TDenseSpace, TReorderer>;

    using SparseMatrix = typename TSparseSpace::MatrixType;

    using Vector = typename TSparseSpace::VectorType;

    using DenseMatrix = typename TDenseSpace::MatrixType;

    AMGCLRawSolver(Parameters rParameters);

    AMGCLRawSolver(AMGCLRawSolver&&) noexcept = default;

    ~AMGCLRawSolver() override;

    /// @copydoc LinearSolver::Solve
    bool Solve(SparseMatrix& rA, Vector& rX, Vector& rB) override;

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrix& rA, DenseMatrix& rX, DenseMatrix& rB) override;

    void  PrintInfo(std::ostream& rOStream) const override;

    void  PrintData(std::ostream& rOStream) const override;

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalData (
        SparseMatrix& rA,
        Vector& rX,
        Vector& rB,
        typename ModelPart::DofsArrayType& rDofSet,
        ModelPart& rModelPart
    ) override;

    static Parameters GetDefaultParameters();

private:
    AMGCLRawSolver(const AMGCLRawSolver& Other) = delete;

    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class AMGCLRawSolver


template<class TSparseSpace, class TDenseSpace,class TReorderer>
std::istream& operator >> (std::istream& rIStream, AMGCLRawSolver< TSparseSpace,
                           TDenseSpace, TReorderer>& rThis)
{
    return rIStream;
}


template<class TSparseSpace, class TDenseSpace, class TReorderer>
std::ostream& operator << (std::ostream& rOStream,
                           const AMGCLRawSolver<TSparseSpace,TDenseSpace,TReorderer>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos
