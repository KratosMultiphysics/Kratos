//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

#pragma once

// System includes

// External includes

// Project includes
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class TrilinosMonotonicityPreservingSolver
 * @ingroup TrilinosApplication
 * @brief A linear solver wrapper (Trilinos compatible) that ensures monotonicity by modifying the system matrix and right-hand side.
 * @details This class wraps another linear solver (the "inner solver").
 * Before calling the inner solver's Solve() method, this class pre-processes the system matrix 'A' and the right-hand side 'B' within the ProvideAdditionalData() method.
 * The pre-processing ensures that the matrix 'A' has non-positive off-diagonal entries.
 * It iterates over all entries A(i,j) of the matrix. If a positive off-diagonal
 * entry is found (A(i,j) > 0 with i != j), this positive value is:
 * 1. Subtracted from the off-diagonal entries A(i,j) and A(j,i).
 * 2. Added to the diagonal entries A(i,i) and A(j,j).
 * To keep the system solution unchanged, the right-hand side vector 'B' is also
 * modified based on the current solution values (dofs_values) associated with
 * the degrees of freedom:
 * B[i] += A(i,j)_original * (dofs_values[j] - dofs_values[i])
 * B[j] += A(i,j)_original * (dofs_values[i] - dofs_values[j])
 * This transformation enforces a condition required for monotonicity (Discrete Maximum Principle) often needed in transport or convection-diffusion problems, without altering the final solution of the linear system. The modified system is then passed to the inner solver.
 * The standard, linear conforming FEM, while robust, may not guarantee the DMP on general or distorted meshes for anisotropic diffusion. The DMP is a desirable qualitative property stating that if the source term is non-negative, the solution should attain its minimum on the boundary and be everywhere non-negative.
 * A key requirement for a scheme to satisfy the DMP is for the LHS matrix A to be an M-matrix. An M-matrix has non-positive off-diagonal entries and a non-negative inverse. For standard FEM, the off-diagonal entries of the LHS matrix are generally not guaranteed to be non-positive, leading to a loss of monotonicity.
 * @author Daniel Diex
 * @tparam TSparseSpaceType The sparse space definition
 * @tparam TDenseSpaceType The dense space definition
 * @tparam TReordererType The reorder considered
 */
template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class TrilinosMonotonicityPreservingSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType,  TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosMonotonicityPreservingSolver
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosMonotonicityPreservingSolver);

    /// Definition of the base type
    using BaseType = LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>;

    /// The definition of the spaces (sparse matrix)
    using SparseMatrixType = typename TSparseSpaceType::MatrixType;

    /// The definition of the spaces (vector)
    using VectorType = typename TSparseSpaceType::VectorType;

    /// The definition of the spaces (dense matrix)
    using DenseMatrixType = typename TDenseSpaceType::MatrixType;

    /// The definition of the linear solver factory type
    using LinearSolverFactoryType = LinearSolverFactory<TSparseSpaceType,TDenseSpaceType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TrilinosMonotonicityPreservingSolver()
    {
    }

    /**
     * @brief Constructor without parameters
     * @param pLinearSolver The linear solver to be scaled
     */
    TrilinosMonotonicityPreservingSolver(
        typename BaseType::Pointer pLinearSolver
        ) : BaseType (),
            mpLinearSolver(pLinearSolver)
    {
    }

    /**
     * @brief Constructor with parameters
     * @param ThisParameters The configuration parameters of the linear solver
     */
    TrilinosMonotonicityPreservingSolver(Parameters ThisParameters)
        : BaseType ()
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
        {
            "solver_type"                    : "monotonicity_preserving",
            "inner_solver_settings"          : {
                "preconditioner_type"            : "amg",
                "solver_type"                    : "amgcl",
                "smoother_type"                  : "ilu0",
                "krylov_type"                    : "lgmres",
                "coarsening_type"                : "aggregation",
                "max_iteration"                  : 100,
                "provide_coordinates"            : false,
                "gmres_krylov_space_dimension"   : 100,
                "verbosity"                      : 1,
                "tolerance"                      : 1e-6,
                "scaling"                        : false,
                "block_size"                     : 1,
                "use_block_matrices_if_possible" : true,
                "coarse_enough"                  : 1000,
                "max_levels"                     : -1,
                "pre_sweeps"                     : 1,
                "post_sweeps"                    : 1,
                "use_gpgpu"                      : false
            }

        }  )" );

        // Now validate agains defaults -- this also ensures no type mismatch
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        mpLinearSolver = LinearSolverFactoryType().Create(ThisParameters["inner_solver_settings"]);

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    TrilinosMonotonicityPreservingSolver(const TrilinosMonotonicityPreservingSolver& Other) : BaseType(Other) {}


    /// Destructor.
    ~TrilinosMonotonicityPreservingSolver() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    TrilinosMonotonicityPreservingSolver& operator=(const TrilinosMonotonicityPreservingSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{
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
    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
    {
        std::vector<int> global_ids(rdof_set.size());
        std::vector<double> dofs_values(rdof_set.size());

        IndexType i = 0;
        for (const auto& dof : rdof_set) {
            const int global_id = dof.EquationId();
            global_ids[i] = global_id;
            dofs_values[i] = dof.GetSolutionStepValue();
            ++i;
        }

        Epetra_Map localmap(-1, global_ids.size(), global_ids.data(), 0, rA.Comm());
        Epetra_Vector dofs_local(Copy, localmap, dofs_values.data());

        Epetra_Import dirichlet_importer(rA.ColMap(), dofs_local.Map());

        // defining a temporary vector to gather all of the values needed
        Epetra_Vector dofs(rA.ColMap());

        // importing in the new temp vector the values
        int ierr = dofs.Import(dofs_local, dirichlet_importer, Insert);
        if (ierr != 0)
            KRATOS_ERROR << "Epetra failure found";

        for (int i = 0; i < rA.NumMyRows(); i++) {
            int numEntries; // number of non-zero entries
            double* vals;   // row non-zero values
            int* cols;      // column indices of row non-zero values
            rA.ExtractMyRowView(i, numEntries, vals, cols);

            for (int j = 0; j < numEntries; j++) {
                const double value = vals[j];
                if (value > 0.0) {
                    const int row_gid = rA.RowMap().GID(i);
                    const int col_gid = rA.ColMap().GID(cols[j]);
                    if (row_gid > col_gid) {
                        Matrix LHS(2, 2);
                        LHS(0, 0) = value;
                        LHS(0, 1) = -value;
                        LHS(1, 0) = -value;
                        LHS(1, 1) = value;
                        Vector dofs_sol(2);
                        int row_lid = localmap.LID(row_gid);
                        dofs_sol[0] = dofs_local[row_lid];
                        dofs_sol[1] = dofs[cols[j]];
                        Vector RHS = -prod(LHS, dofs_sol);
                        Element::EquationIdVectorType equations_ids = {static_cast<std::size_t>(row_gid), static_cast<std::size_t>(col_gid)};
                        TSparseSpaceType::AssembleLHS(rA, LHS, equations_ids);
                        TSparseSpaceType::AssembleRHS(rB, RHS, equations_ids);

                    }
                }
            }
        }

        rA.GlobalAssemble();
        rB.GlobalAssemble();

        if (mpLinearSolver->AdditionalPhysicalDataIsNeeded()) {
            mpLinearSolver->ProvideAdditionalData(rA,rX,rB,rdof_set,r_model_part);
        }
    }

    void InitializeSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mpLinearSolver->InitializeSolutionStep(rA,rX,rB);
    }

    /** This function is designed to be called at the end of the solve step.
     * for example this is the place to remove any data that we do not want to save for later
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void FinalizeSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mpLinearSolver->FinalizeSolutionStep(rA,rX,rB);
    }

    /** This function is designed to clean up all internal data in the solver.
     * Clear is designed to leave the solver object as if newly created.
     * After a clear a new Initialize is needed
     */
    void Clear() override
    {
        mpLinearSolver->Clear();
    }

    /**
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        return mpLinearSolver->Solve(rA,rX,rB);
    }



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Composite Linear Solver. Uses internally the following linear solver " << mpLinearSolver->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer mpLinearSolver;

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class TrilinosMonotonicityPreservingSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::istream& operator >> (std::istream& IStream,
                                  TrilinosMonotonicityPreservingSolver<TSparseSpaceType, TDenseSpaceType,
                                  TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const TrilinosMonotonicityPreservingSolver<TSparseSpaceType, TDenseSpaceType,
                                  TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.