//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//

#pragma once

// System includes
#include <type_traits>
#include <vector>

// External includes
#include <Eigen/SparseLU>

// Project includes
#include "future/linear_solvers/direct_solver.h"

namespace Kratos::Future
{
///@name Kratos Classes
///@{

/**
 * @class EigenSparseLUSolver
 * @ingroup KratosCore
 * @brief Direct linear solver based on the SparseLU factorization of the Eigen library.
 * @details The LHS matrix of the linear system is converted to an Eigen sparse matrix,
 *          factorized with Eigen::SparseLU (with column pivoting) and the solution is
 *          obtained by back-substitution. The factorization is performed inside
 *          PerformSolutionStep (as in SkylineLUFactorizationSolver), so the solver is
 *          self-contained: it does not matter which of the Initialize hooks is called
 *          before the solve. A singular LHS matrix raises a Kratos error.
 * @tparam TLinearAlgebra The version of the linear algebra to be used
 * @author Thomas Oberbichler
 */
template<class TLinearAlgebra>
class EigenSparseLUSolver
    : public DirectSolver<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EigenSparseLUSolver
    KRATOS_CLASS_POINTER_DEFINITION(EigenSparseLUSolver);

    /// Base class definition
    using BaseType = DirectSolver<TLinearAlgebra>;

    /// Linear system type definition
    using LinearSystemType = LinearSystem<TLinearAlgebra>;

    /// Sparse matrix type definition
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Vector type definition
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Type definition for data
    using DataType = typename VectorType::DataType;

    /// Type definition for index
    using IndexType = typename VectorType::IndexType;

    /// Dense vector tag type definition
    using DenseVectorTag = typename LinearSystemTags::DenseVectorTag;

    /// Sparse matrix tag type definition
    using SparseMatrixTag = typename LinearSystemTags::SparseMatrixTag;

    /// Eigen index type (Eigen requires a signed storage index, unlike the CsrMatrix)
    using EigenIndexType = std::make_signed_t<IndexType>;

    /// Eigen sparse matrix type (column major)
    using EigenSparseMatrixType = Eigen::SparseMatrix<DataType, Eigen::ColMajor, EigenIndexType>;

    /// Eigen dense vector type
    using EigenVectorType = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor
     * @param Settings The configuration parameters
     */
    EigenSparseLUSolver(Parameters Settings = Parameters(R"({})"))
        : BaseType()
    {
        // Validate and assign default parameters
        Settings.ValidateAndAssignDefaults(GetDefaultParameters());

        // Assign the validated settings
        this->AssignSettings(Settings);
    }

    /// Copy constructor
    EigenSparseLUSolver(const EigenSparseLUSolver& Other) = delete;

    /// Destructor
    ~EigenSparseLUSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    EigenSparseLUSolver& operator=(const EigenSparseLUSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Factorizes the LHS matrix and solves the linear system.
     * @param rLinearSystem The linear system to be solved
     * @return @p true after a successful factorization and back-substitution
     */
    bool PerformSolutionStep(LinearSystemType& rLinearSystem) override
    {
        if (!rLinearSystem.IsConsistent(SparseMatrixTag::LHS, DenseVectorTag::RHS, DenseVectorTag::Dx)) {
            return false;
        }

        auto& r_lin_op = *rLinearSystem.pGetLinearOperator(SparseMatrixTag::LHS);
        auto& rX = *rLinearSystem.pGetVector(DenseVectorTag::Dx);
        auto& rB = *rLinearSystem.pGetVector(DenseVectorTag::RHS);

        KRATOS_ERROR_IF(r_lin_op.IsMatrixFree()) << "EigenSparseLUSolver cannot be used with matrix-free linear operators." << std::endl;
        const auto& r_A = r_lin_op.GetMatrix();

        // Factorize
        const EigenSparseMatrixType eigen_A = ToEigenSparseMatrix(r_A);
        mEigenSparseLUSolver.compute(eigen_A);
        KRATOS_ERROR_IF(mEigenSparseLUSolver.info() != Eigen::Success)
            << "EigenSparseLUSolver: the LHS matrix is singular and the factorization failed." << std::endl;

        // Back-substitute
        const EigenVectorType eigen_B = ToEigenVector(rB);
        const EigenVectorType eigen_X = mEigenSparseLUSolver.solve(eigen_B);
        KRATOS_ERROR_IF(mEigenSparseLUSolver.info() != Eigen::Success)
            << "EigenSparseLUSolver: the back-substitution failed." << std::endl;

        for (std::size_t i = 0; i < rX.size(); ++i) {
            rX[i] = eigen_X[i];
        }

        return true;
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters(R"({
            "solver_type" : "eigen_sparse_lu"
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
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

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Eigen sparse LU solver finished.";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Eigen::SparseLU<EigenSparseMatrixType> mEigenSparseLUSolver; /// The Eigen factorization

    ///@}
    ///@name Private Operations
    ///@{

    /// Converts a CsrMatrix to an Eigen column-major sparse matrix
    static EigenSparseMatrixType ToEigenSparseMatrix(const MatrixType& rA)
    {
        const auto& r_row_ptr = rA.index1_data();
        const auto& r_col_idx = rA.index2_data();
        const auto& r_vals = rA.value_data();

        std::vector<Eigen::Triplet<DataType, EigenIndexType>> triplets;
        triplets.reserve(r_row_ptr[rA.size1()]);
        for (std::size_t i = 0; i < rA.size1(); ++i) {
            for (IndexType k = r_row_ptr[i]; k < r_row_ptr[i + 1]; ++k) {
                triplets.emplace_back(static_cast<EigenIndexType>(i), static_cast<EigenIndexType>(r_col_idx[k]), r_vals[k]);
            }
        }

        EigenSparseMatrixType result(rA.size1(), rA.size2());
        result.setFromTriplets(triplets.begin(), triplets.end());
        return result;
    }

    /// Converts a SystemVector to an Eigen dense vector
    static EigenVectorType ToEigenVector(const VectorType& rV)
    {
        EigenVectorType result(rV.size());
        for (std::size_t i = 0; i < rV.size(); ++i) {
            result[i] = rV[i];
        }
        return result;
    }

    ///@}
}; // Class EigenSparseLUSolver

///@}

/// input stream function
template<class TLinearAlgebra>
inline std::istream& operator >> (
    std::istream& rIStream,
    EigenSparseLUSolver<TLinearAlgebra>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TLinearAlgebra>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const EigenSparseLUSolver<TLinearAlgebra>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos::Future.
