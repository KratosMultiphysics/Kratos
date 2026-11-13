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
 * @details The LHS matrix of the linear system is exposed to Eigen through a
 *          Eigen::Map: the values are mapped in place (zero copy) and only the index
 *          arrays are copied, because the CsrMatrix index type is unsigned while Eigen
 *          requires a signed storage index. The factorization (Eigen::SparseLU, with
 *          column pivoting) is performed inside InitializeSolutionStep, so
 *          PerformSolutionStep is a cheap back-substitution that can be repeated
 *          (e.g. for multiple right-hand sides) without refactoring. A singular
 *          LHS matrix raises a Kratos error.
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

    /// Eigen sparse matrix type (row major, to match the CsrMatrix layout)
    using EigenSparseMatrixType = Eigen::SparseMatrix<DataType, Eigen::RowMajor, EigenIndexType>;

    /// Read-only Eigen view of a CsrMatrix
    using EigenSparseMapType = Eigen::Map<const EigenSparseMatrixType>;

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
     * @brief Factorizes the LHS matrix of the linear system.
     * @details The factorization is performed here (and not in PerformSolutionStep) so that
     *          the back-substitution can be performed multiple times at a low cost.
     * @param rLinearSystem The linear system to be solved
     */
    void InitializeSolutionStep(LinearSystemType& rLinearSystem) override
    {
        if (!rLinearSystem.IsConsistent(SparseMatrixTag::LHS, DenseVectorTag::RHS, DenseVectorTag::Dx)) {
            return;
        }

        auto& r_lin_op = *rLinearSystem.pGetLinearOperator(SparseMatrixTag::LHS);
        KRATOS_ERROR_IF(r_lin_op.IsMatrixFree()) << "EigenSparseLUSolver cannot be used with matrix-free linear operators." << std::endl;
        const auto& r_A = r_lin_op.GetMatrix();

        // Factorize (the map is read-only: Eigen copies the matrix into its own storage)
        const EigenSparseMapType eigen_A = ToEigenSparseMap(r_A);
        mEigenSparseLUSolver.compute(eigen_A);
        KRATOS_ERROR_IF(mEigenSparseLUSolver.info() != Eigen::Success)
            << "EigenSparseLUSolver: the LHS matrix is singular and the factorization failed." << std::endl;

        mFactorized = true;
    }

    /**
     * @brief Solves the linear system by back-substitution with the factorization
     *        computed in InitializeSolutionStep.
     * @param rLinearSystem The linear system to be solved
     * @return @p true after a successful back-substitution
     */
    bool PerformSolutionStep(LinearSystemType& rLinearSystem) override
    {
        if (!rLinearSystem.IsConsistent(SparseMatrixTag::LHS, DenseVectorTag::RHS, DenseVectorTag::Dx)) {
            return false;
        }

        KRATOS_ERROR_IF(!mFactorized)
            << "EigenSparseLUSolver: InitializeSolutionStep must be called before PerformSolutionStep (the LHS matrix is not factorized)." << std::endl;

        auto& rX = *rLinearSystem.pGetVector(DenseVectorTag::Dx);
        auto& rB = *rLinearSystem.pGetVector(DenseVectorTag::RHS);

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
    bool mFactorized = false;                                    /// True once the LHS matrix has been factorized
    std::vector<EigenIndexType> mIndex1; /// Copy of the CsrMatrix row pointers, in the Eigen index type
    std::vector<EigenIndexType> mIndex2; /// Copy of the CsrMatrix column indices, in the Eigen index type

    ///@}
    ///@name Private Operations
    ///@{

    /// Exposes a CsrMatrix as a read-only Eigen map: the values are mapped in place and only
    /// the index arrays are copied (the CsrMatrix index type is unsigned, Eigen requires signed)
    EigenSparseMapType ToEigenSparseMap(const MatrixType& rA)
    {
        mIndex1.assign(rA.index1_data().begin(), rA.index1_data().end());
        mIndex2.assign(rA.index2_data().begin(), rA.index2_data().end());
        return EigenSparseMapType(rA.size1(), rA.size2(), mIndex2.size(),
                                  mIndex1.data(), mIndex2.data(), rA.value_data().data());
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
