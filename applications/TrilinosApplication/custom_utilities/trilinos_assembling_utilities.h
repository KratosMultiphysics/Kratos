//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes
#include <Epetra_SerialDenseVector.h>
#include <Epetra_FEVector.h>
#include <Epetra_FECrsMatrix.h>
#if (HAVE_TPETRA)
#include <Tpetra_FEMultiVector.hpp>
#include <Tpetra_FECrsMatrix.hpp>
#endif

// Project includes

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
 * @class TrilinosAssemblingUtilities
 * @ingroup TrilinosApplication
 * @brief The Trilinos assembling utilities
 * @details This class provides utility functions for assembling matrices and vectors in Trilinos.
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space type
 */
template<class TSparseSpace>
class TrilinosAssemblingUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosAssemblingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosAssemblingUtilities);

    /// Definition of the matrix type
    using MatrixType = typename TSparseSpace::MatrixType;

    /// Definition of the vector type
    using VectorType = typename TSparseSpace::VectorType;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Definition of the size type
    using SizeType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TrilinosAssemblingUtilities() = delete;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Assembles the relation matrix T of the system with MPC
     * @param rT The T relation matrix
     * @param rTContribution The contribution to the T
     * @param rSlaveEquationId The slave equation ids
     * @param rMasterEquationId The master equation ids
     */
    inline static void AssembleRelationMatrixT(
        MatrixType& rT,
        const Matrix& rTContribution,
        const std::vector<std::size_t>& rSlaveEquationId,
        const std::vector<std::size_t>& rMasterEquationId
        )
    {
        if constexpr (TSparseSpace::LinearAlgebraLibrary() == TrilinosLinearAlgebraLibrary::EPETRA) {
            const unsigned int system_size = rT.NumGlobalRows();

            // Count active indices
            int slave_active_indices = 0;
            for (unsigned int i = 0; i < rSlaveEquationId.size(); i++) {
                if (rSlaveEquationId[i] < system_size) {
                    ++slave_active_indices;
                }
            }
            int master_active_indices = 0;
            for (unsigned int i = 0; i < rMasterEquationId.size(); i++) {
                if (rMasterEquationId[i] < system_size) {
                    ++master_active_indices;
                }
            }

            if (slave_active_indices > 0 && master_active_indices > 0) {
                std::vector<int> indices(master_active_indices);
                std::vector<double> values(master_active_indices);

                // Fill Epetra vectors
                for (unsigned int i = 0; i < rSlaveEquationId.size(); i++) {
                    if (rSlaveEquationId[i] < system_size) {
                        const int current_global_row = rSlaveEquationId[i];

                        unsigned int loc_j = 0;
                        for (unsigned int j = 0; j < rMasterEquationId.size(); j++) {
                            if (rMasterEquationId[j] < system_size) {
                                indices[loc_j] = rMasterEquationId[j];
                                values[loc_j] = rTContribution(i, j);
                                ++loc_j;
                            }
                        }

                        const int ierr = rT.SumIntoGlobalValues(current_global_row, master_active_indices, values.data(), indices.data());
                        KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found" << std::endl;
                    }
                }
            }
        } else if constexpr (TSparseSpace::LinearAlgebraLibrary() == TrilinosLinearAlgebraLibrary::TPETRA) {
        #if (HAVE_TPETRA)
            using GO = typename MatrixType::global_ordinal_type;
            using ST = typename MatrixType::scalar_type;
            const std::size_t system_size = rT.getGlobalNumRows();

            // Open for FE assembly if not already active (mirrors lazy-open in AssembleLHS).
            auto p_fe_T = dynamic_cast<MatrixType*>(&rT);
            if (p_fe_T && !rT.isFillActive()) {
                p_fe_T->beginAssembly();
            }

            for (std::size_t i = 0; i < rSlaveEquationId.size(); ++i) {
                if (rSlaveEquationId[i] < system_size) {
                    const GO global_id = static_cast<GO>(rSlaveEquationId[i]);
                    std::vector<GO> indices;
                    std::vector<ST> values;
                    for (std::size_t j = 0; j < rMasterEquationId.size(); ++j) {
                        if (rMasterEquationId[j] < system_size) {
                            indices.push_back(static_cast<GO>(rMasterEquationId[j]));
                            values.push_back(static_cast<ST>(rTContribution(i, j)));
                        }
                    }
                    if (indices.size() > 0) {
                        rT.sumIntoGlobalValues(global_id, indices.size(), values.data(), indices.data());
                    }
                }
            }
        #else
            KRATOS_ERROR << "You must compile Kratos with TPETRA support" << std::endl;
        #endif
        } else {
            KRATOS_ERROR << "Only EPETRA and TPETRA are supported for now" << std::endl;
        }
    }

    /**
     * @brief Assembles the Constant vector of the system with MPC
     * @param rC The constant vector
     * @param rConstantContribution The RHS contribution
     * @param rEquationId The equation ids
     */
    inline static void AssembleConstantVector(
        VectorType& rC,
        const Vector& rConstantContribution,
        const std::vector<std::size_t>& rSlaveEquationId
        )
    {
        if constexpr (TSparseSpace::LinearAlgebraLibrary() == TrilinosLinearAlgebraLibrary::EPETRA) {
            const unsigned int system_size = rC.GlobalLength();

            // Count active indices
            unsigned int slave_active_indices = 0;
            for (unsigned int i = 0; i < rSlaveEquationId.size(); i++)
                if (rSlaveEquationId[i] < system_size)
                    ++slave_active_indices;

            if (slave_active_indices > 0) {
                // Size Epetra vectors
                Epetra_IntSerialDenseVector indices(slave_active_indices);
                Epetra_SerialDenseVector values(slave_active_indices);

                // Fill Epetra vectors
                unsigned int loc_i = 0;
                for (unsigned int i = 0; i < rSlaveEquationId.size(); i++) {
                    if (rSlaveEquationId[i] < system_size) {
                        indices[loc_i] = rSlaveEquationId[i];
                        values[loc_i] = rConstantContribution[i];
                        ++loc_i;
                    }
                }

                const int ierr = rC.SumIntoGlobalValues(indices, values);
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found" << std::endl;
            }
        } else if constexpr (TSparseSpace::LinearAlgebraLibrary() == TrilinosLinearAlgebraLibrary::TPETRA) {
        #if (HAVE_TPETRA)
            using GO = typename VectorType::global_ordinal_type;
            using ST = typename VectorType::scalar_type;
            const std::size_t system_size = rC.getGlobalLength();

            for (std::size_t i = 0; i < rSlaveEquationId.size(); ++i) {
                if (rSlaveEquationId[i] < system_size) {
                    const GO global_id = static_cast<GO>(rSlaveEquationId[i]);
                    const ST val = static_cast<ST>(rConstantContribution[i]);
                    rC.sumIntoGlobalValue(global_id, size_t(0), val);
                }
            }
        #else
            KRATOS_ERROR << "You must compile Kratos with TPETRA support" << std::endl;
        #endif
        } else {
            KRATOS_ERROR << "Only EPETRA and TPETRA are supported for now" << std::endl;
        }
    }

    /**
     * @brief Sets a value in a vector
     * @param rX The vector considered
     * @param i The index of the value considered
     * @param Value The value considered
     */
    static inline void SetGlobalValue(
        VectorType& rX,
        const IndexType i,
        const double Value
        )
    {
        TSparseSpace::SetGlobalVec(rX, i, Value);
    }

    /**
     * @brief Sets a value in a vector
     * @param rX The vector considered
     * @param i The index of the value considered
     * @param Value The value considered
     */
    static inline void SetGlobalValueWithoutGlobalAssembly(
        VectorType& rX,
        const IndexType i,
        const double Value
        )
    {
        TSparseSpace::SetGlobalVecNoAssemble(rX, i, Value);
    }

    /**
     * @brief Sets a value in a vector (local)
     * @param rX The vector considered
     * @param i The index of the value considered
     * @param Value The value considered
     */
    static inline void SetLocalValue(
        VectorType& rX,
        const IndexType i,
        const double Value
        )
    {
        TSparseSpace::SetLocalVec(rX, i, Value);
    }

    /**
     * @brief Sets a value in a vector (local without global assembly)
     * @param rX The vector considered
     * @param i The index of the value considered
     * @param Value The value considered
     */
    static inline void SetLocalValueWithoutGlobalAssembly(
        VectorType& rX,
        const IndexType i,
        const double Value
        )
    {
        TSparseSpace::SetLocalVecNoAssemble(rX, i, Value);
    }

    /**
     * @brief Sets a value in a matrix
     * @param rA The matrix considered
     * @param i The first index of the value considered
     * @param j The second index of the value considered
     * @param Value The value considered
     */
    static inline void SetGlobalValue(
        MatrixType& rA,
        const IndexType i,
        const IndexType j,
        const double Value
        )
    {
        TSparseSpace::SetGlobalMat(rA, i, j, Value);
    }

    /**
     * @brief Sets a value in a matrix
     * @param rA The matrix considered
     * @param i The first index of the value considered
     * @param j The second index of the value considered
     * @param Value The value considered
     */
    static inline void SetGlobalValueWithoutGlobalAssembly(
        MatrixType& rA,
        const IndexType i,
        const IndexType j,
        const double Value
        )
    {
        TSparseSpace::SetGlobalMatNoAssemble(rA, i, j, Value);
    }

    /**
     * @brief Sets a value in a matrix
     * @param rA The matrix considered
     * @param i The first index of the value considered
     * @param j The second index of the value considered
     * @param Value The value considered
     */
    static inline void SetLocalValue(
        MatrixType& rA,
        const IndexType i,
        const IndexType j,
        const double Value
        )
    {
        TSparseSpace::SetLocalMat(rA, i, j, Value);
    }

    /**
     * @brief Sets a value in a matrix
     * @param rA The matrix considered
     * @param i The first index of the value considered
     * @param j The second index of the value considered
     * @param Value The value considered
     */
    static inline void SetLocalValueWithoutGlobalAssembly(
        MatrixType& rA,
        const IndexType i,
        const IndexType j,
        const double Value
        )
    {
        TSparseSpace::SetLocalMatNoAssemble(rA, i, j, Value);
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

    /**
     * @brief Turn back information as a string.
     * @return Info as a string.
     */
    virtual std::string Info() const
    {
        return "TrilinosAssemblingUtilities";
    }

    /**
     * @brief Print information about this object.
     * @param rOStream The output stream to print on.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TrilinosAssemblingUtilities";
    }

    /**
     * @brief Print object's data.
     * @param rOStream The output stream to print on.
     */
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TrilinosAssemblingUtilities & operator=(TrilinosAssemblingUtilities const& rOther);

    /// Copy constructor.
    TrilinosAssemblingUtilities(TrilinosAssemblingUtilities const& rOther);

    ///@}
}; // Class TrilinosAssemblingUtilities

///@}

} // namespace Kratos.
