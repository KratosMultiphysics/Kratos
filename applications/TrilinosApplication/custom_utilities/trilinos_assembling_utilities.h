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

// Project includes
#include "trilinos_space.h"

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
 * @author Vicente Mataix Ferrandiz
 */
class TrilinosAssemblingUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosAssemblingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosAssemblingUtilities);

    /// Definition of Trilinos space
    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;

    /// Definition of the matrix type
    using MatrixType = TrilinosSparseSpaceType::MatrixType;

    /// Definition of the vector type
    using VectorType = TrilinosSparseSpaceType::VectorType;

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
        const unsigned int system_size = TrilinosSparseSpaceType::Size1(rT);

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
            std::vector<int> indices(slave_active_indices);
            std::vector<double> values(master_active_indices);

            // Fill epetra vectors
            unsigned int loc_i = 0;
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

                    ++loc_i;
                }
            }
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
        const unsigned int system_size = TrilinosSparseSpaceType::Size(rC);

        // Count active indices
        unsigned int slave_active_indices = 0;
        for (unsigned int i = 0; i < rSlaveEquationId.size(); i++)
            if (rSlaveEquationId[i] < system_size)
                ++slave_active_indices;

        if (slave_active_indices > 0) {
            // Size Epetra vectors
            Epetra_IntSerialDenseVector indices(slave_active_indices);
            Epetra_SerialDenseVector values(slave_active_indices);

            // Fill epetra vectors
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
