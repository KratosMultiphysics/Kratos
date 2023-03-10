//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
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
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
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
    TrilinosAssemblingUtilities()
    {
    }

    /// Destructor.
    virtual ~TrilinosAssemblingUtilities()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// TODO: creating the the calculating reaction version
    // 	template<class TOtherMatrixType, class TEquationIdVectorType>

    /**
     * @brief Assembles the LHS of the system
     * @param rA The LHS matrix
     * @param rLHSContribution The contribution to the LHS
     * @param rEquationId The equation ids
     */
    inline static void AssembleLHS(
        MatrixType& rA,
        const Matrix& rLHSContribution,
        const std::vector<std::size_t>& rEquationId
        )
    {
        const unsigned int system_size = TSparseSpace::Size1(rA);

        // Count active indices
        unsigned int active_indices = 0;
        for (unsigned int i = 0; i < rEquationId.size(); i++)
            if (rEquationId[i] < system_size)
                ++active_indices;

        if (active_indices > 0) {
            // Size Epetra vectors
            Epetra_IntSerialDenseVector indices(active_indices);
            Epetra_SerialDenseMatrix values(active_indices, active_indices);

            // Fill epetra vectors
            unsigned int loc_i = 0;
            for (unsigned int i = 0; i < rEquationId.size(); i++) {
                if (rEquationId[i] < system_size) {
                    indices[loc_i] = rEquationId[i];

                    unsigned int loc_j = 0;
                    for (unsigned int j = 0; j < rEquationId.size(); j++) {
                        if (rEquationId[j] < system_size) {
                            values(loc_i, loc_j) = rLHSContribution(i, j);
                            ++loc_j;
                        }
                    }
                    ++loc_i;
                }
            }

            const int ierr = rA.SumIntoGlobalValues(indices, values);
            KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found" << std::endl;
        }
    }

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
        const unsigned int system_size = TSparseSpace::Size1(rT);

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

    //***********************************************************************
    /// TODO: creating the the calculating reaction version
    // 	template<class TOtherVectorType, class TEquationIdVectorType>

    /**
     * @brief Assembles the RHS of the system
     * @param rb The RHS vector
     * @param rRHSContribution The RHS contribution
     * @param rEquationId The equation ids
     */
    inline static void AssembleRHS(
        VectorType& rb,
        const Vector& rRHSContribution,
        const std::vector<std::size_t>& rEquationId
        )
    {
        const unsigned int system_size = TSparseSpace::Size(rb);

        // Count active indices
        unsigned int active_indices = 0;
        for (unsigned int i = 0; i < rEquationId.size(); i++)
            if (rEquationId[i] < system_size)
                ++active_indices;

        if (active_indices > 0) {
            // Size Epetra vectors
            Epetra_IntSerialDenseVector indices(active_indices);
            Epetra_SerialDenseVector values(active_indices);

            // Fill epetra vectors
            unsigned int loc_i = 0;
            for (unsigned int i = 0; i < rEquationId.size(); i++) {
                if (rEquationId[i] < system_size) {
                    indices[loc_i] = rEquationId[i];
                    values[loc_i] = rRHSContribution[i];
                    ++loc_i;
                }
            }

            const int ierr = rb.SumIntoGlobalValues(indices, values);
            KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found" << std::endl;
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
        const unsigned int system_size = TSparseSpace::Size(rC);

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
