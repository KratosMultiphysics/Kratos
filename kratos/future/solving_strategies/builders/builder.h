//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/amgcl_csr_conversion_utilities.h"
#include "utilities/amgcl_csr_spmm_utilities.h"
#include "utilities/builtin_timer.h"
#include "utilities/dof_utilities/dof_array_utilities.h"
#include "utilities/timer.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

//FIXME: We should place this somewhere else
/**
 * @brief Auxiliary container to store the linear system
 * This auxiliary container is intended to store all the arrays requires for the linear system setup
 * @tparam TSparseMatrixType The sparse matrix type
 * @tparam TSystemVectorType The system vector type
 */
template <class TSparseMatrixType, class TSystemVectorType>
struct LinearSystemContainer
{
    typename TSparseMatrixType::Pointer pLhs = nullptr; // Pointer to the LHS matrix

    typename TSystemVectorType::Pointer pRhs = nullptr; // Pointer to the RHS vector

    typename TSystemVectorType::Pointer pDx = nullptr; // Pointer to the solution increment vector

    typename TSparseMatrixType::Pointer pEffectiveLhs = nullptr; // Pointer to the effective LHS matrix (i.e., after applying system constraints)

    typename TSystemVectorType::Pointer pEffectiveRhs = nullptr; // Pointer to the effective RHS vector (i.e., after applying system constraints)

    typename TSystemVectorType::Pointer pEffectiveDx = nullptr; // Pointer to the effective solution increment vector (i.e., after applying system constraints)

    typename TSparseMatrixType::Pointer pEffectiveT = nullptr; // Linear system constraints total relation matrix

    typename TSparseMatrixType::Pointer pConstraintsT = nullptr; // Master-slave constraints relation matrix

    typename TSystemVectorType::Pointer pConstraintsQ = nullptr; // Master-slave constraints constant vector

    typename TSparseMatrixType::Pointer pDirichletT = nullptr; // Dirichlet constraints relation matrix

    typename TSparseMatrixType::Pointer pMassMatrix = nullptr; // Pointer to the mass matrix

    typename TSparseMatrixType::Pointer pDampingMatrix = nullptr; // Pointer to the damping matrix

    void Clear()
    {
        if (pLhs != nullptr) {
            pLhs->Clear();
        }
        if (pRhs != nullptr) {
            pRhs->Clear();
        }
        if (pDx != nullptr) {
            pDx->Clear();
        }
        if (pEffectiveLhs != nullptr) {
            pEffectiveLhs->Clear();
        }
        if (pEffectiveRhs != nullptr) {
            pEffectiveRhs->Clear();
        }
        if (pEffectiveDx != nullptr) {
            pEffectiveDx->Clear();
        }
        if (pEffectiveT != nullptr) {
            pEffectiveT->Clear();
        }
        if (pConstraintsT != nullptr) {
            pConstraintsT->Clear();
        }
        if (pConstraintsQ != nullptr) {
            pConstraintsQ->Clear();
        }
        if (pDirichletT != nullptr) {
            pDirichletT->Clear();
        }
        if (pMassMatrix != nullptr) {
            pMassMatrix->Clear();
        }
        if (pDampingMatrix != nullptr) {
            pDampingMatrix->Clear();
        }
    }
};

/**
 * @class Builder
 * @ingroup KratosCore
 * @brief Base builder class
 * @details This class collects the methods required to handle the building of the sparse sytem matrices
 * This class is thought to never be used, but to serve as basis for all the derived build types (e.g., block and elimination)
 * @author Ruben Zorrilla
 */
template<class TThreadLocalStorage, class TSparseMatrixType, class TSystemVectorType, class TSparseGraphType>
class Builder
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Builder
    KRATOS_CLASS_POINTER_DEFINITION(Builder);

    /// Data type definition from sparse matrix
    using DataType = typename TSparseMatrixType::DataType;

    /// Index type definition from sparse matrix
    using IndexType = typename TSparseMatrixType::IndexType;

    /// DOF type definition
    using DofType = Dof<DataType>;

    /// DOF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    /// DOF pointer vector type definition
    using DofPointerVectorType = typename MasterSlaveConstraint::DofPointerVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Builder() = delete;

    /// Constructor with model part
    Builder(
        const ModelPart& rModelPart,
        Parameters Settings = Parameters(R"({})"))
    : mpModelPart(&rModelPart)
    {
        Parameters default_parameters( R"({
            "name" : "builder",
            "echo_level" : 0
        })");
        Settings.ValidateAndAssignDefaults(default_parameters);

        // Set verbosity level
        mEchoLevel = Settings["echo_level"].GetInt();
    }

    virtual ~Builder() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Allocates the memory for the linear system arrays
     * This method allocates the memory for the linear system arrays
     * Note that the sizes of the resultant arrays depend on the build type
     * @param pDofSet The array of DOFs from elements and conditions
     * @param pEffectiveDofSet The array of DOFs to be solved after the application of constraints
     * @param pLinearSystemContainer Auxiliary container with the linear system arrays
     */
    virtual void AllocateLinearSystemArrays(
        const typename DofsArrayType::Pointer pDofSet,
        const typename DofsArrayType::Pointer pEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer)
        {
            KRATOS_ERROR << "Calling base class 'AllocateLinearSystemArrays'." << std::endl;
        }

    /**
     * @brief Allocates the linear system constraints arrays
     * This method allocates the linear system constraints arrays
     * Note that the sizes of the resultant arrays depend on the build type
     * @param rDofSet The array of DOFs from elements and conditions
     * @param rEffectiveDofSet The effective DOFs array (i.e., those that are not slaves)
     * @param rSlaveToMasterDofsMap The map containing the corresponding master(s) for each slave DOF
     * @param rLinearSystemContainer Auxiliary container with the linear system arrays
     */
    virtual void AllocateLinearSystemConstraints(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        const DofArrayUtilities::SlaveToMasterDofsMap& rSlaveToMasterDofsMap,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class 'AllocateLinearSystemArrays'." << std::endl;
    }

    /**
     * @brief Builds the Dirichlet constraints
     * This method does the building of the Dirichlet constraints according to the build type
     * @param rDofSet The array of DOFs from elements and conditions
     * @param rEffectiveDofSet The effective DOFs array (i.e., those that are not slaves)
     * @param rLinearSystemContainer Auxiliary container with the linear system arrays
     */
    virtual void BuildDirichletConstraints(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class 'BuildDirichletConstraints'." << std::endl;
    }

    /**
     * @brief Set the Up Sparse Matrix Graph
     * This method builds the sparse matrix graph from the connectivities of the
     * elements and conditions in the model part the builder refers to
     * @param rSparseGraph The sparse matrix graph to be filled
     */
    virtual void SetUpSparseMatrixGraph(TSparseGraphType& rSparseGraph)
    {
        // Clear the provided sparse matrix graph
        rSparseGraph.Clear();

        // Add the elements and conditions DOF equation connectivities
        // Note that we add all the DOFs regarless their fixity status
        Element::EquationIdVectorType eq_ids;
        for (auto& r_elem : mpModelPart->Elements()) {
            r_elem.EquationIdVector(eq_ids, mpModelPart->GetProcessInfo());
            rSparseGraph.AddEntries(eq_ids);
        }
        for (auto& r_cond : mpModelPart->Conditions()) {
            r_cond.EquationIdVector(eq_ids, mpModelPart->GetProcessInfo());
            rSparseGraph.AddEntries(eq_ids);
        }
    }

    /**
     * @brief Set the Up Master Slave Constraints Graph object
     * This method builds the sparse graph representing the master-slave constraints DOFs
     * based on the provided slave to master(s) map
     * @param rDofSet The array of DOFs from elements and conditions
     * @param rEffectiveDofSet The effective DOFs array (i.e., those that are not slaves)
     * @param rSlaveToMasterDofsMap The map containing the corresponding master(s) for each slave DOF
     * @param rConstraintsSparseGraph The master-slave constraints graph to be built
     */
    virtual void SetUpMasterSlaveConstraintsGraph(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        const DofArrayUtilities::SlaveToMasterDofsMap& rSlaveToMasterDofsMap,
        TSparseGraphType& rConstraintsSparseGraph)
    {
        // Check if there are constraints to build the constraints sparse graph
        if (!rSlaveToMasterDofsMap.empty()) {
            // Initialize the constraints matrix sparse graph
            // Note that the number of rows is the size of the standard DOFs arrays (the effective one will be the columns)
            rConstraintsSparseGraph = std::move(TSparseGraphType(rDofSet.size()));

            // Loop the elements and conditions DOFs container to add the slave entries to the graph
            for (IndexType i_dof = 0; i_dof < rDofSet.size(); ++i_dof) {
                // Get current DOF
                auto p_dof = *(rDofSet.ptr_begin() + i_dof);
                const IndexType i_dof_eq_id = p_dof->EquationId();

                // Check if current DOF is slave by checking the slaves DOFs map
                // If not present in the slaves DOFs map it should be considered an effective DOF
                // Note that here effective means a masters DOF or a DOF that does not involve any constraint
                auto i_dof_slave_find = rSlaveToMasterDofsMap.find(p_dof);
                if (i_dof_slave_find != rSlaveToMasterDofsMap.end()) { // Slave DOF
                    // Add current slave DOF connectivities to the constraints sparse graph
                    // The slave rows eq ids come from the system ones while the column master ones are the above defined
                    for (auto& rp_master : i_dof_slave_find->second) {
                        auto eff_dof_find = rEffectiveDofSet.find(*rp_master);
                        KRATOS_ERROR_IF(eff_dof_find == rEffectiveDofSet.end()) << "Effective DOF cannot be find." << std::endl;
                        rConstraintsSparseGraph.AddEntry(i_dof_eq_id, eff_dof_find->EffectiveEquationId());
                    }
                } else { // Effective DOF
                    auto eff_dof_find = rEffectiveDofSet.find(*p_dof);
                    KRATOS_ERROR_IF(eff_dof_find == rEffectiveDofSet.end()) << "Effective DOF cannot be find." << std::endl;
                    // mMasterIds.push_back(eff_dof_find->second);
                    rConstraintsSparseGraph.AddEntry(i_dof_eq_id, eff_dof_find->EffectiveEquationId());
                }
            }
        }
    }

    /**
     * @brief Applies the constraints to the linear system of equations
     * This method applies the linear system constraints to the linear system of equations
     * Note that linear system constraints includes both the master-slave and the Dirichlet
     * constraints. The way the constraints are applied will be reimplemented and applied
     * in the derived classes depending on the build type.
     * @param rEffectiveDofSet The effective DOFs array (i.e., those that are not slaves)
     * @param rLinearSystemContainer Auxiliary container with the linear system arrays
     */
    virtual void ApplyLinearSystemConstraints(
        const DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class 'ApplyLinearSystemConstraints'." << std::endl;
    }

    /**
     * @brief Set the solution values vector
     * This method calculates the complete solution vector (x) from the effective DOF set
     * For doing so the effective DOF values are retrieved to then apply the constraints
     * by doing x = T * y + q, being T and q the constraints relation matrix and constant
     * vector respectively and y the auxliary vector containing the effective DOFs values
     * Note that the vector y already includes the Dirichlet constraints as values are
     * retrieved from the DOF (i.e., nodal) database.
     * @param rEffectiveDofSet The effective DOFs array (i.e., those that are not slaves)
     * @param rConstraintsRelationMatrix The constraints relation matrix
     * @param rConstraintsConstantVector The constraints constant vector
     * @param rSolutionVector The solution vector
     */
    void CalculateSolutionVector(
        const DofsArrayType& rEffectiveDofSet,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSystemVectorType& rConstraintsConstantVector,
        TSystemVectorType& rSolutionVector) const
    {
        // Set an auxiliary vector containing the effective solution values
        const std::size_t n_eff_dofs = rEffectiveDofSet.size();
        TSystemVectorType y(n_eff_dofs);
        IndexPartition<IndexType>(rEffectiveDofSet.size()).for_each([&](IndexType Index) {
            // Get effective DOF
            auto p_dof = *(rEffectiveDofSet.ptr_begin() + Index);

            // Get value from DOF and set it in the auxiliary solution values vector
            // Note that the corresponding row is retrieved from the effective DOF id
            y[p_dof->EffectiveEquationId()] = p_dof->GetSolutionStepValue();
        });

        // Check solution vector size
        const std::size_t aux_size = rConstraintsConstantVector.size();
        if (rSolutionVector.size() != aux_size) {
            rSolutionVector = TSystemVectorType(aux_size);
        }

        // Initialize solution vector with the constaints constant vector values
        // Note that we deliberately avoid using the copy constructor as we dont want to overwrite the constraints constant vector
        IndexPartition<IndexType>(rSolutionVector.size()).for_each([&](IndexType i){
            rSolutionVector[i] = rConstraintsConstantVector[i];
        });

        // Compute the solution vector as x = T * y + q
        rConstraintsRelationMatrix.SpMV(1.0, y, 1.0, rSolutionVector); // Note that this performs the operation x = A*y + x
    }

    /**
     * @brief Clear class content
     * This method wipes the content of current builder class
     */
    virtual void Clear()
    {
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Get the Model Part object
     * Returns a reference to the model part the builder is referring to
     * @return const ModelPart& Reference to the builder model part
     */
    const ModelPart& GetModelPart() const
    {
        return *mpModelPart;
    }

    /**
     * @brief Get the Echo Level object
     * Returns the echo level member variable
     * @return const std::size_t Level of information that is output
     */
    const std::size_t GetEchoLevel() const
    {
        return mEchoLevel;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    const ModelPart* mpModelPart = nullptr; /// Pointer to the model part the builder will refer to

    std::size_t mEchoLevel; /// Level of information that is output

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
}; // Class Builder

///@}
///@name Input and output
///@{


///@}
///@} addtogroup block

}  // namespace Kratos.
