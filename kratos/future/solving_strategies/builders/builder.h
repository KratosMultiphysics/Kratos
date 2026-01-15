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
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/builtin_timer.h"
#include "utilities/dof_utilities/dof_array_utilities.h"
#include "utilities/timer.h"

#ifdef KRATOS_USE_FUTURE
#include "future/containers/implicit_strategy_data_container.h"
#endif

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class Builder
 * @ingroup KratosCore
 * @brief Base builder class
 * @details This class collects the methods required to handle the building of the sparse sytem matrices
 * This class is thought to never be used, but to serve as basis for all the derived build types (e.g., block and elimination)
 * @author Ruben Zorrilla
 */
template<class TLinearAlgebra>
class Builder
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Builder
    KRATOS_CLASS_POINTER_DEFINITION(Builder);

    /// Data type definition from sparse matrix
    using DataType = typename TLinearAlgebra::DataType;

    /// Index type definition from sparse matrix
    using IndexType = typename TLinearAlgebra::IndexType;

    /// Matrix type definition
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Vector type definition
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Sparse graph type definition
    using SparseGraphType = TLinearAlgebra::SparseGraphType;

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
     * @param rSparseGraph Reference to the linear system sparse graph
     * @param pLinearSystemContainer Auxiliary container with the linear system arrays
     */
    virtual void AllocateLinearSystem(
        const SparseGraphType& rSparseGraph,
        ImplicitStrategyDataContainer<TLinearAlgebra> &rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class 'AllocateLinearSystem'." << std::endl;
    }

    /**
     * @brief Allocates the memory for the linear system arrays
     * This method calculates the sparse graph and allocates the memory for the linear system arrays
     * Note that the sizes of the resultant arrays depend on the build type
     * @param pLinearSystemContainer Auxiliary container with the linear system arrays
     */
    virtual void AllocateLinearSystem(ImplicitStrategyDataContainer<TLinearAlgebra> &rLinearSystemContainer)
    {
        // Set up the system sparse matrix graph (note that the sparse graph will be destroyed when leaving this scope)
        BuiltinTimer sparse_matrix_graph_time;
        auto p_dof_set = rLinearSystemContainer.pDofSet;
        SparseGraphType sparse_matrix_graph(p_dof_set->size());
        this->SetUpSparseMatrixGraph(sparse_matrix_graph);
        KRATOS_INFO_IF("BlockBuilder", this->GetEchoLevel() > 0) << "Set up sparse matrix graph time: " << sparse_matrix_graph_time << std::endl;

        // Allocate the linear system
        this->AllocateLinearSystem(sparse_matrix_graph, rLinearSystemContainer);
    }

    /**
     * @brief Allocates the linear system constraints arrays
     * This method allocates the linear system constraints arrays
     * Note that the sizes of the resultant arrays depend on the build type
     * @param rLinearSystemContainer Auxiliary container with the linear system arrays
     */
    virtual void AllocateLinearSystemConstraints(ImplicitStrategyDataContainer<TLinearAlgebra>& rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class 'AllocateLinearSystemConstraints'." << std::endl;
    }

    /**
     * @brief Set the Up Sparse Matrix Graph
     * This method builds the sparse matrix graph from the connectivities of the
     * elements and conditions in the model part the builder refers to
     * @param rSparseGraph The sparse matrix graph to be filled
     */
    virtual void SetUpSparseMatrixGraph(SparseGraphType& rSparseGraph)
    {
        // Clear the provided sparse matrix graph
        rSparseGraph.Clear();

        // Add the elements and conditions DOF equation connectivities
        // Note that we add all the DOFs regardless their fixity status
        if constexpr (SparseGraphType::IsThreadSafe) {
            IndexPartition<IndexType>(mpModelPart->NumberOfElements()).for_each([&](IndexType Index) {
                Element::EquationIdVectorType eq_ids; //TODO: we don't use TLS for this (decide what to do once we finish the parallelism discussion)
                auto it_elem = mpModelPart->ElementsBegin() + Index;
                it_elem->EquationIdVector(eq_ids, mpModelPart->GetProcessInfo());
                rSparseGraph.AddEntries(eq_ids);
            });
            IndexPartition<IndexType>(mpModelPart->NumberOfConditions()).for_each([&](IndexType Index) {
                Condition::EquationIdVectorType eq_ids; //TODO: we don't use TLS for this (decide what to do once we finish the parallelism discussion)
                auto it_cond = mpModelPart->ElementsBegin() + Index;
                it_cond->EquationIdVector(eq_ids, mpModelPart->GetProcessInfo());
                rSparseGraph.AddEntries(eq_ids);
            });
            //TODO: Add the constraints in here!
        } else {
            Element::EquationIdVectorType eq_ids;
            for (auto& r_elem : mpModelPart->Elements()) {
                r_elem.EquationIdVector(eq_ids, mpModelPart->GetProcessInfo());
                rSparseGraph.AddEntries(eq_ids);
            }
            for (auto& r_cond : mpModelPart->Conditions()) {
                r_cond.EquationIdVector(eq_ids, mpModelPart->GetProcessInfo());
                rSparseGraph.AddEntries(eq_ids);
            }
            //TODO: Add the constraints in here!
        }
    }

    /**
     * @brief Set the Up Master Slave Constraints Graph object
     * This method builds the sparse graph representing the master-slave constraints DOFs
     * based on the provided slave to master(s) map
     * @param rDofSet The array of DOFs from elements and conditions
     * @param rEffectiveDofSet The effective DOFs array (i.e., those that are not slaves)
     * @param rConstraintsSparseGraph The master-slave constraints graph to be built
     */
    virtual void SetUpMasterSlaveConstraintsGraph(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        SparseGraphType& rConstraintsSparseGraph)
    {
        // Check if there are constraints to build the constraints sparse graph
        const std::size_t n_constraints = mpModelPart->NumberOfMasterSlaveConstraints();
        if (n_constraints) { //TODO: Change the way this is checked (w/o the DOFs map)
            // Initialize the constraints matrix sparse graph
            // Note that the number of rows is the size of the standard DOFs arrays (the effective one will be the columns)
            rConstraintsSparseGraph = std::move(SparseGraphType(rDofSet.size()));

            // Add all effective DOFs
            // Note that here effective means that the DOF is either a master DOF or a DOF that does not involve any constraint
            auto eff_dof_addition_fn = [&](IndexType Index){
                auto p_eff_dof = *(rEffectiveDofSet.ptr_begin() + Index);
                const IndexType i_dof_eq_id = p_eff_dof->EquationId();
                const IndexType i_dof_eff_eq_id = p_eff_dof->EffectiveEquationId();
                rConstraintsSparseGraph.AddEntry(i_dof_eq_id, i_dof_eff_eq_id);
            };

            const std::size_t n_eff_dof = rEffectiveDofSet.size();
            if constexpr (SparseGraphType::IsThreadSafe) {
                IndexPartition<IndexType>(n_eff_dof).for_each(eff_dof_addition_fn);
            } else {
                for (IndexType i_dof = 0; i_dof < n_eff_dof; ++i_dof) {
                    eff_dof_addition_fn(i_dof);
                }
            }

            // Loop the constraints to add the slave DOFs
            // Note that we assume that there constraints are always one to many (not many to many)
            auto slave_dof_addition_fn = [&](IndexType Index){
                auto it_constraint = mpModelPart->MasterSlaveConstraintsBegin() + Index;
                const auto& r_masters = it_constraint->GetMasterDofsVector(); //TODO: we don't use TLS for this (decide what to do once we finish the parallelism discussion)
                const IndexType slave_eq_id = (it_constraint->GetSlaveDofsVector()[0])->EquationId();
                for (const auto& rp_master_dof : r_masters) {
                    const IndexType master_eff_eq_id = rp_master_dof->EffectiveEquationId();
                    rConstraintsSparseGraph.AddEntry(slave_eq_id, master_eff_eq_id);
                }
            };

            if constexpr (SparseGraphType::IsThreadSafe) {
                IndexPartition<IndexType>(n_constraints).for_each(slave_dof_addition_fn);
            } else {
                for (IndexType i_constraint = 0; i_constraint < n_constraints; ++i_constraint) {
                    slave_dof_addition_fn(i_constraint);
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
     * @param rLinearSystemContainer Auxiliary container with the linear system arrays
     */
    virtual void ApplyLinearSystemConstraints(ImplicitStrategyDataContainer<TLinearAlgebra>& rLinearSystemContainer)
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
        const MatrixType& rConstraintsRelationMatrix,
        const VectorType& rConstraintsConstantVector,
        VectorType& rSolutionVector) const
    {
        // Set an auxiliary vector containing the effective solution values
        VectorType y(rEffectiveDofSet.size());
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
            rSolutionVector = VectorType(aux_size);
        }

        // Initialize solution vector with the constaints constant vector values
        // Note that we deliberately avoid using the copy constructor as we dont want to overwrite the constraints constant vector
        IndexPartition<IndexType>(rSolutionVector.size()).for_each([&](IndexType i){
            rSolutionVector[i] = rConstraintsConstantVector[i];
        });

        // Compute the solution vector as x = T * y + q
        rConstraintsRelationMatrix.SpMV(1.0, y, 1.0, rSolutionVector); // Note that this performs the operation x = A*y + q
    }

    /**
     * @brief Clear class content
     * This method wipes the content of current builder class
     */
    virtual void Clear()
    {

        mpModelPart = nullptr;

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
