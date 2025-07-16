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

    typename TSystemVectorType::Pointer pEffectiveQ = nullptr; // Linear system constraints total constant vector

    typename TSparseMatrixType::Pointer pConstraintsT = nullptr; // Master-slave constraints relation matrix

    typename TSystemVectorType::Pointer pConstraintsQ = nullptr; // Master-slave constraints constant vector

    typename TSparseMatrixType::Pointer pDirichletT = nullptr; // Dirichlet constraints relation matrix

    typename TSystemVectorType::Pointer pDirichletQ = nullptr; // Dirichlet constraints constant vector

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
        if (pEffectiveQ != nullptr) {
            pEffectiveQ->Clear();
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
        if (pDirichletQ != nullptr) {
            pDirichletQ->Clear();
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
 * @brief Utility class for handling the build
 * @details This class collects the methods required to
 * handle the building of the sparse sytem matrices
 * @author Ruben Zorrilla
 */
template<class TThreadLocalStorage, class TSparseMatrixType, class TSystemVectorType, class TSparseGraphType>
class Builder
{
public:
    ///@name Type Definitions
    ///@{

    /// Build type definition
    enum class BuildType
    {
        Block,
        Elimination
    };

    /// Scaling type definition
    enum class ScalingType
    {
        NoScaling,
        NormDiagonal,
        MaxDiagonal,
        PrescribedDiagonal
    };

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
            "scaling_type" : "max_diagonal",
            "echo_level" : 0
        })");
        Settings.ValidateAndAssignDefaults(default_parameters);

        // Set build type
        const std::string build_type = Settings["name"].GetString();
        if (build_type == "block_builder") {
            mBuildType = BuildType::Block;
        } else if (build_type == "elimination_builder") {
            mBuildType = BuildType::Elimination;
        } else {
            KRATOS_ERROR << "Provided 'build_type' is '" << build_type << "'. Available options are:\n"
            << "\t- 'block_builder'\n"
            << "\t- 'elimination_builder'" << std::endl;
        }

        // Set scaling type
        const std::string scaling_type = Settings["scaling_type"].GetString();
        if (mBuildType == BuildType::Block) {
            if (scaling_type == "no_scaling") {
                mScalingType = ScalingType::NoScaling;
            } else if (scaling_type == "norm_diagonal") {
                mScalingType = ScalingType::NormDiagonal;
            } else if (scaling_type == "max_diagonal") {
                mScalingType = ScalingType::MaxDiagonal;
            } else if (scaling_type == "prescribed_diagonal") {
                mScalingType = ScalingType::PrescribedDiagonal;
            }
        } else {
            mScalingType = ScalingType::NoScaling; // Note that scaling makes no sense in the elimination build
        }

        // Set verbosity level
        mEchoLevel = Settings["echo_level"].GetInt();
    }

    virtual ~Builder() = default;

    ///@}
    ///@name Operations
    ///@{

    //TODO: In the future, add the one with
    // - LHS alone
    // - MassMatrix
    // - DampingMatrix
    // TODO: To be discussed in the future. It would be great to have one with references. This will require using the move constructors
    virtual void ResizeAndInitializeVectors(
        const TSparseGraphType& rSparseGraph,
        const typename DofsArrayType::Pointer pDofSet,
        const typename DofsArrayType::Pointer pEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class 'ResizeAndInitializeVectors'." << std::endl;
    }

    virtual void ConstructMasterSlaveConstraintsStructure(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer)
    {
        // Clear the provided effective DOFs map
        KRATOS_WARNING_IF("Builder", !rEffectiveDofSet.empty()) << "Provided effective DOFs set is not empty. About to clear it." << std::endl;
        rEffectiveDofSet.clear();

        //FIXME: Do the IsActiveConstraints in here and set a flag that stays "forever"

        // Check if there are constraints to build the effective DOFs map and the corresponding arrays
        const std::size_t n_constraints = rModelPart.NumberOfMasterSlaveConstraints();
        if (n_constraints) {

            // Auxiliary set to store the unordered effective DOFs (masters from constraints and standard ones)
            std::unordered_set<typename DofType::Pointer> effective_dofs_set;

            // Get the master / slave DOFs from the constraints
            std::unordered_map<typename DofType::Pointer, DofPointerVectorType> constraints_slave_dofs;
            const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
            for (IndexType i_const = 0; i_const < n_constraints; ++i_const) {
                // Get current constraint master and slave DOFs
                auto it_const = it_const_begin + i_const;
                const auto& r_slave_dofs = it_const->GetSlaveDofsVector();
                const auto& r_master_dofs = it_const->GetMasterDofsVector();

                // Add the slave DOFs to the slave map
                for (auto& rp_slave : r_slave_dofs) {
                    constraints_slave_dofs.insert(std::make_pair(rp_slave, r_master_dofs));
                }

                // Add the master DOFs to the effective DOFs set
                // Note that we initialize the system ids to zero as these will be overwritten later
                for (auto& rp_master : r_master_dofs) {
                    effective_dofs_set.insert(rp_master);
                }
            }

            // Loop the elements and conditions DOFs container to get the DOFs that are not slave
            for (IndexType i_dof = 0; i_dof < rDofSet.size(); ++i_dof) {
                // Get current DOF
                auto p_dof = *(rDofSet.ptr_begin() + i_dof);

                // Check if current DOF is slave by checking the slaves DOFs map
                // If not present in the slaves DOFs map it should be considered in the resolution of the system
                // Note that this includes masters DOFs or and standard DOFs (those not involved in any constraint)
                if (constraints_slave_dofs.find(p_dof) == constraints_slave_dofs.end()) {
                    // Add current DOF to the effective DOFs set (note that the std::unordered_set guarantees uniqueness)
                    effective_dofs_set.insert(p_dof);
                }
            }

            // Sort the effective DOFs before setting the equation ids
            // Note that we dereference the DOF pointers in order to use the greater operator from dof.h
            std::vector<typename DofType::Pointer> ordered_eff_dofs_vector(effective_dofs_set.begin(), effective_dofs_set.end());
            std::sort(
                ordered_eff_dofs_vector.begin(),
                ordered_eff_dofs_vector.end(),
                [](const typename DofType::Pointer& pA, const typename DofType::Pointer& pB){return *pA > *pB;});

            // Fill the effective DOFs PVS with the sorted effective DOFs container
            rEffectiveDofSet = DofsArrayType(ordered_eff_dofs_vector);

            // Set the effective DOFs equation ids based on the sorted list
            IndexType aux_dof_id = 0;
            for (IndexType i_dof = 0; i_dof < rEffectiveDofSet.size(); ++i_dof) {
                auto p_dof = *(rEffectiveDofSet.ptr_begin() + i_dof);
                p_dof->SetEffectiveEquationId(aux_dof_id);
                ++aux_dof_id;
            }

            // Clear the equation ids vectors
            // mSlaveIds.clear();
            // mMasterIds.clear();

            // Set up constraints matrix sparse graph (note that mEquationSystemSize is the DOF set size in the block build)
            KRATOS_ERROR_IF(this->GetEquationSystemSize() == 0) << "Equation system size is not set yet. Please call 'SetUpSystemIds' before this method." << std::endl;
            TSparseGraphType constraints_sparse_graph(this->GetEquationSystemSize());

            // Loop the elements and conditions DOFs container to add the slave entries to the graph
            for (IndexType i_dof = 0; i_dof < rDofSet.size(); ++i_dof) {
                // Get current DOF
                auto p_dof = *(rDofSet.ptr_begin() + i_dof);
                const IndexType i_dof_eq_id = p_dof->EquationId();

                // Check if current DOF is slave by checking the slaves DOFs map
                // If not present in the slaves DOFs map it should be considered an effective DOF
                // Note that here effective means a masters DOF or a DOF that does not involve any constraint
                auto i_dof_slave_find = constraints_slave_dofs.find(p_dof);
                if (i_dof_slave_find != constraints_slave_dofs.end()) { // Slave DOF
                    // // Add current slave DOF to slave equation ids list
                    // mSlaveIds.push_back(i_dof_eq_id);

                    // Add current slave DOF connectivities to the constraints sparse graph
                    // The slave rows eq ids come from the system ones while the column master ones are the above defined
                    for (auto& rp_master : i_dof_slave_find->second) {
                        auto eff_dof_find = rEffectiveDofSet.find(*rp_master);
                        KRATOS_ERROR_IF(eff_dof_find == rEffectiveDofSet.end()) << "Effective DOF cannot be find." << std::endl;
                        constraints_sparse_graph.AddEntry(i_dof_eq_id, eff_dof_find->EffectiveEquationId());
                    }
                } else { // Effective DOF
                    auto eff_dof_find = rEffectiveDofSet.find(*p_dof);
                    KRATOS_ERROR_IF(eff_dof_find == rEffectiveDofSet.end()) << "Effective DOF cannot be find." << std::endl;
                    // mMasterIds.push_back(eff_dof_find->second);
                    constraints_sparse_graph.AddEntry(i_dof_eq_id, eff_dof_find->EffectiveEquationId());
                }
            }

            // // Loop the effective DOFs container to add the remaining diagonal entries to the graph
            // for (IndexType i_dof = 0; i_dof < rEffectiveDofSet.size(); ++i_dof) {
            //     // Get current effective DOF
            //     auto p_dof = *(rEffectiveDofSet.ptr_begin() + i_dof);
            //     auto eff_dof_find = rEffectiveDofIdMap.find(p_dof);
            //     KRATOS_ERROR_IF(eff_dof_find == rEffectiveDofIdMap.end()) << "Effective DOF cannot be find." << std::endl;
            //     std::cout << "Effective DOF: " << eff_dof_find->second << " - " << eff_dof_find->second << std::endl;
            //     constraints_sparse_graph.AddEntry(eff_dof_find->second, eff_dof_find->second);
            // }

            // Allocate the constraints arrays (note that we are using the move assignment operator in here)
            auto p_aux_q = Kratos::make_shared<TSystemVectorType>(this->GetEquationSystemSize());
            rLinearSystemContainer.pConstraintsQ.swap(p_aux_q);

            auto p_aux_T = Kratos::make_shared<TSparseMatrixType>(constraints_sparse_graph);
            rLinearSystemContainer.pConstraintsT.swap(p_aux_T);
        } else {
            // If there are no constraints the effective DOF set is the standard one
            rEffectiveDofSet = rDofSet;

            // Set the DOFs' effective equation global ids to match the standard ones
            IndexPartition<IndexType>(rEffectiveDofSet.size()).for_each([&](IndexType Index) {
                auto it_dof = rEffectiveDofSet.begin() + Index;
                it_dof->SetEffectiveEquationId(it_dof->EquationId());
            });
        }
    }

    virtual void ConstructDirichletConstraintsStructure(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class 'ConstructDirichletConstraintsStructure'." << std::endl;
    }

    virtual void SetUpSparseMatrixGraph(TSparseGraphType& rSparseGraph)
    {
        KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Current equation system size is 0. Sparse graph cannot be set (call \'SetUpSystemIds\' first)." << std::endl;

        //TODO: This is always thrown as the mEquationSystemSize is already in the graph. I'd use the move (not implemented constructor)
        if (!rSparseGraph.IsEmpty()) {
            KRATOS_WARNING("Builder") << "Provided sparse graph is not empty and will be cleared." << std::endl;
            rSparseGraph.Clear();
        }

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

    virtual std::size_t SetUpSystemIds(const typename DofsArrayType::Pointer pDofSet)
    {
        // Set up the DOFs equation global ids
        IndexPartition<IndexType>(pDofSet->size()).for_each([&](IndexType Index) {
            auto it_dof = pDofSet->begin() + Index;
            it_dof->SetEquationId(Index);
        });

        // Save the size of the system based on the number of DOFs
        mEquationSystemSize = pDofSet->size();

        return mEquationSystemSize;
    }

    virtual void BuildDirichletConstraints(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class 'ApplyLinearSystemConstraints'." << std::endl;
    }

    virtual void ApplyLinearSystemConstraints(
        const DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class 'ApplyLinearSystemConstraints'." << std::endl;
    }

    virtual void ApplyMasterSlaveConstraints(
        typename TSparseMatrixType::Pointer& rpLhs,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSystemVectorType::Pointer& rpRhs,
        TSystemVectorType& rEffectiveRhs,
        TSystemVectorType& rDx,
        TSystemVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSystemVectorType& rConstraintsConstantVector)
    {
        //TODO: Do it as the other assembly functions
        if (mBuildType == BuildType::Block) {
            ApplyMasterSlaveConstraintsImplementation<BuildType::Block>(rpLhs, rpEffectiveLhs, rpRhs, rEffectiveRhs, rDx, rEffectiveDx, rConstraintsRelationMatrix, rConstraintsConstantVector);
        } else if (mBuildType == BuildType::Elimination) {
            ApplyMasterSlaveConstraintsImplementation<BuildType::Elimination>(rpLhs, rpEffectiveLhs, rpRhs, rEffectiveRhs, rDx, rEffectiveDx, rConstraintsRelationMatrix, rConstraintsConstantVector);
        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }
    }

    virtual void ApplyMasterSlaveConstraints(
        typename TSystemVectorType::Pointer& rpRhs,
        TSystemVectorType& rEffectiveRhs,
        TSystemVectorType& rDx,
        TSystemVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSystemVectorType& rConstraintsConstantVector)
    {
        //TODO: Do it as the other assembly functions
        if (mBuildType == BuildType::Block) {
            ApplyMasterSlaveConstraintsImplementation<BuildType::Block>(rpRhs, rEffectiveRhs, rDx, rEffectiveDx, rConstraintsRelationMatrix, rConstraintsConstantVector);
        } else if (mBuildType == BuildType::Elimination) {
            ApplyMasterSlaveConstraintsImplementation<BuildType::Elimination>(rpRhs, rEffectiveRhs, rDx, rEffectiveDx, rConstraintsRelationMatrix, rConstraintsConstantVector);
        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }
    }

    virtual double GetDiagonalScalingFactor(const TSparseMatrixType& rLHS) const
    {
        if (mScalingType == ScalingType::NoScaling) {
            return 1.0;
        } else if (mScalingType == ScalingType::NormDiagonal) {
            return rLHS.NormDiagonal();
            // return rLHS.NormDiagonal() / rLHS.size1(); //TODO: Decide which one
        } else if (mScalingType == ScalingType::MaxDiagonal) {
            return rLHS.MaxDiagonal();
            // return rLHS.MaxDiagonal() / rLHS.size1(); // TODO: Decide which one
        } else if (mScalingType == ScalingType::PrescribedDiagonal) {
            const auto& r_process_info = mpModelPart->GetProcessInfo();
            KRATOS_ERROR_IF_NOT(r_process_info.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined in ProcessInfo container. Please set 'BUILD_SCALE_FACTOR' variable." << std::endl;
            return r_process_info.GetValue(BUILD_SCALE_FACTOR);
        } else {
            KRATOS_ERROR << "Wrong scaling type." << std::endl;
        }
    }

    void CalculateSolutionVector(
        const DofsArrayType& rDofSet,
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
            auto p_dof_find = rEffectiveDofSet.find(*p_dof);
            KRATOS_ERROR_IF(p_dof_find == rEffectiveDofSet.end()) << "DOF cannot be found in effective DOF set." << std::endl;

            // Get value from DOF and set it in the auxiliary solution values vector
            // Note that the corresponding row is retrieved from the effective DOF ids map
            y[p_dof_find->EffectiveEquationId()] = p_dof->GetSolutionStepValue();
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

    virtual void Clear()
    {
        // Clear the system info
        mEquationSystemSize = 0;
    }

    ///@}
    ///@name Input and output
    ///@{

    std::size_t GetEquationSystemSize() const
    {
        return mEquationSystemSize;
    }

    const ModelPart& GetModelPart() const
    {
        return *mpModelPart;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    const ModelPart* mpModelPart = nullptr;

    BuildType mBuildType;

    ScalingType mScalingType;

    double mScalingValue;

    std::size_t mEchoLevel;

    std::size_t mEquationSystemSize = 0; /// Number of degrees of freedom of the problem to be solved

    ///@}
    ///@name Private Operations
    ///@{

    template <BuildType TBuildType>
    void ApplyMasterSlaveConstraintsImplementation(
        typename TSparseMatrixType::Pointer& rpLhs,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSystemVectorType::Pointer& rpRhs,
        TSystemVectorType& rEffectiveRhs,
        TSystemVectorType& rDx,
        TSystemVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSystemVectorType& rConstraintsConstantVector)
    {
        if constexpr (TBuildType == BuildType::Block) {
            // Get the effective size as the number of master DOFs
            // Note that this matches the number of columns in the constraints relation matrix
            const std::size_t n_master = rConstraintsRelationMatrix.size2();

            // Initialize the effective RHS
            if (rEffectiveRhs.size() != n_master) {
                rEffectiveRhs = std::move(TSystemVectorType(n_master));
            }
            rEffectiveRhs.SetValue(0.0);

            // Initialize the effective solution vector
            if (rEffectiveDx.size() != n_master) {
                rEffectiveDx = std::move(TSystemVectorType(n_master));
            }
            rEffectiveDx.SetValue(0.0);

            // Apply constraints to RHS
            rConstraintsRelationMatrix.TransposeSpMV(*rpRhs, rEffectiveRhs);

            // Apply constraints to LHS
            auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(*rpLhs, rConstraintsRelationMatrix);
            auto p_transT = AmgclCSRConversionUtilities::Transpose(rConstraintsRelationMatrix);
            rpEffectiveLhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);

            // // Compute the scale factor value
            // //TODO: think on how to make this user-definable
            // const double scale_factor = rpEffectiveLhs->NormDiagonal();

            // // Apply diagonal values on slave DOFs
            // IndexPartition<IndexType>(mSlaveIds.size()).for_each([&](IndexType Index){
            //     const IndexType slave_eq_id = mSlaveIds[Index];
            //     if (mInactiveSlaveDofs.find(slave_eq_id) == mInactiveSlaveDofs.end()) {
            //         (*rpEffectiveLhs)(slave_eq_id, slave_eq_id) = scale_factor;
            //         (*rpEffectiveRhs)[slave_eq_id] = 0.0;
            //     }
            // });
        } else if (TBuildType == BuildType::Elimination) {

        } else {
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unsupported build type.");
        }
    }

    template <BuildType TBuildType>
    void ApplyMasterSlaveConstraintsImplementation(
        typename TSystemVectorType::Pointer& rpRhs,
        TSystemVectorType& rEffectiveRhs,
        TSystemVectorType& rDx,
        TSystemVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSystemVectorType& rConstraintsConstantVector)
    {
        if constexpr (TBuildType == BuildType::Block) {
            // Get the effective size as the number of master DOFs
            // Note that this matches the number of columns in the constraints relation matrix
            const std::size_t n_master = rConstraintsRelationMatrix.size2();

            // Initialize the effective RHS
            if (rEffectiveRhs.size() != n_master) {
                rEffectiveRhs = std::move(TSystemVectorType(n_master));
            }
            rEffectiveRhs.SetValue(0.0);

            // Initialize the effective solution vector
            if (rEffectiveDx.size() != n_master) {
                rEffectiveDx = std::move(TSystemVectorType(n_master));
            }
            rEffectiveDx.SetValue(0.0);

            // Apply constraints to RHS
            rConstraintsRelationMatrix.TransposeSpMV(*rpRhs, rEffectiveRhs);

            // // Compute the scale factor value
            // //TODO: think on how to make this user-definable
            // const double scale_factor = rpEffectiveLhs->NormDiagonal();

            // // Apply diagonal values on slave DOFs
            // IndexPartition<IndexType>(mSlaveIds.size()).for_each([&](IndexType Index){
            //     const IndexType slave_eq_id = mSlaveIds[Index];
            //     if (mInactiveSlaveDofs.find(slave_eq_id) == mInactiveSlaveDofs.end()) {
            //         rEffectiveRhs[slave_eq_id] = 0.0;
            //     }
            // });
        } else if (TBuildType == BuildType::Elimination) {

        } else {
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unsupported build type.");
        }
    }

    ///@}
}; // Class Builder

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
