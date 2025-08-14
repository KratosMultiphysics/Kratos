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
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "spaces/kratos_space.h"
#include "utilities/builtin_timer.h"
#include "utilities/dof_utilities/dof_array_utilities.h"
#include "utilities/entities_utilities.h"
#include "utilities/openmp_utils.h" //TODO: SOME FILES INCLUDING scheme.h RELY ON THIS. LEAVING AS FUTURE TODO.
#include "utilities/parallel_utilities.h"
#include "utilities/timer.h"

#ifdef KRATOS_USE_FUTURE
#include "future/linear_solvers/amgcl_solver.h"
#include "future/solving_strategies/schemes/implicit_scheme.h"
#include "future/solving_strategies/builders/builder.h"
#endif

namespace Kratos::Future
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
 * @brief Implicit scheme TLS type definition
 * Thread Local Storage container to be used in the parallel assembly of implicit problems
 * @tparam DataType data type of the problem to be solved
 */
template<class TDataType = double >
struct StaticThreadLocalStorage
{
    // Local LHS contribution
    DenseMatrix<TDataType> LocalMatrix;

    // Local RHS constribution
    DenseVector<TDataType> LocalVector;

    // Vector containing the localization in the system of the different terms
    Element::EquationIdVectorType LocalEqIds;

    // Vector containing the slave equation ids
    MasterSlaveConstraint::EquationIdVectorType SlaveEqIds;

    // Vector containing the master equation ids
    MasterSlaveConstraint::EquationIdVectorType MasterEqIds;
};

/**
 * @class StaticScheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of the static scheme
 * @author Ruben Zorrilla
 */
template<class TSparseMatrixType, class TSystemVectorType, class TSparseGraphType>
class StaticScheme : public ImplicitScheme<TSparseMatrixType, TSystemVectorType, TSparseGraphType>
{
public:

    // FIXME: Does not work... ask @Charlie
    // /// Add scheme to Kratos registry
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.KratosMultiphysics", StaticScheme, StaticScheme, TSparseMatrixType, TSystemVectorType)
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.All", StaticScheme, StaticScheme, TSparseMatrixType, TSystemVectorType)

    ///@name Type Definitions
    ///@{

    /// Pointer definition of StaticScheme
    KRATOS_CLASS_POINTER_DEFINITION(StaticScheme);

    /// The definition of the current class
    using BaseType = ImplicitScheme<TSparseMatrixType, TSystemVectorType, TSparseGraphType>;

    /// Index type definition
    using IndexType = typename TSparseMatrixType::IndexType;

    /// Data type definition
    using DataType = typename TSparseMatrixType::DataType;

    /// TLS type
    using TLSType = StaticThreadLocalStorage<DataType>;

    /// DoF type definition
    using DofType = Dof<DataType>;

    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Default constructor
    explicit StaticScheme() = default;

    /// @brief Constructor with parameters
    /// @param rModelPart Reference to the model part
    /// @param ThisParameters Parameters object encapsulating the settings
    explicit StaticScheme(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart, ThisParameters)
    {
        // Validate default parameters
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /// @brief Copy constructor
    /// @param rOther Other StaticScheme
    explicit StaticScheme(StaticScheme& rOther)
    {
        //TODO: Check this... particularly the mpBuilder pointer
    }

    /// @brief Destructor
    virtual ~StaticScheme() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const override
    {
        return Kratos::make_shared<StaticScheme<TSparseMatrixType, TSystemVectorType, TSparseGraphType>>(rModelPart, ThisParameters);
    }

    typename BaseType::Pointer Clone() override
    {
        return Kratos::make_shared<StaticScheme<TSparseMatrixType, TSystemVectorType, TSparseGraphType>>(*this) ;
    }

    void InitializeSolutionStep(
        DofsArrayType::Pointer pDofSet,
        DofsArrayType::Pointer pEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer,
        const bool ReformDofSets = true) override
    {
        KRATOS_TRY

        // Check if the InitializeSolutionStep has been already performed
        if (!this->GetSchemeSolutionStepIsInitialized()) {
            // Set up the system
            if (!(this->GetDofSetIsInitialized()) || ReformDofSets) {
                // Setting up the DOFs list
                BuiltinTimer setup_dofs_time;
                DofArrayUtilities::SlaveToMasterDofsMap slaves_to_master_dofs_map;
                auto [eq_system_size, eff_eq_system_size] = this->SetUpDofArrays(pDofSet, pEffectiveDofSet, slaves_to_master_dofs_map);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Setup DOFs Time: " << setup_dofs_time << std::endl;

                // Set up the equation ids
                BuiltinTimer setup_system_ids_time;
                this->SetUpSystemIds(pDofSet, pEffectiveDofSet);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Set up system time: " << setup_system_ids_time << std::endl;
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Equation system size: " << eq_system_size << std::endl;
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Effective equation system size: " << eff_eq_system_size << std::endl;

                // Allocating the system constraints arrays
                BuiltinTimer constraints_allocation_time;
                this->AllocateLinearSystemConstraints(pDofSet, pEffectiveDofSet, slaves_to_master_dofs_map, rLinearSystemContainer);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Linear system constraints allocation time: " << constraints_allocation_time << std::endl;

                // Allocating the system vectors to their correct sizes
                BuiltinTimer linear_system_allocation_time;
                this->AllocateLinearSystemArrays(pDofSet, pEffectiveDofSet, rLinearSystemContainer);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Linear system allocation time: " << linear_system_allocation_time << std::endl;
            }

            // Initializes solution step for all of the elements, conditions and constraints
            EntitiesUtilities::InitializeSolutionStepAllEntities(this->GetModelPart());

            // Set the flag to avoid calling this twice
            this->SetSchemeSolutionStepIsInitialized(true); // TODO: Discuss with the KTC if these should remain or not
        }

        KRATOS_CATCH("")
    }

    void Predict(
        DofsArrayType::Pointer pDofSet,
        DofsArrayType::Pointer pEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer) override
    {
        KRATOS_TRY

        // Internal solution loop check to avoid repetitions
        KRATOS_ERROR_IF_NOT(this->GetSchemeIsInitialized()) << "Initialize needs to be performed. Call Initialize() once before the solution loop." << std::endl;
        KRATOS_ERROR_IF_NOT(this->GetSchemeSolutionStepIsInitialized()) << "InitializeSolutionStep needs to be performed. Call InitializeSolutionStep() before Predict()." << std::endl;

        // Applying constraints if needed
        const auto& r_model_part = this->GetModelPart();
        const auto& r_comm = r_model_part.GetCommunicator().GetDataCommunicator();
        auto& r_constraints = r_model_part.MasterSlaveConstraints();
        const std::size_t n_constraints_loc = r_constraints.size();
        const std::size_t n_constraints_glob = r_comm.SumAll(n_constraints_loc);

        if (n_constraints_glob != 0) {
            // Assemble constraints constant vector and apply it to the DOF set
            // Note that the constraints constant vector is applied only once in here as we then solve for the solution increment
            auto p_constraints_T = rLinearSystemContainer.pConstraintsT;
            auto p_constraints_Q = rLinearSystemContainer.pConstraintsQ;
            this->BuildMasterSlaveConstraints(*pDofSet, *pEffectiveDofSet, rLinearSystemContainer);

            // Fill the current values vector considering the master-slave constraints
            // Note that this already accounts for the Dirichlet BCs affecting the effective DOF set
            TSystemVectorType x(pDofSet->size());
            (this->GetBuilder()).CalculateSolutionVector(*pEffectiveDofSet, *p_constraints_T, *p_constraints_Q, x);

            // Update DOFs with solution values
            block_for_each(*pDofSet, [&x](DofType& rDof){
                rDof.GetSolutionStepValue() = x[rDof.EquationId()];
            });

            // If the mesh is to be updated, call the MoveMesh() method
            if (this->GetMoveMesh()) {
                this->MoveMesh();
            }
        }

        KRATOS_CATCH("")
    }

    void Update(
        DofsArrayType &rDofSet,
        DofsArrayType &rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer) override
    {
        KRATOS_TRY

        // Get linear system arrays
        auto& r_dx = *(rLinearSystemContainer.pDx);
        auto& r_eff_dx = *(rLinearSystemContainer.pEffectiveDx);

        // First update the constraints loose DOFs with the effective solution vector
        this->UpdateConstraintsLooseDofs(r_eff_dx, rDofSet, rEffectiveDofSet);

        // Get the solution update vector from the effective one
        this->CalculateUpdateVector(rLinearSystemContainer);

        // Update DOFs with solution values (note that we solve for the increments)
        block_for_each(rDofSet, [&r_dx](DofType& rDof){
            if (rDof.IsFree()) {
                rDof.GetSolutionStepValue() += r_dx[rDof.EquationId()];
            }
        });

        // If the mesh is to be updated, call the MoveMesh() method
        if (this->GetMoveMesh()) {
            this->MoveMesh();
        }

        KRATOS_CATCH("")
    }

    void CalculateOutputData(LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer) override
    {
        KRATOS_TRY

        //TODO: Think about creating a datavalue container in here that we can access from outside.

        KRATOS_CATCH("")
    }

    int Check() const override
    {
        KRATOS_TRY

        int check = BaseType::Check();

        return check;

        KRATOS_CATCH("");
    }

    Parameters GetDefaultParameters() const override
    {
        // Current class default parameters
        Parameters default_parameters = Parameters(R"({
            "name" : "static_scheme"
        })");

        // Add base class default parameters
        default_parameters.RecursivelyAddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "static_scheme";
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "StaticScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void AssignSettings(const Parameters ThisParameters) override
    {
        // Assign base scheme settings
        BaseType::AssignSettings(ThisParameters);
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

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
}; // Class Scheme

} // namespace Kratos::Future.
