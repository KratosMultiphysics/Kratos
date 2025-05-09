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
#include "future/solving_strategies/schemes/assembly_helper.h"
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
 * @class StaticScheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of the basic tasks that are needed by the solution strategy.
 * @details It is intended to be the place for tailoring the solution strategies to problem specific tasks.
 * @author Ruben Zorrilla
 */
template<class TSparseMatrixType, class TSparseVectorType, class TSparseGraphType>
class StaticScheme : public ImplicitScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>
{
public:
    // FIXME: Does not work... ask @Charlie
    // /// Add scheme to Kratos registry
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.KratosMultiphysics", StaticScheme, StaticScheme, TSparseMatrixType, TSparseVectorType)
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.All", StaticScheme, StaticScheme, TSparseMatrixType, TSparseVectorType)

    ///@name Type Definitions
    ///@{

    /// Pointer definition of StaticScheme
    KRATOS_CLASS_POINTER_DEFINITION(StaticScheme);

    /// The definition of the current class
    using BaseType = ImplicitScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>;

    /// Size type definition
    using SizeType = std::size_t;

    /// Index type definition
    using IndexType = typename TSparseMatrixType::IndexType;

    /// Data type definition
    using DataType = typename TSparseMatrixType::DataType;

    /// DoF type definition
    using DofType = Dof<DataType>;

    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    /// Effective DOFs map type definition
    using EffectiveDofsMapType = typename BaseType::EffectiveDofsMapType;

    // /// TLS type definition
    // struct ThreadLocalStorage //FIXME: This will be ImplicitThreadLocalStorage when we create the implicit scheme --> Also we need to move them out of here.
    // {
    //     // Local LHS contribution
    //     DenseMatrix<DataType> LocalLhs;

    //     // Local RHS constribution
    //     DenseVector<DataType> LocalRhs;

    //     // Vector containing the localization in the system of the different terms
    //     Element::EquationIdVectorType LocalEqIds;
    // };

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     * @details Initializes the flags
     */
    explicit StaticScheme() = default;

    /**
     * @brief Constructor with Parameters
     */
    explicit StaticScheme(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart, ThisParameters)
    {
        // Validate default parameters
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /** Copy Constructor.
     */
    explicit StaticScheme(StaticScheme& rOther)
    {
        //TODO: Check this... particularly the mpAssemblyHelper pointer
    }

    /** Destructor.
     */
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
        return Kratos::make_shared<StaticScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>>(rModelPart, ThisParameters);
    }

    typename BaseType::Pointer Clone() override
    {
        return Kratos::make_shared<StaticScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>>(*this) ;
    }

    /**
     * @brief Function called once at the beginning of each solution step.
     * @details The basic operations to be carried in there are the following:
     * - managing variables to be kept constant over the time step (for example time-Scheme constants depending on the actual time step)
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */

    /**
     * @brief Function called once at the beginning of each solution step
     * The basic operations to be carried out in here are the following:
     * 1) Set up the DOF array from the element and conditions DOFs (SetUpDofArray)
     * 2) Set up the system ids (i.e., the DOFs equation id) for the element and condition DOF array
     * 3) Construct the master slave constraints structure (this includes the effective DOF array, the effective DOF id map and the relation and constant arrays)
     * 4) Allocate the memory for the system arrays (note that this implies building the sparse matrix graph)
     * 5) Call the InitializeSolutionStep of all entities
     * Further operations might be required depending on the time integration scheme
     * Note that steps from 1 to 4 can be done once if the DOF set does not change (i.e., the mesh and the constraints active/inactive status do not change in time)
     * @param rDofSet The array of DOFs from elements and conditions
     * @param rEffectiveDofSet The array of DOFs to be solved after the application of constraints
     * @param rEffectiveDofIdMap A map relating each effective DOF to its effective id
     * @param rpA The system left hand side matrix
     * @param rpEffectiveLhs The effective left hand side matrix
     * @param rpB The system right hand side vector
     * @param rpEffectiveRhs The effective right hand side vector
     * @param rpDx The solution update vector
     * @param rpEffectiveDx The effective solution update vector
     * @param rConstraintsRelationMatrix The assembled constraints relation matrix (i.e. T)
     * @param rConstraintsConstantVector The assembled constraints constant vector
     * @param ReformDofSet Flag to indicate if the DOFs have changed and need to be updated
     */
    void InitializeSolutionStep(
        DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        typename TSparseMatrixType::Pointer& rpA,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSparseVectorType::Pointer& rpB,
        typename TSparseVectorType::Pointer& rpEffectiveRhs,
        typename TSparseVectorType::Pointer& rpDx,
        typename TSparseVectorType::Pointer& rpEffectiveDx,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector,
        const bool ReformDofSets = true) override
    {
        KRATOS_TRY

        // Check if the InitializeSolutionStep has been already performed
        if (!this->GetSchemeSolutionStepIsInitialized()) {
            // Set up the system
            BuiltinTimer system_construction_time;
            if (!(this->GetDofSetIsInitialized()) || ReformDofSets) {
                // Setting up the DOFs list to be solved
                BuiltinTimer setup_dofs_time;
                this->SetUpDofArray(rDofSet);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Setup DOFs Time: " << setup_dofs_time << std::endl;

                // Set up the equation ids
                BuiltinTimer setup_system_time;
                const SizeType eq_system_size = this->SetUpSystemIds(rDofSet);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Set up system time: " << setup_system_time << std::endl;
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Equation system size: " << eq_system_size << std::endl;

                // Set up the system constraints
                BuiltinTimer constraints_construction_time;
                this->ConstructMasterSlaveConstraintsStructure(rDofSet, rEffectiveDofSet, rEffectiveDofIdMap, rConstraintsRelationMatrix, rConstraintsConstantVector);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Constraints construction time: " << constraints_construction_time << std::endl;

                // Allocating the system vectors to their correct sizes
                BuiltinTimer system_matrix_resize_time;
                this->ResizeAndInitializeVectors(rDofSet, rEffectiveDofSet, rpA, rpEffectiveLhs, rpB, rpEffectiveRhs, rpDx, rpEffectiveDx);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "System matrix resize time: " << system_matrix_resize_time << std::endl;

            } else {
                // Set up the equation ids (note that this needs to be always done as the fixity may have changed and this can affect some build types)
                BuiltinTimer setup_system_time;
                const SizeType eq_system_size = this->SetUpSystemIds(rDofSet);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Set up system time: " << setup_system_time << std::endl;
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Equation system size: " << eq_system_size << std::endl;
            }
            KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "System construction time: " << system_construction_time << std::endl;

            // Initializes solution step for all of the elements, conditions and constraints
            EntitiesUtilities::InitializeSolutionStepAllEntities(this->GetModelPart());

            // Set the flag to avoid calling this twice
            this->SetSchemeSolutionStepIsInitialized(true); // TODO: Discuss with the KTC if these should remain or not
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Performing the prediction of the solution.
     * @warning Must be defined in derived classes
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void Predict(
        DofsArrayType& rDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        TSparseMatrixType& rA,
        TSparseVectorType& rb,
        TSparseVectorType& rDx,
        TSparseVectorType& rEffectiveDx,
        TSparseMatrixType& rConstraintsRelationMatrix) override
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
            const auto& r_process_info = r_model_part.GetProcessInfo();

            block_for_each(r_constraints, [&r_process_info](MasterSlaveConstraint& rConstraint){
                rConstraint.ResetSlaveDofs(r_process_info);
            });

            block_for_each(r_constraints, [&r_process_info](MasterSlaveConstraint& rConstraint){
                rConstraint.Apply(r_process_info);
            });

            // The following is needed since we need to eventually compute time derivatives after applying master-slave relations
            rDx.SetValue(0.0);
            rEffectiveDx.SetValue(0.0);
            this->Update(rDofSet, rEffectiveDofIdMap, rA, rb, rDx, rEffectiveDx, rConstraintsRelationMatrix);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Performing the update of the solution.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void Update(
        DofsArrayType& rDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        TSparseMatrixType& rA,
        TSparseVectorType& rb,
        TSparseVectorType& rDx,
        const TSparseVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix) override
    {
        KRATOS_TRY

        // First update the constraints loose DOFs with the effective solution vector
        this->UpdateConstraintsLooseDofs(rEffectiveDx, rDofSet, rEffectiveDofIdMap);

        // Get the solution update vector from the effective one
        this->CalculateUpdateVector(rConstraintsRelationMatrix, rEffectiveDx, rDx);

        // Update DOFs with solution values (note that we solve for the increments)
        block_for_each(rDofSet, [&rDx](DofType& rDof){
            if (rDof.IsFree()) {
                rDof.GetSolutionStepValue() += rDx[rDof.EquationId()];
            }
        });

        // If the mesh is to be updated, call the MoveMesh() method
        if (this->GetMoveMesh()) {
            this->MoveMesh();
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Functions to be called to prepare the data needed for the output of results.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void CalculateOutputData(
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b) override
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

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        //TODO:
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
