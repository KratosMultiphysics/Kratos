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

    // /// Assembly helper type
    // using AssemblyHelperType = Future::AssemblyHelper<ThreadLocalStorage, TSparseMatrixType, TSparseVectorType, TSparseGraphType>;

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
    void InitializeSolutionStep(
        DofsArrayType& rDofSet,
        typename TSparseMatrixType::Pointer& rpA,
        typename TSparseVectorType::Pointer& rpDx,
        typename TSparseVectorType::Pointer& rpB,
        const bool ReformDofSet = true) override
    {
        KRATOS_TRY

        // Check if the InitializeSolutionStep has been already performed
        if (!this->GetSchemeSolutionStepIsInitialized()) {
            // Set up the system
            BuiltinTimer system_construction_time;
            if (!(this->GetDofSetIsInitialized()) || ReformDofSet) {
                // Setting up the DOFs list to be solved
                BuiltinTimer setup_dofs_time;
                this->SetUpDofArray(rDofSet);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Setup DOFs Time: " << setup_dofs_time << std::endl;

                // Set up the equation ids
                BuiltinTimer setup_system_time;
                const SizeType eq_system_size = this->SetUpSystemIds(rDofSet);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Set up system time: " << setup_system_time << std::endl;
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "Equation system size: " << eq_system_size << std::endl;

                // Allocating the system vectors to their correct sizes
                BuiltinTimer system_matrix_resize_time;
                this->ResizeAndInitializeVectors(rDofSet, rpA, rpDx, rpB);
                KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() > 0) << "System matrix resize time: " << system_matrix_resize_time << std::endl;
            } else {
                // Set up the equation ids (note that this needs to be always done)
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
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b) override
    {
        KRATOS_TRY

        // Internal solution loop check to avoid repetitions
        KRATOS_ERROR_IF_NOT(this->GetSchemeIsInitialized()) << "Initialize needs to be performed. Call Initialize() once before the solution loop." << std::endl;
        KRATOS_ERROR_IF_NOT(this->GetSchemeSolutionStepIsInitialized()) << "InitializeSolutionStep needs to be performed. Call InitializeSolutionStep() before Predict()." << std::endl;

        // If the mesh is to be updated, call the MoveMesh() method
        if (this->GetMoveMesh()) {
            this->MoveMesh();
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
    virtual void Update(
        DofsArrayType& rDofSet,
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b) override
    {
        KRATOS_TRY

        // Update DOFs with solution values (note that we solve for the increments)
        block_for_each(rDofSet, [&Dx](DofType& rDof){
            if (rDof.IsFree()) {
                rDof.GetSolutionStepValue() += Dx[rDof.EquationId()];
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
