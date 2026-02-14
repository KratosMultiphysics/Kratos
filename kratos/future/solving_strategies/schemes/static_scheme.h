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
template<class TLinearAlgebra>
class StaticScheme : public ImplicitScheme<TLinearAlgebra>
{
public:

    // FIXME: Does not work... ask @Charlie
    // /// Add scheme to Kratos registry
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.KratosMultiphysics", StaticScheme, StaticScheme, TLinearAlgebra)
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.All", StaticScheme, StaticScheme, TLinearAlgebra)

    ///@name Type Definitions
    ///@{

    /// Pointer definition of StaticScheme
    KRATOS_CLASS_POINTER_DEFINITION(StaticScheme);

    /// The definition of the current class
    using BaseType = ImplicitScheme<TLinearAlgebra>;

    /// Index type definition
    using IndexType = typename TLinearAlgebra::IndexType;

    /// Data type definition
    using DataType = typename TLinearAlgebra::DataType;

    /// Vector type definition
    using VectorType = typename TLinearAlgebra::VectorType;

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
        return Kratos::make_shared<StaticScheme<TLinearAlgebra>>(rModelPart, ThisParameters);
    }

    typename BaseType::Pointer Clone() override
    {
        return Kratos::make_shared<StaticScheme<TLinearAlgebra>>(*this) ;
    }

    void Predict(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        KRATOS_TRY

        // If needed, reset the DOF sets before applying the constraints and the prediction
        if (this->GetReformDofsAtEachStep()) {
            this->InitializeLinearSystem(rImplicitStrategyData);
        }

        // Applying constraints if needed
        const auto& r_model_part = this->GetModelPart();
        const auto& r_comm = r_model_part.GetCommunicator().GetDataCommunicator();
        auto& r_constraints = r_model_part.MasterSlaveConstraints();
        const unsigned int n_constraints_loc = r_constraints.size();
        const unsigned int n_constraints_glob = r_comm.SumAll(n_constraints_loc);

        if (n_constraints_glob != 0) {
            // Assemble constraints constant vector and apply it to the DOF set
            // Note that the constraints constant vector is applied only once in here as we then solve for the solution increment
            auto p_constraints_T = rImplicitStrategyData.pGetConstraintsT();
            auto p_constraints_Q = rImplicitStrategyData.pGetConstraintsQ();
            this->BuildMasterSlaveConstraints(rImplicitStrategyData);

            // Fill the current values vector considering the master-slave constraints
            // Note that this already accounts for the Dirichlet BCs affecting the effective DOF set
            auto& r_dof_set = *(rImplicitStrategyData.pGetDofSet());
            auto& r_eff_dof_set = *(rImplicitStrategyData.pGetEffectiveDofSet());
            VectorType x(r_dof_set.size());
            (this->GetBuilder()).CalculateSolutionVector(r_eff_dof_set, *p_constraints_T, *p_constraints_Q, x);

            // Update DOFs with solution values
            block_for_each(r_dof_set, [&x](DofType& rDof){
                rDof.GetSolutionStepValue() = x[rDof.EquationId()];
            });

            // If the mesh is to be updated, call the MoveMesh() method
            if (this->GetMoveMesh()) {
                this->MoveMesh();
            }
        }

        KRATOS_CATCH("")
    }

    void Update(ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData) override
    {
        KRATOS_TRY

        // Get linear system arrays
        auto& r_dx = *(rImplicitStrategyData.pGetLinearSystem()->pGetVector(Future::DenseVectorTag::Dx));
        auto& r_eff_dx = *(rImplicitStrategyData.pGetEffectiveLinearSystem()->pGetVector(Future::DenseVectorTag::Dx));
        auto& r_dof_set = *(rImplicitStrategyData.pGetDofSet());
        auto& r_eff_dof_set = *(rImplicitStrategyData.pGetEffectiveDofSet());

        // First update the constraints only DOFs with the effective solution vector
        this->UpdateConstraintsOnlyDofs(r_eff_dx, r_dof_set, r_eff_dof_set);

        // Get the solution update vector from the effective one
        this->CalculateUpdateVector(rImplicitStrategyData);

        // Update DOFs with solution values (note that we solve for the increments)
        block_for_each(r_dof_set, [&r_dx](DofType& rDof){
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
