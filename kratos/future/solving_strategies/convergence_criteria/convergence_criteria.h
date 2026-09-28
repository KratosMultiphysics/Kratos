//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
#include "utilities/reduction_utilities.h"
#include "containers/system_vector.h"
#include "containers/distributed_system_vector.h"
#ifdef KRATOS_USE_FUTURE
#include "future/containers/linear_system_tags.h"
#include "future/solving_strategies/schemes/implicit_scheme.h"
#include "future/solving_strategies/strategies/implicit_strategy_data.h"
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
 * @class ConvergenceCriteria
 * @ingroup KratosCore
 * @brief This is the base class to define the different convergence criterion considered
 * @tparam TLinearAlgebra The linear algebra type
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
*/
template<class TLinearAlgebra>
class ConvergenceCriteria
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION(ConvergenceCriteria);

    /// The definition of the current class
    typedef ConvergenceCriteria< TLinearAlgebra > ClassType;

    /// The data type
    using DataType = typename TLinearAlgebra::DataType;

    /// The index type
    using IndexType = typename TLinearAlgebra::IndexType;

    /// The DOFs array type
    using DofsArrayType = typename ModelPart::DofsArrayType;

    // Scheme pointer type definition
    using SchemePointerType = typename ImplicitScheme<TLinearAlgebra>::Pointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    explicit ConvergenceCriteria()
    {
    }

    /// Constructor with Parameters
    explicit ConvergenceCriteria(
        ModelPart& rModelPart,
        Kratos::Parameters ThisParameters)
        : mpModelPart(&rModelPart)
    {
    }

    /// Copy constructor
    explicit ConvergenceCriteria( ConvergenceCriteria const& rOther) = delete;

    /// Destructor
    virtual ~ConvergenceCriteria() = default;

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates a new instance of the convergence criteria
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     * @return A pointer to the new instance
     */
    virtual typename ClassType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     */
    void SetEchoLevel(const int Level)
    {
        mEchoLevel = Level;
    }

    /**
     * @brief Checks if the solution is converged
     * @param rImplicitStrategyData Data container of the implicit strategy
     * @warning Must be defined on the derived classes
     * @return true if the solution is converged, false otherwise
     */
    virtual bool IsConverged(const ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData)
        {
            KRATOS_ERROR << "Calling the base class IsConverged method. This should be implemented in the derived class." << std::endl;
            return false;
        }

    /**
     * @brief This function initialize the convergence criteria
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void Initialize(const ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
    }

    /**
     * @brief This function initializes the solution step
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void InitializeSolutionStep(const ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
    }

    /**
     * @brief This function initializes the non-linear iteration
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void InitializeNonLinearIteration(const ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
    }

    /**
     * @brief This function finalizes the non-linear iteration
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void FinalizeNonLinearIteration(const ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
    }

    /**
     * @brief This function finalizes the solution step
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void FinalizeSolutionStep(const ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed on the input provided. Checks can be "expensive" as the function is designed to catch user's errors.
     * @return 0 all OK, 1 otherwise
     */
    virtual int Check()
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    virtual Parameters GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "name"       : "convergence_criteria",
            "echo_level" : 1
        })");
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "convergence_criteria";
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the Model Part object
     * Returns a reference to the model part the scheme is referring to
     * @return ModelPart& Reference to the scheme model part
     */
    ModelPart& GetModelPart()
    {
        return *mpModelPart;
    }

    /**
     * @brief Get the Model Part object
     * Returns a reference to the model part the scheme is referring to
     * @return const ModelPart& Reference to the scheme model part
     */
    const ModelPart& GetModelPart() const
    {
        return *mpModelPart;
    }

    /**
     * @brief This returns the level of echo for the solving strategy
     * @return Level of echo for the solving strategy
     */
    int GetEchoLevel() const
    {
        return mEchoLevel;
    }

    ///@}
    ///@name Inquiry
    ///@{

    virtual bool RequiresBuild() const
    {
        KRATOS_ERROR << "Calling base class 'RequiresBuild'." << std::endl;
        return false;
    }

    virtual LinearSystemTags::DenseVectorTag GetConvergenceCheckVectorTag() const
    {
        KRATOS_ERROR << "Calling base class 'GetConvergenceCheckVectorTag'." << std::endl;
        return LinearSystemTags::DenseVectorTag::Dx;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ConvergenceCriteria";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    ModelPart* mpModelPart = nullptr; /// The pointer to the model part

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
    virtual void AssignSettings(const Parameters ThisParameters)
    {
        mEchoLevel = ThisParameters["echo_level"].GetInt();
    }


    std::pair<DataType,std::size_t> CalculateConvergenceVectorNorm(const ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) const
    {
        // Get the effective vector from which convergence will be check and its corresponding effective DOF set
        const auto conv_vect_tag = this->GetConvergenceCheckVectorTag();
        const auto& r_eff_dof_set = *(rImplicitStrategyData.pGetEffectiveDofSet());
        const auto& r_eff_conv_vect = *(rImplicitStrategyData.pGetEffectiveLinearSystem()->pGetVector(conv_vect_tag));

        // Custom reduction to return both the norm and the number of free DOFs
        using CustomReduction = CombinedReduction<SumReduction<DataType>,SumReduction<std::size_t>>;

        // Allocate output variables for reduction
        DataType conv_vect_norm;
        std::size_t n_free_dofs;

        // Loop the effective vector to calculate the norm considering only the values associated to the free DOFs
        // Note that here it is assumed that the effective DOF equation ids are smaller than the problem size (i.e., the size of the effective vector)
        const std::size_t pb_size = r_eff_conv_vect.size();
        const std::size_t n_eff_dofs = r_eff_dof_set.size();
        std::tie(conv_vect_norm, n_free_dofs) = (IndexPartition<IndexType>(n_eff_dofs)).template for_each<CustomReduction>([&](IndexType Index) {
            const auto it_eff_dof = r_eff_dof_set.begin() + Index;
            const std::size_t eff_eq_id = it_eff_dof->EffectiveEquationId();
            if (it_eff_dof->IsFree() && eff_eq_id < pb_size) {
                const DataType value = r_eff_conv_vect[eff_eq_id];
                return std::make_tuple(value * value, 1);
            } else {
                return std::make_tuple(DataType(), 0);
            }
        });

        // Return a pair with the vector norm and the number of free DOFs found during its calculation
        return std::make_pair(sqrt(conv_vect_norm), n_free_dofs);
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

    int mEchoLevel = 0; /// The echo level

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

}; /// Class ConvergenceCriteria
} // namespace Kratos::Future.

