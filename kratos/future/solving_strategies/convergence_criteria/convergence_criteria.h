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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    explicit ConvergenceCriteria()
    {
        SetEchoLevel(0);
    }

    /// Constructor with Parameters
    explicit ConvergenceCriteria(
        ModelPart& rModelPart,
        Kratos::Parameters ThisParameters)
        : mpModelPart(&rModelPart)
    {
        // Validate and assign defaults
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
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
     * @param ThisParameters The configuration parameters
     * @return A pointer to the new instance
     */
    virtual typename ClassType::Pointer Create(Parameters ThisParameters) const
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
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
    virtual bool IsConverged(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData)
    {
        KRATOS_ERROR << "Calling the base class IsConverged method. This should be implemented in the derived class." << std::endl;
        return false;
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rImplicitStrategyData Data container of the implicit strategy
     */
    virtual void Initialize(ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
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

    /**
     * @brief This method computes the norm of the given serial system vector
     * @details Note that only the free DOFs are considered in the norm calculation.
     * @param rDofSet Reference to the container of the problem's DOFs
     * @param rVector The vector whose norm is to be computed
     * @return A pair containing the number of free DOFs and the norm of the vector
     */
    std::pair<std::size_t, DataType> CalculateVectorNorm(
        const DofsArrayType& rDofSet,
        const SystemVector<DataType, IndexType>& rVector)
    {
        // Define custom reduction for parallel computation
        // First item in the reduction tuple: sum of the squared norm of the variation of the DOFs
        // Second item in the reduction tuple: number of free DOFs
        using CustomReductionType = CombinedReduction<SumReduction<DataType>,SumReduction<std::size_t>>;

        // Loop over Dofs and add the contribution of each free DOF to the norm
        DataType vector_norm;
        std::size_t n_free_dofs;
        std::tie(vector_norm, n_free_dofs) = block_for_each<CustomReductionType>(rDofSet, [&rVector](auto& rDof) {
            if (rDof.IsFree()) {
                return std::make_tuple(std::pow(rVector[rDof.EquationId()], 2), 1);
            } else {
                return std::make_tuple(DataType(), 0);
            }
        });

        return std::make_pair(n_free_dofs, std::sqrt(vector_norm));
    }

    /**
     * @brief This method computes the norm of the given distributed systemvector
     * @details Note that only the free DOFs are considered in the norm calculation.
     * @param rDofSet Reference to the container of the problem's DOFs
     * @param rVector The vector whose norm is to be computed
     * @return A pair containing the number of free DOFs and the norm of the vector
     * @note TODO: implementation to be tested when the distributed environment implementation is completed
     */
    std::pair<std::size_t, DataType> CalculateVectorNorm(
        const DofsArrayType& rDofSet,
        const DistributedSystemVector<DataType, IndexType>& rVector)
    {
        // Retrieve the data communicator
        const auto& r_data_communicator = mpModelPart->GetCommunicator().GetDataCommunicator();

        // Define custom reduction for parallel computation
        // First item in the reduction tuple: sum of the squared norm of the variation of the DOFs
        // Second item in the reduction tuple: number of free DOFs
        using CustomReductionType = CombinedReduction<SumReduction<DataType>,SumReduction<std::size_t>>;

        // Loop over Dofs and add the contribution of each free DOF to the norm
        // Note that only local contributions from each rank are considered
        DataType vector_norm;
        std::size_t n_free_dofs;
        std::tie(vector_norm, n_free_dofs) = block_for_each<CustomReductionType>(rDofSet, [&rVector](auto& rDof) {
            IndexType gl_eq_id = rDof.EquationId();
            if (rVector.GetNumbering().IsLocal(gl_eq_id)) {
                if (rDof.IsFree()) {
                    const auto& r_local_data = rVector.GetLocalData();
                    IndexType loc_eq_id = rVector.GetNumbering().LocalId(gl_eq_id);
                    return std::make_tuple(std::pow(r_local_data[loc_eq_id], 2), 1);
                } else {
                    return std::make_tuple(DataType(), 0);
                }
            }
        });

        // Communicator reduction
        n_free_dofs = r_data_communicator.SumAll(n_free_dofs);
        vector_norm = std::sqrt(r_data_communicator.SumAll(vector_norm));

        return std::make_pair(n_free_dofs, vector_norm);
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

    int mEchoLevel; /// The echo level

    ModelPart* mpModelPart = nullptr; /// The pointer to the model part

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

