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
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "future/linear_operators/linear_operator.h"

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
 * @brief Base class for preconditioners
 * @details This class defines the interface for preconditioners to be used in iterative linear solvers.
 * It provides virtual methods for initialization, application, and finalization of the preconditioner.
 * If additional data is required (i.e., reference model part and DOFs container), methods to set it are also provided.
 * The actual implementation of the preconditioner must be done in derived classes as this is a do nothing implementation.
 * @tparam TLinearAlgebra The struct containing the linear algebra types
 */
template<class TLinearAlgebra>
class Preconditioner
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(Preconditioner);

    using DataType = typename TLinearAlgebra::DataType;

    using VectorType = typename TLinearAlgebra::VectorType;

    using DofsArrayType = ModelPart::DofsArrayType;

    using LinearOperatorType = LinearOperator<TLinearAlgebra>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Preconditioner() = default;

    /// Constructor with parameters
    Preconditioner(Parameters Settings)
    {
    }

    /// Copy constructor.
    Preconditioner(const Preconditioner& Other) = default;

    /// Destructor.
    virtual ~Preconditioner() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Preconditioner& operator=(const Preconditioner& Other)
    {
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initialize the preconditioner with the linear operator.
     * This method is called once at the beginning of the solve to set up the preconditioner.
     * Note that a LinearOperator (@see LinearOperator) enables both matrix-free and matrix-based implementations.
     * @param rpLinearOperator Unique pointer to the linear operator representing the system matrix.
     */
    virtual void Initialize(const typename LinearOperatorType::UniquePointer& rpLinearOperator)
    {
    }

    /**
     * @brief Initialize the preconditioner for the current solution step.
     * This method is called every time the coefficients change in the system, for example at the beginning of each solve.
     * For example, if we are implementing a direct solver, this is the place to do the factorization
     * so that then the backward substitution can be performed effectively more than once.
     * @param rpLinearOperator Unique pointer to the linear operator representing the system matrix.
     */
    virtual void InitializeSolutionStep(const typename LinearOperatorType::UniquePointer& rpLinearOperator)
    {
    }

    /**
     * @brief Applies the preconditioner to the given vector.
     * This represents the operation y = M^{-1} x, being M^{-1} an approximation of the linear system matrix inverse.
     * @param rX The vector to be preconditioned.
     * @param rY The vector to store the result.
     */
    virtual void Apply(
        const VectorType& rX,
        VectorType& rY)
    {
        rY = rX;
    }

    /**
     * @brief Applies the transpose of the preconditioner to the given vector.
     * This represents the operation y = M^{-T} x, being M^{-T} an approximation of the linear system matrix transpose inverse.
     * @param rX The vector to be preconditioned.
     * @param rY The vector to store the result.
     */
    virtual void ApplyTranspose(
        const VectorType& rX,
        VectorType& rY)
    {
        rY = rX;
    }

    /**
     * @brief Finalizes the preconditioner for the current solution step.
     * This method is designed to be called at the end of the step.
     * @param rpLinearOperator Unique pointer to the linear operator representing the system matrix.
     */
    virtual void FinalizeSolutionStep(const typename LinearOperatorType::UniquePointer& rpLinearOperator)
    {
    }

    /**
     * @brief Clears the preconditioner.
     * This method is designed to clean up all internal data in the preconditioner.
     * After a clear a new Initialize is needed as well as to eventually set the additional data (@see SetAdditionalData)
     */
    virtual void Clear()
    {
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Set the Additional Data
     * @details Some solvers may require a minimum degree of knowledge of the structure of the matrix.
     * For example, when solving a mixed u-p problem, it is important to identify the row associated with v and p.
     * Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers,
     * which require knowledge of the spatial position of the nodes associated with a given degree of freedom (DOF).
     * This function provides the opportunity to provide such data if needed.
     * @param rModelPart The model part from which the linear system is built
     * @param rDofSet The dofs array of the linear system
     */
    void SetAdditionalData(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet)
    {
        KRATOS_WARNING_IF("Preconditioner", HasAdditionalData()) << "Additional data is already set. Overwriting it" << std::endl;
        mpDofSet = &rDofSet;
        mpModelPart = &rModelPart;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Checks if additional physical data is needed by the preconditioner.
     * @details Some preconditioners may require a minimum degree of knowledge of the structure of the matrix.
     * For instance, when solving a mixed u-p problem, it is important to identify the row associated with v and p.
     * Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers,
     * which require knowledge of the spatial position of the nodes associated with a given degree of freedom (DOF).
     * @return True if additional physical data is needed, false otherwise.
     */
    virtual bool RequiresAdditionalData() const
    {
        return false;
    }

    /**
     * @brief Check if the preconditioner has additional data (model part and dofs)
     * @return true if the preconditioner has additional data
     * @return false if the preconditioner does not have additional data
     */
    bool HasAdditionalData() const
    {
        return mpModelPart != nullptr && mpDofSet != nullptr;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Preconditioner";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Preconditioner";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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

    ModelPart* mpModelPart = nullptr; // Model part of the preconditioner

    typename ModelPart::DofsArrayType* mpDofSet = nullptr; // Dofs of the preconditioner

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

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
}; // Class Preconditioner

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<class TLinearAlgebra>
inline std::istream& operator >> (
    std::istream& IStream,
    Preconditioner<TLinearAlgebra>& rThis)
{
    return IStream;
}

/// output stream function
template<class TLinearAlgebra>
inline std::ostream& operator << (
    std::ostream& OStream,
    const Preconditioner<TLinearAlgebra>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}

///@}

}  // namespace Kratos::Future.
