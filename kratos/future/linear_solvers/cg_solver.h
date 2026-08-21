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
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "future/linear_solvers/iterative_solver.h"

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
 * @class CGSolver
 * @ingroup KratosCore
 * @brief Preconditioned Conjugate Gradient (PCG) linear solver.
 * @details This class implements the Preconditioned Conjugate Gradient (PCG) iterative method for solving linear systems.
 *          It is derived from the IterativeSolver base class.
 * @tparam TLinearAlgebra Type of the linear algebra used (e.g. CSR space, dense space).
 */
template<class TLinearAlgebra>
class CGSolver : public Future::IterativeSolver<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CGSolver
    KRATOS_CLASS_POINTER_DEFINITION(CGSolver);

    /// Base type iterative solver definition
    using BaseType = Future::IterativeSolver<TLinearAlgebra>;

    /// Dense vector type definition from linear algebra traits
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Linear operator type definition
    using LinearOperatorType = LinearOperator<TLinearAlgebra>;

    /// Preconditions pointer type definition
    using PreconditionerPointerType = typename Preconditioner<TLinearAlgebra>::Pointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with parameters
    CGSolver(Parameters Settings = Parameters(R"({})"))
        : BaseType()
    {
        // Validate and assign default parameters
        Settings.ValidateAndAssignDefaults(GetDefaultParameters());

        // Assign the validated settings
        this->AssignSettings(Settings);

        // Check that right preconditioner is not provided
        KRATOS_ERROR_IF(Settings["right_preconditioner_type"].GetString() != "identity") << "Right preconditioner is not supported for CG solver." << std::endl;

        // Create and set left preconditioner
        PreconditionerPointerType p_left_preconditioner = Kratos::make_shared<Preconditioner<TLinearAlgebra>>(); // TODO: implement preconditioner by leveraging the settings and the registry
        this->SetLeftPreconditioner(p_left_preconditioner);
    }

    /// Copy constructor.
    CGSolver(const CGSolver& Other) = delete;

    /// Destructor.
    ~CGSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    CGSolver& operator=(const CGSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters( R"({
            "solver_type" : "cg_solver",
            "right_preconditioner_type" : "identity"
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Conjugate gradient linear solver with left preconditioner " << BaseType::GetLeftPreconditioner()->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

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

    bool IterativeSolve(
        const typename LinearOperatorType::UniquePointer& rpLinearOperator,
        const VectorType& rB,
        VectorType& rX) override
    {
        // Initialize variables
        const int size = rX.size(); // Size of the system
        BaseType::mBNorm = rB.Norm(); // Note that this will be used to compute the relative residual (@see BaseType::IterationNeeded)
        BaseType::mIterationsNumber = 0;

        // r = rB - A * rX
        VectorType r(size);
        r.SetValue(0.0);
        rpLinearOperator->SpMV(rX, r);
        r *= -1.0;
        r += rB;
        BaseType::mResidualNorm = r.Norm();

        // z = M^{-1} * r
        VectorType z(size);
        z.SetValue(0.0);
        this->GetLeftPreconditioner()->Apply(r, z);

        // Initialize iteration variables
        VectorType p(z);
        VectorType q(size);
        q.SetValue(0.0);

        double rho_0 = r.Dot(z);
        double rho_1 = rho_0;

        // Early return condition for zero residual
        if (std::abs(rho_0) < std::numeric_limits<double>::epsilon()) {
            return true;
        }

        // Main loop
        while(BaseType::IterationNeeded() && (std::abs(rho_0) > std::numeric_limits<double>::epsilon())) {

            // q = A * p
            q.SetValue(0.0);
            rpLinearOperator->SpMV(p, q);

            // pq = p^{T} * A * p
            const double pq = p.Dot(q);
            if (std::abs(pq) < std::numeric_limits<double>::epsilon()) {
                break;
            }

            // alpha: step length along the conjugate direction p
            const double alpha = rho_0 / pq;

            // dx = dx + alpha * p
            rX.Add(alpha, p);

            // r = r - alpha * q
            r.Add(-alpha, q);
            BaseType::mResidualNorm = r.Norm(); // Update residual norm for convergence check (@see BaseType::IterationNeeded)

            // z = M^{-1} * r
            z.SetValue(0.0);
            this->GetLeftPreconditioner()->Apply(r, z);

            // rho_1 = r^{T} * z
            rho_1 = r.Dot(z);

            // beta: coefficient to build the next conjugate search direction
            const double beta = (rho_1 / rho_0);

            // p = z + beta * p
            p *= beta;
            p.Add(1.0, z);

            // Next iteration update
            rho_0 = rho_1;
            BaseType::mIterationsNumber++; // Update iterations counter for convergence check (@see BaseType::IterationNeeded)
        }

        return BaseType::IsConverged();
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
}; // Class CGSolver

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
    CGSolver<TLinearAlgebra>& rThis)
{
    return IStream;
}

/// output stream function
template<class TLinearAlgebra>
inline std::ostream& operator << (
    std::ostream& OStream,
    const CGSolver<TLinearAlgebra>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}

///@}

}  // namespace Kratos.



