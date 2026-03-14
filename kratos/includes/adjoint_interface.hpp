//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// --- Kratos Core Includes ---
#include "includes/kratos_export_api.h" // KRATOS_API
#include "containers/variable_data.h" // VariableData
#include "includes/dof.h" // Dof

// --- STL Includes ---
#include <vector> // std::vector
#include <span> // std::span


namespace Kratos {


/// @brief Basic interface for adjoint problems.
class KRATOS_API(KRATOS_CORE) IAdjoint {
public:
    using VARIABLE = VariableData;

    using Scalar = double;

    enum class ResidualTerm {
        Load,
        Mass,
        Damping,
        Stiffness
    }; // enum class Term

    virtual ~IAdjoint() = default;

    /// @name Variable Query
    /// @{

    /// @brief Collect all scalar state variables.
    /// @param[out] rOutput Range of scalar variables this function will populate.
    virtual void GetStateVariables(std::vector<VARIABLE>& rOutput) const;

    /// @brief Collect all scalar variables the entity depends on.
    /// @param[out] rOutput Range of scalar variables this function will populate.
    virtual void GetInfluencingVariables(std::vector<VARIABLE>& rOutput) const;

    /// @}

    static constexpr std::string_view TermName(ResidualTerm Term) noexcept {
        switch (Term) {
            case ResidualTerm::Load: return "Load";             break;
            case ResidualTerm::Mass: return "Mass";             break;
            case ResidualTerm::Damping: return "Damping";       break;
            case ResidualTerm::Stiffness: return "Stiffness";   break;
        }
    }

}; // class IAdjoint


/// @brief Interface for computing the partial derivatives of entity residuals.
/// @details This class is responsible for computing derivatives of the entity's residual.
///          @f[
///             R = F - (M \ddot u + D \dot u + K u)
///          @f]
///          where
///          - @f$R@f$ is the entity's residual,
///          - @f$F@f$ is the external force vector,
///          - @f$u@f$ is the entity's state vector,
///          - @f$M@f$ is the entity's mass matrix,
///          - @f$D@f$ is the entity's damping matrix,
///          - @f$K@f$ is the entity's stiffness matrix.
///          The derivatives are exposed as the derivatives of individual components
///          @f$F@f$, @f$M@f$, @f$D@f$ and @f$K@f$ because different @ref Scheme "time integration schemes"
///          may disregard specific components of the residual.
class KRATOS_API(KRATOS_CORE) IAdjointElement : public IAdjoint {
public:
    /// @name Variable Query Interface
    /// @{

    /// @copydoc IAdjoint::GetInfluencingVariables
    void GetInfluencingVariables(std::vector<IAdjoint::VARIABLE>& rOutput) const final override;

    /// @brief Collect all variables influencing a specific term of the residual.
    /// @tparam Term Term of the residual to be queried.
    /// @param[out] rOutput Array of variables influencing the queried term.
    template <IAdjoint::ResidualTerm Term>
    void GetInfluencingVariables(std::vector<IAdjoint::VARIABLE>& rOutput) const {
        if constexpr (Term == IAdjoint::ResidualTerm::Mass)
            return this->GetMassInfluencingVariables(rOutput);
        else if constexpr (Term == IAdjoint::ResidualTerm::Damping)
            return this->GetDampingInfluencingVariables(rOutput);
        else if constexpr (Term == IAdjoint::ResidualTerm::Stiffness)
            return this->GetStiffnessInfluencingVariables(rOutput);
        else if constexpr (Term == IAdjoint::ResidualTerm::Load)
            return this->GetLoadInfluencingVariables(rOutput);
        else static_assert(Term == IAdjoint::ResidualTerm::Mass, "invalid term");
    }

    /// @}
    /// @name DoF Query
    /// @{

    /// @brief Collect all degrees-of-freedom the entity is defined over.
    /// @param rOutput Array of degrees-of-freedom the entity is defined over.
    virtual void GetDofs(std::vector<const Dof<IAdjoint::Scalar>*>& rOutput) const;

    /// @{
    /// @name Derivative Query Interface
    /// @{

    /// @brief Compute @f$\frac{\partial T}{\partial \xi} \xi@f$ where @f$T@f$ is either the mass, damping, stiffness or load term of the residual.
    /// @tparam Term Term of the residual to compute the derivatives of.
    /// @brief Compute @f$\frac{\partial T}{\partial \xi} \xi@f$.
    /// @details Compute the dot product of the term's derivative
    ///          with respect to a variable, and the components of that variable.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial T}{\partial \xi} \xi         \\
    ///                 \frac{\partial T}{\partial \eta} \eta       \\
    ///                 \vdots                                      \\
    ///                 \frac{\partial T}{\partial \zeta} \zeta
    ///             \end{bmatrix}
    ///          @f]
    ///          where
    ///          - @f$T@f$ is a term of the residual @f$R@f$,
    ///          - @f$\xi@f$ is the first variable the term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the term is differentiated with respect to.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    /// @param Variables
    /// @param rOutput
    template <IAdjoint::ResidualTerm Term>
    void ComputeDerivative(
        std::span<const IAdjoint::VARIABLE> Variables,
        Matrix& rOutput) const {
            if constexpr (Term == IAdjoint::ResidualTerm::Mass)
                this->ComputeMassDerivative(Variables, rOutput);
            else if constexpr (Term == IAdjoint::ResidualTerm::Damping)
                this->ComputeDampingDerivative(Variables, rOutput);
            else if constexpr (Term == IAdjoint::ResidualTerm::Stiffness)
                this->ComputeStiffnessDerivative(Variables, rOutput);
            else if constexpr (Term == IAdjoint::ResidualTerm::Load) {
                this->ComputeLoadDerivative(Variables, rOutput);
            } else static_assert(Term == IAdjoint::ResidualTerm::Load, "invalid term");
    }

    /// @}

protected:
    /// @name Variable Query Implementation
    /// @{

    /// @brief Collect all variables influencing the residual's mass term @f$M@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    virtual void GetMassInfluencingVariables(std::vector<IAdjoint::VARIABLE>& rOutput) const;

    /// @brief Collect all variables influencing the residual's damping term @f$D@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    virtual void GetDampingInfluencingVariables(std::vector<IAdjoint::VARIABLE>& rOutput) const;

    /// @brief Collect all variables influencing the residual's stiffness term @f$K@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    virtual void GetStiffnessInfluencingVariables(std::vector<IAdjoint::VARIABLE>& rOutput) const;

    /// @brief Collect all variables influencing the residual's load term @f$D@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    virtual void GetLoadInfluencingVariables(std::vector<IAdjoint::VARIABLE>& rOutput) const;

    /// @}
    /// @name Derivative Query Implementation
    /// @{

    /// @brief Compute @f$\frac{\partial K}{\partial \xi} \xi@f$.
    /// @details Compute the dot product of the stiffness term's derivative
    ///          with respect to a variable, and the components of that variable.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial K}{\partial \xi} \xi         \\
    ///                 \frac{\partial K}{\partial \eta} \eta       \\
    ///                 \vdots                                      \\
    ///                 \frac{\partial K}{\partial \zeta} \zeta
    ///             \end{bmatrix}
    ///          @f]
    ///          where
    ///          - @f$K = \frac{\partial R}{\partial u}@f$ is the stiffness term,
    ///          - @f$R@f$ is the residual,
    ///          - @f$u@f$ is the vector of state variables,
    ///          - @f$\xi@f$ is the first variable the stiffness term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the stiffness term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the stiffness term is differentiated with respect to.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    virtual void ComputeStiffnessDerivative(
        std::span<const IAdjoint::VARIABLE> Variables,
        Matrix& rOutput) const;

    /// @brief Compute @f$\frac{\partial D}{\partial \xi} \xi@f$.
    /// @details Compute the dot product of the damping term's derivative
    ///          with respect to a variable, and the components of that variable.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial D}{\partial \xi} \xi         \\
    ///                 \frac{\partial D}{\partial \eta} \eta       \\
    ///                 \vdots                                      \\
    ///                 \frac{\partial D}{\partial \zeta} \zeta
    ///             \end{bmatrix}
    ///          @f]
    ///          where
    ///          - @f$D = \frac{\partial R}{\partial \dot u}@f$ is the damping term,
    ///          - @f$R@f$ is the residual,
    ///          - @f$u@f$ is the vector of state variables,
    ///          - @f$\xi@f$ is the first variable the damping term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the damping term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the damping term is differentiated with respect to.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    virtual void ComputeDampingDerivative(
        std::span<const IAdjoint::VARIABLE> Variables,
        Matrix& rOutput) const;

    /// @brief Compute @f$\frac{\partial M}{\partial \xi} \xi@f$.
    /// @details Compute the dot product of the mass term's derivative
    ///          with respect to a variable, and the components of that variable.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial M}{\partial \xi} \xi         \\
    ///                 \frac{\partial M}{\partial \eta} \eta       \\
    ///                 \vdots                                      \\
    ///                 \frac{\partial M}{\partial \zeta} \zeta
    ///             \end{bmatrix}
    ///          @f]
    ///          where
    ///          - @f$M = \frac{\partial R}{\partial \ddot u}@f$ is the mass term,
    ///          - @f$R@f$ is the residual,
    ///          - @f$u@f$ is the vector of state variables,
    ///          - @f$\xi@f$ is the first variable the mass term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the mass term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the mass term is differentiated with respect to.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    virtual void ComputeMassDerivative(
        std::span<const IAdjoint::VARIABLE> Variables,
        Matrix& rOutput) const;

    /// @brief Compute @f$\frac{\partial F}{\partial \xi} \xi@f$.
    /// @details Compute the load term's derivative with respect to a variable.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial F}{\partial \xi}         \\
    ///                 \frac{\partial F}{\partial \eta}        \\
    ///                 \vdots                                  \\
    ///                 \frac{\partial F}{\partial \zeta}
    ///             \end{bmatrix}
    ///          @f]
    ///          where
    ///          - @f$F@f$ is the residual's load term,
    ///          - @f$\xi@f$ is the first variable the load term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the load term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the load term is differentiated with respect to.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[out] rOutput Output matrix containing the derivatives.
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    virtual void ComputeLoadDerivative(
        std::span<const IAdjoint::VARIABLE> Variables,
        Matrix& rOutput) const;

    /// @}

}; // class IAdjointElement


} // namespace Kratos
