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
#include "includes/process_info.h" // ProcessInfo

// --- STL Includes ---
#include <vector> // std::vector
#include <span> // std::span


namespace Kratos {


/// @brief Basic interface for adjoint problems.
class KRATOS_API(KRATOS_CORE) IAdjoint {
public:
    class DynamicVariable : public VariableData {
    public:
        constexpr DynamicVariable() noexcept
            : VariableData(),
              mDynamicIndex(0ul)
        {}

        constexpr DynamicVariable(const VariableData& rRhs) noexcept
            : DynamicVariable(rRhs, 0ul)
        {}

        constexpr DynamicVariable(
            const VariableData& rRhs,
            std::size_t DynamicIndex) noexcept
            : VariableData(rRhs),
              mDynamicIndex(DynamicIndex)
        {}

        constexpr DynamicVariable(
            std::string_view Name,
            std::size_t Size,
            std::size_t DynamicIndex = 0ul) noexcept
            : VariableData(Name, Size),
              mDynamicIndex(DynamicIndex)
        {}

        constexpr DynamicVariable(
            std::string_view Name,
            std::size_t Size,
            const VariableData* pSourceVariable,
            char ComponentIndex,
            std::size_t DynamicIndex = 0ul) noexcept
            : VariableData(Name, Size, pSourceVariable, ComponentIndex),
              mDynamicIndex(DynamicIndex)
        {}

        [[nodiscard]] constexpr std::size_t GetDynamicIndex(std::size_t Index) const noexcept {
            return mDynamicIndex;
        }

        constexpr void SetDynamicIndex(std::size_t Index) noexcept {
            mDynamicIndex = Index;
        }

    private:
        std::size_t mDynamicIndex;
    }; // struct DynamicVariable

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
    virtual void GetStateVariables(std::vector<DynamicVariable>& rOutput) const;

    /// @brief Collect all scalar variables the entity depends on.
    /// @param[out] rOutput Range of scalar variables this function will populate.
    virtual void GetInfluencingVariables(std::vector<DynamicVariable>& rOutput) const;

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
    void GetInfluencingVariables(std::vector<IAdjoint::DynamicVariable>& rOutput) const final override;

    /// @brief Collect all variables influencing a specific term of the residual.
    /// @tparam Term Term of the residual to be queried.
    /// @param[out] rOutput Array of variables influencing the queried term.
    template <IAdjoint::ResidualTerm Term>
    void GetInfluencingVariables(std::vector<IAdjoint::DynamicVariable>& rOutput) const {
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

    /// @brief Compute the partial derivative of a residual term with respect to a set of variables.
    /// @tparam Term Term of the residual to compute the partial derivatives of.
    /// @details The residual term @f$T@f$ can be the load (@f$F@f$), mass (@f$M@f$), damping (@f$D@f$) or stiffness (@f$K@f$) term,
    ///          and this function computes different expressions depending on which term is specified.
    ///          - If @f$T = F@f$, compute the partial derivatives of the load @f$F@f$ with respect to the input variables.
    ///               @f[
    ///                 \begin{bmatrix}
    ///                     \frac{\partial F}{\partial \xi}     \\
    ///                     \frac{\partial F}{\partial \eta}    \\
    ///                     \vdots                              \\
    ///                     \frac{\partial F}{\partial \zeta}
    ///                 \end{bmatrix}
    ///               @f]
    ///          - If @f$T \in \{M,D,K\}@f$, compute the inner product of the term's partial derivative with respect
    ///            to the provided variables, and the vector of state variables' time derivatives corresponding to the term.
    ///             - @f$T = M@f$
    ///               @f[
    ///                 \begin{bmatrix}
    ///                     \frac{\partial M}{\partial \xi} \ddot{u}    \\
    ///                     \frac{\partial M}{\partial \eta} \ddot{u}   \\
    ///                     \vdots                                      \\
    ///                     \frac{\partial M}{\partial \zeta} \ddot{u}
    ///                 \end{bmatrix}
    ///               @f]
    ///             - @f$T = D@f$
    ///               @f[
    ///                 \begin{bmatrix}
    ///                     \frac{\partial D}{\partial \xi} \dot{u}     \\
    ///                     \frac{\partial D}{\partial \eta} \dot{u}    \\
    ///                     \vdots                                      \\
    ///                     \frac{\partial D}{\partial \zeta} \dot{u}
    ///                 \end{bmatrix}
    ///               @f]
    ///             - @f$T = K@f$
    ///               @f[
    ///                 \begin{bmatrix}
    ///                     \frac{\partial K}{\partial \xi} u   \\
    ///                     \frac{\partial K}{\partial \eta} u  \\
    ///                     \vdots                              \\
    ///                     \frac{\partial K}{\partial \zeta} u
    ///                 \end{bmatrix}
    ///               @f]
    ///          where
    ///          - @f$T@f$ is a term of the residual @f$R@f$,
    ///          - @f$\xi@f$ is the first variable the term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the term is differentiated with respect to.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    template <IAdjoint::ResidualTerm Term>
    void ComputeDerivative(
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        Matrix& rOutput) const {
            if constexpr (Term == IAdjoint::ResidualTerm::Mass)
                this->ComputeMassDerivative(Variables, rProcessInfo, rOutput);
            else if constexpr (Term == IAdjoint::ResidualTerm::Damping)
                this->ComputeDampingDerivative(Variables, rProcessInfo, rOutput);
            else if constexpr (Term == IAdjoint::ResidualTerm::Stiffness)
                this->ComputeStiffnessDerivative(Variables, rProcessInfo, rOutput);
            else if constexpr (Term == IAdjoint::ResidualTerm::Load) {
                this->ComputeLoadDerivative(Variables, rProcessInfo, rOutput);
            } else static_assert(Term == IAdjoint::ResidualTerm::Load, "invalid term");
    }

    /// @}

protected:
    /// @name Variable Query Implementation
    /// @{

    /// @brief Collect all variables influencing the residual's mass term @f$M@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    virtual void GetMassInfluencingVariables(std::vector<IAdjoint::DynamicVariable>& rOutput) const;

    /// @brief Collect all variables influencing the residual's damping term @f$D@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    virtual void GetDampingInfluencingVariables(std::vector<IAdjoint::DynamicVariable>& rOutput) const;

    /// @brief Collect all variables influencing the residual's stiffness term @f$K@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    virtual void GetStiffnessInfluencingVariables(std::vector<IAdjoint::DynamicVariable>& rOutput) const;

    /// @brief Collect all variables influencing the residual's load term @f$D@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    virtual void GetLoadInfluencingVariables(std::vector<IAdjoint::DynamicVariable>& rOutput) const;

    /// @}
    /// @name Derivative Query Implementation
    /// @{

    /// @brief Compute @f$\frac{\partial K}{\partial \xi} u@f$.
    /// @details Compute the dot product of the stiffness term's derivative
    ///          with respect to a variable, and the components of that variable.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial K}{\partial \xi} u   \\
    ///                 \frac{\partial K}{\partial \eta} u  \\
    ///                 \vdots                              \\
    ///                 \frac{\partial K}{\partial \zeta} u
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
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @see @ref IAdjointElement::ComputeDerivative "ComputeDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    virtual void ComputeStiffnessDerivative(
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        Matrix& rOutput) const;

    /// @brief Compute @f$\frac{\partial D}{\partial \xi} \dot{u}@f$.
    /// @details Compute the dot product of the damping term's derivative
    ///          with respect to a variable, and the components of that variable.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial D}{\partial \xi} \dot{u}     \\
    ///                 \frac{\partial D}{\partial \eta} \dot{u}    \\
    ///                 \vdots                                      \\
    ///                 \frac{\partial D}{\partial \zeta} \dot{u}
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
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @see @ref IAdjointElement::ComputeDerivative "ComputeDerivative"
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    virtual void ComputeDampingDerivative(
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        Matrix& rOutput) const;

    /// @brief Compute @f$\frac{\partial M}{\partial \xi} \ddot{u}@f$.
    /// @details Compute the dot product of the mass term's derivative
    ///          with respect to a variable, and the components of that variable.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial M}{\partial \xi} \ddot{u}    \\
    ///                 \frac{\partial M}{\partial \eta} \ddot{u}   \\
    ///                 \vdots                                      \\
    ///                 \frac{\partial M}{\partial \zeta} \ddot{u}
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
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @see @ref IAdjointElement::ComputeDerivative "ComputeDerivative"
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    virtual void ComputeMassDerivative(
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        Matrix& rOutput) const;

    /// @brief Compute @f$\frac{\partial F}{\partial \xi}@f$.
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
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[out] rOutput Output matrix containing the derivatives.
    /// @see @ref IAdjointElement::ComputeDerivative "ComputeDerivative"
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    virtual void ComputeLoadDerivative(
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        Matrix& rOutput) const;

    /// @}

}; // class IAdjointElement


} // namespace Kratos
