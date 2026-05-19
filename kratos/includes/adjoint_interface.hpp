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
    /// @brief An extension of @ref VariableData "kratos variables" with a runtime index.
    /// @details @p DynamicVariable augments variables with an index to uniquely identify variables
    ///          that may be contained by multiple objects in an entity. For example, in the context
    ///          of an element, @p DISPLACEMENT_X augmented with a local node index can refer to the
    ///          x-component of the displacement of the corresponding node. Another example would be
    ///          to uniquely identify the @ref YOUNG_MODULUS "young's modulus" of a specific layer
    ///          within a composite material.
    class DynamicVariable : public VariableData {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(DynamicVariable);

        constexpr DynamicVariable() noexcept
            : VariableData(),
              mDynamicIndex(-1)
        {}

        constexpr DynamicVariable(const VariableData& rRhs) noexcept
            : DynamicVariable(rRhs, -1)
        {}

        constexpr DynamicVariable(
            const VariableData& rRhs,
            int DynamicIndex) noexcept
            : VariableData(rRhs),
              mDynamicIndex(DynamicIndex)
        {}

        constexpr DynamicVariable(
            std::string_view Name,
            std::size_t Size,
            int DynamicIndex = -1) noexcept
            : VariableData(Name, Size),
              mDynamicIndex(DynamicIndex)
        {}

        constexpr DynamicVariable(
            std::string_view Name,
            std::size_t Size,
            const VariableData* pSourceVariable,
            char ComponentIndex,
            int DynamicIndex = -1) noexcept
            : VariableData(Name, Size, pSourceVariable, ComponentIndex),
              mDynamicIndex(DynamicIndex)
        {}

        [[nodiscard]] constexpr int GetDynamicIndex() const noexcept {
            return mDynamicIndex;
        }

        constexpr void SetDynamicIndex(int Index) noexcept {
            mDynamicIndex = Index;
        }

        [[nodiscard]] friend bool operator==(
            const DynamicVariable& rLhs,
            const DynamicVariable& rRhs) noexcept {
                return static_cast<const VariableData&>(rLhs) == static_cast<const VariableData&>(rRhs)
                    && rLhs.GetDynamicIndex() == rRhs.GetDynamicIndex();
        }

    private:
        int mDynamicIndex;
    }; // struct DynamicVariable

    using Scalar = double;

    KRATOS_CLASS_POINTER_DEFINITION(IAdjoint);

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
    /// @param[in] rProcessInfo Current state of the computing model part.
    virtual void GetStateVariables(
        std::vector<DynamicVariable>& rOutput,
        const ProcessInfo& rProcessInfo) const;

    /// @brief Collect all scalar variables the entity depends on.
    /// @param[out] rOutput Range of scalar variables this function will populate.
    /// @param[in] rProcessInfo Current state of the computing model part.
    virtual void GetInfluencingVariables(
        std::vector<DynamicVariable>& rOutput,
        const ProcessInfo& rProcessInfo) const;

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
    KRATOS_CLASS_POINTER_DEFINITION(IAdjointElement);

    /// @name Variable Query Interface
    /// @{

    /// @copydoc IAdjoint::GetInfluencingVariables
    void GetInfluencingVariables(
        std::vector<IAdjoint::DynamicVariable>& rOutput,
        const ProcessInfo& rProcessInfo) const final override;

    /// @brief Collect all variables influencing a specific term of the residual.
    /// @tparam Term Term of the residual to be queried.
    /// @param[out] rOutput Array of variables influencing the queried term.
    /// @param[in] rProcessInfo Current state of the computing model part.
    template <IAdjoint::ResidualTerm Term>
    void GetInfluencingVariables(
        std::vector<IAdjoint::DynamicVariable>& rOutput,
        const ProcessInfo& rProcessInfo) const {
            if constexpr (Term == IAdjoint::ResidualTerm::Mass)
                return this->GetMassInfluencingVariables(rOutput, rProcessInfo);
            else if constexpr (Term == IAdjoint::ResidualTerm::Damping)
                return this->GetDampingInfluencingVariables(rOutput, rProcessInfo);
            else if constexpr (Term == IAdjoint::ResidualTerm::Stiffness)
                return this->GetStiffnessInfluencingVariables(rOutput, rProcessInfo);
            else if constexpr (Term == IAdjoint::ResidualTerm::Load)
                return this->GetLoadInfluencingVariables(rOutput, rProcessInfo);
            else static_assert(Term == IAdjoint::ResidualTerm::Mass, "invalid term");
    }

    /// @}
    /// @name Derivative Query Interface
    /// @{

    /// @brief Compute the partial derivative of a residual term with respect to a set of variables.
    /// @tparam Term Term of the residual to compute the partial derivatives of.
    /// @details The residual term @f$T@f$ can be the mass (@f$M@f$), damping (@f$D@f$) or stiffness (@f$K@f$) term,
    ///          (but not the load @f$F@f$). This function the inner product of the term's partial derivative with
    ///          respect to the provided variables, and the provided vector @p rValues.
    ///          @f[
    ///              \begin{bmatrix}
    ///                  \frac{\partial T}{\partial \xi} v    \\
    ///                  \frac{\partial T}{\partial \eta} v   \\
    ///                  \vdots                               \\
    ///                  \frac{\partial T}{\partial \zeta} v
    ///              \end{bmatrix}
    ///          @f]
    ///          where
    ///          - @f$T@f$ is a term of the residual @f$R@f$,
    ///          - @f$\xi@f$ is the first variable the term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the term is differentiated with respect to, and
    ///          - @f$v@f$ is the provided vector @p rValues.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[in] rValues Vector to compute the inner product with.
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[in] iBuffer Relative reverse time step index to compute the derivative for.
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    template <IAdjoint::ResidualTerm Term>
    requires (Term == IAdjoint::ResidualTerm::Stiffness || Term == IAdjoint::ResidualTerm::Damping || Term == IAdjoint::ResidualTerm::Mass)
    void ComputeDerivative(
        Matrix& rOutput,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const Vector& rValues,
        const ProcessInfo& rProcessInfo,
        int iBuffer = 0) const {
            if constexpr (Term == IAdjoint::ResidualTerm::Mass)
                this->ComputeMassDerivative(
                    rOutput,
                    Variables,
                    rValues,
                    rProcessInfo,
                    iBuffer);
            else if constexpr (Term == IAdjoint::ResidualTerm::Damping)
                this->ComputeDampingDerivative(
                    rOutput,
                    Variables,
                    rValues,
                    rProcessInfo,
                    iBuffer);
            else if constexpr (Term == IAdjoint::ResidualTerm::Stiffness)
                this->ComputeStiffnessDerivative(
                    rOutput,
                    Variables,
                    rValues,
                    rProcessInfo,
                    iBuffer);
    }

    /// @brief Compute the partial derivative of the load term with respect to a set of variables.
    /// @tparam Term Term of the residual to compute the partial derivatives of (must be the @ref IAdjoint::ResidualTerm::Load).
    /// @details @f[
    ///              \begin{bmatrix}
    ///                  \frac{\partial F}{\partial \xi}     \\
    ///                  \frac{\partial F}{\partial \eta}    \\
    ///                  \vdots                              \\
    ///                  \frac{\partial F}{\partial \zeta}
    ///              \end{bmatrix}
    ///          @f]
    ///          where
    ///          - @f$T@f$ is a term of the residual @f$R@f$,
    ///          - @f$\xi@f$ is the first variable the term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the term is differentiated with respect to.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[in] iBuffer Relative reverse time step index to compute the derivative for.
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    template <IAdjoint::ResidualTerm Term>
    requires (Term == IAdjoint::ResidualTerm::Load)
    void ComputeDerivative(
        Matrix& rOutput,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer = 0) const {
            return this->ComputeLoadDerivative(
                rOutput,
                Variables,
                rProcessInfo,
                iBuffer);
    }

    /// @}

protected:
    /// @name Variable Query Implementation
    /// @{

    /// @brief Collect all variables influencing the residual's mass term @f$M@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    /// @param[in] rProcessInfo Current state of the computing model part.
    virtual void GetMassInfluencingVariables(
        std::vector<IAdjoint::DynamicVariable>& rOutput,
        const ProcessInfo& rProcessInfo) const;

    /// @brief Collect all variables influencing the residual's damping term @f$D@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    /// @param[in] rProcessInfo Current state of the computing model part.
    virtual void GetDampingInfluencingVariables(
        std::vector<IAdjoint::DynamicVariable>& rOutput,
        const ProcessInfo& rProcessInfo) const;

    /// @brief Collect all variables influencing the residual's stiffness term @f$K@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    /// @param[in] rProcessInfo Current state of the computing model part.
    virtual void GetStiffnessInfluencingVariables(
        std::vector<IAdjoint::DynamicVariable>& rOutput,
        const ProcessInfo& rProcessInfo) const;

    /// @brief Collect all variables influencing the residual's load term @f$D@f$.
    /// @param[out] rOutput Array of variables influencing the queried term.
    /// @param[in] rProcessInfo Current state of the computing model part.
    virtual void GetLoadInfluencingVariables(
        std::vector<IAdjoint::DynamicVariable>& rOutput,
        const ProcessInfo& rProcessInfo) const;

    /// @}
    /// @name Derivative Query Implementation
    /// @{

    /// @brief Compute @f$\frac{\partial K}{\partial \xi} v@f$.
    /// @details Compute the dot product of the stiffness term's derivative
    ///          with respect to a variable, and an input vector @f$v@f$.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial K}{\partial \xi} v   \\
    ///                 \frac{\partial K}{\partial \eta} v  \\
    ///                 \vdots                              \\
    ///                 \frac{\partial K}{\partial \zeta} v
    ///             \end{bmatrix}
    ///          @f]
    ///          where
    ///          - @f$K = \frac{\partial R}{\partial u}@f$ is the stiffness term,
    ///          - @f$R@f$ is the residual,
    ///          - @f$u@f$ is the vector of state variables,
    ///          - @f$\xi@f$ is the first variable the stiffness term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the stiffness term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the stiffness term is differentiated with respect to,
    ///          - @f$v@f$ is the input vector @p rValues.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[in] rValues Vector to compute the inner product with.
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[in] iBuffer Relative reverse time step index to compute the derivative for.
    /// @see @ref IAdjointElement::ComputeDerivative "ComputeDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    virtual void ComputeStiffnessDerivative(
        Matrix& rOutput,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const Vector& rValues,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const;

    /// @brief Compute @f$\frac{\partial D}{\partial \xi} \dot{u}@f$.
    /// @details Compute the dot product of the damping term's derivative
    ///          with respect to a variable, and an input vector @f$v@f$.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial D}{\partial \xi} v     \\
    ///                 \frac{\partial D}{\partial \eta} v    \\
    ///                 \vdots                                \\
    ///                 \frac{\partial D}{\partial \zeta} v
    ///             \end{bmatrix}
    ///          @f]
    ///          where
    ///          - @f$D = \frac{\partial R}{\partial \dot u}@f$ is the damping term,
    ///          - @f$R@f$ is the residual,
    ///          - @f$u@f$ is the vector of state variables,
    ///          - @f$\xi@f$ is the first variable the damping term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the damping term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the damping term is differentiated with respect to,
    ///          - @f$v@f$ is the input vector @p rValues.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[in] rValues Vector to compute the inner product with.
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[in] iBuffer Relative reverse time step index to compute the derivative for.
    /// @see @ref IAdjointElement::ComputeDerivative "ComputeDerivative"
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    virtual void ComputeDampingDerivative(
        Matrix& rOutput,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const Vector& rValues,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const;

    /// @brief Compute @f$\frac{\partial M}{\partial \xi} v@f$.
    /// @details Compute the dot product of the mass term's derivative
    ///          with respect to a variable, and an input vector @f$v@f$.
    ///          @f[
    ///             \begin{bmatrix}
    ///                 \frac{\partial M}{\partial \xi} v    \\
    ///                 \frac{\partial M}{\partial \eta} v   \\
    ///                 \vdots                               \\
    ///                 \frac{\partial M}{\partial \zeta} v
    ///             \end{bmatrix}
    ///          @f]
    ///          where
    ///          - @f$M = \frac{\partial R}{\partial \ddot u}@f$ is the mass term,
    ///          - @f$R@f$ is the residual,
    ///          - @f$u@f$ is the vector of state variables,
    ///          - @f$\xi@f$ is the first variable the stiffness term is differentiated with respect to,
    ///          - @f$\eta@f$ is the second variable the stiffness term is differentiated with respect to,
    ///          - ...
    ///          - @f$\zeta@f$ is the last variable the stiffness term is differentiated with respect to,
    ///          - @f$v@f$ is the input vector @p rValues.
    /// @param[out] rOutput Output matrix containing the dot products of the derivatives
    ///                     and the arrays of variables.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[in] rValues Vector to compute the inner product with.
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[in] iBuffer Relative reverse time step index to compute the derivative for.
    /// @see @ref IAdjointElement::ComputeDerivative "ComputeDerivative"
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeLoadDerivative "ComputeLoadDerivative"
    virtual void ComputeMassDerivative(
        Matrix& rOutput,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const Vector& rValues,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const;

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
    /// @param[out] rOutput Output matrix containing the derivatives.
    /// @param[in] Variables Set of variables to compute the derivative with respect to.
    /// @param[in] rProcessInfo Current state of the computing model part.
    /// @param[in] iBuffer Relative reverse time step index to compute the derivative for.
    /// @see @ref IAdjointElement::ComputeDerivative "ComputeDerivative"
    /// @see @ref IAdjointElement::ComputeStiffnessDerivative "ComputeStiffnessDerivative"
    /// @see @ref IAdjointElement::ComputeDampingDerivative "ComputeDampingDerivative"
    /// @see @ref IAdjointElement::ComputeMassDerivative "ComputeMassDerivative"
    virtual void ComputeLoadDerivative(
        Matrix& rOutput,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const;

    /// @}

}; // class IAdjointElement


} // namespace Kratos
