//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_SENSITIVITY_BUILDER_SCHEME_H_INCLUDED)
#define KRATOS_SENSITIVITY_BUILDER_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "response_functions/adjoint_response_function.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief Scheme used in the Sensitivity Builder
 *
 * This class provides the interface which is used by sensitivity builder.
 * Schemes derrived from this scheme should be able to calculate sensitivities
 * given that, adjoint values are calculated properly.
 *
 * @see SensitivityBuilder
 */
class SensitivityBuilderScheme
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensitivityBuilderScheme);

    using NodeType = ModelPart::NodeType;

    using ConditionType = ModelPart::ConditionType;

    using ElementType = ModelPart::ElementType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SensitivityBuilderScheme()
    {}

    /// Destructor.
    virtual ~SensitivityBuilderScheme() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual int Check(
        const ModelPart& rModelPart,
        ModelPart& rSensitivityModelPart) const
    {
        return 0;
    }

    virtual void Initialize(
        ModelPart& rModelPart,
        ModelPart& rSensitivityModelPart,
        AdjointResponseFunction& rResponseFunction)
    {}

    virtual void InitializeSolutionStep(
        ModelPart& rModelPart,
        ModelPart& rSensitivityModelPart,
        AdjointResponseFunction& rResponseFunction)
    {}

    virtual void FinalizeSolutionStep(
        ModelPart& rModelPart,
        ModelPart& rSensitivityModelPart,
        AdjointResponseFunction& rResponseFunction)
    {}

    virtual void Finalize(
        ModelPart& rModelPart,
        ModelPart& rSensitivityModelPart,
        AdjointResponseFunction& rResponseFunction)
    {}

    virtual void Update(
        ModelPart& rModelPart,
        ModelPart& rSensitivityModelPart,
        AdjointResponseFunction& rResponseFunction)
    {}

    virtual void Clear()
    {}

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on nodal quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentElement           Element being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ElementType& rCurrentElement,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<double>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on condition quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentElement           Element being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ElementType& rCurrentElement,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ConditionType>& rGPSensitivityVector,
        const Variable<double>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on element quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentElement           Element being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ElementType& rCurrentElement,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ElementType>& rGPSensitivityVector,
        const Variable<double>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given condition
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on nodal quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentCondition         Condition being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ConditionType& rCurrentCondition,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<double>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given condition
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on condition quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentCondition         Condition being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ConditionType& rCurrentCondition,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ConditionType>& rGPSensitivityVector,
        const Variable<double>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given condition
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on element quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentCondition         Condition being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ConditionType& rCurrentCondition,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ElementType>& rGPSensitivityVector,
        const Variable<double>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on nodal quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentElement           Element being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ElementType& rCurrentElement,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on condition quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentElement           Element being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ElementType& rCurrentElement,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ConditionType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on element quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentElement           Element being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ElementType& rCurrentElement,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ElementType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given condition
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on nodal quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentCondition         Condition being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ConditionType& rCurrentCondition,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given condition
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on condition quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentCondition         Condition being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ConditionType& rCurrentCondition,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ConditionType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    /**
     * @brief Calculates sensitivity from a given condition
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on element quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentCondition         Condition being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    virtual void CalculateSensitivity(
        ConditionType& rCurrentCondition,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ElementType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "SensitivityBuilderScheme";
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
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class SensitivityBuilderScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_SENSITIVITY_BUILDER_SCHEME_H_INCLUDED defined */
