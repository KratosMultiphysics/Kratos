// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Fusseder Martin,
//                   Armin Geiser,
//                   Daniel Baumgaertner,
//                   Suneth Warnakulasuriya
//

#pragma once

// System include
#include <iostream>
#include <string>

// External include

// Project include
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class EigenfrequencyResponseFunctionUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EigenfrequencyResponseFunctionUtility
    KRATOS_CLASS_POINTER_DEFINITION(EigenfrequencyResponseFunctionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    EigenfrequencyResponseFunctionUtility(
        ModelPart& rModelPart,
        Parameters ResponseSettings);

    EigenfrequencyResponseFunctionUtility(
        ModelPart& rModelPart,
        const double PerturbationSize,
        const std::vector<IndexType>& TracedEigenFrequencyIds,
        const std::vector<double>& WeightingFactors);

    /// Destructor.
    virtual ~EigenfrequencyResponseFunctionUtility() = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() {}

    double CalculateValue();

    void CalculateGradient();

    void CalculateEigenFrequencyMaterialVariableSensitivity(
        const Variable<double>& rMaterialVariable,
        const Variable<double>& rMaterialSensitivityVariable);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "EigenfrequencyResponseFunctionUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "EigenfrequencyResponseFunctionUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void CheckSettingsForGradientAnalysis(Parameters rResponseSettings);

    void DetermineTracedEigenfrequencies(Parameters rResponseSettings);

    bool AreSeveralEigenfrequenciesTraced();

    void CalculateLinearWeights(Parameters rResponseSettings);

    void UseDefaultWeight();

    void CheckIfAllNecessaryEigenvaluesAreComputed();

    void PerformSemiAnalyticSensitivityAnalysis();

    double GetEigenvalue(const IndexType EigenFrequencyId);

    void DetermineEigenvectorOfElement(
        const ModelPart::ElementType& rElement,
        const IndexType EigenFrequencyId,
        Vector& rEigenvectorOfElement,
        const ProcessInfo& rCurrentProcessInfo);

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart &mrModelPart;
    double mDelta;
    std::vector<IndexType> mTracedEigenfrequencyIds;
    std::vector<double> mWeightingFactors;

    ///@}
    ///@name Private static operations
    ///@{

    static void CheckIfEigenValueIsComputed(
        const ModelPart& rModelPart,
        const IndexType EigenFrequencyId);

    ///@}
}; // Class EigenfrequencyResponseFunctionUtility

///@}

} // namespace Kratos.

