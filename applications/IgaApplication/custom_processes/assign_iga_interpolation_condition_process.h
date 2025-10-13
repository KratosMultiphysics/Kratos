//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:    Andrea Gorgi
//  Co-author:      Codex GPT

#pragma once

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"
#include "iga_application_variables.h"

namespace Kratos
{

/**
 * @class AssignIgaInterpolationConditionProcess
 * @ingroup IgaApplication
 * @brief Assigns condition values by interpolating nodal data coming from a reference model part.
 */
class KRATOS_API(IGA_APPLICATION) AssignIgaInterpolationConditionProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{
    using IndexType = std::size_t;
    using SizeType = std::size_t;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    KRATOS_CLASS_POINTER_DEFINITION(AssignIgaInterpolationConditionProcess);
    ///@}

    ///@name Life Cycle
    ///@{
    AssignIgaInterpolationConditionProcess(Model& rModel, Parameters ThisParameters);

    ~AssignIgaInterpolationConditionProcess() override = default;
    ///@}

    ///@name Operations
    ///@{
    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    const Parameters GetDefaultParameters() const override;

    std::string Info() const override
    {
        return "AssignIgaInterpolationConditionProcess";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignIgaInterpolationConditionProcess";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }
    ///@}

private:
    ///@name Internal Operations
    ///@{
    void AssignInterpolatedValuesToConditions(ModelPart& rTargetModelPart, const ModelPart& rReferenceModelPart);

    using VariableListType = std::vector<const Variable<double>*>;
    const VariableListType& GetSupportedVariables() const;
    ///@}

    ///@name Members
    ///@{
    Model* mpModel = nullptr;
    Parameters mParameters;
    ///@}
};

} // namespace Kratos
