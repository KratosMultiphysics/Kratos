//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi, 
//                  Polytimi Zisimoupoulu

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "utilities/function_parser_utility.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/* @class AssignIgaExternalConditionsProcess
 * @ingroup IgaApplication
 **/
class KRATOS_API(IGA_APPLICATION) AssignIgaExternalConditionsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;
    using SizeType = std::size_t;
    typedef Node                                             NodeType;
    typedef Geometry<NodeType>                                  GeometryType;
    typedef GeometryType::Pointer                               GeometryPointerType;

    /// Pointer definition of AssignIgaExternalConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignIgaExternalConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    AssignIgaExternalConditionsProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~AssignIgaExternalConditionsProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
          "echo_level": 0,
          "model_part_name": "IgaModelPart",
          "element_condition_list": []
        })" );

        return default_parameters;
    }

    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "AssignIgaExternalConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignIgaExternalConditionsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Iga functionalities
    ///@{

    Model* mpModel = nullptr;
    Parameters mParameters;
    Parameters mElementConditionList;

    ///@}
    ///@name Iga functionalities
    ///@{

    ///@}
    ///@}

    ///@}
    ///@name Utility
    ///@{

    /**
     * @brief Set the External Condition To Elements And Conditions object
     * 
     * @param sub_model_part 
     * @param variable_name 
     * @param variable_component_name 
     * @param values_string 
     * @param r_process_info 
     * @param time 
     */
    void SetExternalConditionToElementsAndConditions(
        ModelPart& sub_model_part,
        const std::string& variable_name,
        const std::vector<std::string>& variable_component_name,
        const std::vector<std::string>& values_string,
        const ProcessInfo& r_process_info,
        const double time);

    /**
     * @brief Set the Variable Value To Element object
     * 
     * @param rVariableName 
     * @param value 
     * @param p_element 
     */
    void SetVariableValueToElement(
        const std::string& rVariableName,
        const double& value,
        Element::Pointer p_element);

    
    /**
     * @brief Set the variable value to condition
     * 
     * @param rVariableName 
     * @param value
     * @param p_condition 
     */
    void SetVariableValueToCondition(
        const std::string& rVariableName,
        const double& value,
        Condition::Pointer p_condition);

    ///@}
    ///@name Input and output
    ///@{

    ///@}

}; // Class AssignIgaExternalConditionsProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignIgaExternalConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignIgaExternalConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos