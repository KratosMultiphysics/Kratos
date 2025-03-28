//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi, Polytimi Zisimoupoulu

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

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;



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

    void Execute() override;

    void ExecuteInitializeSolutionStep() override {
        Execute();
    };

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
          "echo_level": 0,
          "analysis_model_part_name": "IgaModelPart",
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
    SizeType mEchoLevel;
    Parameters mIgaPhysicsParameters;

    ///@}
    ///@name Iga functionalities
    ///@{

    ///@}
    ///@}

    ///@}
    ///@name Utility
    ///@{

    Parameters ReadParamatersFile(
        const std::string& rDataFileName) const;

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