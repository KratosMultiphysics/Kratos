//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//


#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Delete particle elements and conditions with flag TO_ERASE

class MaterialPointEraseProcess
        : public Process
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodeAndElementEraseProcess
    KRATOS_CLASS_POINTER_DEFINITION(MaterialPointEraseProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MaterialPointEraseProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_TRY;

        const int initial_num_element = mr_model_part.NumberOfElements();
        mr_model_part.RemoveElements( TO_ERASE );
        const int num_removed_elements = initial_num_element - mr_model_part.NumberOfElements();

        KRATOS_WARNING_IF("MaterialPointEraseProcess", num_removed_elements > 0) << num_removed_elements << " particle elements have been erased.\n";

        const int initial_num_condition = mr_model_part.NumberOfConditions();
        mr_model_part.RemoveConditions( TO_ERASE );
        const int num_removed_condition = initial_num_condition - mr_model_part.NumberOfConditions();

        KRATOS_WARNING_IF("MaterialPointEraseProcess", num_removed_condition > 0) << num_removed_condition << " particle conditions have been erased.\n";

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "MaterialPointEraseProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MaterialPointEraseProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:

    ///@name Member Variables
    ///@{

    ModelPart& mr_model_part;

    ///@}

}; // Class NodeAndElementEraseProcess

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MaterialPointEraseProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MaterialPointEraseProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.
