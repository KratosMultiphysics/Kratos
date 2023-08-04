// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{
///@name Type Definitions
///@{


///@}
///@name Kratos Classes
///@{

/**
 * @class DeactivateConditionsOnInactiveElements
 * @ingroup Kratos.GeoMechanicsApplication
 * @brief Deactivate a condition if all elements attached to that are inactive
 * @note Allways neighbour elements should be found before calling this process, otherwise an error will be given.
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) DeactivateConditionsOnInactiveElements
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(DeactivateConditionsOnInactiveElements);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for DeactivateConditionsOnInactiveElements Process
    /**
     * @param rModelPart The model part to check.
     */
    DeactivateConditionsOnInactiveElements( ModelPart& rModelPart ):  Process(),
            mrModelPart(rModelPart)
    {
    }

    DeactivateConditionsOnInactiveElements(const DeactivateConditionsOnInactiveElements&) = delete;
    DeactivateConditionsOnInactiveElements& operator=(DeactivateConditionsOnInactiveElements&) = delete;
    ~DeactivateConditionsOnInactiveElements() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /// deactive conditions which are imposed on inactive elments
    void Execute() override;

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
        return "DeactivateConditionsOnInactiveElements";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  DeactivateConditionsOnInactiveElements& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DeactivateConditionsOnInactiveElements& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}