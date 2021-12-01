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


#ifndef KRATOS_DEACTIVATE_CONDITIONS_ON_INACTIVE_ELEMENTS_PROCESS
#define KRATOS_DEACTIVATE_CONDITIONS_ON_INACTIVE_ELEMENTS_PROCESS

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "includes/kratos_flags.h"


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
class DeactivateConditionsOnInactiveElements : public Process
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

    /// Destructor.
    ~DeactivateConditionsOnInactiveElements() override {}

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


    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override
    {
        KRATOS_TRY

        for (auto itCond = mrModelPart.ConditionsBegin(); itCond != mrModelPart.ConditionsEnd(); ++itCond) {
            const auto &VectorOfNeighbours =  itCond->GetValue(NEIGHBOUR_ELEMENTS);
            KRATOS_ERROR_IF(VectorOfNeighbours.size() == 0) << "Condition without any corresponding element, ID " << itCond->Id() << "\n"
                                                            << "Call a process to find neighbour elements before calling this function."
                                                            << std::endl;

            bool isElementActive = false;
            for (unsigned int i=0; i < VectorOfNeighbours.size(); ++i) {
                if (VectorOfNeighbours[i].IsDefined(ACTIVE)) {
                    if (VectorOfNeighbours[i].Is(ACTIVE)) isElementActive = true;
                } else {
                    isElementActive = true;
                }
            }
            if (!isElementActive) {
                itCond->Set(ACTIVE, false);
            }

        }

        KRATOS_CATCH("")
    }

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

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DeactivateConditionsOnInactiveElements";
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

    /// Assignment operator.
    DeactivateConditionsOnInactiveElements& operator=(DeactivateConditionsOnInactiveElements const& rOther);

    /// Copy constructor.
    DeactivateConditionsOnInactiveElements(DeactivateConditionsOnInactiveElements const& rOther);

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

} // namespace Kratos


#endif // KRATOS_DEACTIVATE_CONDITIONS_ON_INACTIVE_ELEMENTS_PROCESS
