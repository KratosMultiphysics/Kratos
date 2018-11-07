// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//                   Armin Geiser
//


#if !defined(KRATOS_REPLACE_ELEMENTS_AND_CONDITIONS_FOR_ADJOINT_PROBLEM_PROCESS_H_INCLUDED)
#define  KRATOS_REPLACE_ELEMENTS_AND_CONDITIONS_FOR_ADJOINT_PROBLEM_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @brief This class replaces all elements and conditions of a structural
 * model part with their adjoint counterparts. If this process is called on a
 * primal model part, it is transformed to an adjoint model part and vice versa.
 * The variable "IS_ADJOINT" at the model parts process info is set accordingly.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ReplaceElementsAndConditionsForAdjointProblemProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of ReplaceElementsAndConditionsForAdjointProblemProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReplaceElementsAndConditionsForAdjointProblemProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    ReplaceElementsAndConditionsForAdjointProblemProcess(ModelPart& model_part);

    /// Destructor.
    ~ReplaceElementsAndConditionsForAdjointProblemProcess() override {};

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


    /// Execute method is used to execute the ReplaceElementsAndConditionsForAdjointProblemProcess algorithms.
    void Execute()  override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:

    ModelPart& mrModelPart;

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
    void ReplaceToAdjoint();

    void ReplaceToPrimal();

    void UpdateSubModelPart(ModelPart& r_sub_model_part, ModelPart& r_root_model_part);

    bool GetAdjointElementName(const Element& rElement, std::string& rName);

    bool GetAdjointConditionName(const Condition& rCondition, std::string& rName);

    /// Assignment operator.
    ReplaceElementsAndConditionsForAdjointProblemProcess& operator=(ReplaceElementsAndConditionsForAdjointProblemProcess const& rOther) = delete;

    /// Copy constructor.
    //ReplaceElementsAndConditionsForAdjointProblemProcess(ReplaceElementsAndConditionsForAdjointProblemProcess const& rOther);


    ///@}

}; // Class ReplaceElementsAndConditionsForAdjointProblemProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ReplaceElementsAndConditionsForAdjointProblemProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ReplaceElementsAndConditionsForAdjointProblemProcess& rThis);

///@}


}  // namespace Kratos.

#endif // KRATOS_REPLACE_ELEMENTS_AND_CONDITIONS_FOR_ADJOINT_PROBLEM_PROCESS_H_INCLUDED  defined


