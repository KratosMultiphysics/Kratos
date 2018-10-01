// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

#if !defined(KRATOS_LOCAL_STRESS_RESPONSE_H_INCLUDED )
#define  KRATOS_LOCAL_STRESS_RESPONSE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application_variables.h"
#include "stress_response_definitions.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LocalStressResponseFunction {
public:

    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType;

    typedef std::size_t SizeType;
    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of LocalStressResponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(LocalStressResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    LocalStressResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
    virtual ~LocalStressResponseFunction() {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    double CalculateValue();

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
    virtual std::string Info() const {
        return "LocalStressResponseFunction";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {
      rOStream << "LocalStressResponseFunction";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Parameters mResponseSettings;
    unsigned int mIdOfLocation;
    Element::Pointer mpTracedElement;
    StressTreatment mStressTreatment;
    TracedStressType mTracedStressType;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    double CalculateMeanElementStress();

    double CalculateGaussPointStress();

    double CalculateNodeStress();

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.

    /// Copy constructor.
    //LocalStressResponseFunction(LocalStressResponseFunction const& rOther);

    ///@}

}; // Class LocalStressResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_LOCAL_STRESS_RESPONSE_H_INCLUDED  defined
