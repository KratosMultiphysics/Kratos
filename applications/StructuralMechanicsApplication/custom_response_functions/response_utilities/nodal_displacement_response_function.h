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

#if !defined(KRATOS_NODAL_DISPLACEMENT_RESPONSE_H_INCLUDED )
#define  KRATOS_NODAL_DISPLACEMENT_RESPONSE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application_variables.h"

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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) NodalDisplacementResponseFunction {
public:

    ///@name Type Definitions
    ///@{

    typedef Element::DofsVectorType DofsVectorType;
    typedef Node<3>::Pointer PointTypePointer;
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition oNodalDisplacementResesponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(NodalDisplacementResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    NodalDisplacementResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
    virtual ~NodalDisplacementResponseFunction() {
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
        return "NodalDisplacementResponseFunction";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {
      rOStream << "NodalDisplacementResponseFunction";
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
    std::string mTracedDofLabel;
    PointTypePointer mpTracedNode;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
    //NodalDisplacementResponseFunction(NodalDisplacementResponseFunction const& rOther);

    ///@}

}; // Class NodalDisplacementResesponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_NODAL_DISPLACEMENT_RESPONSE_H_INCLUDED  defined
