// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Kevin Braun, https://github.com/MFusseder
//

#ifndef DIRECT_SENSITIVITY_ELEMENT_DATA_VARIABLE_H
#define DIRECT_SENSITIVITY_ELEMENT_DATA_VARIABLE_H

// System includes

// External includes

// Project includes
#include "direct_sensitivity_variable.h"
#include "includes/define.h"
#include "includes/element.h"


namespace Kratos
{

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DirectSensitivityElementDataVariable : public DirectSensitivityVariable
{
public:
    ///@name Type Definitions
    ///@{

    typedef DirectSensitivityVariable BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of AdjointLinearStrainEnergyResponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(DirectSensitivityElementDataVariable);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DirectSensitivityElementDataVariable(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
    ~DirectSensitivityElementDataVariable();
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    void CalculatePseudoLoadVector(Matrix& rPseudoLoadVector, ProcessInfo& rProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    
private:

/*const int mIdTracedElement;
Element mTracedElement;    
const Variable<double> mDesignVariable;
double mDelta;*/
}; // Class AdjointLinearStrainEnergyResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // ADJOINT_LINEAR_STRAIN_ENERGY_RESPONSE_FUNCTION_H_INCLUDED
