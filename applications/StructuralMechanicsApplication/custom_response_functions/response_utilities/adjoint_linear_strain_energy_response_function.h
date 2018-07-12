// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

#ifndef ADJOINT_LINEAR_STRAIN_ENERGY_RESPONSE_FUNCTION_H
#define ADJOINT_LINEAR_STRAIN_ENERGY_RESPONSE_FUNCTION_H

// System includes

// External includes

// Project includes
#include "adjoint_structural_response_function.h"

namespace Kratos
{

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

/** \brief AdjointLinearStrainEnergyResponseFunction
*
* This is a response function which traces the linear strain energy as response.
* It is designed to be used in adjoint sensitivity analysis.
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AdjointLinearStrainEnergyResponseFunction : public AdjointStructuralResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    typedef AdjointStructuralResponseFunction BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of AdjointLinearStrainEnergyResponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(AdjointLinearStrainEnergyResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdjointLinearStrainEnergyResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
    ~AdjointLinearStrainEnergyResponseFunction();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    double CalculateValue(ModelPart& rModelPart) override;

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

    void CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override;


    void CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override;


    void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override;

    void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo) override;

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CheckForBodyForces(Element& rAdjointElem);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

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
