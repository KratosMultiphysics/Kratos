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

#ifndef ADJOINT_POTENTIAL_RESPONSE_FUNCTION_H
#define ADJOINT_POTENTIAL_RESPONSE_FUNCTION_H

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "response_functions/adjoint_response_function.h"


namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

class  AdjointPotentialResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointPotentialResponseFunction);

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Pointer Definitions

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    AdjointPotentialResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
     ~AdjointPotentialResponseFunction() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    void CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo) override;

    double CalculateValue(ModelPart& rModelPart) override;

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    unsigned int mGradientMode;
    double mDelta;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};


} /* namespace Kratos.*/

#endif /*  defined */
