//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nuñez, based on Martin Fusseder work, https://github.com/MFusseder
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

class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) AdjointPotentialResponseFunction : public AdjointResponseFunction
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

    void InitializeSolutionStep() override;

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
