//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//	                Kratos default license: kratos/license.txt
//
//  Main Authors:   Máté Kelemen
//

#ifndef KRATOS_FORCE_AND_TORQUE_UTILS
#define KRATOS_FORCE_AND_TORQUE_UTILS

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"

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

/**
 * 
 */
class KRATOS_API(KRATOS_CORE) ForceAndTorqueUtils
{
public:
    ///@name Type Definitions
    ///@{
    
    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */

    /** Destructor.
     */

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static std::array<array_1d<double,3>,2> SumForceAndTorque(const ModelPart& rModelPart, const array_1d<double,3>& rReferencePoint);

    ///@}
    ///@name Acces
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private Acces
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}
}; /* class ForceAndTorqueUtils */

///@}
} /* namespace Kratos */

#endif /* KRATOS_FORCE_AND_TORQUE_UTILS  defined */