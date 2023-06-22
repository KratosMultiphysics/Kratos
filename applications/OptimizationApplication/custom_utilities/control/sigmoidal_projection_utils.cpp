//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//                   Suneth Warnakulasuriya
//

// System includes
#include <sstream>
#include <type_traits>

// Project includes
#include "expression/variable_expression_io.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "optimization_application_variables.h"

// Include base h
#include "sigmoidal_projection_utils.h"

namespace Kratos
{

///@name Kratos Classes
///@{

template<class TContainerType>
void SigmoidalProjectionUtils::OneLevelForwardProjection(
        const ContainerExpression<TContainerType>& rInputExpression,
        ContainerExpression<TContainerType>& rProjectedExpression)
{
    KRATOS_TRY

    const auto& r_expression_1 = rInputExpression.GetExpression();
    const IndexType local_size_1 = rInputExpression.GetItemComponentCount();
    const IndexType number_of_entities_1 = rInputExpression.GetContainer().size();

    const auto& r_expression_2 = rProjectedExpression.GetExpression();
    const IndexType local_size_2 = rProjectedExpression.GetItemComponentCount();
    const IndexType number_of_entities_2 = rProjectedExpression.GetContainer().size();

    KRATOS_ERROR_IF(local_size_1 != local_size_2)
        << "Data dimension mismatch in OneLevelForwardProjection calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rInputExpression << "\n"
        << "   Container 2: " << rProjectedExpression << "\n";

    KRATOS_ERROR_IF(number_of_entities_1 != number_of_entities_2)
        << "Number of entities mismatch in OneLevelForwardProjection calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rInputExpression << "\n"
        << "   Container 2: " << rProjectedExpression << "\n";

    KRATOS_ERROR_IF(&rInputExpression.GetModelPart() != &rProjectedExpression.GetModelPart())
        << "Model part mismatch in OneLevelForwardProjection calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rInputExpression << "\n"
        << "   Container 2: " << rProjectedExpression << "\n";


    KRATOS_CATCH("");
}

///@}
}