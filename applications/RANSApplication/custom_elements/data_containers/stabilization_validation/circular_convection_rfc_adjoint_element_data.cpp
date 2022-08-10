//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "geometries/geometry_data.h"
#include "includes/process_info.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "circular_convection_rfc_adjoint_element_data.h"

namespace Kratos
{
namespace StabilizationValidationElementData
{

void CircularConvectionRFCAdjointElementData::Check(
    const Element& rElement,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    ScalarEquation::TResidualsDerivatives::Check(rElement, rProcessInfo);

    KRATOS_CATCH("");
}

std::vector<const Variable<double>*> CircularConvectionRFCAdjointElementData::GetDofVariablesList()
{
    return {&RANS_SCALAR_1_ADJOINT_1};
}

std::vector<const Variable<double>*> CircularConvectionRFCAdjointElementData::GetPrimalSecondDerivativeVariablesList()
{
    return {&RANS_AUXILIARY_VARIABLE_1};
}

} // namespace StabilizationValidationElementData

} // namespace Kratos
