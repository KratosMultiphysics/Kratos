// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

// Project includes
#include "check_utilities.h"
#include "includes/exception.h"

namespace Kratos
{

void CheckUtilities::CheckDomainSize(double DomainSize, std::size_t Id, const std::optional<std::string>& PrintName)
{
    constexpr auto min_domain_size = 1.0e-15;
    KRATOS_ERROR_IF(DomainSize < min_domain_size)
        << PrintName.value_or("DomainSize") << " (" << DomainSize << ") is smaller than "
        << min_domain_size << " for element " << Id << std::endl;
}

void CheckUtilities::CheckHasSolutionStepsDataFor(const Geometry<Node>& rGeometry, const VariableData& rVariable)
{
    for (const auto& node : rGeometry) {
        KRATOS_ERROR_IF_NOT(node.SolutionStepsDataHas(rVariable))
            << "Missing variable " << rVariable.Name() << " on node " << node.Id() << std::endl;
    }
}

void CheckUtilities::CheckHasDofsFor(const Geometry<Node>& rGeometry, const Variable<double>& rVariable)
{
    for (const auto& node : rGeometry) {
        KRATOS_ERROR_IF_NOT(node.HasDofFor(rVariable))
            << "Missing degree of freedom for " << rVariable.Name() << " on node " << node.Id()
            << std::endl;
    }
}
} /* namespace Kratos.*/
