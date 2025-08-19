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

void CheckUtilities::CheckNodalVariables(const Geometry<Node>&            rGeometry,
                                         const std::vector<VariableData>& rVariables)
{
    for (std::size_t i = 0; i < rGeometry.PointsNumber(); ++i) {
		for (const auto & r_variable : rVariables){
           if (!rGeometry[i].SolutionStepsDataHas(r_variable))
               KRATOS_ERROR << "Missing variable " << r_variable.Name() << " on node "
                            << rGeometry[i].Id() << std::endl;
		}
	}
}

void CheckUtilities::CheckNodalDof(const Geometry<Node>&            rGeometry,
                                   const std::vector<VariableData>& rVariables)
{
    for (std::size_t i = 0; i < rGeometry.PointsNumber(); ++i) {
        for (const auto& r_variable : rVariables) {
            if (!rGeometry[i].HasDofFor(r_variable))
                KRATOS_ERROR << "Missing the DoF for the variable " << r_variable.Name() << " on node "
                             << rGeometry[i].Id() << std::endl;
        }
    }
}

} /* namespace Kratos.*/
