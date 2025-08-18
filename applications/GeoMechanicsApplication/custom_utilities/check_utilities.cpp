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
                                         const std::vector<VariableData>& rVectorVariables)
{
    for (unsigned int i = 0; i < rGeometry.PointsNumber(); ++i) {
		for (const auto & rVariable : rVectorVariables){
           if (!rGeometry[i].SolutionStepsDataHas(rVariable))
               KRATOS_ERROR << "missing variable " << rVariable.Name() << " on node "
                            << rGeometry[i].Id() << std::endl;
		}
	}
}

} /* namespace Kratos.*/
