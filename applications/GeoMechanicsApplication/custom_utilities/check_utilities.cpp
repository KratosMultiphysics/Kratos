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

void CheckUtilities::CheckProperty(size_t                          Id,
                                   const Properties&               rProperties,
                                   const Kratos::Variable<double>& rVariable,
                                   std::optional<double>           MaxValue)
{
    KRATOS_ERROR_IF_NOT(rProperties.Has(rVariable))
        << rVariable.Name() << " does not exist in the material properties (Id = " << rProperties.Id()
        << ") at element " << Id << std::endl;
    constexpr auto min_value = 0.0;
    if (MaxValue.has_value()) {
        KRATOS_ERROR_IF(rProperties[rVariable] < min_value || rProperties[rVariable] > MaxValue.value())
            << rVariable.Name() << " of material Id = " << rProperties.Id() << " at element " << Id
            << " has an invalid value " << rProperties[rVariable] << " which is outside of the range [ "
            << min_value << ", " << MaxValue.value() << "]" << std::endl;
    } else {
        KRATOS_ERROR_IF(rProperties[rVariable] < min_value)
            << rVariable.Name() << " of material Id = " << rProperties.Id() << " at element " << Id
            << " has an invalid value " << rProperties[rVariable]
            << " which is below the minimum allowed value of " << min_value << std::endl;
    }
}

void CheckUtilities::CheckProperty(size_t                               Id,
                                   const Properties&                    rProperties,
                                   const Kratos::Variable<std::string>& rVariable,
                                   const std::string&                   rName)
{
    KRATOS_ERROR_IF_NOT(rProperties.Has(rVariable))
        << rVariable.Name() << " does not exist in the pressure element's properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rProperties[rVariable] == rName)
        << rVariable.Name() << " has a value of (" << rProperties[rVariable] << ") instead of ("
        << rName << ") at element " << Id << std::endl;
}

void CheckUtilities::CheckForNonZeroZCoordinateIn2D(size_t Dimension, const Geometry<Node>& rGeometry)
{
    if (Dimension == 2) {
        auto pos = std::find_if(rGeometry.begin(), rGeometry.end(),
                                [](const auto& node) { return node.Z() != 0.0; });
        KRATOS_ERROR_IF_NOT(pos == rGeometry.end())
            << " Node with non-zero Z coordinate found. Id: " << pos->Id() << std::endl;
    }
}

} /* namespace Kratos.*/
