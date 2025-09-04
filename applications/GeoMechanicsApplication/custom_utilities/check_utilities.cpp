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
//                   Richard Faasse
//

// Project includes
#include "check_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/exception.h"

#include <sstream>

namespace Kratos
{

void CheckUtilities::CheckDomainSize(double DomainSize, std::size_t Id, const std::optional<std::string>& PrintName)
{
    constexpr auto min_domain_size = 1.0e-15;
    KRATOS_ERROR_IF(DomainSize < min_domain_size)
        << PrintName.value_or("DomainSize") << " (" << DomainSize << ") is smaller than "
        << min_domain_size << " for element " << Id << std::endl;
}

void CheckUtilities::CheckHasNodalSolutionStepData(const Geometry<Node>&             rGeometry,
                                                   const Geo::ConstVariableDataRefs& rVariableRefs)
{
    for (const auto& r_variable_ref : rVariableRefs) {
        std::vector<std::size_t> missing_node_ids;
        for (const auto& node : rGeometry) {
            if (!node.SolutionStepsDataHas(r_variable_ref.get())) {
                missing_node_ids.push_back(node.Id());
            }
        }
        if (!missing_node_ids.empty())
            KRATOS_ERROR << "Missing variable " << r_variable_ref.get().Name() << " on nodes "
                         << PrintVectorContent(missing_node_ids) << std::endl;
    }
}

void CheckUtilities::CheckHasDofs(const Geometry<Node>& rGeometry, const Geo::ConstVariableDataRefs& rVariableRefs)
{
    for (const auto& r_variable_ref : rVariableRefs) {
        std::vector<std::size_t> missing_node_ids;
        for (const auto& node : rGeometry) {
            if (!node.HasDofFor(r_variable_ref.get())) {
                missing_node_ids.push_back(node.Id());
            }
        }

        if (!missing_node_ids.empty())
            KRATOS_ERROR << "Missing the DoF for the variable " << r_variable_ref.get().Name()
                         << " on nodes " << PrintVectorContent(missing_node_ids) << std::endl;
    }
}

std::string CheckUtilities::PrintVectorContent(const std::vector<size_t>& rVector)
{
    std::ostringstream oss;
    for (const auto& r_value : rVector)
        oss << r_value << " ";

    std::string output = oss.str();
    if (!output.empty()) output.pop_back();

    return output;
}

void CheckProperties::RestorePreviousBounds() const
{
    if (!mHistoryOfTypes.empty()) {
        mRangeBoundsType = mHistoryOfTypes.back();
        mHistoryOfTypes.pop_back();
    }
}


void CheckProperties::CheckPermeabilityProperties(size_t Dimension) const
{
    SetNewRangeBounds(Bounds::InclusiveLowerAndExclusiveUpper);
    Check(PERMEABILITY_XX);
    if (Dimension > 1) {
        Check(PERMEABILITY_YY);
        Check(PERMEABILITY_XY);
    }
    if (Dimension > 2) {
        Check(PERMEABILITY_ZZ);
        Check(PERMEABILITY_YZ);
        Check(PERMEABILITY_ZX);
    }
    RestorePreviousBounds();
}

} /* namespace Kratos.*/
