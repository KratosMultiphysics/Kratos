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

void CheckUtilities::CheckHasNodalSolutionStepData(const Geometry<Node>& rGeometry, 
    const Geo::ConstVariableDataRefs& rVariableRefs)
{
    for (const auto& r_variable_ref : rVariableRefs) {
        std::vector<std::size_t> missing_node_ids;
        std::copy_if(rGeometry.begin(), rGeometry.end(), std::back_inserter(missing_node_ids),
                     [&r_variable_ref, &missing_node_ids](const Node& node) {
            if (!node.SolutionStepsDataHas(r_variable_ref.get())) {
                missing_node_ids.push_back(node.Id());
            }
            return false;
        });
        if (!missing_node_ids.empty())
           KRATOS_ERROR << "Missing variable " << r_variable_ref.get().Name() << " on nodes "
                         << PrintVectorContent(missing_node_ids) << std::endl;
    }
}

void CheckUtilities::CheckHasDofs(const Geometry<Node>& rGeometry, const Geo::ConstVariableDataRefs& rVariableRefs)
{
    for (const auto& r_variable_ref : rVariableRefs) {
        std::vector<std::size_t> missing_node_ids;
        std::copy_if(
            rGeometry.begin(), rGeometry.end(), std::back_inserter(missing_node_ids),
            [&r_variable_ref, &missing_node_ids](const Node& node) {
            if (!node.HasDofFor(r_variable_ref.get())) {
                missing_node_ids.push_back(node.Id());
            }
            return false;
        });

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
    output.pop_back();

    return output;
}

} /* namespace Kratos.*/
