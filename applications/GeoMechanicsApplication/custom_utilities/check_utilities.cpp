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

void CheckUtilities::CheckHasNodalSolutionStepData(const Geometry<Node>& rGeometry, 
    const Geo::ConstVariableDataRefs& rVariableRefs)
{
    auto has_all_variables = [&](const auto& node) {
        return std::all_of(rVariableRefs.begin(), rVariableRefs.end(), [&node](const auto& variable_ref) {
            return node.SolutionStepsDataHas(variable_ref.get());
        });
    };

    if (!std::all_of(rGeometry.begin(), rGeometry.end(), has_all_variables)) {
        for (const auto& r_node : rGeometry) {
            for (const auto& r_variable_ref : rVariableRefs) {
                if (!r_node.SolutionStepsDataHas(r_variable_ref.get()))
                    KRATOS_INFO("Missing variable ") << r_variable_ref.get().Name()
                                 << " on node " << r_node.Id() << std::endl;
            }
        }
        KRATOS_ERROR << "Missing variables on nodes" << std::endl;
    }
}

void CheckUtilities::CheckHasDofs(const Geometry<Node>&            rGeometry, 
    const Geo::ConstVariableDataRefs& rVariableRefs)
{
    for (const auto& r_node : rGeometry) {
        for (const auto& r_variable_ref : rVariableRefs) {
            if (!r_node.HasDofFor(r_variable_ref.get()))
                KRATOS_ERROR << "Missing the DoF for the variable " << r_variable_ref.get().Name()
                             << " on node " << r_node.Id() << std::endl;
        }
    }
}

} /* namespace Kratos.*/
