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
#include "check_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/exception.h"

#include <sstream>

namespace Kratos
{

CheckProperties::CheckProperties(const Properties&       rProperties,
                                 const std::string&      rPrintName,
                                 CheckProperties::Bounds RangeBoundsType)
    : mrProperties(rProperties), mrPrintName(rPrintName), mRangeBoundsType(RangeBoundsType)
{
}

CheckProperties::CheckProperties(const Properties&       rProperties,
                                 const std::string&      rPrintName,
                                 std::size_t             ElementId,
                                 CheckProperties::Bounds RangeBoundsType)
    : mrProperties(rProperties), mrPrintName(rPrintName), mElementId(ElementId), mRangeBoundsType(RangeBoundsType)
{
}

CheckProperties CheckProperties::SingleUseBounds(CheckProperties::Bounds RangeBoundsType) const
{
    if (mElementId) {
        return CheckProperties(mrProperties, mrPrintName, *mElementId, RangeBoundsType);
    }
    return CheckProperties(mrProperties, mrPrintName, RangeBoundsType);
}

void CheckProperties::SetNewRangeBounds(CheckProperties::Bounds RangeBoundsType) const
{
    mRangeBoundsType = RangeBoundsType;
}

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

void CheckProperties::CheckPermeabilityProperties(size_t Dimension) const
{
    const auto original_bounds_type = mRangeBoundsType;
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
    mRangeBoundsType = original_bounds_type;
}

void CheckUtilities::CheckForNonZeroZCoordinateIn2D(const Geometry<Node>& rGeometry)
{
    auto pos = std::ranges::find_if(rGeometry, [](const auto& node) {
        // We may want to rethink this local tolerance. Perhaps it should be moved to a more general
        // place. Also, since we compare a length, we may want to account for the length unit.
        constexpr auto absolute_tolerance = 1.0e-12;
        return std::abs(node.Z()) > absolute_tolerance;
    });
    KRATOS_ERROR_IF_NOT(pos == rGeometry.end())
        << "Node with Id: " << pos->Id() << " has non-zero Z coordinate." << std::endl;
}

std::string CheckProperties::double_to_string(double value) const
{
    std::ostringstream oss;
    oss << std::defaultfloat << value;
    return oss.str();
}

std::string CheckProperties::print_property_id() const
{
    std::ostringstream oss;
    oss << " in the " << mrPrintName << " with Id " << mrProperties.Id();
    return oss.str();
}

std::string CheckProperties::print_element_id() const
{
    std::ostringstream oss;
    if (mElementId) oss << " at element with Id " << *mElementId;
    return oss.str();
}
} /* namespace Kratos.*/
