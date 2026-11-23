//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi

// Project includes
#include "custom_conditions/gap_sbm_solid_coupling_condition.h"

namespace Kratos
{

const GapSbmSolidCouplingCondition::GeometryType& GapSbmSolidCouplingCondition::GetGeometryA() const
{
    const auto& r_neighbour_geometries = GetMainNeighbourGeometries();
    KRATOS_ERROR_IF(r_neighbour_geometries.size() < 1)
        << "\"GapSbmSolidCouplingCondition\" #" << Id()
        << " requires NEIGHBOUR_GEOMETRIES[0] for surrogate geometry A." << std::endl;

    return *r_neighbour_geometries[0];
}

const GapSbmSolidCouplingCondition::GeometryType& GapSbmSolidCouplingCondition::GetGeometryB() const
{
    const auto& r_neighbour_geometries = GetMainNeighbourGeometries();
    KRATOS_ERROR_IF(r_neighbour_geometries.size() < 2)
        << "\"GapSbmSolidCouplingCondition\" #" << Id()
        << " requires NEIGHBOUR_GEOMETRIES[1] for surrogate geometry B." << std::endl;
    return *r_neighbour_geometries[1];
}

const GapSbmSolidCouplingCondition::GeometryType& GapSbmSolidCouplingCondition::GetConstitutiveGeometryA() const
{
    return GetTrueGeometry();
}

const GapSbmSolidCouplingCondition::GeometryType& GapSbmSolidCouplingCondition::GetConstitutiveGeometryB() const
{
    return GetTrueGeometry();
}

bool GapSbmSolidCouplingCondition::ShiftGeometryA() const
{
    return true;
}

bool GapSbmSolidCouplingCondition::ShiftGeometryB() const
{
    return true;
}

const std::vector<Geometry<Node>::Pointer>& GapSbmSolidCouplingCondition::GetMainNeighbourGeometries() const
{
    if (this->Has(NEIGHBOUR_GEOMETRIES)) {
        const auto& r_neighbour_geometries = this->GetValue(NEIGHBOUR_GEOMETRIES);
        if (!r_neighbour_geometries.empty()) {
            return r_neighbour_geometries;
        }
    }

    KRATOS_ERROR_IF_NOT(GetGeometry().Has(NEIGHBOUR_GEOMETRIES))
        << "\"GapSbmSolidCouplingCondition\" #" << Id()
        << " has no NEIGHBOUR_GEOMETRIES set on the condition or its main geometry." << std::endl;

    return GetGeometry().GetValue(NEIGHBOUR_GEOMETRIES);
}

} // namespace Kratos
