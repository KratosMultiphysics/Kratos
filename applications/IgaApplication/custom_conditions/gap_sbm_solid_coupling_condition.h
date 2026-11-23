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

#pragma once

// System includes
#include <vector>

// Project includes
#include "custom_conditions/solid_coupling_condition.h"

namespace Kratos
{

/**
 * @brief Gap-SBM solid coupling condition with both sides shifted to the true interface.
 *
 * The condition geometry is the true interface A. The A-side DOFs are taken from
 * NEIGHBOUR_GEOMETRIES[0], while the B-side DOFs are taken from
 * NEIGHBOUR_GEOMETRIES[1].
 */
class KRATOS_API(IGA_APPLICATION) GapSbmSolidCouplingCondition
    : public SolidCouplingCondition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GapSbmSolidCouplingCondition);

    using BaseType = SolidCouplingCondition;
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    GapSbmSolidCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    GapSbmSolidCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    GapSbmSolidCouplingCondition() = default;
    ~GapSbmSolidCouplingCondition() override = default;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GapSbmSolidCouplingCondition>(NewId, pGeom, pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GapSbmSolidCouplingCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"GapSbmSolidCouplingCondition\" #" << Id();
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"GapSbmSolidCouplingCondition\" #" << Id();
    }

protected:
    const GeometryType& GetGeometryA() const override;

    const GeometryType& GetGeometryB() const override;

    const GeometryType& GetConstitutiveGeometryA() const override;

    const GeometryType& GetConstitutiveGeometryB() const override;

    bool ShiftGeometryA() const override;

    bool ShiftGeometryB() const override;

private:
    const std::vector<Geometry<Node>::Pointer>& GetMainNeighbourGeometries() const;
};

} // namespace Kratos
