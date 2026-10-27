// License: BSD License (kratos/license.txt)
#pragma once
#include "custom_conditions/fluid_patch_coupling_condition.h"
#include "custom_conditions/sbm_fluid_gap_interface_condition.h"
#include "integration/line_gauss_legendre_integration_points.h"

namespace Kratos {
/** @brief True-interface Nitsche coupling: extended background / fitted inner Q1.
 * Parent 0 is Cartesian background; parent 1 is the independent fitted patch.
 * Gap volume and connector contributions are required alongside this condition.
 * The pressure-test term is -{q}[u].n for both velocity adjoint choices,
 * matching SymbolicStokes' q div(u) convention. With a positive penalty the
 * velocity block is symmetric and the mixed blocks are transpose negatives.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidPatchGapCouplingCondition2D
    : public FluidPatchCouplingCondition2D
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FluidPatchGapCouplingCondition2D);
    using FluidPatchCouplingCondition2D::FluidPatchCouplingCondition2D;
    Condition::Pointer Create(IndexType Id, GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<FluidPatchGapCouplingCondition2D>(Id, pGeometry, pProperties);
    }
    Condition::Pointer Create(IndexType Id, const NodesArrayType& rNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Create(Id, Kratos::make_shared<Geometry<Node>>(rNodes), pProperties);
    }
    std::string Info() const override { return "FluidPatchGapCouplingCondition2D"; }
    using FluidPatchCouplingCondition2D::CalculateOnIntegrationPoints;
    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rValues, const ProcessInfo& rInfo) override
    {
        if (rVariable != INTEGRATION_COORDINATES) {
            FluidPatchCouplingCondition2D::CalculateOnIntegrationPoints(rVariable, rValues, rInfo);
            return;
        }
        rValues.clear();
        const auto& points = GetValue(COUPLING_INTERFACE_COORDINATES);
        for (const auto& g:LineGaussLegendreIntegrationPoints2::IntegrationPoints()) {
            array_1d<double, 3> point = ZeroVector(3);
            for (std::size_t d = 0; d<3; ++d)
                point[d] = 0.5*((1-g.X())*points(0, d)+(1+g.X())*points(1, d));
            rValues.push_back(point);
        }
    }
protected:
    bool UsesBackgroundExtension() const override { return true; }
private:
    FluidPatchGapCouplingCondition2D() = default;
    friend class Serializer;
    void save(Serializer& rSerializer) const override
    { KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, FluidPatchCouplingCondition2D); }
    void load(Serializer& rSerializer) override
    { KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, FluidPatchCouplingCondition2D); }
};

/** @brief Reuse the existing shared-node gap connector assembly with the
 * multipatch symmetry switch, parent materials and viscosity/h penalty scale.
 * Legacy SBM connector behavior is unchanged.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidPatchGapInterfaceCondition2D
    : public SbmFluidGapInterfaceCondition2D
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FluidPatchGapInterfaceCondition2D);
    using SbmFluidGapInterfaceCondition2D::SbmFluidGapInterfaceCondition2D;
    Condition::Pointer Create(IndexType Id, GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<FluidPatchGapInterfaceCondition2D>(Id, pGeometry, pProperties);
    }
    Condition::Pointer Create(IndexType Id, const NodesArrayType& rNodes,
                              PropertiesType::Pointer pProperties) const override
    { return Create(Id, Kratos::make_shared<Geometry<Node>>(rNodes), pProperties); }
    std::string Info() const override { return "FluidPatchGapInterfaceCondition2D"; }
    using SbmFluidGapInterfaceCondition2D::CalculateOnIntegrationPoints;
    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rValues, const ProcessInfo& rInfo) override
    {
        if (rVariable != INTEGRATION_COORDINATES) {
            SbmFluidGapInterfaceCondition2D::CalculateOnIntegrationPoints(rVariable, rValues, rInfo);
            return;
        }
        rValues.clear();
        const auto& points = GetValue(SURROGATE_BOUNDARY_FACE_COORDINATES);
        for (const auto& g:LineGaussLegendreIntegrationPoints2::IntegrationPoints()) {
            array_1d<double, 3> point = ZeroVector(3);
            for (std::size_t d = 0; d<3; ++d)
                point[d] = 0.5*((1-g.X())*points(0, d)+(1+g.X())*points(1, d));
            rValues.push_back(point);
        }
    }
protected:
    bool UsesPatchCouplingParameters() const override { return true; }
private:
    FluidPatchGapInterfaceCondition2D() = default;
    friend class Serializer;
    void save(Serializer& rSerializer) const override
    { KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SbmFluidGapInterfaceCondition2D); }
    void load(Serializer& rSerializer) override
    { KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SbmFluidGapInterfaceCondition2D); }
};
} // namespace Kratos
