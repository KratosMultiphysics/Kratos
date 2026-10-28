//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   Fluid Dynamics Application
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

#pragma once

// Project includes
#include "includes/condition.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/**
 * @brief Two-sided Nitsche coupling of independent 2D P1 or Q1 fluid patches.
 * @details With [u]=u_A-u_B, n pointing from A to B and arithmetic traction
 * average {sigma n}, the interface form is
 *   -[v].{sigma(u,p)n} - theta {sigma(v,q)n}.[u] + alpha [v].[u].
 * PENALTY_COEFFICIENT>0 selects theta=1 and viscosity-scaled penalty;
 * zero (or the IGA convention -1) selects theta=-1 and alpha=0.
 * The Newtonian deviatoric stress matches Newtonian2DLaw, including its
 * -2 mu div(u)/3 term. Viscosities come from the two parent Properties.
 *
 * NEIGHBOUR_ELEMENTS contains two Triangle2D3 or two Quadrilateral2D4 parents.
 * The generic geometry contains their six or eight independent nodes, in
 * parent order. Both registrations share this geometry-driven implementation.
 * COUPLING_INTERFACE_COORDINATES contains the oriented 2x3 segment; it
 * can be a subsegment of either edge for nonmatching interface partitions.
 * No multiplier, pressure penalty, shifted evaluation or gap is introduced.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidPatchCouplingCondition2D : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FluidPatchCouplingCondition2D);

    FluidPatchCouplingCondition2D(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry) {}

    FluidPatchCouplingCondition2D(IndexType NewId, GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties) {}

    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties) const override;

    Condition::Pointer Create(IndexType NewId, const NodesArrayType& rNodes,
                              PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLHS, VectorType& rRHS,
                              const ProcessInfo& rProcessInfo) override;
    void CalculateLeftHandSide(MatrixType& rLHS, const ProcessInfo& rProcessInfo) override;
    void CalculateRightHandSide(VectorType& rRHS, const ProcessInfo& rProcessInfo) override;
    void CalculateLocalVelocityContribution(MatrixType& rLHS, VectorType& rRHS,
                                             const ProcessInfo& rProcessInfo) override;
    void CalculateDampingMatrix(MatrixType& rLHS, const ProcessInfo& rProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rProcessInfo) const override;
    void GetDofList(DofsVectorType& rDofs, const ProcessInfo& rProcessInfo) const override;
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;
    int Check(const ProcessInfo& rProcessInfo) const override;

    std::string Info() const override { return "FluidPatchCouplingCondition2D"; }

protected:
    /// Only the dedicated Gap-SBM condition may extrapolate the background.
    virtual bool UsesBackgroundExtension() const { return false; }

    FluidPatchCouplingCondition2D() = default;

private:
    void CalculateCouplingMatrix(MatrixType& rLHS, const ProcessInfo& rProcessInfo) const;
    void CalculateResidual(const MatrixType& rLHS, VectorType& rRHS) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }
};

} // namespace Kratos
