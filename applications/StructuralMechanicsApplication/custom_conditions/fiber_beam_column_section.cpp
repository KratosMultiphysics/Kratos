// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//  license:     structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "structural_mechanics_application_variables.h"
#include "custom_conditions/fiber_beam_column_section.hpp"
#include "custom_constitutive/uniaxial_fiber_beam_column_concrete_material_law.hpp"
#include "custom_constitutive/uniaxial_fiber_beam_column_steel_material_law.hpp"

namespace Kratos {

FiberBeamColumnSection::FiberBeamColumnSection(IndexType NewId)
    : mId(NewId) {}

FiberBeamColumnSection::FiberBeamColumnSection(
    IndexType NewId,
    IntegrationPointType integrationPoint,
    PropertiesType::Pointer pProperties)
        : mId(NewId),
          mPosition(integrationPoint(0)),
          mWeight(integrationPoint.Weight())
{
    KRATOS_TRY
    mTolerance = (*pProperties)[ELEMENT_LOOP_TOLERANCE];
    if (std::abs(mTolerance) < std::numeric_limits<double>::epsilon()){
        KRATOS_WARNING("FiberBeamColumnSection") <<
            "WARNING: TOLERANCE FOR ELEMENT LOOP IS NOT SET." << std::endl;
    }
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::Initialize()
{
    CalculateBMatrix();
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        r_fiber.Initialize();
    }
    CalculateLocalFlexibilityMatrix();
}

Matrix FiberBeamColumnSection::GetGlobalFlexibilityMatrix()
{
    KRATOS_TRY
    Matrix aux_matrix  = ZeroMatrix(5, 3);
    Matrix global_flex = ZeroMatrix(5, 5);
    (aux_matrix) += prod(Matrix(trans(mBMatrix)), mLocalFlexibilityMatrix);
    (global_flex) += prod(aux_matrix, mBMatrix);
    return global_flex;
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::CalculateLocalFlexibilityMatrix()
{
    KRATOS_TRY
    // allocate memory
    Matrix local_stiffness = ZeroMatrix(3, 3);
    // get local 3x3 stiffness from fibers' tangential stiffness
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        (local_stiffness) += r_fiber.CreateGlobalFiberStiffnessMatrix();
    }
    // invert stiffness to get 3x3 flexibility
    double det_stiffness = MathUtils<double>::Det(local_stiffness);
    MathUtils<double>::InvertMatrix3(local_stiffness, mLocalFlexibilityMatrix, det_stiffness);
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::CalculateBMatrix()
{
    KRATOS_TRY
    mBMatrix(0, 0) = mPosition/2.0 - 0.5;
    mBMatrix(0, 1) = mPosition/2.0 + 0.5;
    mBMatrix(1, 2) = mPosition/2.0 - 0.5;
    mBMatrix(1, 3) = mPosition/2.0 + 0.5;
    mBMatrix(2, 4) = 1.0;
    KRATOS_CATCH("")
}

bool FiberBeamColumnSection::StateDetermination(const Vector& rElementForceIncrements)
{
    KRATOS_TRY
    Vector force_incr = ZeroVector(3);
    Vector def_incr   = ZeroVector(3);

    force_incr = prod(mBMatrix, rElementForceIncrements);
    mForces += force_incr;

    def_incr = mDeformationResiduals + prod(mLocalFlexibilityMatrix, force_incr);

    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        r_fiber.StateDetermination(def_incr);
    }

    CalculateLocalFlexibilityMatrix();

    KRATOS_WATCH(mLocalFlexibilityMatrix)

    Vector internal_forces = ZeroVector(3);
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        (internal_forces) += r_fiber.CreateGlobalFiberInternalForces();
    }

    mUnbalanceForces = mForces - internal_forces;
    (mDeformationResiduals) = prod(mLocalFlexibilityMatrix, mUnbalanceForces);
    double residual = std::abs(MathUtils<double>::Norm3(mUnbalanceForces));
    return residual < mTolerance;

    KRATOS_CATCH("")
}

void FiberBeamColumnSection::ResetResidual()
{
    KRATOS_TRY
    mDeformationResiduals = ZeroVector(3);
    KRATOS_CATCH("")
}

Vector FiberBeamColumnSection::GetGlobalDeformationResiduals()
{
    KRATOS_TRY
    return prod(Matrix(trans(mBMatrix)), mDeformationResiduals);
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::FinalizeSolutionStep()
{
    KRATOS_TRY
    for(FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        r_fiber.FinalizeSolutionStep();
    }
    KRATOS_CATCH("")
}

std::string FiberBeamColumnSection::Info() const {
    std::stringstream buffer;
    buffer << "FiberBeamColumnSection #" << Id();
    return buffer.str();
}

void FiberBeamColumnSection::PrintInfo(std::ostream& rOStream) const {
    rOStream << Info();
}

void FiberBeamColumnSection::PrintData(std::ostream& rOStream) const {
    rOStream << "    Position : " << mPosition << std::endl;
    rOStream << "    Weight   : " << mWeight;
}

void FiberBeamColumnSection::save(Serializer& rSerializer) const
{
    rSerializer.save("mId", mId);
    rSerializer.save("mPosition", mPosition);
    rSerializer.save("mWeight", mWeight);
    rSerializer.save("mBMatrix", mBMatrix);
    rSerializer.save("mLocalFlexibilityMatrix", mLocalFlexibilityMatrix);
    rSerializer.save("mTolerance", mTolerance);
    rSerializer.save("mForces", mForces);
    rSerializer.save("mUnbalanceForces", mUnbalanceForces);
    rSerializer.save("mDeformationResiduals", mDeformationResiduals);
}
void FiberBeamColumnSection::load(Serializer& rSerializer)
{
    rSerializer.load("mId", mId);
    rSerializer.load("mPosition", mPosition);
    rSerializer.load("mWeight", mWeight);
    rSerializer.load("mBMatrix", mBMatrix);
    rSerializer.load("mLocalFlexibilityMatrix", mLocalFlexibilityMatrix);
    rSerializer.load("mTolerance", mTolerance);
    rSerializer.load("mForces", mForces);
    rSerializer.load("mUnbalanceForces", mUnbalanceForces);
    rSerializer.load("mDeformationResiduals", mDeformationResiduals);
}

} // namespace Kratos.
