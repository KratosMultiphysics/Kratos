//  KRATOS  ____       _               _
//         / ___|  ___(_)___ _ __ ___ (_) ___
//         \___ \ / _ \ / __| '_ ` _ \| |/ __|
//          ___) |  __/ \__ \ | | | | | | (__
//         |____/ \___|_|___/_| |_| |_|_|\___|
//
//  License:     BSD License
//  license:     structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//    Co authors: Long Chen
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "seismic_application_variables.h"
#include "custom_classes/fiber_beam_column_section.h"

namespace Kratos {

FiberBeamColumnSection::FiberBeamColumnSection(IndexType NewId)
    : BaseType(NewId) {}

FiberBeamColumnSection::FiberBeamColumnSection(
    IndexType NewId,
    IntegrationPointType integrationPoint)
        : BaseType(NewId),
          mPosition(integrationPoint(0)),
          mWeight(integrationPoint.Weight())
{
}

void FiberBeamColumnSection::Initialize()
{
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        r_fiber.Initialize();
    }
    UpdateLocalFlexibilityMatrix();
}

void FiberBeamColumnSection::GetGlobalFlexibilityMatrix(Matrix& rGlobalFlexibilityMatrix)
{
    KRATOS_TRY

    if (rGlobalFlexibilityMatrix.size1() != 5 || rGlobalFlexibilityMatrix.size2() != 5) {
        rGlobalFlexibilityMatrix.resize(5, 5, false);
    }
    noalias(rGlobalFlexibilityMatrix) = ZeroMatrix(5, 5);

    Matrix aux_matrix  = ZeroMatrix(5, 3);
    Matrix b_matrix;
    CalculateBMatrix(b_matrix);
    noalias(aux_matrix) += prod(Matrix(trans(b_matrix)), mLocalFlexibilityMatrix);
    noalias(rGlobalFlexibilityMatrix) += prod(aux_matrix, b_matrix);

    KRATOS_CATCH("")
}

void FiberBeamColumnSection::UpdateLocalFlexibilityMatrix()
{
    KRATOS_TRY
    // allocate memory
    Matrix local_stiffness = ZeroMatrix(3, 3);
    // get local 3x3 stiffness from fibers' tangential stiffness
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        Matrix fiber_global_stiffness;
        r_fiber.GetGlobalStiffnessMatrix(fiber_global_stiffness);
        local_stiffness += fiber_global_stiffness;
    }
    // invert stiffness to get 3x3 flexibility
    double det_stiffness = MathUtils<double>::Det(local_stiffness);
    MathUtils<double>::InvertMatrix3(local_stiffness, mLocalFlexibilityMatrix, det_stiffness);
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::CalculateBMatrix(Matrix& rBMatrix)
{
    KRATOS_TRY

    if (rBMatrix.size1() != 3 || rBMatrix.size2() != 5) {
        rBMatrix.resize(3, 5, false);
    }
    noalias(rBMatrix) = ZeroMatrix(3, 5);

    rBMatrix(0, 0) = mPosition/2.0 - 0.5;
    rBMatrix(0, 1) = mPosition/2.0 + 0.5;
    rBMatrix(1, 2) = mPosition/2.0 - 0.5;
    rBMatrix(1, 3) = mPosition/2.0 + 0.5;
    rBMatrix(2, 4) = 1.0;

    KRATOS_CATCH("")
}

bool FiberBeamColumnSection::StateDetermination(const Vector& rElementForceIncrements, const double Tolerance)
{
    KRATOS_TRY
    Vector force_incr = ZeroVector(3);
    Vector def_incr   = ZeroVector(3);

    Matrix b_matrix;
    CalculateBMatrix(b_matrix);
    noalias(force_incr) = prod(b_matrix, rElementForceIncrements);
    mForces += force_incr;

    noalias(def_incr) = mDeformationResiduals + prod(mLocalFlexibilityMatrix, force_incr);

    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        r_fiber.StateDetermination(def_incr);
    }

    UpdateLocalFlexibilityMatrix();

    Vector internal_forces = ZeroVector(3);
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        Vector fiber_global_forces;
        r_fiber.GetGlobalInternalForces(fiber_global_forces);
        internal_forces += fiber_global_forces;
    }

    mUnbalanceForces = mForces - internal_forces;
    noalias(mDeformationResiduals) = prod(mLocalFlexibilityMatrix, mUnbalanceForces);
    double residual = std::abs(MathUtils<double>::Norm(mUnbalanceForces));
    return residual < Tolerance;

    KRATOS_CATCH("")
}

void FiberBeamColumnSection::ResetResidual()
{
    KRATOS_TRY
    noalias(mDeformationResiduals) = ZeroVector(3);
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::GetGlobalDeformationResiduals(Vector& rGlobalResiduals)
{
    KRATOS_TRY

    if (rGlobalResiduals.size() != 5) {
        rGlobalResiduals.resize(5, false);
    }

    Matrix b_matrix;
    CalculateBMatrix(b_matrix);
    noalias(rGlobalResiduals) = prod(Matrix(trans(b_matrix)), mDeformationResiduals);
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
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    rSerializer.save("mPosition", mPosition);
    rSerializer.save("mWeight", mWeight);
    rSerializer.save("mLocalFlexibilityMatrix", mLocalFlexibilityMatrix);
    rSerializer.save("mForces", mForces);
    rSerializer.save("mUnbalanceForces", mUnbalanceForces);
    rSerializer.save("mDeformationResiduals", mDeformationResiduals);
}
void FiberBeamColumnSection::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    rSerializer.load("mPosition", mPosition);
    rSerializer.load("mWeight", mWeight);
    rSerializer.load("mLocalFlexibilityMatrix", mLocalFlexibilityMatrix);
    rSerializer.load("mForces", mForces);
    rSerializer.load("mUnbalanceForces", mUnbalanceForces);
    rSerializer.load("mDeformationResiduals", mDeformationResiduals);
}

} // namespace Kratos.
