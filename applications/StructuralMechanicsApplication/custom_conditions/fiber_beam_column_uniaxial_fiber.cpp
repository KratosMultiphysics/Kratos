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
#include<iostream>

// External includes

// Project includes
#include "custom_conditions/fiber_beam_column_uniaxial_fiber.hpp"

namespace Kratos {

FiberBeamColumnUniaxialFiber::FiberBeamColumnUniaxialFiber(IndexType NewId)
    : mId(NewId) {}

FiberBeamColumnUniaxialFiber::FiberBeamColumnUniaxialFiber(
    IndexType NewId, double Y, double Z, double Area, UniaxialFiberBeamColumnMaterialLaw::Pointer pMaterial)
    : mId(NewId)
{
    mTransformationVector[0] = -Y;
    mTransformationVector[1] =  Z;
    mTransformationVector[2] = 1.0;
    mArea = Area;
    mpMaterial = pMaterial;
}

void FiberBeamColumnUniaxialFiber::Initialize()
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

Matrix FiberBeamColumnUniaxialFiber::CreateGlobalFiberStiffnessMatrix(){
    KRATOS_TRY
    double ea = mpMaterial->GetTangentModulus() * mArea;
    Matrix global_stiffness = ZeroMatrix(3);
    CalculateTransformationMatrix(global_stiffness);
    global_stiffness *= ea;
    return global_stiffness;
    KRATOS_CATCH("")
}

void FiberBeamColumnUniaxialFiber::CalculateTransformationMatrix(Matrix& rTransformationMatrix){
    KRATOS_TRY
    rTransformationMatrix = outer_prod(mTransformationVector, mTransformationVector);
    KRATOS_CATCH("")
}

void FiberBeamColumnUniaxialFiber::SetStrainIncrements(Vector ChangeSectionDeformationIncrements)
{
    KRATOS_TRY
    mChangeStrainIncrement = inner_prod(mTransformationVector, ChangeSectionDeformationIncrements);
    mStrainIncrement += mChangeStrainIncrement;
    mStrain = mConvergedStrain + mStrainIncrement;
    mpMaterial->SetStrain(mStrain);
    mpMaterial->CalculateMaterialResponse();
    mStress = mpMaterial->GetStress();
    KRATOS_CATCH("")
}

Vector FiberBeamColumnUniaxialFiber::CreateGlobalFiberResistingForces()
{
    KRATOS_TRY
    return mTransformationVector * mArea * mStress;
    KRATOS_CATCH("")
}

void FiberBeamColumnUniaxialFiber::FinalizeSolutionStep()
{
    KRATOS_TRY
    mConvergedStrain = mStrain;
    mStrainIncrement = 0.0;
    mpMaterial->FinalizeMaterialResponse();
    KRATOS_CATCH("")
}

void FiberBeamColumnUniaxialFiber::PrintInfo(std::ostream& rOStream) const {
    rOStream << "FiberBeamColumnUniaxialFiber #" << Id();
}

void FiberBeamColumnUniaxialFiber::PrintData(std::ostream& rOStream) const {
    rOStream << "    Position : [" << -1.0 * mTransformationVector[0] << "," << mTransformationVector[1] << "]" << std::endl;
    rOStream << "    Area     : " << mArea << std::endl;
    rOStream << "    Material : " << *mpMaterial << std::endl;
}


void FiberBeamColumnUniaxialFiber::save(Serializer& rSerializer) const
{
    rSerializer.save("mId", mId);
    rSerializer.save("mTransformationVector", mTransformationVector);
    rSerializer.save("mArea", mArea);
    rSerializer.save("mpMaterial", mpMaterial);
    rSerializer.save("mChangeStrainIncrement", mChangeStrainIncrement);
    rSerializer.save("mStrainIncrement", mStrainIncrement);
    rSerializer.save("mStrain", mStrain);
    rSerializer.save("mConvergedStrain", mConvergedStrain);
    rSerializer.save("mStress", mStress);
    rSerializer.save("mConvergedStress", mConvergedStress);
}
void FiberBeamColumnUniaxialFiber::load(Serializer& rSerializer)
{
    rSerializer.load("mId", mId);
    rSerializer.load("mTransformationVector", mTransformationVector);
    rSerializer.load("mArea", mArea);
    rSerializer.load("mpMaterial", mpMaterial);
    rSerializer.load("mChangeStrainIncrement", mChangeStrainIncrement);
    rSerializer.load("mStrainIncrement", mStrainIncrement);
    rSerializer.load("mStrain", mStrain);
    rSerializer.load("mConvergedStrain", mConvergedStrain);
    rSerializer.load("mStress", mStress);
    rSerializer.load("mConvergedStress", mConvergedStress);
}

} // namespace Kratos.
