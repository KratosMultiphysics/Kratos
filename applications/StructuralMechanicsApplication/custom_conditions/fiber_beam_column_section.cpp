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
    : mId(NewId), mPosition(integrationPoint(0)), mWeight(integrationPoint.Weight())
{
    KRATOS_TRY

    const unsigned int number_of_fibers_y = (*pProperties)[NUMBER_FIBERS_Y];
    const unsigned int number_of_fibers_z = (*pProperties)[NUMBER_FIBERS_Z];
    const double width = (*pProperties)[BEAM_WIDTH];
    const double height = (*pProperties)[BEAM_HEIGHT];

    // check if all properties are given
    if ((number_of_fibers_y == 0) || (number_of_fibers_z == 0) || (width == 0) || (height == 0)) {
        KRATOS_ERROR << "Error in setting the number of fibers or beam dimensions." << std::endl;
    }

    if (mFibers.size() != number_of_fibers_y * number_of_fibers_z) {
        mFibers.resize(number_of_fibers_y * number_of_fibers_z);
    }

    IndexType id = 1 + (mId-1)*number_of_fibers_y*number_of_fibers_z;
    double y = 0.0;
    double z = 0.0;
    double area = (width / number_of_fibers_y) * (height / number_of_fibers_z);
    for (unsigned int i = 0; i < number_of_fibers_y; ++i){
        y = width / number_of_fibers_y * (i + 0.5);
        for (unsigned int j = 0; j < number_of_fibers_z; j++)
        {
            z = height / number_of_fibers_z * (j + 0.5);
            auto p_material = Kratos::make_shared<UniaxialFiberBeamColumnSteelMaterialLaw>(pProperties);
            mFibers[i*number_of_fibers_y + j] = FiberBeamColumnUniaxialFiber(
                id, y, z, area, p_material );
            // // Steel fibers
            // if (((i == 1) || (i == number_of_fibers_y - 2)) &&
            //     ((j == 1) || (j == number_of_fibers_y - 2))) {
            //     auto p_material = Kratos::make_shared<UniaxialFiberBeamColumnSteelMaterialLaw>(pProperties);
            //     mFibers[i*number_of_fibers_y + j] = FiberBeamColumnUniaxialFiber(
            //         id, y, z, area, p_material );
            // // Concrete fibers
            // } else {
            //     auto p_material = Kratos::make_shared<UniaxialFiberBeamColumnConcreteMaterialLaw>(pProperties);
            //     mFibers[i*number_of_fibers_y + j] = FiberBeamColumnUniaxialFiber(
            //         id, y, z, area, p_material );
            // }
            id++;
        }

    }

    KRATOS_CATCH("")
}

// FiberBeamColumnSection::FiberBeamColumnSection(FiberBeamColumnSection const& rOther)
//     : mId(rOther.mId), mPosition(rOther.mPosition), mWeight(rOther.mWeight) {}

// FiberBeamColumnSection::~FiberBeamColumnSection() {}

// FiberBeamColumnSection& FiberBeamColumnSection::operator=(FiberBeamColumnSection const& rOther)
// {
//     return *this;
// }

void FiberBeamColumnSection::Initialize()
{
    CalculateLocalFlexibilityMatrix(mLocalFlexibilityMatrix);
    CalculateBMatrix(mBMatrix);
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        r_fiber.Initialize();
    }
}

void FiberBeamColumnSection::CalculateLocalStiffnessMatrix(Matrix& rLocalStiffnessMatrix) {
    KRATOS_TRY
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        rLocalStiffnessMatrix += r_fiber.CreateGlobalFiberStiffnessMatrix();
    }
    KRATOS_CATCH("")
}

Matrix FiberBeamColumnSection::GetGlobalFlexibilityMatrix()
{
    KRATOS_TRY

    // allocate memory
    Matrix aux_matrix         = ZeroMatrix(5, 3);
    Matrix global_flexibility = ZeroMatrix(5, 5);

    aux_matrix += prod(trans(mBMatrix), mLocalFlexibilityMatrix);
    global_flexibility += prod(aux_matrix, mBMatrix);

    return global_flexibility;
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::CalculateLocalFlexibilityMatrix(Matrix& rLocalFlexibilityMatrix)
{
    KRATOS_TRY

    // allocate memory
    Matrix local_stiffness = ZeroMatrix(3);
    // get local 3x3 stiffness from fibers' tangential stiffness
    CalculateLocalStiffnessMatrix(local_stiffness);
    // invert stiffness to get 3x3 flexibility
    double det_stiffness = MathUtils<double>::Det(local_stiffness);
    MathUtils<double>::InvertMatrix3(local_stiffness, rLocalFlexibilityMatrix, det_stiffness);

    KRATOS_CATCH("")
}

void FiberBeamColumnSection::CalculateBMatrix(Matrix& rBMatrix)
{
    KRATOS_TRY
    rBMatrix = ZeroMatrix(3, 5);
    rBMatrix(0, 0) = mPosition/2.0 - 0.5;
    rBMatrix(0, 1) = mPosition/2.0 + 0.5;
    rBMatrix(1, 2) = mPosition/2.0 - 0.5;
    rBMatrix(1, 3) = mPosition/2.0 + 0.5;
    rBMatrix(2, 4) = 1.0;
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::SetSectionForces(Vector ChangeElementForceIncrements)
{
    KRATOS_TRY

    noalias(mChangeForceIncrements) = prod(mBMatrix, ChangeElementForceIncrements); //FIXME: not mem
    mForceIncrements += mChangeForceIncrements;
    mForces = mConvergedForces + mForceIncrements;

    KRATOS_CATCH("")
}

void FiberBeamColumnSection::CalculateDeformationIncrements()
{
    KRATOS_TRY
    mChangeDeformationIncrements = mDeformationResiduals + prod(mLocalFlexibilityMatrix, mChangeForceIncrements);
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::CalculateFiberState()
{
    KRATOS_TRY
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers){
        r_fiber.SetStrainIncrements(mChangeDeformationIncrements);
    }
    CalculateLocalFlexibilityMatrix(mLocalFlexibilityMatrix);
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::CalculateResiduals()
{
    KRATOS_TRY
    Vector resisting_forces = ZeroVector(3);
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        resisting_forces += r_fiber.CreateGlobalFiberResistingForces();
    }
    mUnbalanceForces = mForces - resisting_forces;
    mDeformationResiduals = prod(mLocalFlexibilityMatrix, mUnbalanceForces);
    KRATOS_CATCH("")
}

Vector FiberBeamColumnSection::GetGlobalDeformationResiduals()
{
    KRATOS_TRY
    return prod(trans(mBMatrix), mDeformationResiduals);
    KRATOS_CATCH("")
}

bool FiberBeamColumnSection::CheckConvergence()
{
    KRATOS_TRY
    return MathUtils<double>::Norm(mUnbalanceForces) < mTolerance;
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::FinalizeSolutionStep()
{
    KRATOS_TRY

    std::fill(mDeformationResiduals.begin(), mDeformationResiduals.end(), 0.00);

    mConvergedForces = mForces;
    std::fill(mForceIncrements.begin(), mForceIncrements.end(), 0.00);

    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers){
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
    rSerializer.save("mChangeForceIncrements", mChangeForceIncrements);
    rSerializer.save("mForceIncrements", mForceIncrements);
    rSerializer.save("mForces", mForces);
    rSerializer.save("mConvergedForces", mConvergedForces);
    rSerializer.save("mUnbalanceForces", mUnbalanceForces);
    rSerializer.save("mChangeDeformationIncrements", mChangeDeformationIncrements);
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
    rSerializer.load("mChangeForceIncrements", mChangeForceIncrements);
    rSerializer.load("mForceIncrements", mForceIncrements);
    rSerializer.load("mForces", mForces);
    rSerializer.load("mConvergedForces", mConvergedForces);
    rSerializer.load("mUnbalanceForces", mUnbalanceForces);
    rSerializer.load("mChangeDeformationIncrements", mChangeDeformationIncrements);
    rSerializer.load("mDeformationResiduals", mDeformationResiduals);
}

} // namespace Kratos.
