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
#include "custom_constitutive/uniaxial_fiber_beam_column_concrete_material_law2.hpp"
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

    // const unsigned int number_of_fibers_y = (*pProperties)[NUMBER_FIBERS_Y];
    // const unsigned int number_of_fibers_z = (*pProperties)[NUMBER_FIBERS_Z];
    // const double width = (*pProperties)[BEAM_WIDTH];
    // const double height = (*pProperties)[BEAM_HEIGHT];

    // const double w = width / number_of_fibers_y;
    // const double h = height / number_of_fibers_z;

    // // check if all properties are given
    // if ((number_of_fibers_y == 0) || (number_of_fibers_z == 0) || (width == 0) || (height == 0)) {
    //     KRATOS_ERROR << "Error in setting the number of fibers or beam dimensions." << std::endl;
    // }

    // if (mFibers.size() != number_of_fibers_y * number_of_fibers_z) {
    //     mFibers.resize(number_of_fibers_y * number_of_fibers_z);
    // }

    // IndexType id = 1 + (mId-1)*number_of_fibers_y*number_of_fibers_z;
    // unsigned int counter = 0;
    // double y = 0.0;
    // double z = 0.0;
    // double area = (width / number_of_fibers_y) * (height / number_of_fibers_z);

    // for (unsigned int i = 0; i < number_of_fibers_y; ++i){
    //     // y = width / number_of_fibers_y * (i + 0.5);
    //     y = 0.5*(w-width) + i*w;

    //     for (unsigned int j = 0; j < number_of_fibers_z; ++j)
    //     {
    //         // z = height / number_of_fibers_z * (j + 0.5);
    //         z = 0.5*(h-height) + j*h;

    //         // Steel fibers
    //         if (((i == 1) || (i == number_of_fibers_y - 2)) &&
    //             ((j == 1) || (j == number_of_fibers_y - 2))) {
    //             UniaxialFiberBeamColumnSteelMaterialLaw::Pointer p_material =
    //                 Kratos::make_shared<UniaxialFiberBeamColumnSteelMaterialLaw>(pProperties);
    //             mFibers[counter] = FiberBeamColumnUniaxialFiber(
    //                 id, y, z, area, p_material );

    //         // // Unconfined Concrete fibers
    //         // } else if ( (i==0) || (j==0) || (i==number_of_fibers_y-1) || (j==number_of_fibers_z-1) ){
    //         //     std::cout << "i: " << i << " , j: " << j << std::endl;
    //         //     UniaxialFiberBeamColumnConcreteMaterialLaw::Pointer p_material =
    //         //     Kratos::make_shared<UniaxialFiberBeamColumnConcreteMaterialLaw>(pProperties);
    //         //     mFibers[counter] = FiberBeamColumnUniaxialFiber(
    //         //         id, y, z, area, p_material );

    //         // Confined Concrete fibers
    //         } else {
    //             UniaxialFiberBeamColumnConcreteMaterialLaw::Pointer p_material =
    //             Kratos::make_shared<UniaxialFiberBeamColumnConcreteMaterialLaw>(pProperties);
    //             p_material->Confine();
    //             mFibers[counter] = FiberBeamColumnUniaxialFiber(
    //                 id, y, z, area, p_material );
    //         }
    //         id++;
    //         counter++;
    //     }
    // }

    // Matrix concrete_fibers_data = (*pProperties)[CONCRETE_FIBERS_DATA];
    // Matrix steel_fibers_data = (*pProperties)[STEEL_FIBERS_DATA];
    // mFibers.reserve(concrete_fibers_data.size1() + steel_fibers_data.size1());

    // for (unsigned int i = 0; i < concrete_fibers_data.size1(); ++i) {
    //     auto p_material = Kratos::make_shared<UniaxialFiberBeamColumnConcreteMaterialLaw>(pProperties);
    //     p_material->Confine();
    //     mFibers.push_back( FiberBeamColumnUniaxialFiber( 0,
    //         concrete_fibers_data(i, 0),
    //         concrete_fibers_data(i, 1),
    //         concrete_fibers_data(i, 2),
    //         p_material) );
    // }

    // for (unsigned int i = 0; i < steel_fibers_data.size1(); ++i) {
    //     auto p_material = Kratos::make_shared<UniaxialFiberBeamColumnSteelMaterialLaw>(pProperties);
    //     mFibers.push_back( FiberBeamColumnUniaxialFiber( 0,
    //         steel_fibers_data(i, 0),
    //         steel_fibers_data(i, 1),
    //         steel_fibers_data(i, 2),
    //         p_material) );
    // }

    KRATOS_CATCH("")
}

void FiberBeamColumnSection::Initialize()
{
    CalculateBMatrix(mBMatrix);
    CalculateLocalFlexibilityMatrix();
    for (unsigned int i = 0; i < mFibers.size(); ++i) {
        mFibers[i].Initialize();
    }
}

// void FiberBeamColumnSection::CalculateLocalStiffnessMatrix(Matrix& rLocalStiffnessMatrix) {
//     KRATOS_TRY
//     for (unsigned int i = 0; i < mFibers.size(); ++i) {
//         rLocalStiffnessMatrix += mFibers[i].CreateGlobalFiberStiffnessMatrix();
//     }
//     KRATOS_CATCH("")
// }

Matrix FiberBeamColumnSection::GetGlobalFlexibilityMatrix()
{
    KRATOS_TRY

    // allocate memory
    Matrix aux_matrix         = ZeroMatrix(5, 3);
    Matrix global_flexibility = ZeroMatrix(5, 5);

    noalias(aux_matrix) += prod(Matrix(trans(mBMatrix)), mLocalFlexibilityMatrix);
    noalias(global_flexibility) += prod(aux_matrix, mBMatrix);
    // FIXME try .. you know

    return global_flexibility;
    KRATOS_CATCH("")
}

void FiberBeamColumnSection::CalculateLocalFlexibilityMatrix()
{
    KRATOS_TRY

    // allocate memory
    Matrix local_stiffness = ZeroMatrix(3, 3);
    // get local 3x3 stiffness from fibers' tangential stiffness
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        noalias(local_stiffness) += r_fiber.CreateGlobalFiberStiffnessMatrix();
    }
    // invert stiffness to get 3x3 flexibility
    double det_stiffness = MathUtils<double>::Det(local_stiffness);
    MathUtils<double>::InvertMatrix3(local_stiffness, mLocalFlexibilityMatrix, det_stiffness);

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

bool FiberBeamColumnSection::StateDetermination(Vector ChangeElementForceIncrements)
{
    KRATOS_TRY
    Vector chng_force_incr = ZeroVector(3);
    Vector chng_def_incr = ZeroVector(3);

    noalias(chng_force_incr) = prod(mBMatrix, ChangeElementForceIncrements);
    mForceIncrements += chng_force_incr;
    mForces = mConvergedForces + mForceIncrements;

    noalias(chng_def_incr) = mDeformationResiduals + prod(mLocalFlexibilityMatrix, chng_force_incr);
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        r_fiber.StateDetermination(chng_def_incr);
    }

    CalculateLocalFlexibilityMatrix();

    Vector resisting_forces = ZeroVector(3);
    for (FiberBeamColumnUniaxialFiber& r_fiber : mFibers) {
        noalias(resisting_forces) += r_fiber.CreateGlobalFiberResistingForces();
    }

    mUnbalanceForces = mForces - resisting_forces;
    noalias(mDeformationResiduals) = prod(mLocalFlexibilityMatrix, mUnbalanceForces);
    return std::abs(MathUtils<double>::Norm3(mUnbalanceForces)) < mTolerance;

    KRATOS_CATCH("")
}

void FiberBeamColumnSection::ResetResidual()
{
    KRATOS_TRY
    mDeformationResiduals = ZeroVector(3);
    KRATOS_CATCH("")
}


// void FiberBeamColumnSection::SetSectionForces(Vector ChangeElementForceIncrements)
// {
//     KRATOS_TRY

//     noalias(mChangeForceIncrements) = prod(mBMatrix, ChangeElementForceIncrements); //FIXME: not mem
//     mForceIncrements += mChangeForceIncrements;
//     mForces = mConvergedForces + mForceIncrements;

//     KRATOS_CATCH("")
// }

// void FiberBeamColumnSection::CalculateDeformationIncrements()
// {
//     KRATOS_TRY
//     mChangeDeformationIncrements = mDeformationResiduals + prod(mLocalFlexibilityMatrix, mChangeForceIncrements);
//     KRATOS_CATCH("")
// }

// void FiberBeamColumnSection::CalculateFiberState()
// {
//     KRATOS_TRY
//     for (unsigned int i = 0; i < mFibers.size(); ++i) {
//         mFibers[i].SetStrainIncrements(mChangeDeformationIncrements);
//     }
//     CalculateLocalFlexibilityMatrix(mLocalFlexibilityMatrix);
//     KRATOS_CATCH("")
// }

// void FiberBeamColumnSection::CalculateResiduals()
// {
//     KRATOS_TRY
//     Vector resisting_forces = ZeroVector(3);
//     for (unsigned int i = 0; i < mFibers.size(); ++i) {
//         resisting_forces += mFibers[i].CreateGlobalFiberResistingForces();
//     }
//     mUnbalanceForces = mForces - resisting_forces;
//     mDeformationResiduals = prod(mLocalFlexibilityMatrix, mUnbalanceForces);
//     KRATOS_CATCH("")
// }

Vector FiberBeamColumnSection::GetGlobalDeformationResiduals()
{
    KRATOS_TRY
    return prod(Matrix(trans(mBMatrix)), mDeformationResiduals);
    KRATOS_CATCH("")
}

// bool FiberBeamColumnSection::CheckConvergence()
// {
//     KRATOS_TRY
//     return MathUtils<double>::Norm(mUnbalanceForces) < mTolerance;
//     KRATOS_CATCH("")
// }

void FiberBeamColumnSection::FinalizeSolutionStep()
{
    KRATOS_TRY

    // std::fill(mDeformationResiduals.begin(), mDeformationResiduals.end(), 0.00);

    mConvergedForces = mForces;
    mForceIncrements = ZeroVector(3);

    for(unsigned int i = 0; i < mFibers.size(); ++i) {
        mFibers[i].FinalizeSolutionStep();
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
