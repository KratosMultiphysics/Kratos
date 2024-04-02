// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// External includes
#include "custom_constitutive/small_strain_udsm_2D_plane_strain_law.hpp"

namespace Kratos
{

ConstitutiveLaw::Pointer SmallStrainUDSM2DPlaneStrainLaw::Clone() const
{
    KRATOS_TRY

    return Kratos::make_shared<SmallStrainUDSM2DPlaneStrainLaw>(*this);

    KRATOS_CATCH("")
}

void SmallStrainUDSM2DPlaneStrainLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& rStrainVector = rValues.GetStrainVector();

    for (unsigned int i = 0; i < VoigtSize; ++i) {
        mDeltaStrainVector[i] = rStrainVector(i) - mStrainVectorFinalized[i];
    }
}

void SmallStrainUDSM2DPlaneStrainLaw::SetExternalStressVector(Vector& rStressVector)
{
    KRATOS_TRY
    for (unsigned int i = 0; i < VoigtSize; ++i) {
        rStressVector(i) = mStressVector[i];
    }
    KRATOS_CATCH("")
}

void SmallStrainUDSM2DPlaneStrainLaw::SetInternalStressVector(const Vector& rStressVector)
{
    KRATOS_TRY
    for (unsigned int i = 0; i < VoigtSize; ++i) {
        mStressVectorFinalized[i] = rStressVector(i);
    }
    KRATOS_CATCH("")
}

void SmallStrainUDSM2DPlaneStrainLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
    KRATOS_TRY
    for (unsigned int i = 0; i < VoigtSize; ++i) {
        mStrainVectorFinalized[i] = rStrainVector(i);
    }
    KRATOS_CATCH("")
}

void SmallStrainUDSM2DPlaneStrainLaw::CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues,
                                                             Matrix& rConstitutiveMatrix)
{
    KRATOS_TRY

    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer Fortran style matrix to C++ style
        for (unsigned int i = 0; i < VoigtSize; ++i) {
            for (unsigned int j = 0; j < VoigtSize; ++j) {
                rConstitutiveMatrix(i, j) = mMatrixD[j][i];
            }
        }
    } else {
        for (unsigned int i = 0; i < VoigtSize; ++i) {
            for (unsigned int j = 0; j < VoigtSize; ++j) {
                rConstitutiveMatrix(i, j) = mMatrixD[i][j];
            }
        }
    }

    KRATOS_CATCH("")
}

void SmallStrainUDSM2DPlaneStrainLaw::CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues,
                                                                 Vector& rStrainVector)
{
    // 1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // for shells/membranes in case the DeformationGradient is of size 3x3
    BoundedMatrix<double, 2, 2> F2x2;
    for (unsigned int i = 0; i < 2; ++i)
        for (unsigned int j = 0; j < 2; ++j)
            F2x2(i, j) = F(i, j);

    Matrix E_tensor = prod(trans(F2x2), F2x2);

    for (unsigned int i = 0; i < 2; ++i)
        E_tensor(i, i) -= 1.0;

    E_tensor *= 0.5;
    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

Vector& SmallStrainUDSM2DPlaneStrainLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == STATE_VARIABLES) {
        SmallStrainUDSM3DLaw::GetValue(rThisVariable, rValue);
    } else if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        if (rValue.size() != VoigtSize) rValue.resize(VoigtSize);
        for (unsigned int i = 0; i < VoigtSize; ++i) {
            rValue[i] = mStressVectorFinalized[i];
        }
    }
    return rValue;
}

void SmallStrainUDSM2DPlaneStrainLaw::SetValue(const Variable<Vector>& rThisVariable,
                                               const Vector&           rValue,
                                               const ProcessInfo&      rCurrentProcessInfo)
{
    if (rThisVariable == STATE_VARIABLES) {
        SmallStrainUDSM3DLaw::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    } else if ((rThisVariable == CAUCHY_STRESS_VECTOR) && (rValue.size() == VoigtSize)) {
        this->SetInternalStressVector(rValue);
    }
}

} // namespace Kratos