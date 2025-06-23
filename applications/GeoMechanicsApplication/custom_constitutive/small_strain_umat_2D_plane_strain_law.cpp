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

#include "custom_constitutive/small_strain_umat_2D_plane_strain_law.hpp"

namespace Kratos
{

ConstitutiveLaw::Pointer SmallStrainUMAT2DPlaneStrainLaw::Clone() const
{
    KRATOS_TRY

    return Kratos::make_shared<SmallStrainUMAT2DPlaneStrainLaw>(*this);

    KRATOS_CATCH("")
}

void SmallStrainUMAT2DPlaneStrainLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& rStrainVector = rValues.GetStrainVector();

    for (unsigned int i = 0; i < VoigtSize; ++i) {
        mDeltaStrainVector[i] = rStrainVector(i) - mStrainVectorFinalized[i];
    }
}

void SmallStrainUMAT2DPlaneStrainLaw::SetExternalStressVector(Vector& rStressVector)
{
    KRATOS_TRY
    std::copy_n(mStressVector.begin(), VoigtSize, rStressVector.begin());
    KRATOS_CATCH("")
}

void SmallStrainUMAT2DPlaneStrainLaw::SetInternalStressVector(const Vector& rStressVector)
{
    KRATOS_TRY
    std::copy_n(rStressVector.begin(), VoigtSize, mStressVectorFinalized.begin());
    KRATOS_CATCH("")
}

void SmallStrainUMAT2DPlaneStrainLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
    KRATOS_TRY
    std::copy_n(rStrainVector.begin(), VoigtSize, mStrainVectorFinalized.begin());
    KRATOS_CATCH("")
}

void SmallStrainUMAT2DPlaneStrainLaw::CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues,
                                                             Matrix& rConstitutiveMatrix)
{
    KRATOS_TRY

    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer fortran style matrix to C++ style
        for (unsigned int i = 0; i < VoigtSize; i++) {
            for (unsigned int j = 0; j < VoigtSize; j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[j][i];
            }
        }
    } else {
        for (unsigned int i = 0; i < VoigtSize; i++) {
            for (unsigned int j = 0; j < VoigtSize; j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[i][j];
            }
        }
    }

    KRATOS_CATCH("")
}

Vector& SmallStrainUMAT2DPlaneStrainLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == STATE_VARIABLES) {
        SmallStrainUMAT3DLaw::GetValue(rThisVariable, rValue);
    } else if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        if (rValue.size() != VoigtSize) rValue.resize(VoigtSize);
        for (unsigned int i = 0; i < VoigtSize; ++i) {
            rValue[i] = mStressVectorFinalized[i];
        }
    }
    return rValue;
}

void SmallStrainUMAT2DPlaneStrainLaw::SetValue(const Variable<Vector>& rVariable,
                                               const Vector&           rValue,
                                               const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == STATE_VARIABLES) {
        SmallStrainUMAT3DLaw::SetValue(rVariable, rValue, rCurrentProcessInfo);
    } else if ((rVariable == CAUCHY_STRESS_VECTOR) && (rValue.size() == VoigtSize)) {
        this->SetInternalStressVector(rValue);
    }
}

} // namespace Kratos