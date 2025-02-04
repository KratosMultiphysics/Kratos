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

// System includes

// External includes

#include "custom_constitutive/small_strain_udsm_3D_interface_law.hpp"

namespace Kratos
{

ConstitutiveLaw::Pointer SmallStrainUDSM3DInterfaceLaw::Clone() const
{
    KRATOS_TRY

    return Kratos::make_shared<SmallStrainUDSM3DInterfaceLaw>(*this);

    KRATOS_CATCH("")
}

void SmallStrainUDSM3DInterfaceLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& rStrainVector = rValues.GetStrainVector();

    mDeltaStrainVector[INDEX_3D_ZZ] = rStrainVector(INDEX_3D_INTERFACE_ZZ) - mStrainVectorFinalized[INDEX_3D_ZZ];
    mDeltaStrainVector[INDEX_3D_YZ] = rStrainVector(INDEX_3D_INTERFACE_YZ) - mStrainVectorFinalized[INDEX_3D_YZ];
    mDeltaStrainVector[INDEX_3D_XZ] = rStrainVector(INDEX_3D_INTERFACE_XZ) - mStrainVectorFinalized[INDEX_3D_XZ];
}

void SmallStrainUDSM3DInterfaceLaw::SetExternalStressVector(Vector& rStressVector)
{
    rStressVector(INDEX_3D_INTERFACE_ZZ) = mStressVector[INDEX_3D_ZZ];
    rStressVector(INDEX_3D_INTERFACE_YZ) = mStressVector[INDEX_3D_YZ];
    rStressVector(INDEX_3D_INTERFACE_XZ) = mStressVector[INDEX_3D_XZ];
}

void SmallStrainUDSM3DInterfaceLaw::SetInternalStressVector(const Vector& rStressVector)
{
    KRATOS_TRY
    std::fill(mStressVectorFinalized.begin(), mStressVectorFinalized.end(), 0.0);

    mStressVectorFinalized[INDEX_3D_ZZ] = rStressVector(INDEX_3D_INTERFACE_ZZ);
    mStressVectorFinalized[INDEX_3D_YZ] = rStressVector(INDEX_3D_INTERFACE_YZ);
    mStressVectorFinalized[INDEX_3D_XZ] = rStressVector(INDEX_3D_INTERFACE_XZ);

    KRATOS_CATCH("")
}

void SmallStrainUDSM3DInterfaceLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
    std::fill(mStrainVectorFinalized.begin(), mStrainVectorFinalized.end(), 0.0);

    mStrainVectorFinalized[INDEX_3D_ZZ] = rStrainVector(INDEX_3D_INTERFACE_ZZ);
    mStrainVectorFinalized[INDEX_3D_YZ] = rStrainVector(INDEX_3D_INTERFACE_YZ);
    mStrainVectorFinalized[INDEX_3D_XZ] = rStrainVector(INDEX_3D_INTERFACE_XZ);
}

void SmallStrainUDSM3DInterfaceLaw::CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer fortran style matrix to C++ style
        for (unsigned int i = 0; i < VoigtSize; i++) {
            for (unsigned int j = 0; j < VoigtSize; j++) {
                rConstitutiveMatrix(i, j) =
                    mMatrixD[getIndex3D(static_cast<indexStress3DInterface>(j))]
                            [getIndex3D(static_cast<indexStress3DInterface>(i))];
            }
        }
    } else {
        for (unsigned int i = 0; i < VoigtSize; i++) {
            for (unsigned int j = 0; j < VoigtSize; j++) {
                rConstitutiveMatrix(i, j) =
                    mMatrixD[getIndex3D(static_cast<indexStress3DInterface>(i))]
                            [getIndex3D(static_cast<indexStress3DInterface>(j))];
            }
        }
    }
}

indexStress3D SmallStrainUDSM3DInterfaceLaw::getIndex3D(const indexStress3DInterface index3D) const
{
    switch (index3D) {
    case INDEX_3D_INTERFACE_ZZ:
        return INDEX_3D_ZZ;
    case INDEX_3D_INTERFACE_YZ:
        return INDEX_3D_YZ;
    case INDEX_3D_INTERFACE_XZ:
        return INDEX_3D_XZ;
    default:
        KRATOS_ERROR << "invalid index: " << index3D << std::endl;
    }
}

Vector& SmallStrainUDSM3DInterfaceLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == STATE_VARIABLES) {
        SmallStrainUDSM3DLaw::GetValue(rThisVariable, rValue);
    } else if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        if (rValue.size() != VoigtSize) rValue.resize(VoigtSize);

        rValue[INDEX_3D_INTERFACE_ZZ] = mStressVectorFinalized[INDEX_3D_ZZ];
        rValue[INDEX_3D_INTERFACE_YZ] = mStressVectorFinalized[INDEX_3D_YZ];
        rValue[INDEX_3D_INTERFACE_XZ] = mStressVectorFinalized[INDEX_3D_XZ];
    }
    return rValue;
}

void SmallStrainUDSM3DInterfaceLaw::SetValue(const Variable<Vector>& rThisVariable,
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