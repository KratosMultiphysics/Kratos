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

#include "custom_constitutive/small_strain_udsm_2D_interface_law.hpp"

namespace Kratos
{

ConstitutiveLaw::Pointer SmallStrainUDSM2DInterfaceLaw::Clone() const
{
    KRATOS_TRY

    return Kratos::make_shared<SmallStrainUDSM2DInterfaceLaw>(*this);

    KRATOS_CATCH("")
}

void SmallStrainUDSM2DInterfaceLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues)
{
    const auto& r_strain_vector = rValues.GetStrainVector();

    mDeltaStrainVector[INDEX_3D_ZZ] =
        r_strain_vector[INDEX_2D_INTERFACE_ZZ] - mStrainVectorFinalized[INDEX_3D_ZZ];
    mDeltaStrainVector[INDEX_3D_XZ] =
        r_strain_vector[INDEX_2D_INTERFACE_XZ] - mStrainVectorFinalized[INDEX_3D_XZ];
}

void SmallStrainUDSM2DInterfaceLaw::SetExternalStressVector(Vector& rStressVector)
{
    rStressVector(INDEX_2D_INTERFACE_ZZ) = mStressVector[INDEX_3D_ZZ];
    rStressVector(INDEX_2D_INTERFACE_XZ) = mStressVector[INDEX_3D_XZ];
}

void SmallStrainUDSM2DInterfaceLaw::SetInternalStressVector(const Vector& rStressVector)
{
    auto& r_sig0 = GetSig0();

    std::fill_n(r_sig0.begin(), StressVectorSize, 0.0);

    r_sig0[INDEX_3D_ZZ] = rStressVector[INDEX_2D_INTERFACE_ZZ];
    r_sig0[INDEX_3D_XZ] = rStressVector[INDEX_2D_INTERFACE_XZ];
}

void SmallStrainUDSM2DInterfaceLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
    std::fill(mStrainVectorFinalized.begin(), mStrainVectorFinalized.end(), 0.0);

    mStrainVectorFinalized[INDEX_3D_ZZ] = rStrainVector[INDEX_2D_INTERFACE_ZZ];
    mStrainVectorFinalized[INDEX_3D_XZ] = rStrainVector[INDEX_2D_INTERFACE_XZ];
}

void SmallStrainUDSM2DInterfaceLaw::CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer fortran style matrix to C++ style
        for (unsigned int i = 0; i < GetStrainSize(); i++) {
            for (unsigned int j = 0; j < GetStrainSize(); j++) {
                rConstitutiveMatrix(i, j) =
                    mMatrixD[getIndex3D(static_cast<indexStress2DInterface>(j))]
                            [getIndex3D(static_cast<indexStress2DInterface>(i))];
            }
        }
    } else {
        for (unsigned int i = 0; i < GetStrainSize(); i++) {
            for (unsigned int j = 0; j < GetStrainSize(); j++) {
                rConstitutiveMatrix(i, j) =
                    mMatrixD[getIndex3D(static_cast<indexStress2DInterface>(i))]
                            [getIndex3D(static_cast<indexStress2DInterface>(j))];
            }
        }
    }
}

indexStress3D SmallStrainUDSM2DInterfaceLaw::getIndex3D(const indexStress2DInterface index2D) const
{
    switch (index2D) {
    case INDEX_2D_INTERFACE_ZZ:
        return INDEX_3D_ZZ;
    case INDEX_2D_INTERFACE_XZ:
        return INDEX_3D_XZ;
    default:
        KRATOS_ERROR << "invalid index: " << index2D << std::endl;
    }
}

Vector& SmallStrainUDSM2DInterfaceLaw::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == STATE_VARIABLES) {
        SmallStrainUDSM3DLaw::GetValue(rVariable, rValue);
    } else if (rVariable == CAUCHY_STRESS_VECTOR) {
        rValue.resize(GetStrainSize());

        auto& r_sig0                  = GetSig0();
        rValue[INDEX_2D_INTERFACE_ZZ] = r_sig0[INDEX_3D_ZZ];
        rValue[INDEX_2D_INTERFACE_XZ] = r_sig0[INDEX_3D_XZ];
    }
    return rValue;
}

void SmallStrainUDSM2DInterfaceLaw::SetValue(const Variable<Vector>& rVariable,
                                             const Vector&           rValue,
                                             const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == STATE_VARIABLES) {
        SmallStrainUDSM3DLaw::SetValue(rVariable, rValue, rCurrentProcessInfo);
    } else if ((rVariable == CAUCHY_STRESS_VECTOR) && (rValue.size() == GetStrainSize())) {
        this->SetInternalStressVector(rValue);
    }
}

SizeType SmallStrainUDSM2DInterfaceLaw::GetStrainSize() const { return VOIGT_SIZE_2D_INTERFACE; }

} // Namespace Kratos
