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

#include "custom_constitutive/small_strain_udsm_3D_interface_law.hpp"

#include <type_traits>

namespace Kratos
{

ConstitutiveLaw::Pointer SmallStrainUDSM3DInterfaceLaw::Clone() const
{
    auto pResult = std::make_shared<SmallStrainUDSM3DInterfaceLaw>();
    CloneDataMembersTo(*pResult);
    return pResult;
}

void SmallStrainUDSM3DInterfaceLaw::UpdateInternalDeltaStrainVector(Parameters& rValues)
{
    const auto& r_strain_vector = rValues.GetStrainVector();

    mDeltaStrainVector[INDEX_3D_ZZ] =
        r_strain_vector[INDEX_3D_INTERFACE_ZZ] - mStrainVectorFinalized[INDEX_3D_ZZ];
    mDeltaStrainVector[INDEX_3D_YZ] =
        r_strain_vector[INDEX_3D_INTERFACE_YZ] - mStrainVectorFinalized[INDEX_3D_YZ];
    mDeltaStrainVector[INDEX_3D_XZ] =
        r_strain_vector[INDEX_3D_INTERFACE_XZ] - mStrainVectorFinalized[INDEX_3D_XZ];
}

void SmallStrainUDSM3DInterfaceLaw::SetExternalStressVector(Vector& rStressVector)
{
    rStressVector[INDEX_3D_INTERFACE_ZZ] = mStressVector[INDEX_3D_ZZ];
    rStressVector[INDEX_3D_INTERFACE_YZ] = mStressVector[INDEX_3D_YZ];
    rStressVector[INDEX_3D_INTERFACE_XZ] = mStressVector[INDEX_3D_XZ];
}

void SmallStrainUDSM3DInterfaceLaw::SetInternalStressVector(const Vector& rStressVector)
{
    auto& r_sig0 = GetSig0();

    std::fill_n(r_sig0.begin(), StressVectorSize, 0.0);

    r_sig0[INDEX_3D_ZZ] = rStressVector[INDEX_3D_INTERFACE_ZZ];
    r_sig0[INDEX_3D_YZ] = rStressVector[INDEX_3D_INTERFACE_YZ];
    r_sig0[INDEX_3D_XZ] = rStressVector[INDEX_3D_INTERFACE_XZ];
}

void SmallStrainUDSM3DInterfaceLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
    std::fill(mStrainVectorFinalized.begin(), mStrainVectorFinalized.end(), 0.0);

    mStrainVectorFinalized[INDEX_3D_ZZ] = rStrainVector[INDEX_3D_INTERFACE_ZZ];
    mStrainVectorFinalized[INDEX_3D_YZ] = rStrainVector[INDEX_3D_INTERFACE_YZ];
    mStrainVectorFinalized[INDEX_3D_XZ] = rStrainVector[INDEX_3D_INTERFACE_XZ];
}

void SmallStrainUDSM3DInterfaceLaw::CopyConstitutiveMatrix(Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer fortran style matrix to C++ style
        for (unsigned int i = 0; i < GetStrainSize(); i++) {
            for (unsigned int j = 0; j < GetStrainSize(); j++) {
                rConstitutiveMatrix(i, j) =
                    mMatrixD[getIndex3D(static_cast<indexStress3DInterface>(j))]
                            [getIndex3D(static_cast<indexStress3DInterface>(i))];
            }
        }
    } else {
        for (unsigned int i = 0; i < GetStrainSize(); i++) {
            for (unsigned int j = 0; j < GetStrainSize(); j++) {
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

void SmallStrainUDSM3DInterfaceLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainUDSMLaw)
}

void SmallStrainUDSM3DInterfaceLaw::load(Serializer& rSerializer){
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainUDSMLaw)}

Vector& SmallStrainUDSM3DInterfaceLaw::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == STATE_VARIABLES) {
        SmallStrainUDSMLaw::GetValue(rVariable, rValue);
    } else if (rVariable == CAUCHY_STRESS_VECTOR) {
        rValue.resize(GetStrainSize());

        auto& r_sig0                  = GetSig0();
        rValue[INDEX_3D_INTERFACE_ZZ] = r_sig0[INDEX_3D_ZZ];
        rValue[INDEX_3D_INTERFACE_YZ] = r_sig0[INDEX_3D_YZ];
        rValue[INDEX_3D_INTERFACE_XZ] = r_sig0[INDEX_3D_XZ];
    }
    return rValue;
}

void SmallStrainUDSM3DInterfaceLaw::SetValue(const Variable<Vector>& rVariable,
                                             const Vector&           rValue,
                                             const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == STATE_VARIABLES) {
        SmallStrainUDSMLaw::SetValue(rVariable, rValue, rCurrentProcessInfo);
    } else if ((rVariable == CAUCHY_STRESS_VECTOR) && (rValue.size() == GetStrainSize())) {
        this->SetInternalStressVector(rValue);
    }
}

SizeType SmallStrainUDSM3DInterfaceLaw::WorkingSpaceDimension() { return N_DIM_3D; }

SizeType SmallStrainUDSM3DInterfaceLaw::GetStrainSize() const { return VOIGT_SIZE_3D_INTERFACE; }

std::string SmallStrainUDSM3DInterfaceLaw::Info() const { return "SmallStrainUDSM3DInterfaceLaw"; }

void SmallStrainUDSM3DInterfaceLaw::PrintInfo(std::ostream& rOStream) const { rOStream << Info(); }

void SmallStrainUDSM3DInterfaceLaw::PrintData(std::ostream& rOStream) const
{
    rOStream << "SmallStrainUDSM3DInterfaceLaw Data";
}

// Instances of this class cannot be copied, but they can be moved. Check that at compile time.
static_assert(!std::is_copy_constructible_v<SmallStrainUDSM3DInterfaceLaw>);
static_assert(!std::is_copy_assignable_v<SmallStrainUDSM3DInterfaceLaw>);
static_assert(std::is_move_constructible_v<SmallStrainUDSM3DInterfaceLaw>);
static_assert(std::is_move_assignable_v<SmallStrainUDSM3DInterfaceLaw>);

} // namespace Kratos