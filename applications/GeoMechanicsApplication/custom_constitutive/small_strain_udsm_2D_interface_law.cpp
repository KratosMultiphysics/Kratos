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

#include "custom_constitutive/small_strain_udsm_2D_interface_law.h"
#include "custom_constitutive/constitutive_law_dimension.h"

#include <type_traits>

namespace
{

using namespace Kratos;

indexStress3D GetIndex3D(const indexStress2DInterface Index2D)
{
    switch (Index2D) {
    case INDEX_2D_INTERFACE_ZZ:
        return INDEX_3D_ZZ;
    case INDEX_2D_INTERFACE_XZ:
        return INDEX_3D_XZ;
    default:
        KRATOS_ERROR << "invalid index: " << Index2D << std::endl;
    }
}

} // namespace

namespace Kratos
{

ConstitutiveLaw::Pointer SmallStrainUDSM2DInterfaceLaw::Clone() const
{
    auto pResult = std::make_shared<SmallStrainUDSM2DInterfaceLaw>();
    CloneDataMembersTo(*pResult);
    return pResult;
}

void SmallStrainUDSM2DInterfaceLaw::UpdateInternalDeltaStrainVector(Parameters& rValues)
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

void SmallStrainUDSM2DInterfaceLaw::CopyConstitutiveMatrix(Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer fortran style matrix to C++ style
        for (unsigned int i = 0; i < GetStrainSize(); i++) {
            for (unsigned int j = 0; j < GetStrainSize(); j++) {
                rConstitutiveMatrix(i, j) =
                    mMatrixD[GetIndex3D(static_cast<indexStress2DInterface>(j))]
                            [GetIndex3D(static_cast<indexStress2DInterface>(i))];
            }
        }
    } else {
        for (unsigned int i = 0; i < GetStrainSize(); i++) {
            for (unsigned int j = 0; j < GetStrainSize(); j++) {
                rConstitutiveMatrix(i, j) =
                    mMatrixD[GetIndex3D(static_cast<indexStress2DInterface>(i))]
                            [GetIndex3D(static_cast<indexStress2DInterface>(j))];
            }
        }
    }
}

void SmallStrainUDSM2DInterfaceLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainUDSMLaw)
}

void SmallStrainUDSM2DInterfaceLaw::load(Serializer& rSerializer){
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainUDSMLaw)}

Vector& SmallStrainUDSM2DInterfaceLaw::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == STATE_VARIABLES) {
        SmallStrainUDSMLaw::GetValue(rVariable, rValue);
    } else if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
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
        SmallStrainUDSMLaw::SetValue(rVariable, rValue, rCurrentProcessInfo);
    } else if ((rVariable == CAUCHY_STRESS_VECTOR || rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) &&
               rValue.size() == GetStrainSize()) {
        this->SetInternalStressVector(rValue);
    }
}

SizeType SmallStrainUDSM2DInterfaceLaw::WorkingSpaceDimension() { return N_DIM_2D; }

SizeType SmallStrainUDSM2DInterfaceLaw::GetStrainSize() const { return VOIGT_SIZE_2D_INTERFACE; }

std::string SmallStrainUDSM2DInterfaceLaw::Info() const { return "SmallStrainUDSM2DInterfaceLaw"; }

void SmallStrainUDSM2DInterfaceLaw::PrintData(std::ostream& rOStream) const
{
    rOStream << "SmallStrainUDSM2DInterfaceLaw Data";
}

// Instances of this class can neither be copied nor moved. Check that at compile time.
static_assert(!std::is_copy_constructible_v<SmallStrainUDSM2DInterfaceLaw>);
static_assert(!std::is_copy_assignable_v<SmallStrainUDSM2DInterfaceLaw>);
static_assert(!std::is_move_constructible_v<SmallStrainUDSM2DInterfaceLaw>);
static_assert(!std::is_move_assignable_v<SmallStrainUDSM2DInterfaceLaw>);

} // Namespace Kratos
