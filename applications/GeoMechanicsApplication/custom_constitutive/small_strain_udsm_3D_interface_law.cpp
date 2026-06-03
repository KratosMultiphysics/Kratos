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

#include "custom_constitutive/small_strain_udsm_3D_interface_law.h"
#include "custom_constitutive/constitutive_law_dimension.h"

#include <type_traits>

namespace
{

using namespace Kratos;
using namespace std::string_literals;

std::size_t GetIndex3D(const std::size_t Index3D)
{
    switch (Index3D) {
    case 2:
        return 2;
    case 1:
        return 4;
    case 0:
        return 5;
    default:
        KRATOS_ERROR << "invalid index: " << Index3D << std::endl;
    }
}

} // namespace

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

    mDeltaStrainVector[2] = r_strain_vector[2] - mStrainVectorFinalized[2];
    mDeltaStrainVector[4] = r_strain_vector[1] - mStrainVectorFinalized[4];
    mDeltaStrainVector[5] = r_strain_vector[0] - mStrainVectorFinalized[5];
}

void SmallStrainUDSM3DInterfaceLaw::SetExternalStressVector(Vector& rStressVector)
{
    rStressVector[2] = mStressVector[2];
    rStressVector[1] = mStressVector[4];
    rStressVector[0] = mStressVector[5];
}

void SmallStrainUDSM3DInterfaceLaw::SetInternalStressVector(const Vector& rStressVector)
{
    auto& r_sig0 = GetSig0();

    std::fill_n(r_sig0.begin(), StressVectorSize, 0.0);

    r_sig0[2] = rStressVector[2];
    r_sig0[4] = rStressVector[1];
    r_sig0[5] = rStressVector[0];
}

void SmallStrainUDSM3DInterfaceLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
    std::fill(mStrainVectorFinalized.begin(), mStrainVectorFinalized.end(), 0.0);

    mStrainVectorFinalized[2] = rStrainVector[2];
    mStrainVectorFinalized[4] = rStrainVector[1];
    mStrainVectorFinalized[5] = rStrainVector[0];
}

void SmallStrainUDSM3DInterfaceLaw::CopyConstitutiveMatrix(Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer fortran style matrix to C++ style
        for (unsigned int i = 0; i < GetStrainSize(); i++) {
            for (unsigned int j = 0; j < GetStrainSize(); j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[GetIndex3D(j)][GetIndex3D(i)];
            }
        }
    } else {
        for (unsigned int i = 0; i < GetStrainSize(); i++) {
            for (unsigned int j = 0; j < GetStrainSize(); j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[GetIndex3D(i)][GetIndex3D(j)];
            }
        }
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
    } else if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
        rValue.resize(GetStrainSize());

        auto& r_sig0 = GetSig0();
        rValue[2]    = r_sig0[2];
        rValue[1]    = r_sig0[4];
        rValue[0]    = r_sig0[5];
    }
    return rValue;
}

void SmallStrainUDSM3DInterfaceLaw::SetValue(const Variable<Vector>& rVariable,
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

SizeType SmallStrainUDSM3DInterfaceLaw::WorkingSpaceDimension() { return N_DIM_3D; }

SizeType SmallStrainUDSM3DInterfaceLaw::GetStrainSize() const { return VOIGT_SIZE_3D_INTERFACE; }

std::string SmallStrainUDSM3DInterfaceLaw::Info() const { return "SmallStrainUDSM3DInterfaceLaw"s; }

void SmallStrainUDSM3DInterfaceLaw::PrintData(std::ostream& rOStream) const
{
    rOStream << "SmallStrainUDSM3DInterfaceLaw Data";
}

// Instances of this class can neither be copied nor moved. Check that at compile time.
static_assert(!std::is_copy_constructible_v<SmallStrainUDSM3DInterfaceLaw>);
static_assert(!std::is_copy_assignable_v<SmallStrainUDSM3DInterfaceLaw>);
static_assert(!std::is_move_constructible_v<SmallStrainUDSM3DInterfaceLaw>);
static_assert(!std::is_move_assignable_v<SmallStrainUDSM3DInterfaceLaw>);

} // namespace Kratos