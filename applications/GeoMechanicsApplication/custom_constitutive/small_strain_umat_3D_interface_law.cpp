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

#include "custom_constitutive/small_strain_umat_3D_interface_law.h"
#include "constitutive_law_dimension.h"

namespace Kratos
{
using namespace std::string_literals;

SmallStrainUMAT3DInterfaceLaw::SmallStrainUMAT3DInterfaceLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : SmallStrainUMATLaw<VOIGT_SIZE_3D>(std::move(pConstitutiveDimension))
{
}

ConstitutiveLaw::Pointer SmallStrainUMAT3DInterfaceLaw::Clone() const
{
    KRATOS_TRY
    return Kratos::make_shared<SmallStrainUMAT3DInterfaceLaw>(*this);
    KRATOS_CATCH("")
}

void SmallStrainUMAT3DInterfaceLaw::UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& rStrainVector = rValues.GetStrainVector();

    mDeltaStrainVector[2] = rStrainVector(2) - mStrainVectorFinalized[2];
    mDeltaStrainVector[4] = rStrainVector(1) - mStrainVectorFinalized[4];
    mDeltaStrainVector[5] = rStrainVector(0) - mStrainVectorFinalized[5];
}

void SmallStrainUMAT3DInterfaceLaw::SetExternalStressVector(Vector& rStressVector)
{
    rStressVector[2] = mStressVector[2];
    rStressVector[1] = mStressVector[4];
    rStressVector[0] = mStressVector[5];
}

void SmallStrainUMAT3DInterfaceLaw::SetInternalStressVector(const Vector& rStressVector)
{
    KRATOS_TRY
    std::fill(mStressVectorFinalized.begin(), mStressVectorFinalized.end(), 0.0);

    mStressVectorFinalized[2] = rStressVector[2];
    mStressVectorFinalized[4] = rStressVector[1];
    mStressVectorFinalized[5] = rStressVector[0];
    KRATOS_CATCH("")
}

void SmallStrainUMAT3DInterfaceLaw::SetInternalStrainVector(const Vector& rStrainVector)
{
    std::fill(mStrainVectorFinalized.begin(), mStrainVectorFinalized.end(), 0.0);

    mStrainVectorFinalized[2] = rStrainVector(2);
    mStrainVectorFinalized[4] = rStrainVector(1);
    mStrainVectorFinalized[5] = rStrainVector(0);
}

void SmallStrainUMAT3DInterfaceLaw::CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    if (rValues.GetMaterialProperties()[IS_FORTRAN_UDSM]) {
        // transfer fortran style matrix to C++ style
        for (unsigned int i = 0; i < VoigtSize; i++) {
            for (unsigned int j = 0; j < VoigtSize; j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[getIndex3D(j)][getIndex3D(i)];
            }
        }
    } else {
        for (unsigned int i = 0; i < VoigtSize; i++) {
            for (unsigned int j = 0; j < VoigtSize; j++) {
                rConstitutiveMatrix(i, j) = mMatrixD[getIndex3D(i)][getIndex3D(j)];
            }
        }
    }
}

std::size_t SmallStrainUMAT3DInterfaceLaw::getIndex3D(std::size_t index3D)
{
    switch (index3D) {
    case 0:
        return 5;
    case 1:
        return 4;
    case 2:
        return 2;
    default:
        KRATOS_ERROR << "invalid index: " << index3D << std::endl;
    }
}

void SmallStrainUMAT3DInterfaceLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainUMATLaw)
}

void SmallStrainUMAT3DInterfaceLaw::load(Serializer& rSerializer){
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainUMATLaw)}

Vector& SmallStrainUMAT3DInterfaceLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == STATE_VARIABLES) {
        SmallStrainUMATLaw::GetValue(rThisVariable, rValue);
    } else if (rThisVariable == CAUCHY_STRESS_VECTOR || rThisVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
        if (rValue.size() != VoigtSize) rValue.resize(VoigtSize);

        rValue[2] = mStressVectorFinalized[2];
        rValue[1] = mStressVectorFinalized[4];
        rValue[0] = mStressVectorFinalized[5];
    }
    return rValue;
}

void SmallStrainUMAT3DInterfaceLaw::SetValue(const Variable<Vector>& rVariable,
                                             const Vector&           rValue,
                                             const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == STATE_VARIABLES) {
        SmallStrainUMATLaw::SetValue(rVariable, rValue, rCurrentProcessInfo);
    } else if ((rVariable == CAUCHY_STRESS_VECTOR || rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) &&
               rValue.size() == VoigtSize) {
        this->SetInternalStressVector(rValue);
    }
}

SizeType SmallStrainUMAT3DInterfaceLaw::WorkingSpaceDimension() { return Dimension; }

SizeType SmallStrainUMAT3DInterfaceLaw::GetStrainSize() const { return VoigtSize; }

ConstitutiveLaw::StrainMeasure SmallStrainUMAT3DInterfaceLaw::GetStrainMeasure()
{
    return StrainMeasure_Infinitesimal;
}

ConstitutiveLaw::StressMeasure SmallStrainUMAT3DInterfaceLaw::GetStressMeasure()
{
    return StressMeasure_Cauchy;
}

std::string SmallStrainUMAT3DInterfaceLaw::Info() const { return "SmallStrainUMAT3DInterfaceLaw"s; }

void SmallStrainUMAT3DInterfaceLaw::PrintInfo(std::ostream& rOStream) const { rOStream << Info(); }

void SmallStrainUMAT3DInterfaceLaw::PrintData(std::ostream& rOStream) const
{
    rOStream << "SmallStrainUMAT3DInterfaceLaw Data";
}
} // namespace Kratos