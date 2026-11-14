// License: BSD License
// Kratos default license: kratos/license.txt

#pragma once

#include "generic_small_strain_isotropic_plasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"
#include "custom_constitutive/auxiliary_files/hardening_laws/swift_hardening_law.h"

namespace Kratos
{

/** @brief Public Swift law backed by generic associated J2 plasticity.
 * @details History and material-response lifecycle belong to the generic law;
 * the generic integrator uses SwiftHardeningLaw for sigma_y(p) and its derivative.
 * This wrapper preserves the registered name and Python/C++ construction API.
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) SmallStrainSwiftJ2Plasticity3D
    : public GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<
        VonMisesYieldSurface<VonMisesPlasticPotential<6>>, SwiftHardeningLaw>>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainSwiftJ2Plasticity3D);
    using BaseType = GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<
        VonMisesYieldSurface<VonMisesPlasticPotential<6>>, SwiftHardeningLaw>>;
    using BaseType::Has;
    using BaseType::GetValue;
    using BaseType::SetValue;

    SmallStrainSwiftJ2Plasticity3D() = default;
    SmallStrainSwiftJ2Plasticity3D(const SmallStrainSwiftJ2Plasticity3D&) = default;
    ~SmallStrainSwiftJ2Plasticity3D() override = default;

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<SmallStrainSwiftJ2Plasticity3D>(*this);
    }

    std::string Info() const override { return "SmallStrainSwiftJ2Plasticity3DLaw"; }

    /// The legacy state vector omits p; use the two explicit history variables.
    bool Has(const Variable<Vector>& rVariable) override
    {
        return rVariable != INTERNAL_VARIABLES && BaseType::Has(rVariable);
    }

    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override
    {
        KRATOS_ERROR_IF(rVariable == INTERNAL_VARIABLES)
            << "Swift does not support INTERNAL_VARIABLES; use PLASTIC_STRAIN_VECTOR and ACCUMULATED_PLASTIC_STRAIN";
        return BaseType::GetValue(rVariable, rValue);
    }

    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue,
        const ProcessInfo& rProcessInfo) override
    {
        KRATOS_ERROR_IF(rVariable == INTERNAL_VARIABLES)
            << "Swift does not support INTERNAL_VARIABLES; use PLASTIC_STRAIN_VECTOR and ACCUMULATED_PLASTIC_STRAIN";
        BaseType::SetValue(rVariable, rValue, rProcessInfo);
    }

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }
};

} // namespace Kratos
