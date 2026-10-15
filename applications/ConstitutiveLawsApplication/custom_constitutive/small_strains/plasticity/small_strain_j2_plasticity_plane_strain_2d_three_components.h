// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <array>
#include "small_strain_j2_plasticity_plane_strain_2d.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{
/** @brief Plane-strain J2 plasticity with the external ordering [xx, yy, xy].
 * The inherited integrator retains [xx, yy, zz, xy] plastic history internally.
 * Strains must be supplied by the element, with engineering shear strain.
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) SmallStrainJ2PlasticityPlaneStrain2DThreeComponents
    : public SmallStrainJ2PlasticityPlaneStrain2D
{
public:
    using BaseType = SmallStrainJ2PlasticityPlaneStrain2D;
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainJ2PlasticityPlaneStrain2DThreeComponents);
    using BaseType::Has;
    using BaseType::GetValue;
    using BaseType::CalculateValue;

    SmallStrainJ2PlasticityPlaneStrain2DThreeComponents()
    {
        mPlasticStrain = ZeroVector(4);
        mAccumulatedPlasticStrain = 0.0;
    }

    SmallStrainJ2PlasticityPlaneStrain2DThreeComponents(
        const SmallStrainJ2PlasticityPlaneStrain2DThreeComponents& rOther)
        : BaseType(rOther)
    {
        // The inherited copy constructor does not copy the plastic history.
        mPlasticStrain = rOther.mPlasticStrain;
        mAccumulatedPlasticStrain = rOther.mAccumulatedPlasticStrain;
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<SmallStrainJ2PlasticityPlaneStrain2DThreeComponents>(*this);
    }

    SizeType GetStrainSize() const override { return 3; }

    void GetLawFeatures(Features& rFeatures) override
    {
        BaseType::GetLawFeatures(rFeatures);
        rFeatures.mStrainSize = 3;
    }

    void InitializeMaterial(const Properties&, const GeometryType&, const Vector&) override
    {
        mPlasticStrain = ZeroVector(4);
        mAccumulatedPlasticStrain = 0.0;
    }

    bool Has(const Variable<double>& rVariable) override
    {
        return rVariable == EQUIVALENT_PLASTIC_STRAIN || BaseType::Has(rVariable);
    }

    double& GetValue(const Variable<double>& rVariable, double& rValue) override
    {
        if (rVariable == EQUIVALENT_PLASTIC_STRAIN) {
            rValue = mAccumulatedPlasticStrain;
            return rValue;
        }
        return BaseType::GetValue(rVariable, rValue);
    }

    double& CalculateValue(Parameters& rValues, const Variable<double>& rVariable, double& rValue) override
    {
        if (rVariable == EQUIVALENT_PLASTIC_STRAIN || rVariable == ACCUMULATED_PLASTIC_STRAIN) {
            Vector plastic_strain;
            double accumulated_plastic_strain = 0.0;
            CalculateStressResponse(rValues, plastic_strain, accumulated_plastic_strain);
            rValue = accumulated_plastic_strain;
        } else if (rVariable == STRAIN_ENERGY) {
            Vector strain = ExpandStrain(rValues.GetStrainVector());
            ProcessInfo process_info = ExpandedProcessInfo(rValues.GetProcessInfo());
            Parameters values = rValues;
            values.SetStrainVector(strain);
            values.SetProcessInfo(process_info);
            BaseType::CalculateValue(values, rVariable, rValue);
        } else {
            BaseType::CalculateValue(rValues, rVariable, rValue);
        }
        return rValue;
    }

    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Plane-strain J2 plasticity, three external strain components";
    }

protected:
    void CalculateStressResponse(Parameters& rValues, Vector& rPlasticStrain,
                                 double& rAccumulatedPlasticStrain) override
    {
        KRATOS_TRY
        KRATOS_ERROR_IF_NOT(rValues.GetOptions().Is(USE_ELEMENT_PROVIDED_STRAIN))
            << "This law requires element-provided infinitesimal strains." << std::endl;
        Vector strain = ExpandStrain(rValues.GetStrainVector());
        Vector stress = ZeroVector(4);
        Matrix tangent = ZeroMatrix(4, 4);
        ProcessInfo process_info = ExpandedProcessInfo(rValues.GetProcessInfo());
        Parameters values = rValues;
        values.SetStrainVector(strain);
        values.SetStressVector(stress);
        values.SetConstitutiveMatrix(tangent);
        values.SetProcessInfo(process_info);
        BaseType::CalculateStressResponse(values, rPlasticStrain, rAccumulatedPlasticStrain);

        constexpr std::array<SizeType, 3> indices{0, 1, 3};
        if (rValues.GetOptions().Is(COMPUTE_STRESS)) {
            auto& r_stress = rValues.GetStressVector();
            r_stress.resize(3, false);
            for (SizeType i = 0; i < 3; ++i) {
                r_stress[i] = stress[indices[i]];
            }
        }
        if (rValues.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
            auto& r_tangent = rValues.GetConstitutiveMatrix();
            r_tangent.resize(3, 3, false);
            for (SizeType i = 0; i < 3; ++i) {
                for (SizeType j = 0; j < 3; ++j) {
                    r_tangent(i, j) = tangent(indices[i], indices[j]);
                }
            }
        }
        KRATOS_CATCH("")
    }

private:
    static Vector ExpandStrain(const Vector& rStrain)
    {
        KRATOS_ERROR_IF(rStrain.size() != 3)
            << "Expected [epsilon_xx, epsilon_yy, gamma_xy]." << std::endl;
        Vector strain = ZeroVector(4);
        strain[0] = rStrain[0];
        strain[1] = rStrain[1];
        strain[3] = rStrain[2];
        return strain;
    }

    static ProcessInfo ExpandedProcessInfo(const ProcessInfo& rProcessInfo)
    {
        ProcessInfo process_info = rProcessInfo;
        if (process_info.Has(INITIAL_STRAIN)) {
            process_info[INITIAL_STRAIN] = ExpandStrain(rProcessInfo[INITIAL_STRAIN]);
        }
        return process_info;
    }

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
