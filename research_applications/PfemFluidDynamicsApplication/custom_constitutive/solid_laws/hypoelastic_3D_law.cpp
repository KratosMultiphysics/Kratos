//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Alessandro Franci
//  Collaborators:
//
//-------------------------------------------------------------
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/solid_laws/hypoelastic_3D_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    Hypoelastic3DLaw::Hypoelastic3DLaw() : PfemSolidConstitutiveLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    Hypoelastic3DLaw::Hypoelastic3DLaw(const Hypoelastic3DLaw &rOther) : PfemSolidConstitutiveLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer Hypoelastic3DLaw::Clone() const { return Kratos::make_shared<Hypoelastic3DLaw>(*this); }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    Hypoelastic3DLaw::~Hypoelastic3DLaw() {}

    ConstitutiveLaw::SizeType Hypoelastic3DLaw::WorkingSpaceDimension() { return 3; }

    ConstitutiveLaw::SizeType Hypoelastic3DLaw::GetStrainSize() const { return 6; }

    void Hypoelastic3DLaw::CalculateMaterialResponseCauchy(Parameters &rParameters)
    {

        const Vector &r_strain_vector = rParameters.GetStrainVector();
        Vector &r_stress_vector = rParameters.GetStressVector();

        const double young_modulus = this->GetEffectiveMaterialParameter(rParameters, YOUNG_MODULUS);
        const double poisson_ratio = this->GetEffectiveMaterialParameter(rParameters, POISSON_RATIO);
        const double time_step = rParameters.GetProcessInfo()[DELTA_TIME];

        const double second_lame = time_step * young_modulus / (2.0 * (1.0 + poisson_ratio));

        const double strain_trace = r_strain_vector[0] + r_strain_vector[1] + r_strain_vector[2];

        r_stress_vector[0] += 2.0 * second_lame * (r_strain_vector[0] - strain_trace / 3.0);
        r_stress_vector[1] += 2.0 * second_lame * (r_strain_vector[1] - strain_trace / 3.0);
        r_stress_vector[2] += 2.0 * second_lame * (r_strain_vector[2] - strain_trace / 3.0);
        r_stress_vector[3] += 2.0 * second_lame * r_strain_vector[3];
        r_stress_vector[4] += 2.0 * second_lame * r_strain_vector[4];
        r_stress_vector[5] += 2.0 * second_lame * r_strain_vector[5];
    }

    int Hypoelastic3DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                const ProcessInfo &rCurrentProcessInfo) const
    {

        KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] < 0.0)
            << "Incorrect or missing YOUNG_MODULUS provided in process info for Hypoelastic3DLaw: "
            << rMaterialProperties[YOUNG_MODULUS] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[POISSON_RATIO] < 0.0)
            << "Incorrect or missing POISSON_RATIO provided in process info for Hypoelastic3DLaw: "
            << rMaterialProperties[POISSON_RATIO] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0)
            << "Incorrect or missing DENSITY provided in process info for Hypoelastic3DLaw: "
            << rMaterialProperties[DENSITY] << std::endl;

        return 0;
    }

    std::string Hypoelastic3DLaw::Info() const { return "Hypoelastic3DLaw"; }

    double Hypoelastic3DLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
    {
        return rParameters.GetMaterialProperties()[rVariable];
    }

    void Hypoelastic3DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PfemSolidConstitutiveLaw)
    }

    void Hypoelastic3DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PfemSolidConstitutiveLaw)
    }

} // Namespace Kratos
