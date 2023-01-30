//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Ruben Zorilla
//  Collaborators:  Massimiliano Zecchetto
//
//-------------------------------------------------------------
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/fluid_laws/newtonian_3D_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    Newtonian3DLaw::Newtonian3DLaw() : PfemFluidConstitutiveLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    Newtonian3DLaw::Newtonian3DLaw(const Newtonian3DLaw &rOther) : PfemFluidConstitutiveLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer Newtonian3DLaw::Clone() const { return Kratos::make_shared<Newtonian3DLaw>(*this); }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    Newtonian3DLaw::~Newtonian3DLaw() {}

    ConstitutiveLaw::SizeType Newtonian3DLaw::WorkingSpaceDimension() { return 3; }

    ConstitutiveLaw::SizeType Newtonian3DLaw::GetStrainSize() const { return 6; }

    void Newtonian3DLaw::CalculateMaterialResponseCauchy(Parameters &rValues)
    {

        const Flags &r_options = rValues.GetOptions();
        const Vector &r_strain_vector = rValues.GetStrainVector();
        Vector &r_stress_vector = rValues.GetStressVector();

        const double effective_dynamic_viscosity = rValues.GetMaterialProperties()[DYNAMIC_VISCOSITY];

        const double strain_trace = r_strain_vector[0] + r_strain_vector[1] + r_strain_vector[2];

        r_stress_vector[0] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[0] - strain_trace / 3.0);
        r_stress_vector[1] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[1] - strain_trace / 3.0);
        r_stress_vector[2] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[2] - strain_trace / 3.0);
        r_stress_vector[3] = 2.0 * effective_dynamic_viscosity * r_strain_vector[3];
        r_stress_vector[4] = 2.0 * effective_dynamic_viscosity * r_strain_vector[4];
        r_stress_vector[5] = 2.0 * effective_dynamic_viscosity * r_strain_vector[5];

        if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
        {
            this->EffectiveViscousConstitutiveMatrix3D(effective_dynamic_viscosity, rValues.GetConstitutiveMatrix());
        }
    }

    int Newtonian3DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                              const ProcessInfo &rCurrentProcessInfo) const
    {

        KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] < 0.0)
            << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for Newtonian3DLaw: "
            << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[BULK_MODULUS] < 0.0)
            << "Incorrect or missing BULK_MODULUS provided in process info for Newtonian3DLaw: "
            << rMaterialProperties[BULK_MODULUS] << std::endl;

        return 0;
    }

    std::string Newtonian3DLaw::Info() const { return "Newtonian3DLaw"; }

    double Newtonian3DLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
    {
        return rParameters.GetMaterialProperties()[rVariable];
    }

    void Newtonian3DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
    }

    void Newtonian3DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
    }

} // Namespace Kratos
