//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Riccardo Rossi, Alessandro Franci
//  Collaborators:  Massimiliano Zecchetto
//
//-------------------------------------------------------------
//

// System includes
#include <iostream>

// External includes
#include <cmath>

// Project includes
#include "custom_constitutive/fluid_laws/frictional_viscoplastic_2D_law.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    FrictionalViscoplastic2DLaw::FrictionalViscoplastic2DLaw() : PfemFluidConstitutiveLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    FrictionalViscoplastic2DLaw::FrictionalViscoplastic2DLaw(const FrictionalViscoplastic2DLaw &rOther) : PfemFluidConstitutiveLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer FrictionalViscoplastic2DLaw::Clone() const { return Kratos::make_shared<FrictionalViscoplastic2DLaw>(*this); }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    FrictionalViscoplastic2DLaw::~FrictionalViscoplastic2DLaw() {}

    ConstitutiveLaw::SizeType FrictionalViscoplastic2DLaw::WorkingSpaceDimension() { return 2; }

    ConstitutiveLaw::SizeType FrictionalViscoplastic2DLaw::GetStrainSize() { return 3; }

    void FrictionalViscoplastic2DLaw::CalculateMaterialResponseCauchy(Parameters &rValues)
    {

        Flags &r_options = rValues.GetOptions();

        const Properties &r_properties = rValues.GetMaterialProperties();

        Vector &r_strain_vector = rValues.GetStrainVector();
        Vector &r_stress_vector = rValues.GetStressVector();

        const double dynamic_viscosity = this->GetEffectiveDynamicViscosity(rValues);
        const double friction_angle = r_properties[INTERNAL_FRICTION_ANGLE];
        const double cohesion = r_properties[COHESION];
        const double adaptive_exponent = r_properties[ADAPTIVE_EXPONENT];
        double effective_dynamic_viscosity = 0;

        const double old_pressure = this->CalculateInGaussPoint(PRESSURE, rValues, 1);
        const double new_pressure = this->CalculateInGaussPoint(PRESSURE, rValues, 0);

        const double theta_momentum = this->GetThetaMomentumForPressureIntegration();
        double mean_pressure = (1.0 - theta_momentum) * old_pressure + theta_momentum * new_pressure;

        if (mean_pressure > 0.0) // cutoff for tractions
        {
            mean_pressure = 0;
        }

        const double equivalent_strain_rate =
            std::sqrt(2.0 * r_strain_vector[0] * r_strain_vector[0] + 2.0 * r_strain_vector[1] * r_strain_vector[1] +
                      4.0 * r_strain_vector[2] * r_strain_vector[2]);

        // Ensuring that the case of equivalent_strain_rate = 0 is not problematic
        const double tolerance = 1e-12;
        if (equivalent_strain_rate < tolerance)
        {
            effective_dynamic_viscosity = dynamic_viscosity;
        }
        else
        {
            const double friction_angle_rad = friction_angle * Globals::Pi / 180.0;
            const double tanFi = std::tan(friction_angle_rad);
            double regularization = 1.0 - std::exp(-adaptive_exponent * equivalent_strain_rate);
            effective_dynamic_viscosity = dynamic_viscosity + regularization * (cohesion + tanFi * fabs(mean_pressure) / equivalent_strain_rate);
        }

        const double strain_trace = r_strain_vector[0] + r_strain_vector[1];

        r_stress_vector[0] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[0] - strain_trace / 3.0);
        r_stress_vector[1] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[1] - strain_trace / 3.0);
        r_stress_vector[2] = 2.0 * effective_dynamic_viscosity * r_strain_vector[2];

        if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
        {
            this->EffectiveViscousConstitutiveMatrix2D(effective_dynamic_viscosity, rValues.GetConstitutiveMatrix());
        }
    }

    std::string FrictionalViscoplastic2DLaw::Info() const { return "FrictionalViscoplastic2DLaw"; }

    //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW******************
    //*****************************************************************************

    int FrictionalViscoplastic2DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                           const ProcessInfo &rCurrentProcessInfo)
    {

        if (rMaterialProperties[DYNAMIC_VISCOSITY] < 0.0)
        {
            KRATOS_ERROR << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for FrictionalViscoplastic2DLaw: "
                         << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;
        }

        if (rMaterialProperties[INTERNAL_FRICTION_ANGLE] < 0.0)
        {
            KRATOS_ERROR << "Incorrect or missing INTERNAL_FRICTION_ANGLE provided in process info for FrictionalViscoplastic2DLaw: "
                         << rMaterialProperties[INTERNAL_FRICTION_ANGLE] << std::endl;
        }

        if (rMaterialProperties[COHESION] < 0.0)
        {
            KRATOS_ERROR << "Incorrect or missing COHESION provided in process info for FrictionalViscoplastic2DLaw: "
                         << rMaterialProperties[COHESION] << std::endl;
        }

        if (rMaterialProperties[ADAPTIVE_EXPONENT] < 0.0)
        {
            KRATOS_ERROR << "Incorrect or missing ADAPTIVE_EXPONENT provided in process info for FrictionalViscoplastic2DLaw: "
                         << rMaterialProperties[ADAPTIVE_EXPONENT] << std::endl;
        }

        if (rMaterialProperties[BULK_MODULUS] <= 0.0)
        {
            KRATOS_ERROR << "Incorrect or missing BULK_MODULUS provided in process info for FrictionalViscoplastic2DLaw: "
                         << rMaterialProperties[BULK_MODULUS] << std::endl;
        }

        return 0;
    }

    double FrictionalViscoplastic2DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters &rParameters) const
    {
        return rParameters.GetConstitutiveMatrix()(2, 2);
    }

    double FrictionalViscoplastic2DLaw::GetEffectiveDensity(ConstitutiveLaw::Parameters &rParameters) const
    {
        return rParameters.GetMaterialProperties()[DENSITY];
    }

    double FrictionalViscoplastic2DLaw::GetEffectiveDynamicViscosity(ConstitutiveLaw::Parameters &rParameters) const
    {
        return rParameters.GetMaterialProperties()[DYNAMIC_VISCOSITY];
    }

    double FrictionalViscoplastic2DLaw::GetEffectiveFrictionAngle(ConstitutiveLaw::Parameters &rParameters) const
    {
        return rParameters.GetMaterialProperties()[INTERNAL_FRICTION_ANGLE];
    }

    double FrictionalViscoplastic2DLaw::GetEffectiveCohesion(ConstitutiveLaw::Parameters &rParameters) const
    {
        return rParameters.GetMaterialProperties()[COHESION];
    }

    void FrictionalViscoplastic2DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
    }

    void FrictionalViscoplastic2DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
    }

} // Namespace Kratos
