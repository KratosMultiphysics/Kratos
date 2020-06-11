//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Alessandro Franci
//  Collaborators:  Massimiliano Zecchetto
//
//-------------------------------------------------------------
//

// System includes
#include <iostream>

// External includes
#include <cmath>

// Project includes
#include "custom_constitutive/fluid_laws/bercovier_mu_I_rheology_3D_law.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos {

//********************************CONSTRUCTOR*********************************
//****************************************************************************

BercovierMuIRheology3DLaw::BercovierMuIRheology3DLaw() : PfemFluidConstitutiveLaw() {}

//******************************COPY CONSTRUCTOR******************************
//****************************************************************************

BercovierMuIRheology3DLaw::BercovierMuIRheology3DLaw(const BercovierMuIRheology3DLaw& rOther)
    : PfemFluidConstitutiveLaw(rOther) {}

//***********************************CLONE************************************
//****************************************************************************

ConstitutiveLaw::Pointer BercovierMuIRheology3DLaw::Clone() const {
    return Kratos::make_shared<BercovierMuIRheology3DLaw>(*this);
}

//*********************************DESTRUCTOR*********************************
//****************************************************************************

BercovierMuIRheology3DLaw::~BercovierMuIRheology3DLaw() {}

ConstitutiveLaw::SizeType BercovierMuIRheology3DLaw::WorkingSpaceDimension() { return 3; }

ConstitutiveLaw::SizeType BercovierMuIRheology3DLaw::GetStrainSize() { return 6; }

void BercovierMuIRheology3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues) {

    Flags& r_options = rValues.GetOptions();

    const Properties& r_properties = rValues.GetMaterialProperties();

    Vector& r_strain_vector = rValues.GetStrainVector();
    Vector& r_stress_vector = rValues.GetStressVector();

    const double static_friction = r_properties[STATIC_FRICTION];
    const double dynamic_friction = r_properties[DYNAMIC_FRICTION];
    const double delta_friction = dynamic_friction - static_friction;
    const double inertial_number_zero = r_properties[INERTIAL_NUMBER_ZERO];
    const double grain_diameter = r_properties[GRAIN_DIAMETER];
    const double grain_density = r_properties[GRAIN_DENSITY];
    const double regularization_coeff = r_properties[REGULARIZATION_COEFFICIENT];
    double inertial_number = 0;
    double effective_dynamic_viscosity = 0;

    const double old_pressure = this->CalculateInGaussPoint(PRESSURE, rValues, 1);
    const double new_pressure = this->CalculateInGaussPoint(PRESSURE, rValues, 0);
    const GeometryType& r_geometry = rValues.GetElementGeometry();

    const double theta_momentum = r_geometry[0].GetValue(THETA_MOMENTUM);
    double mean_pressure = (1.0 - theta_momentum) * old_pressure + theta_momentum * new_pressure;
    if (mean_pressure > 0.0) {
        mean_pressure = 0.0000001;
    }

    const double equivalent_strain_rate =
        std::sqrt(2.0 * r_strain_vector[0] * r_strain_vector[0] + 2.0 * r_strain_vector[1] * r_strain_vector[1] +
                  2.0 * r_strain_vector[2] * r_strain_vector[2] + 4.0 * r_strain_vector[3] * r_strain_vector[3] +
                  4.0 * r_strain_vector[4] * r_strain_vector[4] + 4.0 * r_strain_vector[5] * r_strain_vector[5]);

    if (mean_pressure != 0) {
        inertial_number = equivalent_strain_rate * grain_diameter / std::sqrt(std::fabs(mean_pressure) / grain_density);
    }

    if (equivalent_strain_rate != 0 && std::fabs(mean_pressure) != 0) {
        const double first_viscous_term =
            static_friction / std::sqrt(std::pow(equivalent_strain_rate, 2) + std::pow(regularization_coeff, 2));
        const double second_viscous_term =
            delta_friction * inertial_number / ((inertial_number_zero + inertial_number) * equivalent_strain_rate);
        effective_dynamic_viscosity = (first_viscous_term + second_viscous_term) * std::fabs(mean_pressure);
    } else {
        effective_dynamic_viscosity = 1.0;
    }

    const double strain_trace = r_strain_vector[0] + r_strain_vector[1] + r_strain_vector[2];

    r_stress_vector[0] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[0] - strain_trace / 3.0);
    r_stress_vector[1] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[1] - strain_trace / 3.0);
    r_stress_vector[2] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[2] - strain_trace / 3.0);
    r_stress_vector[3] = 2.0 * effective_dynamic_viscosity * r_strain_vector[3];
    r_stress_vector[4] = 2.0 * effective_dynamic_viscosity * r_strain_vector[4];
    r_stress_vector[5] = 2.0 * effective_dynamic_viscosity * r_strain_vector[5];

    if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        this->EffectiveViscousConstitutiveMatrix3D(effective_dynamic_viscosity, rValues.GetConstitutiveMatrix());
    }
}

std::string BercovierMuIRheology3DLaw::Info() const { return "BercovierMuIRheology3DLaw"; }

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW******************
//*****************************************************************************

int BercovierMuIRheology3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                                     const ProcessInfo& rCurrentProcessInfo) {
    KRATOS_CHECK_VARIABLE_KEY(STATIC_FRICTION);
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_FRICTION);
    KRATOS_CHECK_VARIABLE_KEY(INERTIAL_NUMBER_ZERO);
    KRATOS_CHECK_VARIABLE_KEY(GRAIN_DIAMETER);
    KRATOS_CHECK_VARIABLE_KEY(GRAIN_DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(REGULARIZATION_COEFFICIENT);
    KRATOS_CHECK_VARIABLE_KEY(BULK_MODULUS);

    if (rMaterialProperties[STATIC_FRICTION] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing STATIC_FRICTION provided in process info for BercovierMuIRheology3DLaw: "
                     << rMaterialProperties[STATIC_FRICTION] << std::endl;
    }

    if (rMaterialProperties[DYNAMIC_FRICTION] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing DYNAMIC_FRICTION provided in process info for BercovierMuIRheology3DLaw: "
                     << rMaterialProperties[DYNAMIC_FRICTION] << std::endl;
    }

    if (rMaterialProperties[INERTIAL_NUMBER_ZERO] < 0.0) {
        KRATOS_ERROR
            << "Incorrect or missing INERTIAL_NUMBER_ZERO provided in process info for BercovierMuIRheology3DLaw: "
            << rMaterialProperties[INERTIAL_NUMBER_ZERO] << std::endl;
    }

    if (rMaterialProperties[GRAIN_DIAMETER] <= 0.0) {
        KRATOS_ERROR << "Incorrect or missing GRAIN_DIAMETER provided in process info for BercovierMuIRheology3DLaw: "
                     << rMaterialProperties[GRAIN_DIAMETER] << std::endl;
    }

    if (rMaterialProperties[GRAIN_DENSITY] <= 0.0) {
        KRATOS_ERROR << "Incorrect or missing GRAIN_DENSITY provided in process info for BercovierMuIRheology3DLaw: "
                     << rMaterialProperties[GRAIN_DENSITY] << std::endl;
    }

    if (rMaterialProperties[REGULARIZATION_COEFFICIENT] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing REGULARIZATION_COEFFICIENT provided in process info for "
                        "BercovierMuIRheology3DLaw: "
                     << rMaterialProperties[REGULARIZATION_COEFFICIENT] << std::endl;
    }

    if (rMaterialProperties[BULK_MODULUS] <= 0.0) {
        KRATOS_ERROR << "Incorrect or missing BULK_MODULUS provided in process info for "
                        "BercovierMuIRheology3DLaw: "
                     << rMaterialProperties[BULK_MODULUS] << std::endl;
    }

    return 0;
}

double BercovierMuIRheology3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetConstitutiveMatrix()(5, 5);
}

double BercovierMuIRheology3DLaw::GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetMaterialProperties()[DENSITY];
}

void BercovierMuIRheology3DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
}

void BercovierMuIRheology3DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
}

}  // Namespace Kratos
