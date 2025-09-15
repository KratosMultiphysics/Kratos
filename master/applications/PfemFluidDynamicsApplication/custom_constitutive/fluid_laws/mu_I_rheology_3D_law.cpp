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
#include "custom_constitutive/fluid_laws/mu_I_rheology_3D_law.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    MuIRheology3DLaw::MuIRheology3DLaw() : PfemFluidConstitutiveLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    MuIRheology3DLaw::MuIRheology3DLaw(const MuIRheology3DLaw &rOther)
        : PfemFluidConstitutiveLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer MuIRheology3DLaw::Clone() const
    {
        return Kratos::make_shared<MuIRheology3DLaw>(*this);
    }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    MuIRheology3DLaw::~MuIRheology3DLaw() {}

    ConstitutiveLaw::SizeType MuIRheology3DLaw::WorkingSpaceDimension() { return 3; }

    ConstitutiveLaw::SizeType MuIRheology3DLaw::GetStrainSize() const { return 6; }

    void MuIRheology3DLaw::CalculateMaterialResponseCauchy(Parameters &rValues)
    {

        Flags &r_options = rValues.GetOptions();

        const Properties &r_properties = rValues.GetMaterialProperties();

        Vector &r_strain_vector = rValues.GetStrainVector();
        Vector &r_stress_vector = rValues.GetStressVector();

        const double static_friction = r_properties[STATIC_FRICTION];
        const double dynamic_friction = r_properties[DYNAMIC_FRICTION];
        const double delta_friction = dynamic_friction - static_friction;
        const double inertial_number_zero = r_properties[INERTIAL_NUMBER_ZERO];
        const double grain_diameter = r_properties[GRAIN_DIAMETER];
        const double grain_density = r_properties[GRAIN_DENSITY];
        const double regularization_coeff = r_properties[REGULARIZATION_COEFFICIENT];
        double effective_dynamic_viscosity = 0;

        const double equivalent_strain_rate =
            std::sqrt(2.0 * r_strain_vector[0] * r_strain_vector[0] + 2.0 * r_strain_vector[1] * r_strain_vector[1] +
                      2.0 * r_strain_vector[2] * r_strain_vector[2] + 4.0 * r_strain_vector[3] * r_strain_vector[3] +
                      4.0 * r_strain_vector[4] * r_strain_vector[4] + 4.0 * r_strain_vector[5] * r_strain_vector[5]);

        const double old_pressure = this->CalculateInGaussPoint(PRESSURE, rValues, 1);
        const double new_pressure = this->CalculateInGaussPoint(PRESSURE, rValues, 0);
        const double theta_momentum = GetThetaMomentumForPressureIntegration();
        double mean_pressure = (1.0 - theta_momentum) * old_pressure + theta_momentum * new_pressure;
        const double pressure_tolerance = -1.0e-07;
        if (mean_pressure > pressure_tolerance)
        {
            mean_pressure = pressure_tolerance;
        }

        const double exponent = -equivalent_strain_rate / regularization_coeff;

        const double second_viscous_term = delta_friction * grain_diameter / (inertial_number_zero * std::sqrt(std::fabs(mean_pressure) / grain_density) + equivalent_strain_rate * grain_diameter);
        if (equivalent_strain_rate != 0)
        {
            const double first_viscous_term = static_friction * (1 - std::exp(exponent)) / equivalent_strain_rate;
            effective_dynamic_viscosity = (first_viscous_term + second_viscous_term) * std::fabs(mean_pressure);
        }
        else
        {
            effective_dynamic_viscosity = 1.0; // this is for the first iteration and first time step
        }

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

    std::string MuIRheology3DLaw::Info() const { return "MuIRheology3DLaw"; }

    //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW******************
    //*****************************************************************************

    int MuIRheology3DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                const ProcessInfo &rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF(rMaterialProperties[STATIC_FRICTION] < 0.0)
            << "Incorrect or missing STATIC_FRICTION provided in process info for MuIRheology3DLaw: "
            << rMaterialProperties[STATIC_FRICTION] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_FRICTION] < 0.0)
            << "Incorrect or missing DYNAMIC_FRICTION provided in process info for MuIRheology3DLaw: "
            << rMaterialProperties[DYNAMIC_FRICTION] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[INERTIAL_NUMBER_ZERO] < 0.0)
            << "Incorrect or missing INERTIAL_NUMBER_ZERO provided in process info for MuIRheology3DLaw: "
            << rMaterialProperties[INERTIAL_NUMBER_ZERO] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[GRAIN_DIAMETER] < 0.0)
            << "Incorrect or missing GRAIN_DIAMETER provided in process info for MuIRheology3DLaw: "
            << rMaterialProperties[GRAIN_DIAMETER] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[GRAIN_DENSITY] < 0.0)
            << "Incorrect or missing GRAIN_DENSITY provided in process info for MuIRheology3DLaw: "
            << rMaterialProperties[GRAIN_DENSITY] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[REGULARIZATION_COEFFICIENT] < 0.0)
            << "Incorrect or missing REGULARIZATION_COEFFICIENT provided in process info for MuIRheology3DLaw: "
            << rMaterialProperties[REGULARIZATION_COEFFICIENT] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[BULK_MODULUS] < 0.0)
            << "Incorrect or missing BULK_MODULUS provided in process info for MuIRheology3DLaw: "
            << rMaterialProperties[BULK_MODULUS] << std::endl;

        return 0;
    }

    double MuIRheology3DLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
    {
        return rParameters.GetMaterialProperties()[rVariable];
    }

    void MuIRheology3DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
    }

    void MuIRheology3DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
    }

} // Namespace Kratos
