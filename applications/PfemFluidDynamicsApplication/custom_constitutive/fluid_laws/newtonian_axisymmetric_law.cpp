//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Massimiliano Zecchetto
//  Collaborators:
//
//-------------------------------------------------------------
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/fluid_laws/newtonian_axisymmetric_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos {

//********************************CONSTRUCTOR*********************************
//****************************************************************************

NewtonianAxisymmetricLaw::NewtonianAxisymmetricLaw() : PfemFluidConstitutiveLaw() {}

//******************************COPY CONSTRUCTOR******************************
//****************************************************************************

NewtonianAxisymmetricLaw::NewtonianAxisymmetricLaw(const NewtonianAxisymmetricLaw& rOther) : PfemFluidConstitutiveLaw(rOther) {}

//***********************************CLONE************************************
//****************************************************************************

ConstitutiveLaw::Pointer NewtonianAxisymmetricLaw::Clone() const { return Kratos::make_shared<NewtonianAxisymmetricLaw>(*this); }

//*********************************DESTRUCTOR*********************************
//****************************************************************************

NewtonianAxisymmetricLaw::~NewtonianAxisymmetricLaw() {}

ConstitutiveLaw::SizeType NewtonianAxisymmetricLaw::WorkingSpaceDimension() { return 2; }

ConstitutiveLaw::SizeType NewtonianAxisymmetricLaw::GetStrainSize() { return 4; }

void NewtonianAxisymmetricLaw::CalculateMaterialResponseCauchy(Parameters& rValues) {

    const Flags& r_options = rValues.GetOptions();
    const Vector& r_strain_vector = rValues.GetStrainVector();
    Vector& r_stress_vector = rValues.GetStressVector();

    double effective_dynamic_viscosity = this->GetEffectiveViscosity(rValues); //const TO BE ADDED LATER

	// temporary workaround for slip conditions - start
	const GeometryType& r_geometry = rValues.GetElementGeometry();
	for (SizeType i = 0; i < r_geometry.PointsNumber(); i++) {
		if (r_geometry[i].Is(RIGID) && r_geometry[i].Coordinates()[0] == 0.0) {
			//effective_dynamic_viscosity = 0.0;
			break;
		}
	}
	// temporary workaround for slip conditions - end

	const double strain_trace = r_strain_vector[0] + r_strain_vector[1] + r_strain_vector[2];

	r_stress_vector[0] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[0] - strain_trace / 3.0);
    r_stress_vector[1] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[1] - strain_trace / 3.0);
    r_stress_vector[2] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[2] - strain_trace / 3.0);
    r_stress_vector[3] = 2.0 * effective_dynamic_viscosity * r_strain_vector[3];

    if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        this->EffectiveViscousConstitutiveMatrixAxisymmetric(effective_dynamic_viscosity, rValues.GetConstitutiveMatrix());
    }
}

int NewtonianAxisymmetricLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                          const ProcessInfo& rCurrentProcessInfo) {
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(BULK_MODULUS);

    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for NewtonianAxisymmetricLaw: "
        << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[BULK_MODULUS] <= 0.0)
        << "Incorrect or missing BULK_MODULUS provided in process info for NewtonianAxisymmetricLaw: "
        << rMaterialProperties[BULK_MODULUS] << std::endl;

    return 0;
}

std::string NewtonianAxisymmetricLaw::Info() const { return "NewtonianAxisymmetricLaw"; }

double NewtonianAxisymmetricLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    const double effective_viscosity = r_properties[DYNAMIC_VISCOSITY];
    return effective_viscosity;
}

double NewtonianAxisymmetricLaw::GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    const double effective_density = r_properties[DENSITY];
    return effective_density;
}

void NewtonianAxisymmetricLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
}

void NewtonianAxisymmetricLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
}

}  // Namespace Kratos
