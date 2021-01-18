// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_soft_torque_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_soft_torque::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_soft_torque(*this));
        return p_clone;
    }

    void DEM_KDEM_soft_torque::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_soft_torque to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_KDEM_soft_torque::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                    SphericContinuumParticle* neighbor,
                                                    double equiv_young,
                                                    double distance,
                                                    double calculation_area,
                                                    double LocalCoordSystem[3][3],
                                                    double ElasticLocalRotationalMoment[3],
                                                    double ViscoLocalRotationalMoment[3],
                                                    double equiv_poisson,
                                                    double indentation) {

        KRATOS_TRY

        double LocalDeltaRotatedAngle[3]    = {0.0};

        array_1d<double, 3> GlobalDeltaRotatedAngle;
        noalias(GlobalDeltaRotatedAngle) = element->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);

        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaRotatedAngle, LocalDeltaRotatedAngle);

        const double equivalent_radius = sqrt(calculation_area / Globals::Pi);
        const double Inertia_I = 0.25 * Globals::Pi * equivalent_radius * equivalent_radius * equivalent_radius * equivalent_radius;
        const double Inertia_J = 2.0 * Inertia_I; // This is the polar inertia
        double rotational_moment_coeff = element->GetProperties()[ROTATIONAL_MOMENT_COEFFICIENT];

        ElasticLocalRotationalMoment[0] = -rotational_moment_coeff * equiv_young * Inertia_I * LocalDeltaRotatedAngle[0] / distance;
        ElasticLocalRotationalMoment[1] = -rotational_moment_coeff * equiv_young * Inertia_I * LocalDeltaRotatedAngle[1] / distance;
        ElasticLocalRotationalMoment[2] = -rotational_moment_coeff * equiv_young * Inertia_J * LocalDeltaRotatedAngle[2] / distance;

        KRATOS_CATCH("")
    }//ComputeParticleRotationalMoments

} // namespace Kratos
