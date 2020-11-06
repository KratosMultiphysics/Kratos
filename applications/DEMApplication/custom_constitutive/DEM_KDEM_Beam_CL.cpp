// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_Beam_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_Beam::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_Beam(*this));
        return p_clone;
    }

    void DEM_KDEM_Beam::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_Beam to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_KDEM_Beam::Check(Properties::Pointer pProp) const {
        DEMContinuumConstitutiveLaw::Check(pProp);

        if(!pProp->Has(BEAM_CROSS_SECTION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_CROSS_SECTION should be present in the properties when using DEM_KDEM_Beam. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_CROSS_SECTION) = 1.0;
        }

        if(!pProp->Has(BEAM_LENGTH)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_LENGTH should be present in the properties when using DEM_KDEM_Beam. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_LENGTH) = 1.0;
        }

        if(!pProp->Has(BEAM_ELEMENTS_DISTANCE)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_ELEMENTS_DISTANCE should be present in the properties when using DEM_KDEM_Beam. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_ELEMENTS_DISTANCE) = 0.0;
        }

        if(!pProp->Has(BEAM_PLANAR_MOMENT_OF_INERTIA_XX)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_PLANAR_MOMENT_OF_INERTIA_XX should be present in the properties when using DEM_KDEM_Beam. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_PLANAR_MOMENT_OF_INERTIA_XX) = 1.0;
        }

        if(!pProp->Has(BEAM_PLANAR_MOMENT_OF_INERTIA_YY)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_PLANAR_MOMENT_OF_INERTIA_YY should be present in the properties when using DEM_KDEM_Beam. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_PLANAR_MOMENT_OF_INERTIA_YY) = 1.0;
        }

        if(!pProp->Has(BEAM_MOMENT_OF_INERTIA_PER_METER_X)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_MOMENT_OF_INERTIA_PER_METER_X should be present in the properties when using DEM_KDEM_Beam. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_MOMENT_OF_INERTIA_PER_METER_X) = 0.0;
        }

        if(!pProp->Has(BEAM_MOMENT_OF_INERTIA_PER_METER_Y)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_MOMENT_OF_INERTIA_PER_METER_Y should be present in the properties when using DEM_KDEM_Beam. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_MOMENT_OF_INERTIA_PER_METER_Y) = 1.0;
        }

        if(!pProp->Has(BEAM_MOMENT_OF_INERTIA_PER_METER_Z)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_MOMENT_OF_INERTIA_PER_METER_Z should be present in the properties when using DEM_KDEM_Beam. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_MOMENT_OF_INERTIA_PER_METER_Z) = 1.0;
        }
    }

    void DEM_KDEM_Beam::CalculateElasticConstants(double& kn_el,
                                                  double& kt_el_0,
                                                  double& kt_el_1,
                                                  double initial_dist,
                                                  double equiv_young,
                                                  double equiv_poisson,
                                                  double calculation_area,
                                                  SphericContinuumParticle* element1,
                                                  SphericContinuumParticle* element2) {

        KRATOS_TRY

        kn_el = equiv_young * calculation_area / initial_dist;

        const double Inertia_Ixx = 0.5 * (element1->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_XX] + element2->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_XX]);
        const double Inertia_Iyy = 0.5 * (element1->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_YY] + element2->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_YY]);

        kt_el_0 = 3.0 * equiv_young * Inertia_Iyy / (calculation_area * initial_dist);
        kt_el_1 = 3.0 * equiv_young * Inertia_Ixx / (calculation_area * initial_dist);

        KRATOS_CATCH("")
    }

    void DEM_KDEM_Beam::CalculateViscoDampingCoeff(double& equiv_visco_damp_coeff_normal,
                                                   double& equiv_visco_damp_coeff_tangential_0,
                                                   double& equiv_visco_damp_coeff_tangential_1,
                                                   SphericContinuumParticle* element1,
                                                   SphericContinuumParticle* element2,
                                                   const double kn_el,
                                                   const double kt_el_0,
                                                   const double kt_el_1) {

        KRATOS_TRY

        const double equiv_mass  = 0.5 * (element1->GetMass() + element2->GetMass());

        const double beam_total_mass = element1->GetProperties()[BEAM_LENGTH] * element1->GetProperties()[BEAM_CROSS_SECTION] * element1->GetDensity();
        const double aux_mass = beam_total_mass / equiv_mass;

        const double equiv_gamma = 0.5 * (element1->GetProperties()[DAMPING_GAMMA] + element2->GetProperties()[DAMPING_GAMMA]);

        equiv_visco_damp_coeff_normal       = equiv_gamma * aux_mass * sqrt(equiv_mass * kn_el  );
        equiv_visco_damp_coeff_tangential_0 = equiv_gamma * aux_mass * sqrt(equiv_mass * kt_el_0);
        equiv_visco_damp_coeff_tangential_1 = equiv_gamma * aux_mass * sqrt(equiv_mass * kt_el_1);

        KRATOS_CATCH("")
    }

    void DEM_KDEM_Beam::CalculateForces(const ProcessInfo& r_process_info,
                                        double OldLocalElasticContactForce[3],
                                        double LocalElasticContactForce[3],
                                        double LocalElasticExtraContactForce[3],
                                        double LocalCoordSystem[3][3],
                                        double LocalDeltDisp[3],
                                        const double kn_el,
                                        const double kt_el_0,
                                        const double kt_el_1,
                                        double& contact_sigma,
                                        double& contact_tau,
                                        double& failure_criterion_state,
                                        double equiv_young,
                                        double equiv_shear,
                                        double indentation,
                                        double calculation_area,
                                        double& acumulated_damage,
                                        SphericContinuumParticle* element1,
                                        SphericContinuumParticle* element2,
                                        int i_neighbour_count,
                                        int time_steps,
                                        bool& sliding,
                                        double &equiv_visco_damp_coeff_normal,
                                        double &equiv_visco_damp_coeff_tangential_0,
                                        double &equiv_visco_damp_coeff_tangential_1,
                                        double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3]) {

        KRATOS_TRY

        CalculateNormalForces(LocalElasticContactForce,
                              kn_el,
                              indentation);

        CalculateTangentialForces(OldLocalElasticContactForce,
                                  LocalElasticContactForce,
                                  LocalDeltDisp,
                                  kt_el_0,
                                  kt_el_1);

        CalculateViscoDampingCoeff(equiv_visco_damp_coeff_normal,
                                   equiv_visco_damp_coeff_tangential_0,
                                   equiv_visco_damp_coeff_tangential_1,
                                   element1,
                                   element2,
                                   kn_el,
                                   kt_el_0,
                                   kt_el_1);

        CalculateViscoDamping(LocalRelVel,
                              ViscoDampingLocalContactForce,
                              equiv_visco_damp_coeff_normal,
                              equiv_visco_damp_coeff_tangential_0,
                              equiv_visco_damp_coeff_tangential_1);

        KRATOS_CATCH("")
    }

    void DEM_KDEM_Beam::CalculateNormalForces(double LocalElasticContactForce[3],
                                              const double kn_el,
                                              double indentation) {

        KRATOS_TRY

            LocalElasticContactForce[2] = kn_el * indentation;

        KRATOS_CATCH("")
    }

    void DEM_KDEM_Beam::CalculateTangentialForces(double OldLocalElasticContactForce[3],
                                                  double LocalElasticContactForce[3],
                                                  double LocalDeltDisp[3],
                                                  const double kt_el_0,
                                                  const double kt_el_1) {

        KRATOS_TRY

        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_el_0 * LocalDeltDisp[0]; // 0: first tangential
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_el_1 * LocalDeltDisp[1]; // 1: second tangential

        KRATOS_CATCH("")
    }

    void DEM_KDEM_Beam::CalculateViscoDamping(double LocalRelVel[3],
                                              double ViscoDampingLocalContactForce[3],
                                              double equiv_visco_damp_coeff_normal,
                                              double equiv_visco_damp_coeff_tangential_0,
                                              double equiv_visco_damp_coeff_tangential_1) {

        KRATOS_TRY

        ViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal       * LocalRelVel[2];
        ViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential_0 * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential_1 * LocalRelVel[1];

        KRATOS_CATCH("")
    }

    void DEM_KDEM_Beam::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
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
        double LocalDeltaAngularVelocity[3] = {0.0};

        array_1d<double, 3> GlobalDeltaRotatedAngle;
        noalias(GlobalDeltaRotatedAngle) = element->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
        array_1d<double, 3> GlobalDeltaAngularVelocity;
        noalias(GlobalDeltaAngularVelocity) = element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaRotatedAngle, LocalDeltaRotatedAngle);
        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaAngularVelocity, LocalDeltaAngularVelocity);

        const double norm_distance = (element->GetRadius() + neighbor->GetRadius()) / distance; // If spheres are not tangent the Damping coefficient, DeltaRotatedAngle and DeltaAngularVelocity have to be normalized
        const double norm_length = element->GetProperties()[BEAM_LENGTH] / distance; // Total beam lenght / distance

        ////////////////////////////////////////////////////// ElasticLocalRotationalMoment //////////////////////////////////////////////////////

        const double equiv_shear   = equiv_young / (2.0 * (1 + equiv_poisson));

        const double Inertia_Ixx = 0.5 * (element->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_XX] + neighbor->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_XX]);
        const double Inertia_Iyy = 0.5 * (element->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_YY] + neighbor->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_YY]);
        const double Inertia_J  = Inertia_Ixx + Inertia_Iyy;

        const double k_rot_x = equiv_young * Inertia_Ixx * norm_distance / distance;
        const double k_rot_y = equiv_young * Inertia_Iyy * norm_distance / distance;
        const double k_tor   = equiv_shear * Inertia_J  / distance;

        ElasticLocalRotationalMoment[0] = -k_rot_x * LocalDeltaRotatedAngle[0];
        ElasticLocalRotationalMoment[1] = -k_rot_y * LocalDeltaRotatedAngle[1];
        ElasticLocalRotationalMoment[2] = -k_tor   * LocalDeltaRotatedAngle[2];

        ////////////////////////////////////////////////////// ViscoLocalRotationalMoment //////////////////////////////////////////////////////

        const double equiv_gamma = 0.5 * (element->GetProperties()[DAMPING_GAMMA] + neighbor->GetProperties()[DAMPING_GAMMA]);

        const double a = std::sqrt(12.0 * element->GetProperties()[BEAM_MOMENT_OF_INERTIA_PER_METER_Y] - 1.0);
        const double b = std::sqrt(12.0 * element->GetProperties()[BEAM_MOMENT_OF_INERTIA_PER_METER_Z] - 1.0);

        const double equiv_mass  = 0.5 * (element->GetMass() + neighbor->GetMass());
        const double beam_total_mass = element->GetProperties()[BEAM_LENGTH] * element->GetProperties()[BEAM_CROSS_SECTION] * element->GetDensity();
        const double aux_mass = beam_total_mass / equiv_mass;

        const double equiv_moment_of_inertiaX = 0.083333333 * (a * a + distance * distance) * equiv_mass;
        const double equiv_moment_of_inertiaY = 0.083333333 * (b * b + distance * distance) * equiv_mass;
        const double equiv_moment_of_inertiaZ = element->GetProperties()[BEAM_MOMENT_OF_INERTIA_PER_METER_X] * equiv_mass;

        const double visc_param_rot_x = equiv_gamma * aux_mass * norm_length * sqrt(equiv_moment_of_inertiaX * k_rot_x);
        const double visc_param_rot_y = equiv_gamma * aux_mass * norm_length * sqrt(equiv_moment_of_inertiaY * k_rot_y);
        const double visc_param_tor = equiv_gamma * aux_mass * sqrt(equiv_moment_of_inertiaZ * k_tor);

        ViscoLocalRotationalMoment[0] = -visc_param_rot_x * LocalDeltaAngularVelocity[0];
        ViscoLocalRotationalMoment[1] = -visc_param_rot_y * LocalDeltaAngularVelocity[1];
        ViscoLocalRotationalMoment[2] = -visc_param_tor * LocalDeltaAngularVelocity[2];

        KRATOS_CATCH("")
    }//ComputeParticleRotationalMoments

} // namespace Kratos
