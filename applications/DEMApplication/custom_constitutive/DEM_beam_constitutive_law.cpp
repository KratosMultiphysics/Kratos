// Project includes
#include "DEM_beam_constitutive_law.h"
#include "custom_elements/spheric_continuum_particle.h"
#include "dem_contact.h"

namespace Kratos {

    DEMBeamConstitutiveLaw::DEMBeamConstitutiveLaw() {} // Class DEMBeamConstitutiveLaw

    void DEMBeamConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEMBeamConstitutiveLaw to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_BEAM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEMBeamConstitutiveLaw::Initialize(SphericContinuumParticle* element1, SphericContinuumParticle* element2, Properties::Pointer pProps) {
        mpProperties = pProps;
    }

    void DEMBeamConstitutiveLaw::SetConstitutiveLawInPropertiesWithParameters(Properties::Pointer pProp, const Parameters& parameters, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEMBeamConstitutiveLaw to Properties " << pProp->Id() <<" with parameters"<< std::endl;
        pProp->SetValue(DEM_BEAM_CONSTITUTIVE_LAW_POINTER, this->Clone());

        this->Check(pProp);
    }

    void DEMBeamConstitutiveLaw::Check(Properties::Pointer pProp) const {

        if(!pProp->Has(STATIC_FRICTION)) {
            if(!pProp->Has(FRICTION)) { //deprecated since April 6th, 2020
                KRATOS_WARNING("DEM")<<std::endl;
                KRATOS_WARNING("DEM")<<"WARNING: Variable STATIC_FRICTION or FRICTION should be present in the properties when using DEMContinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
                KRATOS_WARNING("DEM")<<std::endl;
                pProp->GetValue(STATIC_FRICTION) = 0.0;
            }
            else {
                pProp->GetValue(STATIC_FRICTION) = pProp->GetValue(FRICTION);
            }
        }
        if(!pProp->Has(DYNAMIC_FRICTION)) {
            if(!pProp->Has(FRICTION)) { //deprecated since April 6th, 2020
                KRATOS_WARNING("DEM")<<std::endl;
                KRATOS_WARNING("DEM")<<"WARNING: Variable DYNAMIC_FRICTION or FRICTION should be present in the properties when using DEMBeamConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
                KRATOS_WARNING("DEM")<<std::endl;
                pProp->GetValue(DYNAMIC_FRICTION) = 0.0;
            }
            else {
                pProp->GetValue(DYNAMIC_FRICTION) = pProp->GetValue(FRICTION);
            }
        }
        if(!pProp->Has(FRICTION_DECAY)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable FRICTION_DECAY should be present in the properties when using DEMBeamConstitutiveLaw. 500.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(FRICTION_DECAY) = 500.0;
        }
        if(!pProp->Has(YOUNG_MODULUS)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable YOUNG_MODULUS should be present in the properties when using DEMBeamConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(YOUNG_MODULUS) = 0.0;
        }
        if(!pProp->Has(POISSON_RATIO)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable POISSON_RATIO should be present in the properties when using DEMBeamConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(POISSON_RATIO) = 0.0;
        }
        if(!pProp->Has(COEFFICIENT_OF_RESTITUTION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable COEFFICIENT_OF_RESTITUTION should be present in the properties when using DEMBeamConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(COEFFICIENT_OF_RESTITUTION) = 0.0;
        }

        if(!pProp->Has(CROSS_AREA)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable CROSS_AREA should be present in the properties when using DEMBeamConstitutiveLaw. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(CROSS_AREA) = 1.0;
        }

        if(!pProp->Has(BEAM_LENGTH)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_LENGTH should be present in the properties when using DEMBeamConstitutiveLaw. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_LENGTH) = 1.0;
        }

        if(!pProp->Has(BEAM_PARTICLES_DISTANCE)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_PARTICLES_DISTANCE should be present in the properties when using DEMBeamConstitutiveLaw. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_PARTICLES_DISTANCE) = 0.0;
        }

        if(!pProp->Has(I22)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable I22 should be present in the properties when using DEMBeamConstitutiveLaw. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(I22) = 1.0;
        }

        if(!pProp->Has(I33)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable I33 should be present in the properties when using DEMBeamConstitutiveLaw. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(I33) = 1.0;
        }

        if(!pProp->Has(BEAM_INERTIA_ROT_UNIT_LENGHT_X)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_INERTIA_ROT_UNIT_LENGHT_X should be present in the properties when using DEMBeamConstitutiveLaw. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_INERTIA_ROT_UNIT_LENGHT_X) = 0.0;
        }

        if(!pProp->Has(BEAM_INERTIA_ROT_UNIT_LENGHT_Y)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_INERTIA_ROT_UNIT_LENGHT_Y should be present in the properties when using DEMBeamConstitutiveLaw. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_INERTIA_ROT_UNIT_LENGHT_Y) = 1.0;
        }

        if(!pProp->Has(BEAM_INERTIA_ROT_UNIT_LENGHT_Z)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_INERTIA_ROT_UNIT_LENGHT_Z should be present in the properties when using DEMBeamConstitutiveLaw. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_INERTIA_ROT_UNIT_LENGHT_Z) = 1.0;
        }

        if(!pProp->Has(DEM_BEAM_CONSTITUTIVE_LAW_POINTER)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable BEAM_INERTIA_ROT_UNIT_LENGHT_Z should be present in the properties when using DEMBeamConstitutiveLaw. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(BEAM_INERTIA_ROT_UNIT_LENGHT_Z) = 1.0;
        }
    }

    std::string DEMBeamConstitutiveLaw::GetTypeOfLaw() {
        std::string type_of_law = "Beam";
        return type_of_law;
    }

    DEMBeamConstitutiveLaw::Pointer DEMBeamConstitutiveLaw::Clone() const {
        DEMBeamConstitutiveLaw::Pointer p_clone(new DEMBeamConstitutiveLaw(*this));
        return p_clone;
    }

    DEMBeamConstitutiveLaw::~DEMBeamConstitutiveLaw() {}

    void DEMBeamConstitutiveLaw::CalculateElasticConstants(double& kn_el,
                                                           double& kt_el_0,
                                                           double& kt_el_1,
                                                           double initial_dist,
                                                           double equiv_young,
                                                           double equiv_poisson,
                                                           double calculation_area,
                                                           SphericContinuumParticle* element1,
                                                           SphericContinuumParticle* element2, double indentation) {

        KRATOS_TRY

        kn_el = equiv_young * calculation_area / initial_dist;

        const double& Inertia_Iyy = (*mpProperties)[I22];
        const double& Inertia_Izz = (*mpProperties)[I33];

        kt_el_0 = 3.0 * equiv_young * Inertia_Izz / (calculation_area * initial_dist);
        kt_el_1 = 3.0 * equiv_young * Inertia_Iyy / (calculation_area * initial_dist);

        KRATOS_CATCH("")
    }

    void DEMBeamConstitutiveLaw::CalculateViscoDampingCoeff(double& equiv_visco_damp_coeff_normal,
                                                            double& equiv_visco_damp_coeff_tangential_0,
                                                            double& equiv_visco_damp_coeff_tangential_1,
                                                            SphericContinuumParticle* element1,
                                                            SphericContinuumParticle* element2,
                                                            const double kn_el,
                                                            const double kt_el_0,
                                                            const double kt_el_1) {

        KRATOS_TRY

        const double equiv_mass  = 0.5 * (element1->GetMass() + element2->GetMass());

        const double beam_total_mass = (*mpProperties)[BEAM_LENGTH] * (*mpProperties)[CROSS_AREA] * element1->GetDensity();
        const double aux_mass = beam_total_mass / equiv_mass;

        const double& equiv_gamma = (*mpProperties)[DAMPING_GAMMA];

        equiv_visco_damp_coeff_normal       = equiv_gamma * aux_mass * sqrt(equiv_mass * kn_el  );
        equiv_visco_damp_coeff_tangential_0 = equiv_gamma * aux_mass * sqrt(equiv_mass * kt_el_0);
        equiv_visco_damp_coeff_tangential_1 = equiv_gamma * aux_mass * sqrt(equiv_mass * kt_el_1);

        KRATOS_CATCH("")
    }

    void DEMBeamConstitutiveLaw::CalculateForces(const ProcessInfo& r_process_info,
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
                                  LocalRelVel,
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

    void DEMBeamConstitutiveLaw::CalculateNormalForces(double LocalElasticContactForce[3],
                                                       const double kn_el,
                                                       double indentation) {

        KRATOS_TRY

            LocalElasticContactForce[2] = kn_el * indentation;

        KRATOS_CATCH("")
    }

    void DEMBeamConstitutiveLaw::CalculateTangentialForces(double OldLocalElasticContactForce[3],
                                                           double LocalElasticContactForce[3],
                                                           double LocalDeltDisp[3],
                                                           double LocalRelVel[3],
                                                           const double kt_el_0,
                                                           const double kt_el_1) {

        KRATOS_TRY

        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_el_0 * LocalDeltDisp[0]; // 0: first tangential
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_el_1 * LocalDeltDisp[1]; // 1: second tangential

        KRATOS_CATCH("")
    }

    void DEMBeamConstitutiveLaw::CalculateViscoDamping(double LocalRelVel[3],
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

    void DEMBeamConstitutiveLaw::CalculateMoments(SphericContinuumParticle* element, 
                    SphericContinuumParticle* neighbor, 
                    double equiv_young, 
                    double distance, 
                    double calculation_area,
                    double LocalCoordSystem[3][3], 
                    double ElasticLocalRotationalMoment[3], 
                    double ViscoLocalRotationalMoment[3], 
                    double equiv_poisson, 
                    double indentation, 
                    double normalLocalContactForce,
                    double GlobalContactForce[3],
                    double LocalCoordSystem_2[3],
                    const int i_neighbor_count) 
    {
        KRATOS_TRY

        int failure_type = element->mIniNeighbourFailureId[i_neighbor_count];
        //int continuum_ini_neighbors_size = element->mContinuumInitialNeighborsSize;

        if (failure_type == 0) {
                ComputeParticleRotationalMoments(element, 
                                        neighbor, 
                                        equiv_young, 
                                        distance, 
                                        calculation_area,
                                        LocalCoordSystem, 
                                        ElasticLocalRotationalMoment, 
                                        ViscoLocalRotationalMoment, 
                                        equiv_poisson, 
                                        indentation);
        }              

        DemContact::ComputeParticleContactMoments(normalLocalContactForce,
                                                GlobalContactForce,
                                                LocalCoordSystem_2,
                                                element,
                                                neighbor,
                                                indentation,
                                                i_neighbor_count);

        KRATOS_CATCH("")
    }

    void DEMBeamConstitutiveLaw::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
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
        const double norm_length = (*mpProperties)[BEAM_LENGTH] / distance; // Total beam lenght / distance

        ////////////////////////////////////////////////////// ElasticLocalRotationalMoment //////////////////////////////////////////////////////

        const double equiv_shear   = equiv_young / (2.0 * (1 + equiv_poisson));

        const double& Inertia_Iyy = (*mpProperties)[I22];
        const double& Inertia_Izz = (*mpProperties)[I33];
        const double Inertia_J = Inertia_Iyy + Inertia_Izz;

        const double k_rot_x = equiv_young * Inertia_Iyy * norm_distance / distance;
        const double k_rot_y = equiv_young * Inertia_Izz * norm_distance / distance;
        const double k_tor   = equiv_shear * Inertia_J  / distance;

        ElasticLocalRotationalMoment[0] = -k_rot_x * LocalDeltaRotatedAngle[0];
        ElasticLocalRotationalMoment[1] = -k_rot_y * LocalDeltaRotatedAngle[1];
        ElasticLocalRotationalMoment[2] = -k_tor   * LocalDeltaRotatedAngle[2];

        ////////////////////////////////////////////////////// ViscoLocalRotationalMoment //////////////////////////////////////////////////////

        const double& equiv_gamma = (*mpProperties)[DAMPING_GAMMA];

        const double a = std::sqrt(12.0 * (*mpProperties)[BEAM_INERTIA_ROT_UNIT_LENGHT_Y] - 1.0);
        const double b = std::sqrt(12.0 * (*mpProperties)[BEAM_INERTIA_ROT_UNIT_LENGHT_Z] - 1.0);

        const double equiv_mass  = 0.5 * (element->GetMass() + neighbor->GetMass());
        const double beam_total_mass = (*mpProperties)[BEAM_LENGTH] * (*mpProperties)[CROSS_AREA] * element->GetDensity();
        const double aux_mass = beam_total_mass / equiv_mass;

        const double equiv_moment_of_inertiaX = 0.083333333 * (a * a + distance * distance) * equiv_mass;
        const double equiv_moment_of_inertiaY = 0.083333333 * (b * b + distance * distance) * equiv_mass;
        const double equiv_moment_of_inertiaZ = (*mpProperties)[BEAM_INERTIA_ROT_UNIT_LENGHT_X] * equiv_mass;

        const double visc_param_rot_x = equiv_gamma * aux_mass * norm_length * sqrt(equiv_moment_of_inertiaX * k_rot_x);
        const double visc_param_rot_y = equiv_gamma * aux_mass * norm_length * sqrt(equiv_moment_of_inertiaY * k_rot_y);
        const double visc_param_tor = equiv_gamma * aux_mass * sqrt(equiv_moment_of_inertiaZ * k_tor);

        ViscoLocalRotationalMoment[0] = -visc_param_rot_x * LocalDeltaAngularVelocity[0];
        ViscoLocalRotationalMoment[1] = -visc_param_rot_y * LocalDeltaAngularVelocity[1];
        ViscoLocalRotationalMoment[2] = -visc_param_tor * LocalDeltaAngularVelocity[2];

        KRATOS_CATCH("")
    }//ComputeParticleRotationalMoments

    bool DEMBeamConstitutiveLaw::CheckRequirementsOfStressTensor() { return false; }

} // namespace Kratos
