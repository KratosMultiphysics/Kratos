// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_Dempack_CL.h"
//#include "custom_elements/spheric_particle.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    void DEM_Dempack::Initialize(SphericContinuumParticle* owner_sphere) {

    KRATOS_TRY
        mHistoryMaxInd              = 0.0; //maximum indentation achieved
        mHistoryMaxForce            = 0.0; //maximum force achieved
        mHistoryDamage              = 0.0; //cumulated_damage
        mHistoryDegradation         = 1.0; //degradation factor for G reducing in Dempack;
        mHistoryDisp                = 0.0; //displacement;
        mHistoryShearFlag           = 0.0; //shear limit achived;
    KRATOS_CATCH("")
    }

    DEMContinuumConstitutiveLaw::Pointer DEM_Dempack::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_Dempack(*this));
        return p_clone;
    }

    void DEM_Dempack::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_Dempack to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_Dempack::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

        double rmin = radius;
        if (other_radius < radius) rmin = other_radius;
        calculation_area = Globals::Pi * rmin*rmin;
    }

    void DEM_Dempack::GetContactArea(const double radius,
                                     const double other_radius,
                                     const Vector& vector_of_initial_areas,
                                     const int neighbour_position,
                                     double& calculation_area) {

        CalculateContactArea(radius, other_radius, calculation_area);
    }

    double DEM_Dempack::LocalMaxSearchDistance(const int i,
                                               SphericContinuumParticle* element1,
                                               SphericContinuumParticle* element2) {

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();
        double mDamageMaxDisplacementFactor;
        double mTensionLimit;

        // calculation of equivalent young modulus
        double myYoung = element1->GetYoung();
        double other_young = element2->GetYoung();
        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);

        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();

        double calculation_area = 0;
        CalculateContactArea(my_radius, other_radius, calculation_area);

        double radius_sum = my_radius + other_radius;
        double initial_delta = element1->GetInitialDelta(i);
        double initial_dist = radius_sum - initial_delta;

        // calculation of elastic constants
        double kn_el = equiv_young * calculation_area / initial_dist;

        if (&element1_props == &element2_props) {
             mDamageMaxDisplacementFactor = element1_props[DAMAGE_FACTOR];
             mTensionLimit = element1_props[CONTACT_SIGMA_MIN];
        } else {
            mDamageMaxDisplacementFactor = 0.5*(element1_props[DAMAGE_FACTOR] + element2_props[DAMAGE_FACTOR]);
            mTensionLimit = 0.5*(element1_props[CONTACT_SIGMA_MIN] + element2_props[CONTACT_SIGMA_MIN]);
        }

        const double Ntstr_el = mTensionLimit * calculation_area;
        double u1 = Ntstr_el / kn_el;
        double u2 = u1 * (1 + mDamageMaxDisplacementFactor)*10;
        return u2;

    }



    void DEM_Dempack::CalculateElasticConstants(double &kn_el,
                                                double &kt_el,
                                                double initial_dist,
                                                double equiv_young,
                                                double equiv_poisson,
                                                double calculation_area,
                                                SphericContinuumParticle* element1,
                                                SphericContinuumParticle* element2) {

        KRATOS_TRY
        double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));
        kn_el = equiv_young * calculation_area / initial_dist;
        kt_el = equiv_shear * calculation_area / initial_dist;
        KRATOS_CATCH("")
    }

    void DEM_Dempack::CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                                                double &equiv_visco_damp_coeff_tangential,
                                                SphericContinuumParticle* element1,
                                                SphericContinuumParticle* element2,
                                                const double kn_el,
                                                const double kt_el) {

        KRATOS_TRY
        double aux_norm_to_tang = 0.0;               // sqrt(kt_el / kn_el);
        const double mRealMass = element1->GetMass();
        const double &other_real_mass = element2->GetMass();
        const double mCoefficientOfRestitution = element1->GetProperties()[COEFFICIENT_OF_RESTITUTION];
        const double mOtherCoefficientOfRestitution = element2->GetProperties()[COEFFICIENT_OF_RESTITUTION];
        const double equiv_coefficientOfRestitution = 0.5 * (mCoefficientOfRestitution + mOtherCoefficientOfRestitution);

        equiv_visco_damp_coeff_normal = (1-equiv_coefficientOfRestitution) * 2.0 * sqrt(kn_el / (mRealMass + other_real_mass)) * (sqrt(mRealMass * other_real_mass)); // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )
        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; // not used in Dempack
        KRATOS_CATCH("")
    }

    void DEM_Dempack::CalculateForces(const ProcessInfo& r_process_info,
                                    double OldLocalElasticContactForce[3],
                                    double LocalElasticContactForce[3],
                                    double LocalElasticExtraContactForce[3],
                                    double LocalCoordSystem[3][3],
                                    double LocalDeltDisp[3],
                                    const double kn_el,
                                    const double kt_el,
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
                                    int search_control,
                                    DenseVector<int>& search_control_vector,
                                    double &equiv_visco_damp_coeff_normal,
                                    double &equiv_visco_damp_coeff_tangential,
                                    double LocalRelVel[3],
                                    double ViscoDampingLocalContactForce[3]) {

        KRATOS_TRY

        CalculateNormalForces(LocalElasticContactForce,
                kn_el,
                equiv_young,
                indentation,
                calculation_area,
                acumulated_damage,
                element1,
                element2,
                i_neighbour_count,
                time_steps,
                r_process_info);

        CalculateTangentialForces(OldLocalElasticContactForce,
                LocalElasticContactForce,
                LocalElasticExtraContactForce,
                LocalCoordSystem,
                LocalDeltDisp,
                kt_el,
                equiv_shear,
                contact_sigma,
                contact_tau,
                indentation,
                calculation_area,
                failure_criterion_state,
                element1,
                element2,
                i_neighbour_count,
                sliding,
                search_control,
                search_control_vector,
                r_process_info);

        CalculateViscoDampingCoeff(equiv_visco_damp_coeff_normal,
                                   equiv_visco_damp_coeff_tangential,
                                   element1,
                                   element2,
                                   kn_el,
                                   kt_el);

        CalculateViscoDamping(LocalRelVel,
                              ViscoDampingLocalContactForce,
                              indentation,
                              equiv_visco_damp_coeff_normal,
                              equiv_visco_damp_coeff_tangential,
                              sliding,
                              element1->mIniNeighbourFailureId[i_neighbour_count]);

        KRATOS_CATCH("")
    }

    void DEM_Dempack::CalculateNormalForces(double LocalElasticContactForce[3],
            const double kn_el,
            double equiv_young,
            double indentation,
            double calculation_area,
            double& acumulated_damage,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            int time_steps,
            const ProcessInfo& r_process_info) {

        KRATOS_TRY

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();
        double mN1;
        double mN2;
        double mN3;
        double mC1;
        double mC2;
        double mC3;
        double mYoungPlastic;
        double mPlasticityLimit;
        double mDamageMaxDisplacementFactor;
        double mTensionLimit;


        if (&element1_props == &element2_props) {

             mN1 = element1_props[SLOPE_FRACTION_N1];
             mN2 = element1_props[SLOPE_FRACTION_N2];
             mN3 = element1_props[SLOPE_FRACTION_N3];
             mC1 = element1_props[SLOPE_LIMIT_COEFF_C1];
             mC2 = element1_props[SLOPE_LIMIT_COEFF_C2];
             mC3 = element1_props[SLOPE_LIMIT_COEFF_C3];
             mYoungPlastic = element1_props[YOUNG_MODULUS_PLASTIC];
             mPlasticityLimit = element1_props[PLASTIC_YIELD_STRESS];
             mDamageMaxDisplacementFactor = element1_props[DAMAGE_FACTOR];
             mTensionLimit = element1_props[CONTACT_SIGMA_MIN]; //N/m2
        } else {

            mN1 = 0.5*(element1_props[SLOPE_FRACTION_N1] + element2_props[SLOPE_FRACTION_N1] );
            mN2 = 0.5*(element1_props[SLOPE_FRACTION_N2] + element2_props[SLOPE_FRACTION_N2] );
            mN3 = 0.5*(element1_props[SLOPE_FRACTION_N3] + element2_props[SLOPE_FRACTION_N3] );
            mC1 = 0.5*(element1_props[SLOPE_LIMIT_COEFF_C1] + element2_props[SLOPE_LIMIT_COEFF_C1]);
            mC2 = 0.5*(element1_props[SLOPE_LIMIT_COEFF_C2] + element2_props[SLOPE_LIMIT_COEFF_C2]);
            mC3 = 0.5*(element1_props[SLOPE_LIMIT_COEFF_C3] + element2_props[SLOPE_LIMIT_COEFF_C3]);
            mYoungPlastic = 0.5*(element1_props[YOUNG_MODULUS_PLASTIC] + element2_props[YOUNG_MODULUS_PLASTIC]);
            mPlasticityLimit = 0.5*(element1_props[PLASTIC_YIELD_STRESS] + element2_props[PLASTIC_YIELD_STRESS]);
            mDamageMaxDisplacementFactor = 0.5*(element1_props[DAMAGE_FACTOR] + element2_props[DAMAGE_FACTOR]);
            mTensionLimit = 0.5*(element1_props[CONTACT_SIGMA_MIN] + element2_props[CONTACT_SIGMA_MIN]); //N/m2
        }

        const double kn_b = kn_el / mN1;
        const double kn_c = kn_el / mN2;
        const double kn_d = kn_el / mN3;
        const double kp_el = mYoungPlastic / equiv_young * kn_el;
        const double Yields_el = mPlasticityLimit * calculation_area;

        const double Ncstr1_el = mC1 * calculation_area;
        const double Ncstr2_el = mC2 * calculation_area;
        const double Ncstr3_el = mC3 * calculation_area;
        const double Ntstr_el = mTensionLimit * calculation_area;
        double u_max = mHistoryMaxInd;

        double& fn = LocalElasticContactForce[2]; //[2] means 'normal' contact force

        if (indentation >= 0.0) { //COMPRESSION

            fn = kn_el * indentation;
            double u_ela1 = Ncstr1_el / kn_el;
            double u_ela2 = u_ela1 + (Ncstr2_el - Ncstr1_el) / (kn_b);
            double u_ela3 = u_ela2 + (Ncstr3_el - Ncstr2_el) / (kn_c);

            if ((indentation > u_max) || (time_steps <= 1)) {   //maximum historical indentation OR first step

                mHistoryMaxInd = indentation;                   //Saves maximum historical indentation

                if (indentation > u_ela3) {                     //4th stage
                    fn = Ncstr3_el + (indentation - u_ela3) * kn_d;
                    mHistoryDegradation = kn_d / kn_el;
                } else if (indentation > u_ela2) {              //3rd stage
                    fn = Ncstr2_el + (indentation - u_ela2) * kn_c;
                    mHistoryDegradation = kn_c / kn_el;
                } else {
                    if (indentation > u_ela1) {                 //2nd stage
                        fn = Ncstr1_el + (indentation - u_ela1) * kn_b;
                        mHistoryDegradation = kn_b / kn_el;
                    }
                }
                mHistoryMaxForce = fn;                          //Update max compressive force.
            } else {                                            //Below max.
                if (mHistoryMaxForce > 0.0) {                   //max compressive force

                    double u_plas;
                    if (Yields_el <= Ncstr1_el) {               //Plastic point located in the first stage.

                        u_plas = Yields_el / kn_el;
                    } else {
                        if (Yields_el <= Ncstr2_el) {           //Plastic point located in the second stage
                            u_plas = u_ela1 + (Yields_el - Ncstr1_el) / (kn_b);
                        } else if (Yields_el <= Ncstr3_el) {    //Plastic point located in the third stage
                            u_plas = u_ela2 + (Yields_el - Ncstr2_el) / (kn_c);
                        } else {                                //Plastic point located in the fourth stage
                            u_plas = u_ela3 + (Yields_el - Ncstr3_el) / (kn_d);
                        }
                    }
                    if (u_plas < u_max) {                       //below max BUT already in plastic zone
                        fn = mHistoryMaxForce - kp_el * (u_max - indentation); // Plastic unloading zone
                        mHistoryDegradation = kp_el / kn_el;
                    } else {                                    // Elastic unloading zone
                        if (indentation > u_ela3) {             // 4th stage
                            fn = Ncstr3_el + (indentation - u_ela3) * kn_d;
                        } else if (indentation > u_ela2) {      // 3rd stage
                            fn = Ncstr2_el + (indentation - u_ela2) * kn_c;
                        } else {
                            if (indentation > u_ela1) {         // 2nd stage
                                fn = Ncstr1_el + (indentation - u_ela1) * kn_b;
                            }
                        }
                    }
                }       //Compressive preloading.
            }           //Below max value.
        }               //Compression
        else {          //Tensile
            fn = kn_el * indentation;
            double u1 = Ntstr_el / kn_el;
            double u2 = u1 * (1 + mDamageMaxDisplacementFactor);

            if (failure_type == 0) {

                if (fabs(indentation) > u2) {                   // FULL DAMAGE
                    failure_type = 4;                           //tension failure
                    acumulated_damage = 1.0;
                    fn = 0.0;
                } else {
                    if (fabs(indentation) > u1) {
                        double u_frac = (fabs(indentation) - u1) / (u2 - u1);
                        //failure_criterion_state = fabs(indentation)/u2;
                        acumulated_damage = u_frac;

                        if (u_frac > mHistoryDamage) {
                            mHistoryDamage = u_frac;
                        }
                        double kn_damage = u1 / (fabs(indentation)) * kn_el * (1.0 - mHistoryDamage); // normal adhesive force (gap +)
                        fn = kn_damage * indentation;
                        //fn = indentation * kn_el*(1.0 -  mHistory[mapping_new_cont][2]);  // normal adhesive force (gap +)
                    }
                }
            }
            else {fn = 0.0;}
        } //Tension

    KRATOS_CATCH("")
    }



      /////// check old tangential forces calculations with unstable rotations




//     /*
    void DEM_Dempack::CalculateTangentialForces(double OldLocalElasticContactForce[3],
                                                double LocalElasticContactForce[3],
						double LocalElasticExtraContactForce[3],
                                                double LocalCoordSystem[3][3],
                                                double LocalDeltDisp[3],
                                                const double kt_el,
                                                const double equiv_shear,
                                                double& contact_sigma,
                                                double& contact_tau,
                                                double indentation,
                                                double calculation_area,
                                                double& failure_criterion_state,
                                                SphericContinuumParticle* element1,
                                                SphericContinuumParticle* element2,
                                                int i_neighbour_count,
                                                bool& sliding,
                                                int search_control,
                                                DenseVector<int>& search_control_vector,
                                                const ProcessInfo& r_process_info) {


        KRATOS_TRY

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();

        double mTensionLimit;
        double mTauZero;
        double mInternalFriccion;
        double mShearEnergyCoef;

        if (&element1_props == &element2_props) {

            mTensionLimit = element1_props[CONTACT_SIGMA_MIN];
            mTauZero = element1_props[CONTACT_TAU_ZERO];
            mInternalFriccion = element1_props[CONTACT_INTERNAL_FRICC];
            mShearEnergyCoef = element1_props[SHEAR_ENERGY_COEF];
        }
        else{

            mTensionLimit = 0.5*(element1_props[CONTACT_SIGMA_MIN] + element2_props[CONTACT_SIGMA_MIN]);
            mTauZero = 0.5*(element1_props[CONTACT_TAU_ZERO] + element2_props[CONTACT_TAU_ZERO]);
            mInternalFriccion = 0.5*(element1_props[CONTACT_INTERNAL_FRICC] + element2_props[CONTACT_INTERNAL_FRICC]);
            mShearEnergyCoef = 0.5*(element1_props[SHEAR_ENERGY_COEF] + element2_props[SHEAR_ENERGY_COEF]);
        }

        double degradation = 1.0;               //Tangential. With degradation:

        if (i_neighbour_count < int(element1->mContinuumInitialNeighborsSize)) {
            if (indentation >= 0.0) {           //COMPRESSION
                degradation = mHistoryDegradation;
            } else {
                degradation = (1.0 - mHistoryDamage);
            }
        }

        if (failure_type == 0) {                //This means it has not broken

            if (mHistoryShearFlag == 0.0) {

                LocalElasticContactForce[0] += -degradation * kt_el * LocalDeltDisp[0]; // 0: first tangential
                LocalElasticContactForce[1] += -degradation * kt_el * LocalDeltDisp[1]; // 1: second tangential
            }

            double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                    + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

            contact_tau = ShearForceNow / calculation_area;
            contact_sigma = LocalElasticContactForce[2] / calculation_area;

            double tau_strength = mTauZero;

            if (contact_sigma >= 0) {
                tau_strength = mTauZero + mInternalFriccion * contact_sigma;
            }

            if (contact_tau > tau_strength) {
                mHistoryShearFlag = 1.0;
            }

            if (mHistoryShearFlag != 0.0) {

                double increment_disp = sqrt(LocalDeltDisp[0] * LocalDeltDisp[0] + LocalDeltDisp[1] * LocalDeltDisp[1]);
                mHistoryDisp += increment_disp;

                double u1_tau = tau_strength * calculation_area / kt_el;
                double damage_tau = 1.0;

                if (mShearEnergyCoef != 0.0) {
                    damage_tau = mHistoryDisp / (u1_tau * (mShearEnergyCoef));
                }

                double aux = (1.0 - damage_tau) * (tau_strength / contact_tau);
                LocalElasticContactForce[0] = aux * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = aux * LocalElasticContactForce[1];

                //failure_criterion_state = 1.0;                        // *under revision
                failure_criterion_state = (1+mShearEnergyCoef*damage_tau) /(1+mShearEnergyCoef);  // *under revision
                if(contact_sigma<0){ failure_criterion_state = GeometryFunctions::max( failure_criterion_state, -contact_sigma / mTensionLimit ); }   // *revisar

                if (damage_tau >= 1.0) {
                    failure_type = 2;                   // shear
                    failure_criterion_state = 1.0;      // *under revision
                    sliding = true;
                }
            } else {

                if (contact_sigma < 0) {
                    failure_criterion_state = GeometryFunctions::max(contact_tau / tau_strength, -contact_sigma / mTensionLimit);
                } else failure_criterion_state = contact_tau / tau_strength;
                if (failure_criterion_state > 1.0) failure_criterion_state = 1.0;
            }
        }
        if (search_control == 0) {
            if (failure_type != 0) {
                search_control_vector[OpenMPUtils::ThisThread()] = 1;
            }
        }
    KRATOS_CATCH("")
    }

    void DEM_Dempack::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                    SphericContinuumParticle* neighbor,
                                                    double equiv_young,
                                                    double distance,
                                                    double calculation_area,
                                                    double LocalCoordSystem[3][3],
                                                    double ElasticLocalRotationalMoment[3],
                                                    double ViscoLocalRotationalMoment[3],
                                                    double equiv_poisson,
                                                    double indentation) {}  //ComputeParticleRotationalMoments





} /* namespace Kratos.*/
