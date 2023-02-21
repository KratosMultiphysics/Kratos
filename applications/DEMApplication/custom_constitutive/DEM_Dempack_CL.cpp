// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_Dempack_CL.h"
//#include "custom_elements/spheric_particle.h"
#include "custom_elements/spheric_continuum_particle.h"
#include "dem_contact.h"

namespace Kratos {

    void DEM_Dempack::Initialize(SphericContinuumParticle* element1, SphericContinuumParticle* element2, Properties::Pointer pProps) {

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


    void DEM_Dempack::TransferParametersToProperties(const Parameters& parameters, Properties::Pointer pProp)  {
        BaseClassType::TransferParametersToProperties(parameters, pProp);
        pProp->SetValue(SLOPE_FRACTION_N1, parameters["SLOPE_FRACTION_N1"].GetDouble());
        pProp->SetValue(SLOPE_FRACTION_N2, parameters["SLOPE_FRACTION_N2"].GetBool());
        pProp->SetValue(SLOPE_FRACTION_N3, parameters["SLOPE_FRACTION_N3"].GetDouble());
        pProp->SetValue(SLOPE_LIMIT_COEFF_C1, parameters["SLOPE_LIMIT_COEFF_C1"].GetDouble());
        pProp->SetValue(SLOPE_LIMIT_COEFF_C2, parameters["SLOPE_LIMIT_COEFF_C2"].GetDouble());
        pProp->SetValue(SLOPE_LIMIT_COEFF_C3, parameters["SLOPE_LIMIT_COEFF_C3"].GetDouble());
        pProp->SetValue(SLOPE_LIMIT_COEFF_C3, parameters["YOUNG_MODULUS_PLASTIC"].GetDouble());
        pProp->SetValue(SLOPE_LIMIT_COEFF_C3, parameters["PLASTIC_YIELD_STRESS"].GetDouble());
        pProp->SetValue(SLOPE_LIMIT_COEFF_C3, parameters["DAMAGE_FACTOR"].GetDouble());
        pProp->SetValue(SLOPE_LIMIT_COEFF_C3, parameters["CONTACT_SIGMA_MIN"].GetDouble());
        pProp->SetValue(SLOPE_LIMIT_COEFF_C3, parameters["CONTACT_TAU_ZERO"].GetDouble());
        pProp->SetValue(SLOPE_LIMIT_COEFF_C3, parameters["CONTACT_INTERNAL_FRICC"].GetDouble());
        pProp->SetValue(SLOPE_LIMIT_COEFF_C3, parameters["SHEAR_ENERGY_COEF"].GetDouble());
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

        mDamageMaxDisplacementFactor = (*mpProperties)[DAMAGE_FACTOR];
        mTensionLimit = (*mpProperties)[CONTACT_SIGMA_MIN];

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
                                                SphericContinuumParticle* element2, double indentation) {

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
        const double& equiv_coefficientOfRestitution = (*mpProperties)[COEFFICIENT_OF_RESTITUTION];

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
                ViscoDampingLocalContactForce,
                LocalCoordSystem,
                LocalDeltDisp,
                LocalRelVel,
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

        mN1 = (*mpProperties)[SLOPE_FRACTION_N1];
        mN2 = (*mpProperties)[SLOPE_FRACTION_N2];
        mN3 = (*mpProperties)[SLOPE_FRACTION_N3];
        mC1 = (*mpProperties)[SLOPE_LIMIT_COEFF_C1];
        mC2 = (*mpProperties)[SLOPE_LIMIT_COEFF_C2];
        mC3 = (*mpProperties)[SLOPE_LIMIT_COEFF_C3];
        mYoungPlastic = (*mpProperties)[YOUNG_MODULUS_PLASTIC];
        mPlasticityLimit = (*mpProperties)[PLASTIC_YIELD_STRESS];
        mDamageMaxDisplacementFactor = (*mpProperties)[DAMAGE_FACTOR];
        mTensionLimit = (*mpProperties)[CONTACT_SIGMA_MIN]; //N/m2

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


    void DEM_Dempack::Check(Properties::Pointer pProp) const {
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
                KRATOS_WARNING("DEM")<<"WARNING: Variable DYNAMIC_FRICTION or FRICTION should be present in the properties when using DEMContinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
                KRATOS_WARNING("DEM")<<std::endl;
                pProp->GetValue(DYNAMIC_FRICTION) = 0.0;
            }
            else {
                pProp->GetValue(DYNAMIC_FRICTION) = pProp->GetValue(FRICTION);
            }
        }
        if(!pProp->Has(FRICTION_DECAY)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable FRICTION_DECAY should be present in the properties when using DEMContinuumConstitutiveLaw. 500.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(FRICTION_DECAY) = 500.0;
        }
        if(!pProp->Has(COEFFICIENT_OF_RESTITUTION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable COEFFICIENT_OF_RESTITUTION should be present in the properties when using DEMContinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(COEFFICIENT_OF_RESTITUTION) = 0.0;
        }
    }


    void DEM_Dempack::CalculateTangentialForces(double OldLocalElasticContactForce[3],
                                                double LocalElasticContactForce[3],
						                        double LocalElasticExtraContactForce[3],
                                                double ViscoDampingLocalContactForce[3],
                                                double LocalCoordSystem[3][3],
                                                double LocalDeltDisp[3],
                                                double LocalRelVel[3],
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
                                                const ProcessInfo& r_process_info) {


        KRATOS_TRY

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

        double mTensionLimit;
        double mTauZero;
        double mInternalFriccion;
        double mShearEnergyCoef;

        mTensionLimit = (*mpProperties)[CONTACT_SIGMA_MIN];
        mTauZero = (*mpProperties)[CONTACT_TAU_ZERO];
        mInternalFriccion = (*mpProperties)[CONTACT_INTERNAL_FRICC];
        mShearEnergyCoef = (*mpProperties)[SHEAR_ENERGY_COEF];

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

    KRATOS_CATCH("")
    }

    void DEM_Dempack::CalculateMoments(SphericContinuumParticle* element, 
                    SphericContinuumParticle* neighbor, 
                    double equiv_young, 
                    double distance, 
                    double calculation_area,
                    double LocalCoordSystem[3][3], 
                    double ElasticLocalRotationalMoment[3], 
                    double ViscoLocalRotationalMoment[3], 
                    double equiv_poisson, 
                    double indentation, 
                    double LocalElasticContactForce[3],
                    double normalLocalContactForce,
                    double GlobalElasticContactForces[3],
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
                                        indentation, 
                                        LocalElasticContactForce);
        }             

        DemContact::ComputeParticleContactMoments(normalLocalContactForce,
                                                GlobalElasticContactForces,
                                                LocalCoordSystem_2,
                                                element,
                                                neighbor,
                                                indentation,
                                                i_neighbor_count);

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
                                                    double indentation,
                                                    double LocalElasticContactForce[3]) {}  //ComputeParticleRotationalMoments





} /* namespace Kratos.*/
