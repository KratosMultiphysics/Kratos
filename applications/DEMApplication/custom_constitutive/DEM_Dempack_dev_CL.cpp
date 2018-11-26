// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_Dempack_dev_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_Dempack_dev::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_Dempack_dev(*this));
        return p_clone;
    }

    void DEM_Dempack_dev::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_Dempack_dev to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }


    void DEM_Dempack_dev::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

        double rmin = radius;
        if (other_radius < radius) rmin = other_radius;
        calculation_area = Globals::Pi * rmin*rmin;
    }
    double DEM_Dempack_dev::CalculateContactArea(double radius, double other_radius, Vector& v) {
        double a = 0.0;
        CalculateContactArea(radius, other_radius, a);
        unsigned int old_size = v.size();
        v.resize(old_size + 1);
        v[old_size]=a;
        return a;
    }

    void DEM_Dempack_dev::GetContactArea(const double radius,
                                     const double other_radius,
                                     const Vector& vector_of_initial_areas,
                                     const int neighbour_position,
                                     double& calculation_area)
    {
        if (vector_of_initial_areas.size()) calculation_area = vector_of_initial_areas[neighbour_position];
        else CalculateContactArea(radius, other_radius, calculation_area);
        //CalculateContactArea(radius, other_radius, calculation_area);
    }

    void DEM_Dempack_dev::CalculateElasticConstants(double &kn_el,
                                                double &kt_el,
                                                double initial_dist,
                                                double equiv_young,
                                                double equiv_poisson,
                                                double calculation_area,
                                                SphericContinuumParticle* element1,
                                                SphericContinuumParticle* element2) {

        KRATOS_TRY

//        double rmin = element1->GetRadius();    // test rebalance solo resistencia
//        const double other_radius = element2->GetRadius();
//        if (other_radius < rmin) rmin = other_radius;
//        double effective_calculation_area = Globals::Pi * rmin*rmin;

        double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));
        kn_el = equiv_young * calculation_area / initial_dist;
        kt_el = equiv_shear * calculation_area / initial_dist;

//        std::ofstream outputfile("knkt.txt", std::ios_base::out | std::ios_base::app);
//        outputfile << kn_el << " " << kt_el << "\n";
//        outputfile.close();

        KRATOS_CATCH("")
    }



    void DEM_Dempack_dev::CalculateNormalForces(double LocalElasticContactForce[3],
            const double kn_el,
            double equiv_young,
            double indentation,
            double calculation_area,
            double& acumulated_damage,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            int time_steps) {

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


//        double rmin = element1->GetRadius();
//        const double other_radius = element2->GetRadius();
//        if (other_radius < rmin) rmin = other_radius;
//        double effective_calculation_area = Globals::Pi * rmin*rmin;
          double effective_calculation_area = calculation_area;

        if (&element1_props == &element2_props) {

             mN1 = element1_props[SLOPE_FRACTION_N1];
             mN2 = element1_props[SLOPE_FRACTION_N2];
             mN3 = element1_props[SLOPE_FRACTION_N3];
             mC1 = element1_props[SLOPE_LIMIT_COEFF_C1]*1e6;
             mC2 = element1_props[SLOPE_LIMIT_COEFF_C2]*1e6;
             mC3 = element1_props[SLOPE_LIMIT_COEFF_C3]*1e6;
             mYoungPlastic = element1_props[YOUNG_MODULUS_PLASTIC];
             mPlasticityLimit = element1_props[PLASTIC_YIELD_STRESS]*1e6;
             mDamageMaxDisplacementFactor = element1_props[DAMAGE_FACTOR];
             mTensionLimit = element1_props[CONTACT_SIGMA_MIN]*1e6; //N/m2
        } else {

            mN1 = 0.5*(element1_props[SLOPE_FRACTION_N1] + element2_props[SLOPE_FRACTION_N1] );
            mN2 = 0.5*(element1_props[SLOPE_FRACTION_N2] + element2_props[SLOPE_FRACTION_N2] );
            mN3 = 0.5*(element1_props[SLOPE_FRACTION_N3] + element2_props[SLOPE_FRACTION_N3] );
            mC1 = 0.5*1e6*(element1_props[SLOPE_LIMIT_COEFF_C1] + element2_props[SLOPE_LIMIT_COEFF_C1]);
            mC2 = 0.5*1e6*(element1_props[SLOPE_LIMIT_COEFF_C2] + element2_props[SLOPE_LIMIT_COEFF_C2]);
            mC3 = 0.5*1e6*(element1_props[SLOPE_LIMIT_COEFF_C3] + element2_props[SLOPE_LIMIT_COEFF_C3]);
            mYoungPlastic = 0.5*(element1_props[YOUNG_MODULUS_PLASTIC] + element2_props[YOUNG_MODULUS_PLASTIC]);
            mPlasticityLimit = 0.5*1e6*(element1_props[PLASTIC_YIELD_STRESS] + element2_props[PLASTIC_YIELD_STRESS]);
            mDamageMaxDisplacementFactor = 0.5*(element1_props[DAMAGE_FACTOR] + element2_props[DAMAGE_FACTOR]);
            mTensionLimit = 0.5*1e6*(element1_props[CONTACT_SIGMA_MIN] + element2_props[CONTACT_SIGMA_MIN]); //N/m2
        }

        const double kn_b = kn_el / mN1;
        const double kn_c = kn_el / mN2;
        const double kn_d = kn_el / mN3;
        const double kp_el = mYoungPlastic / equiv_young * kn_el;
        const double Yields_el = mPlasticityLimit * effective_calculation_area;

        const double Ncstr1_el = mC1 * effective_calculation_area;
        const double Ncstr2_el = mC2 * effective_calculation_area;
        const double Ncstr3_el = mC3 * effective_calculation_area;
        const double Ntstr_el = mTensionLimit * effective_calculation_area;
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
/*
    void DEM_Dempack_dev::CalculateTangentialForces(double OldLocalElasticContactForce[3],
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

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

//        double rmin = element1->GetRadius();
//        const double other_radius = element2->GetRadius();
//        if (other_radius < rmin) rmin = other_radius;
//        double effective_calculation_area = Globals::Pi * rmin*rmin;
        double effective_calculation_area = calculation_area;

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();

        double mTensionLimit;
        double mTauZero;
        double mInternalFriccion;
        double mShearEnergyCoef;

        if (&element1_props == &element2_props) {

            mTensionLimit = element1_props[CONTACT_SIGMA_MIN]*1e6;
            mTauZero = element1_props[CONTACT_TAU_ZERO]*1e6;
            mInternalFriccion = element1_props[CONTACT_INTERNAL_FRICC];
            mShearEnergyCoef = element1_props[SHEAR_ENERGY_COEF];
        }
        else{

            mTensionLimit = 0.5*1e6*(element1_props[CONTACT_SIGMA_MIN] + element2_props[CONTACT_SIGMA_MIN]);
            mTauZero = 0.5*1e6*(element1_props[CONTACT_TAU_ZERO] + element2_props[CONTACT_TAU_ZERO]);
            mInternalFriccion = 0.5*(element1_props[CONTACT_INTERNAL_FRICC] + element2_props[CONTACT_INTERNAL_FRICC]);
            mShearEnergyCoef = 0.5*(element1_props[SHEAR_ENERGY_COEF] + element2_props[SHEAR_ENERGY_COEF]);
        }

        double degradation = 1.0; //Tangential. With degradation:

        if (i_neighbour_count < int(element1->mContinuumInitialNeighborsSize)) {
            if (indentation >= 0.0) { //COMPRESSION
                degradation = mHistoryDegradation;
            } else {
                degradation = (1.0 - mHistoryDamage);
            }
        }

        if (failure_type == 0) { //This means it has not broken
            if (r_process_info[SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION]) {
                AddContributionOfShearStrainParallelToBond(OldLocalElasticContactForce, LocalElasticExtraContactForce, element1->mNeighbourElasticExtraContactForces[i_neighbour_count], LocalCoordSystem, kt_el, calculation_area,  element1, element2);
            }

            if (mHistoryShearFlag == 0.0) {

                LocalElasticContactForce[0] += -degradation * kt_el * LocalDeltDisp[0]; // 0: first tangential
                LocalElasticContactForce[1] += -degradation * kt_el * LocalDeltDisp[1]; // 1: second tangential
            }

            double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                    + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

            contact_tau = ShearForceNow / effective_calculation_area;
            contact_sigma = LocalElasticContactForce[2] / effective_calculation_area;

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

                double u1_tau = tau_strength * effective_calculation_area / kt_el;
                double damage_tau = 1.0;

                if (mShearEnergyCoef != 0.0) {
                    damage_tau = mHistoryDisp / (u1_tau * (mShearEnergyCoef));
                }

                double aux = (1.0 - damage_tau) * (tau_strength / contact_tau);

                LocalElasticContactForce[0] = aux * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = aux * LocalElasticContactForce[1];

                failure_criterion_state = 1.0;

                //failure_criterion_state = (1+mShearEnergyCoef*damage_tau) /(1+mShearEnergyCoef);
                //if(contact_sigma<0){ failure_criterion_state = GeometryFunctions::max( failure_criterion_state, -contact_sigma / mTensionLimit ); }

                if (damage_tau >= 1.0) {
                    failure_type = 2; // shear
                    //failure_criterion_state = 1.0;    //
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
    }
*/



    void DEM_Dempack_dev::CalculateTangentialForces(double OldLocalElasticContactForce[3],
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
        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_el * LocalDeltDisp[0]; // 0: first tangential
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_el * LocalDeltDisp[1]; // 1: second tangential

        double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                  + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

//        double rmin = element1->GetRadius();  // test de rebalance de areas solo en resistencia
//        const double other_radius = element2->GetRadius();
//        if (other_radius < rmin) rmin = other_radius;
//        double effective_calculation_area = Globals::Pi * rmin*rmin;
        double effective_calculation_area = calculation_area;


        if (failure_type == 0) { // This means it has not broken
            if (r_process_info[SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION]) {
                AddContributionOfShearStrainParallelToBond(OldLocalElasticContactForce, LocalElasticExtraContactForce, element1->mNeighbourElasticExtraContactForces[i_neighbour_count], LocalCoordSystem, kt_el, calculation_area,  element1, element2);
            }
            Properties& element1_props = element1->GetProperties();
            Properties& element2_props = element2->GetProperties();
            const double mTauZero = 0.5 * 1e6 * (element1_props[CONTACT_TAU_ZERO] + element2_props[CONTACT_TAU_ZERO]);
            const double mInternalFriction = 0.5 * (element1_props[CONTACT_INTERNAL_FRICC] + element2_props[CONTACT_INTERNAL_FRICC]);

            contact_tau = ShearForceNow / effective_calculation_area;
            contact_sigma = LocalElasticContactForce[2] / effective_calculation_area;

            double tau_strength = mTauZero;

            if (contact_sigma >= 0) {
                tau_strength = mTauZero + mInternalFriction * contact_sigma;
            }

            if (contact_tau > tau_strength) {
                failure_type = 2; // shear
            }
        }
        else {
            const double equiv_tg_of_fri_ang = 0.5 * (element1->GetTgOfFrictionAngle() + element2->GetTgOfFrictionAngle());
            double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];
            if (Frictional_ShearForceMax < 0.0) {
                Frictional_ShearForceMax = 0.0;
            }

            if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)) {
                LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[1];
                sliding = true;
            }
        }

        KRATOS_CATCH("")
    }


    void DEM_Dempack_dev::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
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
        //GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, mContactMoment, LocalRotationalMoment);

        const double equivalent_radius = sqrt(calculation_area / Globals::Pi);
        const double Inertia_I = 0.25 * Globals::Pi * equivalent_radius * equivalent_radius * equivalent_radius * equivalent_radius;
        const double Inertia_J = 2.0 * Inertia_I; // This is the polar inertia
        const double debugging_rotational_factor = 5.0; //1.0; // Hardcoded only for testing purposes. Obviously, this parameter should be always 1.0

        const double element_mass  = element->GetMass();
        const double neighbor_mass = neighbor->GetMass();
        const double equiv_mass    = element_mass * neighbor_mass / (element_mass + neighbor_mass);

        // Viscous parameter taken from J.S.Marshall, 'Discrete-element modeling of particle aerosol flows', section 4.3. Twisting resistance
        const double alpha = 0.9; // TODO: Hardcoded only for testing purposes. This value depends on the restitution coefficient and goes from 0.1 to 1.0
        const double visc_param = 0.5 * equivalent_radius * equivalent_radius * alpha * sqrt(1.33333333333333333 * equiv_mass * equiv_young * equivalent_radius);

        //equiv_young or G in torsor (LocalRotationalMoment[2]) ///////// TODO
        ElasticLocalRotationalMoment[0] = -debugging_rotational_factor * equiv_young * Inertia_I * LocalDeltaRotatedAngle[0] / distance;
        ///- debugging_rotational_factor * equiv_shear * (calculation_area / distance) * (OtherWeightedRadius * MyLocalDeltaDisplacement[0] - MyWeightedRadius * OtherLocalDeltaDisplacement[0]);

        ElasticLocalRotationalMoment[1] = -debugging_rotational_factor * equiv_young * Inertia_I * LocalDeltaRotatedAngle[1] / distance;
        ///- debugging_rotational_factor * equiv_shear * (calculation_area / distance) * (OtherWeightedRadius * MyLocalDeltaDisplacement[1] - MyWeightedRadius * OtherLocalDeltaDisplacement[1]);

        ElasticLocalRotationalMoment[2] = -debugging_rotational_factor * equiv_young * Inertia_J * LocalDeltaRotatedAngle[2] / distance;

        ViscoLocalRotationalMoment[0] = -visc_param * LocalDeltaAngularVelocity[0];
        ViscoLocalRotationalMoment[1] = -visc_param * LocalDeltaAngularVelocity[1];
        ViscoLocalRotationalMoment[2] = -visc_param * LocalDeltaAngularVelocity[2];

        // TODO: Judge if the rotation spring is broken or not
        /*
        double ForceN  = LocalElasticContactForce[2];
        double ForceS  = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
        double MomentS = sqrt(LocalRotaSpringMoment[0] * LocalRotaSpringMoment[0] + LocalRotaSpringMoment[1] * LocalRotaSpringMoment[1]);
        double MomentN = LocalRotaSpringMoment[2];
        // bending stress and axial stress add together, use edge of the bar will failure first
        double TensiMax = -ForceN / calculation_area + MomentS / Inertia_I * equiv_radius;
        double ShearMax =  ForceS / calculation_area + fabs(MomentN) / Inertia_J * equiv_radius;
        if (TensiMax > equiv_tension || ShearMax > equiv_cohesion) {
            mRotaSpringFailureType[i_neighbor_count] = 1;
            LocalRotaSpringMoment[0] = LocalRotaSpringMoment[1] = LocalRotaSpringMoment[2] = 0.0;
            //LocalRotaSpringMoment[1] = 0.0;
            //LocalRotaSpringMoment[2] = 0.0;
        }
        */
        //GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRotationalMoment, mContactMoment);
        KRATOS_CATCH("")
    }//ComputeParticleRotationalMoments


    void DEM_Dempack_dev::AddPoissonContribution(const double equiv_poisson, double LocalCoordSystem[3][3], double& normal_force,
                                          double calculation_area, Matrix* mSymmStressTensor, SphericContinuumParticle* element1,
                                          SphericContinuumParticle* element2, const ProcessInfo& r_process_info, const int i_neighbor_count, const double indentation) {

        if (!r_process_info[POISSON_EFFECT_OPTION]) return;
        if (element1->mIniNeighbourFailureId[i_neighbor_count] > 0  &&  indentation < 0.0) return;

        double force[3];
        Matrix average_stress_tensor = ZeroMatrix(3,3);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                average_stress_tensor(i,j) = 0.5 * ((*mSymmStressTensor)(i,j) + (*(element2->mSymmStressTensor))(i,j));
            }
        }

        for (int i = 0; i < 3; i++) {

            force[i] = (average_stress_tensor)(i,0) * LocalCoordSystem[0][0] +
                       (average_stress_tensor)(i,1) * LocalCoordSystem[0][1] +
                       (average_stress_tensor)(i,2) * LocalCoordSystem[0][2]; // StressTensor*unitaryNormal0
        }

        double sigma_x = force[0] * LocalCoordSystem[0][0] +
                         force[1] * LocalCoordSystem[0][1] +
                         force[2] * LocalCoordSystem[0][2]; // projection to normal to obtain value of the normal stress

        for (int i = 0; i < 3; i++) {

            force[i] = (average_stress_tensor)(i,0) * LocalCoordSystem[1][0] +
                       (average_stress_tensor)(i,1) * LocalCoordSystem[1][1] +
                       (average_stress_tensor)(i,2) * LocalCoordSystem[1][2]; // StressTensor*unitaryNormal1
        }

        double sigma_y = force[0] * LocalCoordSystem[1][0] +
                         force[1] * LocalCoordSystem[1][1] +
                         force[2] * LocalCoordSystem[1][2]; // projection to normal to obtain value of the normal stress

        double poisson_force = calculation_area * equiv_poisson * (sigma_x + sigma_y);

        normal_force -= poisson_force;

    } //AddPoissonContribution

    void DEM_Dempack_dev::AddContributionOfShearStrainParallelToBond(double OldLocalElasticContactForce[3],
                                                              double LocalElasticExtraContactForce[3],
                                                              array_1d<double, 3>& OldElasticExtraContactForce,
                                                              double LocalCoordSystem[3][3],
                                                              const double kt_el,
                                                              const double calculation_area,
                                                              SphericContinuumParticle* element1,
                                                              SphericContinuumParticle* element2) {

        if (element1->mSymmStressTensor == NULL) return;
        //if(element1->IsSkin() || element2->IsSkin()) return;

        double average_stress_tensor[3][3];

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                average_stress_tensor[i][j]     = 0.5 * ((*(element1->mSymmStressTensor   ))(i,j) + (*(element2->mSymmStressTensor   ))(i,j));
            }
        }

        double current_sigma_local[3][3];
        GeometryFunctions::TensorGlobal2Local(LocalCoordSystem, average_stress_tensor, current_sigma_local);

        array_1d<double, 3> OldLocalElasticExtraContactForce;
        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, OldElasticExtraContactForce, OldLocalElasticExtraContactForce);

        double force_due_to_stress0 = calculation_area * current_sigma_local[0][2];
        double force_due_to_stress1 = calculation_area * current_sigma_local[1][2];

        LocalElasticExtraContactForce[0] = -OldLocalElasticContactForce[0] - force_due_to_stress0;
        LocalElasticExtraContactForce[1] = -OldLocalElasticContactForce[1] - force_due_to_stress1;

        if(fabs(LocalElasticExtraContactForce[0]) > fabs(force_due_to_stress0) ) {
            LocalElasticExtraContactForce[0] = LocalElasticExtraContactForce[0] / fabs(LocalElasticExtraContactForce[0]) * fabs(force_due_to_stress0);
        }
        if(fabs(LocalElasticExtraContactForce[1]) > fabs(force_due_to_stress1) ) {
            LocalElasticExtraContactForce[1] = LocalElasticExtraContactForce[1] / fabs(LocalElasticExtraContactForce[1]) * fabs(force_due_to_stress1);
        }
    }

} /* namespace Kratos.*/
