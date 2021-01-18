// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_ExponentialHC_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    void DEM_ExponentialHC::Initialize(SphericContinuumParticle* owner_sphere) {
        KRATOS_TRY
        mHistoryMaxInd              = 0.0; //maximum indentation achieved
        mHistoryMaxForce            = 0.0; //maximum force achieved
        mHistoryDamage              = 0.0; //cumulated_damage
        mHistoryDegradation         = 1.0; //degradation factor for G reducing in Dempack;

        mGamma1                     = 0.0;
        mGamma2                     = 0.0;
        mGamma3                     = 0.0;
        mMaxDef                     = 0.0;

        KRATOS_CATCH("")
    }

    DEMContinuumConstitutiveLaw::Pointer DEM_ExponentialHC::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_ExponentialHC(*this));
        return p_clone;
    }

    void DEM_ExponentialHC::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_ExponentialHC to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_ExponentialHC::CalculateNormalForces(double LocalElasticContactForce[3],
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

        int &mNeighbourFailureId_count = element1->mIniNeighbourFailureId[i_neighbour_count];

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();
        double mDamageMaxDisplacementFactor;
        double mTensionLimit;

        if(&element1_props == &element2_props ){


             mDamageMaxDisplacementFactor = element1_props[DAMAGE_FACTOR];
             mTensionLimit = element1_props[CONTACT_SIGMA_MIN]; //N/m2
        }

        else{


            mDamageMaxDisplacementFactor = 0.5*(element1_props[DAMAGE_FACTOR] + element2_props[DAMAGE_FACTOR]);
            mTensionLimit = 0.5*(element1_props[CONTACT_SIGMA_MIN] + element2_props[CONTACT_SIGMA_MIN]); //N/m2
        }

        mGamma1 = 0.2;
        mGamma2 = 16;
        mGamma3 = 0.275;
        mMaxDef = 0.002;

        double &mNeighbourDelta_count = element1->mIniNeighbourDelta[i_neighbour_count];

        const double &other_radius = element2->GetRadius();
        const double my_radius = element1->GetRadius();
        const double initial_delta = mNeighbourDelta_count;
        const double initial_dist = (other_radius + my_radius - initial_delta);
        double current_def = indentation / initial_dist;
        double u_max = mHistoryMaxInd;
        double kn_plas = kn_el;          // double kn_plas = mYoungPlastic / equiv_young * kn_el;
        double Ntstr_el = mTensionLimit * calculation_area;
        double &fn = LocalElasticContactForce[2];
        double kn_exp = kn_el * mGamma3 + kn_el * mGamma1 * exp(mGamma2 * (current_def - mMaxDef));

        if (kn_exp > kn_el) {
            kn_exp = kn_el;}

        if (indentation >= 0.0) {       // COMPRESSION STAGE
            fn = kn_el * indentation;   //  fuerza en parte lineal
            if ((indentation > u_max) || (time_steps <= 1)) { //maximum historical vemos si esta en carga, comparando la indentation actual con la maxima historica
                {
                    mHistoryMaxInd = indentation; // Guarda el threshold de la maxima indentation
                    if (indentation > initial_dist * mMaxDef) //// dempack_indentation >= C1*Area/kn_el   se supera el limite para el cambio de pendiente.
                    {
                        fn = kn_el * initial_dist * mMaxDef + kn_exp * (indentation - initial_dist * mMaxDef);
                    }
                }
                mHistoryMaxForce = fn;          //actualitzar la força màxima a compressió.
            } else {                            ////Per sota del màxim. esta en descarga, la distancia entre particulas aumenta respecto la historica, current_dist > mHistDist
                if (mHistoryMaxForce > 0.0) {   //Màxim en compressió.

                    double u_plas;              //  por ahora coincide con el cambio de pendiente
                    if (indentation <= initial_dist * mMaxDef) { // //si el punt de plastificació està en la primera rama elastica..
                        u_plas = indentation;
                    } else {u_plas = initial_dist * mMaxDef + fn/kn_exp;
                    }
                    if (u_plas < u_max) {      //si nosaltres estem per sota del maxim pero ja estem plastificant
                        fn = mHistoryMaxForce - kn_plas * (u_max - indentation); // Esta en zona de descarga plastica (pot estar en carga/descarga)
                        mHistoryDegradation = kn_plas / kn_el; // used in tangential degradation
                    } else {
                        if (indentation > initial_dist * mMaxDef) {  // por encima de la rama elastica
                            fn = kn_el * initial_dist * mMaxDef + kn_exp * (indentation - initial_dist * mMaxDef);
                        }
                    }
                }
            }
        } else {                                //tension   same as dempack atm
            fn = kn_el * indentation;
            double u1 = Ntstr_el / kn_el;
            double u2 = u1 * (1 + mDamageMaxDisplacementFactor);

            if (fabs(indentation) > u2) {       // FULL DAMAGE
                mNeighbourFailureId_count = 4;  //tension failure
                acumulated_damage = 1.0;
                fn = 0.0;
            } else {
                if (fabs(indentation) > u1) {
                    double u_frac = (fabs(indentation) - u1) / (u2 - u1);
                    acumulated_damage = u_frac;

                    if (u_frac > mHistoryDamage) {
                        mHistoryDamage = u_frac;
                    }
                    double kn_damage = u1 / (fabs(indentation)) * kn_el * (1.0 - mHistoryDamage);
                    fn = kn_damage * indentation;
                }
            }
        } //Tension
    }

} /* namespace Kratos.*/









///////////////////////////////////method scheme how to//////////////////////////////////////////
//void SphericContinuumParticle::NonlinearNormalForceCalculation(double LocalElasticContactForce[3],
//        double kn_el,
//        double distance,
//        double max_dist,
//        double initial_dist) {
//
//    //mGamma1 = rCurrentProcessInfo[DONZE_G1];
//    //mGamma2 = rCurrentProcessInfo[DONZE_G2];
//    //mGamma3 = rCurrentProcessInfo[DONZE_G3];
//    //mMaxDef = rCurrentProcessInfo[DONZE_MAX_DEF];
//
//    //sabemos que double indentation = initial_dist - distance;
//
//    double current_def = (initial_dist - distance) / initial_dist;
//    double kn1 = kn_el;
//    double kn2 = kn1 * mGamma3 + kn1 * mGamma1 * exp(mGamma2 * (current_def - mMaxDef));
//    if (kn2 > kn1) {
//        kn2 = kn1;
//    }
//    double max_dist = initial_dist * (1 - mMaxDef);
//
//    //initial_dist ???
//    double kn_plas = kn1; // modificable en el futuro con un input
//
//    double &fn = LocalElasticContactForce[2] q inicialmente entra como 0.0
//
//    if indentation >= 0 {
//        fn = kn1 * (initial_dist - distance); // = kn1 * indentation; fuerza en parte lineal
//
//        if (distance < mHistDist) //// indentation >= u_max   vemos si esta en carga, comparando la distancia actual entre particulas con la maxima historica
//        {
//            mHistDist = distance;
//
//            if (distance < max_dist) //// indentation >= C1*Area/kn_el   se supera el limite para el cambio de pendiente.
//            {
//                fn = kn1 * (initial_dist - max_dist) + kn2 * (max_dist - distance);
//            }
//        }
//    }
//    double mHistoryMaxForce = fn; //guardamos el maximo historico de fuerza fn
//
//    else // esta en descarga, la distancia entre particulas aumenta respecto la historica, current_dist > mHistDist
//    {
//        if (hist_fn > 0); //  fuerza normal esta en el rango de compresion de la curva
//        {
//            plast_dist = max_dist; // initial_dist*(1-plast_def) distancia associada al valor de plast_def impuesto, por ahora coincide con el cambio de pendiente
//            if (plast_dist > mHistDist) // mientras se este por encima de la maxima historica, estamos en plasticidad.
//            {
//                fn = hist_fn + kn_el_plas * (mHistDist - distance); // en descarga: 500 - kn_plas(10 - 12). en carga 500 - kn_el(10 - 11) pero con distance > mHistDistance
//            } else // esta en descarga pero no en la zona plastica, descarga por la linea elastica
//            {
//                if (distance < max_dist) // se supera el limite para le cambio de pendiente, mientras plast_dist=max_dist nunca pasara, nunca descargara por kn2
//                {
//                    NonlinearNormalForceCalculation(LocalElasticContactForce, kn1, kn2, distance, max_dist, initial_dist);
//                } else // descarga por la primera rama elastica
//                {
//                    LocalElasticContactForce[2] = kn1 * (initial_dist - distance); // fuerza en parte lineal
//                }
//            }
//        }
//    }
//}












