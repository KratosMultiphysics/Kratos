#include "translational_RK4_scheme.h"

namespace Kratos {

    void TranslationalRungeKuttaScheme::SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void TranslationalRungeKuttaScheme::SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void TranslationalRungeKuttaScheme::UpdateTranslationalVariables(
            int StepFlag,
            Node < 3 >& i,
            array_1d<double, 3 >& coor,
            array_1d<double, 3 >& displ,
            array_1d<double, 3 >& delta_displ,
            array_1d<double, 3 >& vel,
            const array_1d<double, 3 >& initial_coor,
            const array_1d<double, 3 >& force,
            const double force_reduction_factor,
            const double mass,
            const double delta_t,
            const bool Fix_vel[3]) {

        double mass_inv = 1.0 / mass;

        if(StepFlag == 0) //Init
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    mInitialDispl[k] = displ[k];

                    mRungeKuttaL(StepFlag,k) = delta_t * force[k] * mass_inv;
                    mRungeKuttaK(StepFlag,k) = delta_t * vel[k];

                    KRATOS_WATCH(mRungeKuttaL(StepFlag,k))
                    KRATOS_WATCH(mRungeKuttaK(StepFlag,k))

                    delta_displ[k] = mRungeKuttaK(StepFlag,k) * 0.5 ;
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];

                } else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            }
        }
        else if(StepFlag == 1) //Step1
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    mRungeKuttaL(StepFlag,k) = delta_t * force[k] * mass_inv;
                    mRungeKuttaK(StepFlag,k) = delta_t * vel[k] +  delta_t * mRungeKuttaL(StepFlag-1,k) * 0.5;

                    KRATOS_WATCH(mRungeKuttaL(StepFlag,k))
                    KRATOS_WATCH(mRungeKuttaK(StepFlag,k))

                    delta_displ[k] = mRungeKuttaK(StepFlag,k) * 0.5 ;
                    displ[k] = mInitialDispl[k] + delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            }
        }
        else if(StepFlag == 2) //Step2
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    mRungeKuttaL(StepFlag,k) = delta_t * force[k] * mass_inv;
                    mRungeKuttaK(StepFlag,k) = delta_t * vel[k] +  delta_t * mRungeKuttaL(StepFlag-1,k) * 0.5;

                    KRATOS_WATCH(mRungeKuttaL(StepFlag,k))
                    KRATOS_WATCH(mRungeKuttaK(StepFlag,k))

                    delta_displ[k] = mRungeKuttaK(StepFlag,k);
                    displ[k] = mInitialDispl[k] + delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            }
        }
        else if(StepFlag == 3) //Step3
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    mRungeKuttaL(StepFlag,k) = delta_t * force[k] * mass_inv;
                    mRungeKuttaK(StepFlag,k) = delta_t * vel[k] +  delta_t * mRungeKuttaL(StepFlag-1,k);

                    delta_displ[k] = 1.0/6.0 * (mRungeKuttaK(0,k) + 2.0 * mRungeKuttaK(1,k) + 2.0 * mRungeKuttaK(2,k) + mRungeKuttaK(3,k));
                    KRATOS_WATCH(delta_displ[k])

                    displ[k] = mInitialDispl[k] + delta_displ[k];
                    KRATOS_WATCH(displ[k])

                    coor[k] = initial_coor[k] + displ[k];
                    KRATOS_WATCH(coor[k])

                    vel[k] += 1.0/6.0 * (mRungeKuttaL(0,k) + 2.0 * mRungeKuttaL(1,k) + 2.0 * mRungeKuttaL(2,k) + mRungeKuttaL(3,k));
                    KRATOS_WATCH(vel[k])

                }
            }
        }
    }

}
