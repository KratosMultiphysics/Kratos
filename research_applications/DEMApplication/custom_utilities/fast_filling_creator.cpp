/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu
// Date: Sep 2023
/////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <random>
#include <functional>

#include "fast_filling_creator.h"
#include "create_and_destroy.h"
#include "custom_elements/spheric_continuum_particle.h"
#include "custom_elements/cluster3D.h"
#include "custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "custom_constitutive/DEM_continuum_constitutive_law.h"
#include "dem_fem_utilities.h"
#include "GeometryFunctions.h"


namespace Kratos {

    /// Constructor

    Fast_Filling_Creator::Fast_Filling_Creator(const int seed):
                        Fast_Filling_Creator(Parameters(R"({})"), seed){}

    Fast_Filling_Creator::Fast_Filling_Creator(const Parameters& r_fast_filling_creator_settings, const int seed):
     mFastFillingCreatorSettings(Parameters(r_fast_filling_creator_settings))
        {
        std::mt19937 gen(seed);
        mGenerator = gen;

        mTotalNumberOfParticlesInjected = 0;
        mTotalMassInjected = 0.0;
    }
    
    double Fast_Filling_Creator::GetRandomParticleRadius(ParticleCreatorDestructor& creator) {

        KRATOS_TRY

        if (mFastFillingCreatorSettings["PROBABILITY_DISTRIBUTION"].GetString() == "piecewise_linear" || mFastFillingCreatorSettings["PROBABILITY_DISTRIBUTION"].GetString() == "discrete"){

            const Parameters& rv_settings = mFastFillingCreatorSettings["random_variable_settings"];
            int seed = rv_settings["seed"].GetInt();
            if (!rv_settings["do_use_seed"].GetBool()){
                seed = std::random_device{}();
            }
            if (mFastFillingCreatorSettings["PROBABILITY_DISTRIBUTION"].GetString() == "piecewise_linear"){

                mInletsRandomVariables[mFastFillingCreatorSettings["NAME"].GetString()] = std::unique_ptr<PiecewiseLinearRandomVariable>(new PiecewiseLinearRandomVariable(rv_settings, seed));
            }

            else if (mFastFillingCreatorSettings["PROBABILITY_DISTRIBUTION"].GetString() == "discrete"){
                mInletsRandomVariables[mFastFillingCreatorSettings["NAME"].GetString()] = std::unique_ptr<DiscreteRandomVariable>(new DiscreteRandomVariable(rv_settings, seed));
            }

            else {
                KRATOS_ERROR << "Unknown Fast Filling Creator random variable: " << mFastFillingCreatorSettings["PROBABILITY_DISTRIBUTION"].GetString() << ".";
            }
        }

        double radius = creator.SelectRadius(mFastFillingCreatorSettings, mInletsRandomVariables);

        return radius;

        KRATOS_CATCH("")
    }

    bool Fast_Filling_Creator::CheckHasIndentationOrNot(const double x_1, const double y_1, const double z_1, const double r_1,
                                              const double x_2, const double y_2, const double z_2, const double r_2){
        
        KRATOS_TRY
        
        double distance = std::sqrt(std::pow(x_1 - x_2, 2) + std::pow(y_1 - y_2, 2) + std::pow(z_1 - z_2, 2));

        return distance <= (r_1 + r_2); 

        KRATOS_CATCH("")                                    
    }

} // namespace Kratos
