/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu
// Date: Sep 2023
/////////////////////////////////////////////////

#if !defined(FAST_FILLING_PARTICLES_H)
#define FAST_FILLING_PARTICLES_H

// System includes
#include <string>
#include <iostream>
#include <random>

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "includes/indexed_object.h"
#include "containers/global_pointers_vector.h"
#include "includes/constitutive_law.h"
#include "includes/condition.h"
#include "custom_elements/discrete_element.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_utilities/random_variable.h"
#include "custom_utilities/piecewise_linear_random_variable.h"
#include "custom_utilities/discrete_random_variable.h"
#include "custom_utilities/properties_proxies.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    class ParticleCreatorDestructor;

    class  KRATOS_API(DEM_APPLICATION) Fast_Filling_Creator
    {
    public:

        typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;
        typedef GlobalPointersVector<Element> ParticleWeakVectorType;
        typedef ModelPart::ElementsContainerType ElementsArrayType;

        KRATOS_CLASS_POINTER_DEFINITION(Fast_Filling_Creator);

        /// Constructor:
        Fast_Filling_Creator(const int seed=42);
        Fast_Filling_Creator(const Parameters& r_fast_filling_creator_settings, const int seed=42);

        /// Destructor.
        virtual ~Fast_Filling_Creator(){}

        Fast_Filling_Creator(const Fast_Filling_Creator&) = delete;
        Fast_Filling_Creator& operator=(const Fast_Filling_Creator&) = delete;

        virtual double GetRandomParticleRadius(ParticleCreatorDestructor& creator);
        virtual bool CheckHasIndentationOrNot(const double x_1, const double y_1, const double z_1, const double r_1,
                                              const double x_2, const double y_2, const double z_2, const double r_2);

    protected:

    private:

        int  mTotalNumberOfParticlesInjected;
        double mTotalMassInjected;
        //int mSeed;
        std::vector<double> mMassInjected;
        std::mt19937 mGenerator;

        std::map<std::string, std::unique_ptr<RandomVariable>> mInletsRandomVariables;
        std::map<std::string, Parameters> mInletsRandomSettings;
        Parameters mFastFillingCreatorSettings;
    };
}// namespace Kratos.

#endif // FAST_FILLING_PARTICLES_H defined

