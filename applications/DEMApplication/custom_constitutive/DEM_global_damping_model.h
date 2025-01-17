//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#pragma once

// Includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/AuxiliaryFunctions.h"

namespace Kratos
{
    class SphericParticle; // Forward declaration of spheric particle

    class KRATOS_API(DEM_APPLICATION) DEMGlobalDampingModel
    {
        public:
            // Pointer definition
            KRATOS_CLASS_POINTER_DEFINITION(DEMGlobalDampingModel);

            // Constructor / destructor methods
            DEMGlobalDampingModel() {}
            virtual ~DEMGlobalDampingModel() {}

            // Clone
            virtual DEMGlobalDampingModel::Pointer Clone() const;
            virtual std::unique_ptr<DEMGlobalDampingModel> CloneUnique();

            // Public methods
            virtual void SetGlobalDampingModelInProperties(Properties::Pointer pProp, bool verbose = true);
            virtual void AddGlobalDampingForceAndMoment(SphericParticle* p_element, array_1d<double,3>& total_forces, array_1d<double,3>& total_moment) {}

            // Public attributes
            double mGlobalDamping = 0.0;

        private:
            // Serializer
            friend class Serializer;
            virtual void save(Serializer& rSerializer) const {}
            virtual void load(Serializer& rSerializer) {}
    };

    // This definition is done here to avoid recursive inclusion of header files
    KRATOS_DEFINE_APPLICATION_VARIABLE(DEM_APPLICATION, DEMGlobalDampingModel::Pointer, DEM_GLOBAL_DAMPING_MODEL_POINTER)
}
