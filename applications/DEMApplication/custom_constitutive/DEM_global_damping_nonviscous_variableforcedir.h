//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#pragma once

// Project includes
#include "DEM_global_damping.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos
{
    class KRATOS_API(DEM_APPLICATION) DEMGlobalDampingNonViscousVarForceDir : public DEMGlobalDampingModel
    {
        public:
            // Pointer definition
            KRATOS_CLASS_POINTER_DEFINITION(DEMGlobalDampingNonViscousVarForceDir);

            // Constructor / Destructor
            DEMGlobalDampingNonViscousVarForceDir() {}
            virtual ~DEMGlobalDampingNonViscousVarForceDir() {}

            // Clone
            DEMGlobalDampingModel::Pointer Clone() const override;
            std::unique_ptr<DEMGlobalDampingModel> CloneUnique() override;

            // Public methods
            void AddGlobalDampingForceAndMoment(SphericParticle* p_element, array_1d<double,3>& total_forces, array_1d<double,3>& total_moment) override;

        private:
            // Serializer
            friend class Serializer;
            void save(Serializer & rSerializer) const override {}
            void load(Serializer & rSerializer) override {}
    };
}
