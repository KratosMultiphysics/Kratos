//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Salva Latorre
//                  Ignasi de Pouplana
//

// Application includes
#include "custom_elements/spheric_continuum_particle.h"
#include "sand_production_utilities.hpp"

#include "dem_structures_coupling_application_variables.h"


namespace Kratos
{

/// Default constructor.

SandProductionUtilities::SandProductionUtilities(){}

/// Destructor.

SandProductionUtilities::~SandProductionUtilities(){}


void SandProductionUtilities::MarkSandProductionParticlesForErasing(ModelPart& r_model_part) {

    KRATOS_TRY

    Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    const double sigma_3_average = fabs(r_process_info[SIGMA_3_AVERAGE]);
    const double percentage_of_stress = 0.1; // 10%, obv

    #pragma omp parallel for
    for (int k = 0; k < (int)rElements.size(); k++){
        Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin() + k;

        if ((*particle_pointer_it)->Is(DEMFlags::BELONGS_TO_A_CLUSTER)) continue;
        if ((*particle_pointer_it)->Is(BLOCKED)) continue;
        if ((*particle_pointer_it)->Is(DEMFlags::STICKY)) continue;

        Element* p_element = particle_pointer_it->get();
        SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*>(p_element);
        BoundedMatrix<double, 3, 3> stress_tensor = (*(pDemElem->mSymmStressTensor));
        Vector principal_stresses(3);
        noalias(principal_stresses) = AuxiliaryFunctions::EigenValuesDirectMethod(stress_tensor);
        const double particle_sigma_3 = fabs(*std::min_element(principal_stresses.begin(), principal_stresses.end()));

        bool the_particle_plays_a_structural_role = (particle_sigma_3 > percentage_of_stress * sigma_3_average);

        // TODO: Do this correctly using the latest Ignasi's json
        the_particle_plays_a_structural_role = false;
        //

        // TODO: We need to think about this thoroughly:
        //       1. If we first erase a sphere (because of different reasons), we should then count it as SP,
        //       but we may erase structural spheres (spheres that play a role in the resistance of the specimen)
        //       2. If we count SP without erasing all counted spheres, we may have counting errors,
        //       but we avoid erasing structural spheres.
        if ((*particle_pointer_it)->Is(DEMFlags::IS_SAND_PRODUCTION) && !the_particle_plays_a_structural_role) {
            (*particle_pointer_it)->GetGeometry()[0].Set(TO_ERASE);
            (*particle_pointer_it)->Set(TO_ERASE);
        }
    }

    KRATOS_CATCH("")
}

std::string SandProductionUtilities::Info() const
{
    return "";
}

/// Print information about this object.

void SandProductionUtilities::PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.

void SandProductionUtilities::PrintData(std::ostream& rOStream) const
{
}


}  // namespace Kratos.
