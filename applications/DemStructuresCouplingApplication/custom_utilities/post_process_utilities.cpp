//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Guillermo Casas
//                  Ignasi de Pouplana
//

// Application includes
#include "custom_elements/spheric_continuum_particle.h"
#include "post_process_utilities.hpp"

#include "dem_structures_coupling_application_variables.h"


namespace Kratos
{

PostProcessUtilities::PostProcessUtilities(ModelPart& rModelPart) : mrModelPart(rModelPart) {
    KRATOS_TRY

    KRATOS_CATCH("")
}

/// Destructor.

PostProcessUtilities::~PostProcessUtilities(){}

void PostProcessUtilities::GetStickyStatus(pybind11::list& is_sticky_list)
{
    ModelPart::ElementsContainerType& rElements = mrModelPart.GetCommunicator().LocalMesh().Elements();
    const int num_elements = (int)rElements.size();

    this->Clear(is_sticky_list);

    // #pragma omp parallel for
    for(int i = 0; i < num_elements; i++)
    {
        ModelPart::ElementsContainerType::ptr_iterator ptr_itElem = rElements.ptr_begin() + i;

        Element* p_element = ptr_itElem->get();
        SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*> (p_element);

        is_sticky_list.append(pDemElem->Is(DEMFlags::STICKY));
    }

}


void PostProcessUtilities::GetInitialContinuumBonds(pybind11::list& initial_continuum_bonds_list) {

    ModelPart::ElementsContainerType& rElements = mrModelPart.GetCommunicator().LocalMesh().Elements();
    const int num_elements = (int)rElements.size();

    this->Clear(initial_continuum_bonds_list);

    // #pragma omp parallel for
    for(int i = 0; i < num_elements; i++)
    {
        ModelPart::ElementsContainerType::ptr_iterator ptr_itElem = rElements.ptr_begin() + i;

        Element* p_element = ptr_itElem->get();
        SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*> (p_element);

        initial_continuum_bonds_list.append(pDemElem->mContinuumInitialNeighborsSize);
    }
}


void PostProcessUtilities::GetCurrentContinuumBonds(pybind11::list& current_continuum_bonds_list) {

    ModelPart::ElementsContainerType& rElements = mrModelPart.GetCommunicator().LocalMesh().Elements();
    const int num_elements = (int)rElements.size();

    this->Clear(current_continuum_bonds_list);

    // #pragma omp parallel for
    for(int i = 0; i < num_elements; i++)
    {
        ModelPart::ElementsContainerType::ptr_iterator ptr_itElem = rElements.ptr_begin() + i;

        Element* p_element = ptr_itElem->get();
        SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*> (p_element);

        const unsigned int continuum_initial_neighbors_size = pDemElem->mContinuumInitialNeighborsSize;
        unsigned int broken_bonds_counter = 0.0;

        for (unsigned int j = 0; j < pDemElem->mNeighbourElements.size(); j++) {

            if (pDemElem->mNeighbourElements[j] == NULL) continue;

            if ((j < continuum_initial_neighbors_size) && pDemElem->mIniNeighbourFailureId[j] > 0) {
                ++broken_bonds_counter;
            }
        }
        const unsigned int intact_continuum_bonds = continuum_initial_neighbors_size - broken_bonds_counter;
        current_continuum_bonds_list.append(intact_continuum_bonds);
    }
}

void PostProcessUtilities::Clear(pybind11::list& my_list){

    for(unsigned int i = 0; i < my_list.size(); i++)
    {
        my_list.attr("pop")();
    }

}

std::string PostProcessUtilities::Info() const
{
    return "";
}

/// Print information about this object.

void PostProcessUtilities::PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.

void PostProcessUtilities::PrintData(std::ostream& rOStream) const
{
}

}  // namespace Kratos.
