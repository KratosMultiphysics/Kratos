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

void PostProcessUtilities::Clear(pybind11::list& my_list){

    for(unsigned int i = 0; i < my_list.size(); i++)
    {
        my_list.attr("pop")();
    }

}

}  // namespace Python.
