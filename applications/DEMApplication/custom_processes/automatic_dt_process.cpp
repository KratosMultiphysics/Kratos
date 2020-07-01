//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// Application includes
#include "custom_processes/automatic_dt_process.hpp"

namespace Kratos
{

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void AutomaticDTProcess::ExecuteInitialize() override
    {
        KRATOS_TRY;

        // const int NNodes = static_cast<int>(mrModelPart.Nodes().size());
        // ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
        // #pragma omp parallel for
        // for(int i = 0; i<NNodes; i++) {
        //     ModelPart::NodesContainerType::iterator it = it_begin + i;
        //     it->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
        // }
        // ModelPart::ElementsContainerType& rElements = mrModelPart.GetCommunicator().LocalMesh().Elements();

        const int NElems = static_cast<int>(mrModelPart.GetCommunicator().LocalMesh().Elements().size());
        ModelPart::ElementsContainerType::ptr_iterator ptr_itElem_begin = rElements.ptr_begin();

        // Find smallest particle
        double min_radius = std::numeric_limits<double>::infinity();
        SphericContinuumParticle* pSmallestDemElem;
        #pragma omp parallel for private(min_radius,pSmallestDemElem)
        for (int i = 0; i < NElems; i++) {
            ModelPart::ElementsContainerType::ptr_iterator ptr_itElem = ptr_itElem_begin + i;
            Element* p_element = ptr_itElem->get();
            SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*>(p_element);
            const double radius = pDemElem->GetRadius();
            if (radius < min_radius) {
                pSmallestDemElem = pDemElem;
                min_radius = radius;
            }
        }

        // Calculate stiffness of the smallest particle with itself
        double initial_dist = 2.0*min_radius;
        double myYoung = pSmallestDemElem->GetYoung();
        double myPoisson = pSmallestDemElem->GetPoisson();
        Vector& cont_ini_neigh_area = pSmallestDemElem->GetValue(NEIGHBOURS_CONTACT_AREAS);
        double calculation_area = 0.0;
        double kn_el = 0.0;
        double kt_el = 0.0;
        for (int i = 0; i < pSmallestDemElem->mContinuumInitialNeighborsSize; i++) {
            pSmallestDemElem->mContinuumConstitutiveLawArray[i]->GetContactArea(min_radius, min_radius, cont_ini_neigh_area, i, calculation_area);
            pSmallestDemElem->mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(kn_el, kt_el, initial_dist, myYoung, myPoisson, calculation_area, pSmallestDemElem, pSmallestDemElem);
        }

        // Calculate mass of the smallest particle
        const double particle_density = pSmallestDemElem->GetDensity();
        const double particle_volume = pSmallestDemElem->CalculateVolume();
        const double min_mass = particle_volume * particle_density;

        // Calculate critical delta time
        const double critical_delta_time = std::sqrt(min_mass/kn_el);
        mrModelPart.GetProcessInfo().SetValue(DELTA_TIME,mCorrectionFactor*critical_delta_time);

        KRATOS_WATCH(critical_delta_time)

        KRATOS_CATCH("");
    }


    ///------------------------------------------------------------------------------------

} // namespace Kratos.
