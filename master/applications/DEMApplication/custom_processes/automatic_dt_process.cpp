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
    void AutomaticDTProcess::ExecuteBeforeSolutionLoop()
    {
        KRATOS_TRY;

        ModelPart::ElementsContainerType& rElements = mrModelPart.GetCommunicator().LocalMesh().Elements();
        const int NElems = static_cast<int>(rElements.size());
        ModelPart::ElementsContainerType::ptr_iterator ptr_itElem_begin = rElements.ptr_begin();

        // Find smallest particle
        double min_radius = std::numeric_limits<double>::infinity();
        SphericContinuumParticle* pSmallestDemElem = dynamic_cast<SphericContinuumParticle*>(ptr_itElem_begin->get());

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
        double calculation_area = 0.0;
        double kn_el = 0.0;
        double kt_el = 0.0;
        double indentation = 0.0;
        DEMContinuumConstitutiveLaw::Pointer NewContinuumConstitutiveLaw = pSmallestDemElem->GetProperties()[DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER]->Clone();
        NewContinuumConstitutiveLaw->CalculateContactArea(min_radius, min_radius, calculation_area);
        NewContinuumConstitutiveLaw->CalculateElasticConstants(kn_el, kt_el, initial_dist, myYoung, myPoisson, calculation_area, pSmallestDemElem, pSmallestDemElem, indentation);

        // Calculate mass of the smallest particle
        const double particle_density = pSmallestDemElem->GetDensity();
        const double particle_volume = pSmallestDemElem->CalculateVolume();
        const double min_mass = particle_volume * particle_density;

        // Calculate critical delta time
        const double critical_delta_time = std::sqrt(min_mass/kn_el);
        mrModelPart.GetProcessInfo().SetValue(DELTA_TIME,mCorrectionFactor*critical_delta_time);

        KRATOS_INFO("Automatic DT process") << "Calculated critical time step: " << critical_delta_time << " seconds." << std::endl;
        KRATOS_INFO("Automatic DT process") << "Using a correction factor of: " << mCorrectionFactor << ", the resulting time step is: " << mCorrectionFactor*critical_delta_time << " seconds." << std::endl;

        KRATOS_CATCH("");
    }


    ///------------------------------------------------------------------------------------

} // namespace Kratos.
