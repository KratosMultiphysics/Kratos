//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andi Makarim Katili
//


#ifndef KRATOS_MPM_VOLUME_SUM_UTILITY
#define KRATOS_MPM_VOLUME_SUM_UTILITY

// System includes

// External includes

// Project includes
#include "mpm_application_variables.h"
#include "containers/model.h"
#include "includes/element.h"

namespace Kratos
{
namespace MPMVolumeSumUtility
{
    void AddModelPartMPMVolumeIntoGrid(const ModelPart& rModelPart)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i) {
            auto pMPMElement = (rModelPart.ElementsBegin() + i);
            KRATOS_ERROR_IF(!pMPMElement->GetGeometry().GetGeometryParent(0).Has(TOTAL_MP_VOLUME)) << "The TOTAL_MP_VOLUME variable is not defined in the parent geometry." << std::endl;
            
            std::vector<double> current_mp_volume;
            pMPMElement->CalculateOnIntegrationPoints(MP_VOLUME, current_mp_volume, rModelPart.GetProcessInfo());
            
            double& rTOTAL_MP_VOLUME = pMPMElement->GetGeometry().GetGeometryParent(0).GetValue(TOTAL_MP_VOLUME);
            AtomicAdd(rTOTAL_MP_VOLUME, current_mp_volume[0]);
        }
    }

} // end namespace MPMVolumeSumUtility
} // end namespace Kratos
#endif // KRATOS_MPM_VOLUME_SUM_UTILITY


