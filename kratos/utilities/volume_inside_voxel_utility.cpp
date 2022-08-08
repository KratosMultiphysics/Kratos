//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//

//Project includes
#include "volume_inside_voxel_utility.h"

namespace Kratos {

    double VolumeInsideVoxelUtility::NodesApproximation(
        const GeometryType& rVoxel        
    ) {
        double volume = 0;
        PointsArrayType nodes = rVoxel.Points();
        for (int i = 0; i < nodes.size(); i++) {
            if (nodes[i].GetSolutionStepValue(DISTANCE) > 0) {
                volume+=(1.0/nodes.size()); 
            } 
        }
        return volume;
    }

    /***********************************************************************************
     **********************************************************************************/
}