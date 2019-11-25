//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah,
//                   Ruben Zorilla,
//                   Eduardo Soudah
//

#include "includes/variables.h"
#include "utilities/variable_redistribution_utility.h"
#include "wss_statistics_utilities.h"
#include "fluid_dynamics_biomedical_application_variables.h"

namespace Kratos {

void WssStatisticsUtilities::CalculateWSS(ModelPart &rModelPart)
{
    // Distribute the REACTION as a surface load
    const double tolerance = 1.0e-5;
    const unsigned int max_it = 100;
    VariableRedistributionUtility::DistributePointValues(rModelPart, REACTION, FACE_LOAD, tolerance, max_it);

    // Allocate auxiliary arrays
    array_1d<double,3> normal_load;
    array_1d<double,3> tangential_load;

    // Loop the WSS model part conditions
    const unsigned int n_nodes = rModelPart.NumberOfNodes();
    #pragma omp parallel for
    for (int i_node = 0; i_node < n_nodes; ++i_node) {
        auto it_node = rModelPart.NodesBegin() + i_node;

        // Normalize nodal normal
        array_1d<double,3> normal = it_node->FastGetSolutionStepValue(NORMAL);
        const double normal_norm = norm_2(normal);
        normal /= normal_norm;

        // Calculate the FACE_LOAD (distributed REACTION) projections
        const auto& r_face_load = it_node->FastGetSolutionStepValue(FACE_LOAD);
        const double projection = inner_prod(r_face_load, normal);

        normal_load[0] = projection * normal[0];
        normal_load[1] = projection * normal[1];
        normal_load[2] = projection * normal[2];

        tangential_load[0] = r_face_load[0] - normal_load[0];
        tangential_load[1] = r_face_load[1] - normal_load[1];
        tangential_load[2] = r_face_load[2] - normal_load[2];

        const double wss = norm_2(tangential_load);

        // Save computed magnitudes
        it_node->FastGetSolutionStepValue(WSS) = wss;
        it_node->FastGetSolutionStepValue(WSS_NORMAL_STRESS) = normal_load;
        it_node->FastGetSolutionStepValue(WSS_TANGENTIAL_STRESS) = tangential_load;
    }
}

void WssStatisticsUtilities::CalculateTWSS(ModelPart &rModelPart)
{
    // Allocate auxiliary arrays
    array_1d<double,3> previous_tangential;
    array_1d<double,3> tangential;

    // Get the step counter
    const unsigned int step = rModelPart.GetProcessInfo()[STEP];

    //TWSS :Time-averaged wall shear stress (save in the last value)
    double twss = 0.0;
    const unsigned int n_nodes = rModelPart.NumberOfNodes();
    #pragma omp parallel for
    for (int i_node = 0; i_node < n_nodes; ++i_node)
    {
        auto it_node = rModelPart.NodesBegin() + i_node;
        //Calculate the sum of the vector components of WSS for all times step
        previous_tangential=it_node->FastGetSolutionStepValue(TEMPORAL_OSI,0);
        tangential=it_node->FastGetSolutionStepValue(TANGENTIAL_STRESS,0);
        previous_tangential[0]=previous_tangential[0]+((tangential[0]-previous_tangential[0])/step);
        previous_tangential[1]=previous_tangential[1]+((tangential[1]-previous_tangential[1])/step);
        previous_tangential[2]=previous_tangential[2]+((tangential[2]-previous_tangential[2])/step);
        it_node->FastGetSolutionStepValue(TEMPORAL_OSI)=previous_tangential;
        //Calculates the sum of the WSS magnitudes for all time steps
        twss= it_node->FastGetSolutionStepValue(TAWSS,0) + (norm_2(tangential) - (it_node->FastGetSolutionStepValue(TAWSS,0)/(step)) );
        it_node->FastGetSolutionStepValue(TAWSS)=twss;
    }
}

void WssStatisticsUtilities::CalculateOSI(ModelPart &rModelPart)
{
    // Allocate auxiliary arrays
    double aux_ECAP = 0.0;
    double aux_RRT = 0.0;
    double aux_OSI = 0.0;
    double aux_TWSS = 0.0;
    double aux_WSS = 0.0;
    array_1d<double,3> SumWSS;

    // Get the step counter
    const unsigned int step = rModelPart.GetProcessInfo()[STEP];

    const unsigned int n_nodes = rModelPart.NumberOfNodes();
    #pragma omp parallel for
    for (int i_node = 0; i_node < n_nodes; ++i_node)
    {
        auto it_node = rModelPart.NodesBegin() + i_node;
        //Calculate the sum of the vector components of WSS for all times step
        SumWSS=it_node->FastGetSolutionStepValue(TEMPORAL_OSI,0);
        SumWSS *= (1/step);
        aux_WSS=norm_2(SumWSS);
        //Calculates the magnitude of the time-averaged WSS vector
        aux_TWSS=(it_node->FastGetSolutionStepValue(TWSS,0)/step);
        aux_OSI=(0.5* (1.0-(aux_WSS/aux_TWSS)));
        aux_RRT=1/((1-2*OSI)*aux_WSS);
        aux_ECAP= (OSI/aux_WSS);
        it_node->FastGetSolutionStepValue(ECAP)=aux_ECAP;
        it_node->FastGetSolutionStepValue(RRT)=aux_RRT;
        it_node->FastGetSolutionStepValue(OSI)=aux_OSI;
        it_node->FastGetSolutionStepValue(TWSS)=aux_TWSS;
    }
}



}