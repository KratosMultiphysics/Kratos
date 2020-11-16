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
//                   Ruben Zorilla
//

#include "includes/variables.h"
#include "utilities/variable_redistribution_utility.h"
#include "utilities/body_normal_calculation_utils.h"
#include "parabolic_profile_utilities.h"

namespace Kratos {

void CalculateParabolicProfile::Main(ModelPart &rModelPart) 
{
    const double max_dist = ComputeMaxDist(rModelPart);
    ImposeParabolic(rModelPart, max_dist);
}



double CalculateParabolicProfile::ComputeMaxDist(ModelPart &rModelPart) 
//, ModelPart &rSkinModelPart)
{
    // Initialize parabolic profile variables
    double dist_min=100000000000000000000.0;
    double dist_max=0.0;
    double dist_aux=0.0;
     // Get the step counter
    const unsigned int step = rModelPart.GetProcessInfo()[STEP];
    const unsigned int buffer_size = rModelPart.GetBufferSize();    
    if (step > buffer_size) {
        const unsigned int n_nodes = rModelPart.NumberOfNodes();
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node)
        {
            auto it_node = rModelPart.NodesBegin() + i_node;
            //Calculate distance for each node
            dist_aux=it_node->GetValue(DISTANCE);
            if(dist_aux<dist_min){
                dist_min=dist_aux;}
            if(dist_aux>dist_min){
                dist_max=dist_aux;}    
        }
    return dist_max;
}

void CalculateParabolicProfile::ImposeParabolic(
    ModelPart &rModelPart,
    const double MaxDist,
    const double ValueIn) 
//ModelPart &rSkinModelPart)
{
    // INPUT:Skin Inlet, value_in(Velocity or Flow)
    // OUTPUT: Value_out (velocity or Flow)   
     // Get the step counter
    const unsigned int step = rModelPart.GetProcessInfo()[STEP];
    const unsigned int buffer_size = rModelPart.GetBufferSize();    
    if (step > buffer_size) {
        const unsigned int n_nodes = rModelPart.NumberOfNodes();
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node)
        {
            array_1d<double,3> value_out =ZeroVector(3);  
            auto it_node = rModelPart.NodesBegin() + i_node;
            //Calculate distance for each node
            const double dist_aux = it_node->FastGetSolutionStepValue(DISTANCE);
            const auto& normals_in = it_node->GetValue(NORMAL);
            if (dist_aux < 0.0):
                dist_aux=0.0
            if (dist_aux > maxdist):
                dist_aux=MaxDist

            value_out = ValueIn * (1-((MaxDist-dist_aux)**2/(MaxDist**2)))
            const double nnorm = norm_2(normals_in);
            value_final *=-value_out/nnorm;
            it_node->FastGetSolutionStepValue(VELOCITY) = value_final;
            KRATOS_WATCH(value_final)
        }
    
}
