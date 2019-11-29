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
#include "utilities/body_normal_calculation_utils.h"
#include "wss_statistics_utilities.h"
#include "fluid_dynamics_biomedical_application_variables.h"

namespace Kratos {

void WssStatisticsUtilities::CalculateWSS(ModelPart &rModelPart) 
//, ModelPart &rSkinModelPart)
{
    // Initialize WSS variables
    array_1d<double,3> aux_zero = ZeroVector(3);
    const unsigned int n_nodes = rModelPart.NumberOfNodes();
    for (int i_node = 0; i_node < n_nodes; ++i_node) {
        auto it_node = rModelPart.NodesBegin() + i_node;
        it_node->SetValue(WSS, 0.0);
        it_node->SetValue(ECAP, 0.0);
        it_node->SetValue(RRT, 0.0);
        it_node->SetValue(TWSS, 0.0);
        it_node->SetValue(TAWSS, 0.0);
        it_node->SetValue(WSS_NORMAL_STRESS, aux_zero);
        it_node->SetValue(WSS_TANGENTIAL_STRESS, aux_zero);
        it_node->SetValue(TEMPORAL_OSI, aux_zero);
    }
    
    // Check if buffer size is filled (we need REACTION to be already computed)
    const unsigned int buffer_size = rModelPart.GetBufferSize();
    const unsigned int step = rModelPart.GetProcessInfo()[STEP];
    if (step > buffer_size) {
        //Distribute the REACTION as a surface load

        const double tolerance = 1.0e-5;
        const unsigned int max_it = 100;
        
        VariableRedistributionUtility::DistributePointValues(rModelPart, REACTION, FACE_LOAD, tolerance, max_it);
        
        // Allocate auxiliary arrays
        array_1d<double,3> normal_load;
        array_1d<double,3> tangential_load;

        // Loop the WSS model part conditions
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = rModelPart.NodesBegin() + i_node;

            // Normalize nodal normal
            array_1d<double,3> normal = it_node->FastGetSolutionStepValue(NORMAL);
            const double normal_norm = norm_2(normal);
            if (normal_norm > 1.0e-12) {
                normal /= normal_norm;
            } else {
                std::cout << "Node " << str(it_node->Id()) << " has normal_norm " << str(normal_norm) << std::endl;
            }

            // Calculate the FACE_LOAD (distributed REACTION) projections
            const auto& r_face_load = it_node->FastGetSolutionStepValue(FACE_LOAD);
            //const double projection = r_face_load[0]*normal[0] + r_face_load[1]*normal[1] + r_face_load[2]*normal[2];
            const double projection = inner_prod(r_face_load, normal);
            //KRATOS_WATCH(it_node->FastGetSolutionStepValue(FACE_LOAD));
            //KRATOS_WATCH(normal_norm);
            //KRATOS_WATCH(r_face_load);
            //KRATOS_WATCH(projection);
            normal_load[0] = projection * normal[0];
            normal_load[1] = projection * normal[1];
            normal_load[2] = projection * normal[2];
            //normal_load = projection * normal;
            tangential_load[0] = r_face_load[0] - normal_load[0];
            tangential_load[1] = r_face_load[1] - normal_load[1];
            tangential_load[2] = r_face_load[2] - normal_load[2];
            //tangential_load = r_face_load - normal_load;

            const double wss = norm_2(tangential_load);
            // KRATOS_WATCH(wss);
            // Save computed magnitudes
            it_node->GetValue(WSS) = wss;
            it_node->GetValue(WSS_NORMAL_STRESS) = normal_load;
            it_node->GetValue(WSS_TANGENTIAL_STRESS) = tangential_load;
            //KRATOS_WATCH(it_node->GetValue(WSS_TANGENTIAL_STRESS))
        }
    }
}

void WssStatisticsUtilities::CalculateTWSS(ModelPart &rModelPart)
{ 

    // Allocate auxiliary arrays
    array_1d<double,3> previous_tangential;
    array_1d<double,3> tangential;
    // Allocate auxiliary arrays
    double aux_ECAP = 0.0;
    double aux_RRT = 0.0;
    double aux_OSI = 0.0;
    double aux_TWSS = 0.0;
    double aux_WSS = 0.0;
    double factor = 0.0;
    array_1d<double,3> sum_WSS;

    // Get the step counter
    const unsigned int step = rModelPart.GetProcessInfo()[STEP];
    const unsigned int buffer_size = rModelPart.GetBufferSize();    
    if (step > buffer_size) {
        //TWSS :Time-averaged wall shear stress (save in the last value)
        double twss = 0.0;
        const unsigned int n_nodes = rModelPart.NumberOfNodes();
        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node)
        {
            auto it_node = rModelPart.NodesBegin() + i_node;
            //Calculate the sum of the vector components of WSS for all times step
            previous_tangential=it_node->GetValue(TEMPORAL_OSI);
            //KRATOS_WATCH(previous_tangential)
            tangential=it_node->GetValue(WSS_TANGENTIAL_STRESS);
            //KRATOS_WATCH(tangential)
            previous_tangential[0]=previous_tangential[0]+((tangential[0]-previous_tangential[0])/step);
            previous_tangential[1]=previous_tangential[1]+((tangential[1]-previous_tangential[1])/step);
            previous_tangential[2]=previous_tangential[2]+((tangential[2]-previous_tangential[2])/step);
            it_node->GetValue(TEMPORAL_OSI)=previous_tangential;
            //KRATOS_WATCH(previous_tangential)
            //Calculates the sum of the WSS magnitudes for all time steps
            twss= it_node->GetValue(TAWSS,0) + (norm_2(tangential) - (it_node->GetValue(TAWSS,0)/(step)) );
            it_node->GetValue(TWSS)=twss;
            sum_WSS=it_node->GetValue(TEMPORAL_OSI);
            //KRATOS_WATCH(sum_WSS)
            //KRATOS_WATCH(step)            
            factor = 1.0/static_cast<double>(step);
            // KRATOS_WATCH(factor)
            sum_WSS *= factor;
            // KRATOS_WATCH(sum_WSS)  
            aux_WSS=norm_2(sum_WSS);
            // KRATOS_WATCH(aux_WSS)
            //Calculates the magnitude of the time-averaged WSS vector
            aux_TWSS=(it_node->GetValue(TWSS,0)/step);
            // KRATOS_WATCH(aux_WSS)
            if((aux_WSS/aux_TWSS)>1){
                    aux_OSI = 0.0
                else{
                    aux_OSI=(0.5* (1.0-(aux_WSS/aux_TWSS)));
                }
            }
            if (aux_WSS > 1.0e-12) {
                if (aux_OSI==0.5){
                    aux_RRT=0.0;
                }
                else {
                    aux_RRT=1/((1-2*aux_OSI)*aux_WSS);
                }
                aux_ECAP= (aux_OSI/aux_WSS);  
            }      
            it_node->GetValue(ECAP)=aux_ECAP;
            it_node->GetValue(RRT)=aux_RRT;
            it_node->GetValue(OSI)=aux_OSI;
            it_node->GetValue(TWSS)=aux_TWSS;
        }
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
    array_1d<double,3> sum_WSS;

    // // Get the step counter
    // const unsigned int step = rModelPart.GetProcessInfo()[STEP];
    // const unsigned int buffer_size = rModelPart.GetBufferSize();    
    // if (step > buffer_size) {
    //     const unsigned int n_nodes = rModelPart.NumberOfNodes();
    //     #pragma omp parallel for
    //     for (int i_node = 0; i_node < n_nodes; ++i_node)
    //     {
    //         auto it_node = rModelPart.NodesBegin() + i_node;
    //         //Calculate the sum of the vector components of WSS for all times step
    //         sum_WSS=it_node->GetValue(TEMPORAL_OSI);        
    //         sum_WSS *= (1/step);
    //         aux_WSS=norm_2(sum_WSS);
    //         KRATOS_WATCH(aux_WSS)
    //         //Calculates the magnitude of the time-averaged WSS vector
    //         aux_TWSS=(it_node->GetValue(TWSS,0)/step);
    //         // KRATOS_WATCH(aux_WSS)
    //         aux_OSI=(0.5* (1.0-(aux_WSS/aux_TWSS)));
    //         aux_RRT=1/((1-2*OSI)*aux_WSS);
    //         aux_ECAP= (OSI/aux_WSS);
    //         it_node->GetValue(ECAP)=aux_ECAP;
    //         it_node->GetValue(RRT)=aux_RRT;
    //         it_node->GetValue(OSI)=aux_OSI;
    //         it_node->GetValue(TWSS)=aux_TWSS;
    //     }
    // }
}

// // // CONDITIONS 
// void WssStatisticsUtilities::CalculateWSSGauss(ModelPart &rModelPart) {
    
//     KRATOS_WATCH ("Computing CONDITIONS WSS")

//     // Normal Calculation
//     BodyNormalCalculationUtils bncu;
//     bncu.CalculateBodyNormals(rModelPart, 3);

//     // Distribute the REACTION as a surface load
//     const double tolerance = 1.0e-5;
//     const unsigned int max_it = 100;
//     VariableRedistributionUtility::DistributePointValues(rModelPart, REACTION, FACE_LOAD, tolerance, max_it);

//     // Allocate auxiliary arrays
//     array_1d<double,3> normal_load;
//     array_1d<double,3> tangential_load;
//     array_1d<double,3> aux_face_load;

//     ModelPart::ConditionIterator CondBegin;
//     ModelPart::ConditionIterator CondEnd;

//     for (ModelPart::ConditionIterator itCond = CondBegin; itCond != CondEnd; ++itCond)
//     {
//         Condition::ConditionType& rNodes = itCond->GetNodes();
//         // Normalize nodal normal
//         array_1d<double,3> normal = itCond->GetValue(NORMAL);
//         const double normal_norm = norm_2(normal);
//         normal /= normal_norm;
//         for (unsigned int i = 0; i < rNodes.PointsNumber(); ++i)
//             {
//             node_load=i->GetValue(FACE_LOAD);
//             aux_face_load[0]= aux_face_load[0] + node_load[0];
//             aux_face_load[1]= aux_face_load[1] + node_load[1];
//             aux_face_load[2]= aux_face_load[2] + node_load[2];
//             }
//         aux_face_load *= 0.33333333333;
//         // Calculate the FACE_LOAD (distributed REACTION) projections
//         const double projection = inner_prod(aux_face_load, normal);

//         normal_load[0] = projection * normal[0];
//         normal_load[1] = projection * normal[1];
//         normal_load[2] = projection * normal[2];

//         tangential_load[0] = aux_face_load[0] - normal_load[0];
//         tangential_load[1] = aux_face_load[1] - normal_load[1];
//         tangential_load[2] = aux_face_load[2] - normal_load[2];

//         const double wss = norm_2(tangential_load);

//         // Save computed magnitudes
//         cond_it->GetValue(WSS) = wss;
//         cond_it->GetValue(WSS_NORMAL_STRESS) = normal_load;
//         cond_it->GetValue(WSS_TANGENTIAL_STRESS) = tangential_load;
//     }
// }




}