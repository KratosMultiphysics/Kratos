#pragma once 

#include "custom_utilities/neighbours_search_utilities.h"

namespace Kratos
{

double NeighboursSearchUtilities::ComputeSmoothingLength(const ModelPart& rThisModelPart, double Coeff)
{
    double h = Coeff * ComputeInterparticleMinDist(rThisModelPart);
    return h;
}

double NeighboursSearchUtilities::ComputeInterparticleMinDist(const ModelPart& rThisModelPart)
{
    // This function should be updated with a more efficient version using bins/tree
    const auto& rnodes = rThisModelPart.Nodes();
    const std::size_t domain_size = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    double min_dist = std::numeric_limits<double>::max();

    for (auto node1 = rnodes.begin(); node1 != rnodes.end(); ++node1){
       
        const auto& coords1 = node1->Coordinates();
        auto node2 = node1;
        ++node2;
        
        for (; node2 != rnodes.end(); ++node2){
            const auto& coords2 = node2->Coordinates();

            const double dx = std::abs(coords1[0] - coords2[0]);
            const double dy = std::abs(coords1[1] - coords2[1]);
            const double dz = std::abs(coords1[2] - coords2[2]);

            const double dist = dx * dx + dy * dy + dz * dz;

            if (dist < min_dist) min_dist = dist;
        }
    }

    return std::sqrt(domain_size * min_dist);
}

}