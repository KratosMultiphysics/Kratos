// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#include "custom_processes/neighbours_search_process.h"

namespace Kratos
{
    void NeighboursSearchProcess::Execute()
    {
        /*See if it is possible to use bins to search the neighbours.
        Especially for large scale simulation and or eulerian framework a more efficient 
        version of this function could be taken into account*/
        KRATOS_TRY
        const std::size_t domain_size = mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
        auto& relements = mrThisModelPart.Elements();
        const double h = NeighboursSearchUtilities::ComputeSmoothingLength(mrThisModelPart, mrThisParameters["coefficient"].GetDouble());
        const double radius = 2 * h;
        mrThisModelPart.GetProcessInfo().SetValue(SMOOTHING_LENGTH, h);

        // Allocating memory in each particle element for the list of neighbours 
        for (auto elem = relements.begin(); elem != relements.end(); ++elem ){
            elem->SetValue(NEIGHBOURS, std::vector<Element::Pointer>());
        }
        
        for (auto elem1 = relements.begin(); elem1 != relements.end(); ++elem1){
            auto& elem1_neigh = elem1->GetValue(NEIGHBOURS);
            elem1_neigh.push_back(elem1->shared_from_this());
            
            const auto& coords1 = elem1->GetGeometry()[0].Coordinates();
            auto elem2 = elem1;
            ++elem2;

            int counter = 0;

            for (; elem2 != relements.end(); ++elem2){

                const auto& coords2 = elem2->GetGeometry()[0].Coordinates();

                double dx = std::abs(coords1[0] - coords2[0]);
                double dy = std::abs(coords1[1] - coords2[1]);
                double dz = (domain_size == 3) ? std::abs(coords1[2] - coords2[2]) : 0.0;

                double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                if (dist <= radius){
                    auto& elem2_neigh = elem2->GetValue(NEIGHBOURS);
                    elem2_neigh.push_back(elem1->shared_from_this());
                    elem1_neigh.push_back(elem2->shared_from_this());
                }
            }
        } 
        KRATOS_CATCH("")
    }

    void NeighboursSearchProcess::ExecuteInitialize()
    {
        this->Execute();
        KRATOS_INFO("NeighboursSearch") << "The neighbours search process was executed" << std::endl;
    }
}