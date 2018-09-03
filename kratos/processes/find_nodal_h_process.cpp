//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Vicente Mataix Ferrandiz
//                    
//

// System includes
#include <limits>

// External includes

// Project includes
#include "processes/find_nodal_h_process.h"

namespace Kratos
{
void FindNodalHProcess::Execute()
{
    KRATOS_TRY
    
    // Check if variables are available       
    KRATOS_ERROR_IF_NOT(mrModelPart.NodesBegin()->SolutionStepsDataHas( NODAL_H )) << "Variable NODAL_H not in the model part!";
    
    #pragma omp parallel for 
    for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = mrModelPart.NodesBegin() + i;
        it_node->GetSolutionStepValue(NODAL_H, 0) = std::numeric_limits<double>::max();
    }
    
    for(unsigned int i=0; i<mrModelPart.Elements().size(); ++i) {
        auto it_element = mrModelPart.ElementsBegin() + i;
        auto& geom = it_element->GetGeometry();
        
        for(unsigned int k=0; k<geom.size()-1; ++k) {
            double& h1 = geom[k].FastGetSolutionStepValue(NODAL_H);
            for(unsigned int l=k+1; l<geom.size(); ++l) {
                double hedge = norm_2(geom[l].Coordinates() - geom[k].Coordinates());
                double& h2 = geom[l].FastGetSolutionStepValue(NODAL_H);
                
                // Get minimum between the existent value and the considered edge length 
                geom[k].FastGetSolutionStepValue(NODAL_H) = std::min(h1, hedge);
                geom[l].FastGetSolutionStepValue(NODAL_H) = std::min(h2, hedge);
            }
        }
    }
    
    mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(NODAL_H);

    KRATOS_CATCH("")
} // class FindNodalHProcess
} // namespace Kratos
