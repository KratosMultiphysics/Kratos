//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "containers/system_vector.h"
#include "containers/distributed_system_vector.h"
#include "utilities/parallel_utilities.h"
#include "interface_vec_container.h"

namespace Kratos
{

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
template<class TVector>
void InterfaceVecContainer<TVector>::UpdateSystemVectorFromModelPart(
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions)
{
    const int num_nodes_local = mrModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = mrModelPart.GetCommunicator().LocalMesh().NodesBegin();

    IndexPartition<>(num_nodes_local).for_each(
        [&, this](std::size_t i){
            (*mpInterfaceVector)[i] = (nodes_begin + i)->FastGetSolutionStepValue(rVariable);
        }
    );
}

template<class TVector>
void InterfaceVecContainer<TVector>::UpdateModelPartFromSystemVector(
    const Variable<double>& rVariable,
    const Kratos::Flags& rMappingOptions)
{
    const int num_nodes_local = mrModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = mrModelPart.GetCommunicator().LocalMesh().NodesBegin();

    IndexPartition<>(num_nodes_local).for_each(
        [&, this](std::size_t i){
            (nodes_begin + i)->FastGetSolutionStepValue(rVariable) = (*mpInterfaceVector)[i];
        }
    );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class InterfaceVecContainer< SystemVector<> >;
template class InterfaceVecContainer< DistributedSystemVector<> >;

}  // namespace Kratos.
