//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


// System includes


// External includes


// Project includes
#include "utilities/parallel_utilities.h"
#include "mapping_matrix_builder.h"


namespace Kratos
{

template<bool IsDistributed>
void MappingMatrixBuilder<IsDistributed>::BuildMappingMatrix(
    MapperLocalSystemPointerVector& rLocalSystems,
    MappingMatrixPointerType& rpMappingMatrix,
    InterfaceVectorPointerType& rpInterfaceVectorOrigin,
    InterfaceVectorPointerType& rpInterfaceVectorDestination)
{
    GraphType matrix_graph;
    GraphType vec_origin_graph;
    GraphType vec_destination_graph;

    MapperLocalSystem::EquationIdVectorType origin_ids;
    MapperLocalSystem::EquationIdVectorType destination_ids;

    for (auto& r_local_sys : rLocalSystems) {
        r_local_sys->EquationIdVectors(origin_ids, destination_ids);
        matrix_graph.AddEntries(origin_ids, destination_ids);
        vec_origin_graph.AddEntries(origin_ids);
        vec_destination_graph.AddEntries(destination_ids);
    }

    MappingMatrixPointerType p_matrix(Kratos::make_unique<MappingMatrixType>(matrix_graph));
    rpMappingMatrix.swap(p_matrix);

    InterfaceVectorPointerType p_vec_o(Kratos::make_unique<InterfaceVectorType>(vec_origin_graph));
    rpInterfaceVectorOrigin.swap(p_vec_o);

    InterfaceVectorPointerType p_vec_d(Kratos::make_unique<InterfaceVectorType>(vec_destination_graph));
    rpInterfaceVectorDestination.swap(p_vec_d);
}



///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MappingMatrixBuilder< false >; // non-distributed version
// template class MappingMatrixBuilder< true >; // distributed version

}  // namespace Kratos.


