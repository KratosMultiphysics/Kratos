//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "processes/process.h"
#include "python/add_global_pointers_to_python.h"
#include "containers/global_pointers_vector.h"
#include "utilities/pointer_communicator.h"
#include "processes/find_global_nodal_neighbours_process.h"

namespace Kratos {
namespace Python {

std::unordered_map<int, std::vector<int> > GetNeighbourIds(
            FindGlobalNodalNeighboursProcess& self, 
            const DataCommunicator& rComm,
            ModelPart::NodesContainerType& rNodes
            )
{
    std::unordered_map<int, std::vector<int> > output;

    GlobalPointersVector< Node<3> > gp_list;
    for(auto& node : rNodes)
        for(auto& gp : node.GetValue(NEIGHBOUR_NODES))
            gp_list.push_back(gp);
    gp_list.Unique();

    GlobalPointerCommunicator<Node<3>> pointer_comm(rComm, gp_list);
    auto result_proxy = pointer_comm.Apply<int>(
            [](GlobalPointer<Node<3>>& gp){return gp->Id();}
    );

    for(auto& node : rNodes)
    {
        auto& neighbours = node.GetValue(NEIGHBOUR_NODES);
        std::vector<int> tmp(neighbours.size());
        for(unsigned int i=0; i<neighbours.size(); ++i)
        {
            tmp[i] = result_proxy.Get(neighbours[i]);
        }
        output[node.Id()] = tmp;   
    }

    return output;
}


void AddGlobalPointersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // NOTE: this function is special in that it accepts a "pyObject" - this is the reason for which it is defined in this same file
    py::class_< GlobalPointer<Node<3>>  >(m,"GlobalNodePointer");
    py::class_< GlobalPointer<Element>  >(m,"GlobalElementPointer");
    py::class_< GlobalPointer<Condition>  >(m,"GlobalConditionPointer");

    py::class_< GlobalPointersVector<Node<3>> >(m,"GlobalNodePointersVector");

    py::class_<FindGlobalNodalNeighboursProcess, FindGlobalNodalNeighboursProcess::Pointer, Process>
        (m,"FindGlobalNodalNeighboursProcess")
            .def(py::init<const DataCommunicator&, ModelPart&>())
    .def("ClearNeighbours",&FindGlobalNodalNeighboursProcess::ClearNeighbours)
    .def("GetNeighbourIds",GetNeighbourIds)
    ;

}

} // namespace Python.
} // Namespace Kratos
