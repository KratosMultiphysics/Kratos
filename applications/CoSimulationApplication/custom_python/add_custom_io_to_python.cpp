// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "custom_python/add_custom_io_to_python.h"
#include "includes/model_part.h"

// IO
#include "custom_io/co_sim_EMPIRE_API.h"

namespace Kratos {
namespace Python {

template<bool TIsDataField>
void Wrapper_SendArray(char* name, int sizeOfArray, std::vector<double> signal)
{
    // Wrapper is needed bcs pybind cannot do the conversion to raw-ptr automatically
    if (TIsDataField) {
        EMPIRE_API_sendDataField(name, sizeOfArray, &signal[0]);
    } else {
        EMPIRE_API_sendSignal_double(name, sizeOfArray, &signal[0]);
    }
}

template<bool TIsDataField>
void Wrapper_ReceiveArray(char* name, int sizeOfArray, pybind11::list signal)
{
    KRATOS_ERROR_IF(static_cast<int>(signal.size()) != sizeOfArray) << "The size of the list has to be specified before, expected size of " << sizeOfArray << ", current size: " << signal.size() << std::endl;

    // Wrapper is needed bcs pybind cannot do the conversion to raw-ptr automatically
    // also the list can only be modified in place otherwise the references are not working
    std::vector<double> vec_signal(sizeOfArray);
    if (TIsDataField) {
        EMPIRE_API_recvDataField(name, sizeOfArray, &vec_signal[0]);
    } else {
        EMPIRE_API_recvSignal_double(name, sizeOfArray, &vec_signal[0]);
    }

    // copy back the received values
    for (int i=0; i<sizeOfArray; ++i) {
        signal[i] = vec_signal[i];
    }
}

void Wrapper_EMPIRE_API_sendMesh(const ModelPart& rModelPart)
{
    // extract information from ModelPart
    const int numNodes = rModelPart.NumberOfNodes();
    const int numElems = rModelPart.NumberOfElements();

    std::vector<double> nodes(numNodes*3);
    std::vector<int> nodeIDs(numNodes);
    std::size_t node_counter = 0;
    for (const auto& r_node : rModelPart.Nodes()) {
        const auto& r_coords = r_node.GetInitialPosition();
        nodes[node_counter*3]   = r_coords[0];
        nodes[node_counter*3+1] = r_coords[1];
        nodes[node_counter*3+2] = r_coords[2];
        nodeIDs[node_counter++] = r_node.Id();
    }

    std::vector<int> numNodesPerElem(numElems);
    std::vector<int> elems;
    elems.reserve(numElems*3);
    std::size_t elem_counter = 0;
    for (const auto& r_elem : rModelPart.Elements()) {
        const auto& r_geom = r_elem.GetGeometry();
        numNodesPerElem[elem_counter++] = r_geom.PointsNumber();
        for (const auto& r_node : r_geom) {
            elems.push_back(r_node.Id());
        }
    }

    EMPIRE_API_sendMesh(const_cast<char*>(rModelPart.Name().c_str()), numNodes, numElems, &nodes[0], &nodeIDs[0], &numNodesPerElem[0], &elems[0]);
}

void Wrapper_EMPIRE_API_recvMesh(ModelPart& rModelPart)
{
    KRATOS_ERROR_IF(rModelPart.NumberOfNodes() > 0) << "ModelPart is not empty, it has nodes!" << std::endl;
    KRATOS_ERROR_IF(rModelPart.NumberOfProperties() > 0) << "ModelPart is not empty, it has properties!" << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed() > 0) << "ModelPart cannot be distributed!" << std::endl;

    int numNodes;
    int numElems;
    double** nodes = new double*;
    int** nodeIDs = new int*;
    int** numNodesPerElem = new int*;
    int** elem = new int*;

    EMPIRE_API_recvMesh(const_cast<char*>(rModelPart.Name().c_str()), &numNodes, &numElems, nodes, nodeIDs, numNodesPerElem, elem);

    const std::unordered_map<int, std::string> element_name_map = {
        // {1 : "Element2D1N"}, // does not yet exist
        {2 , "Element2D2N"},
        {3 , "Element2D3N"},
        {4 , "Element2D4N"}
    };

    // fill ModelPart with received entities
    for (int i=0; i<numNodes; ++i) {
        rModelPart.CreateNewNode((*nodeIDs)[i], (*nodes)[i*3], (*nodes)[i*3+1], (*nodes)[i*3+2]);
    }

    auto p_props = rModelPart.CreateNewProperties(0);

    int counter=0;
    for (int i=0; i<numElems; ++i) {
        const int num_nodes_elem = (*numNodesPerElem)[i];
        std::vector<ModelPart::IndexType> elem_node_ids(num_nodes_elem);
        for (int j=0; j<numElems; ++j) {
            elem_node_ids[j] = (*elem)[counter++];
        }
        rModelPart.CreateNewElement(element_name_map.at(num_nodes_elem), i+1, elem_node_ids, p_props);
    }

    // TODO deallocate memory!
}

void  AddCustomIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto mEMPIREAPI = m.def_submodule("EMPIRE_API");

    mEMPIREAPI.def("EMPIRE_API_Connect", EMPIRE_API_Connect);
    mEMPIREAPI.def("EMPIRE_API_Disconnect", EMPIRE_API_Disconnect);

    mEMPIREAPI.def("EMPIRE_API_getUserDefinedText", EMPIRE_API_getUserDefinedText);

    mEMPIREAPI.def("EMPIRE_API_sendMesh", Wrapper_EMPIRE_API_sendMesh);
    mEMPIREAPI.def("EMPIRE_API_recvMesh", Wrapper_EMPIRE_API_recvMesh);

    mEMPIREAPI.def("EMPIRE_API_sendDataField", Wrapper_SendArray<true>);
    mEMPIREAPI.def("EMPIRE_API_recvDataField", Wrapper_ReceiveArray<true>);

    mEMPIREAPI.def("EMPIRE_API_sendSignal_double", Wrapper_SendArray<false>);
    mEMPIREAPI.def("EMPIRE_API_recvSignal_double", Wrapper_ReceiveArray<false>);

    mEMPIREAPI.def("EMPIRE_API_recvConvergenceSignal", EMPIRE_API_recvConvergenceSignal);
    mEMPIREAPI.def("EMPIRE_API_sendConvergenceSignal", EMPIRE_API_sendConvergenceSignal);
}

}  // namespace Python.
} // Namespace Kratos

