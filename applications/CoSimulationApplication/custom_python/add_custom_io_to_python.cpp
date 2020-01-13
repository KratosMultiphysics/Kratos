// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_python/add_custom_io_to_python.h"

// IO
#include "custom_io/co_sim_EMPIRE_API.h"

namespace Kratos {
namespace Python {

namespace EMPIRE_API_Wrappers { // helpers namespace

template<bool TIsDataField>
void SendArray(const std::string& rName, const int sizeOfArray, const std::vector<double>& signal)
{
    // Wrapper is needed bcs pybind cannot do the conversion to raw-ptr automatically
    if (TIsDataField) {
        EMPIRE_API_sendDataField(rName.c_str(), sizeOfArray, signal.data());
    } else {
        EMPIRE_API_sendSignal_double(rName.c_str(), sizeOfArray, signal.data());
    }
}

template<bool TIsDataField>
void ReceiveArray(const std::string& rName, const int sizeOfArray, pybind11::list signal)
{
    KRATOS_ERROR_IF(static_cast<int>(signal.size()) != sizeOfArray) << "The size of the list has to be specified before, expected size of " << sizeOfArray << ", current size: " << signal.size() << std::endl;

    // Wrapper is needed bcs pybind cannot do the conversion to raw-ptr automatically
    // also the list can only be modified in place otherwise the references are not working
    std::vector<double> vec_signal(sizeOfArray);
    if (TIsDataField) {
        EMPIRE_API_recvDataField(rName.c_str(), sizeOfArray, vec_signal.data());
    } else {
        EMPIRE_API_recvSignal_double(rName.c_str(), sizeOfArray, vec_signal.data());
    }

    // copy back the received values
    for (int i=0; i<sizeOfArray; ++i) {
        signal[i] = vec_signal[i];
    }
}

void sendDataField_scalar(const ModelPart& rModelPart, const std::string& rName, const Variable<double>& rVariable)
{
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rVariable)) << "Missing nodal solutionstepvariable: " << rVariable.Name() << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    const int size = rModelPart.NumberOfNodes();
    std::vector<double> values(size);

    std::size_t counter=0;
    for (const auto& r_node : rModelPart.Nodes()) {
        values[counter++] = r_node.FastGetSolutionStepValue(rVariable);
    }

    EMPIRE_API_sendDataField(rName.c_str(), size, values.data());
}

void sendDataField_scalar_DefaultName(ModelPart& rModelPart, const Variable<double>& rVariable)
{
    sendDataField_scalar(rModelPart, rVariable.Name(), rVariable);
}

void recvDataField_scalar(ModelPart& rModelPart, const std::string rName, const Variable<double>& rVariable)
{
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rVariable)) << "Missing nodal solutionstepvariable: " << rVariable.Name() << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    const int size = rModelPart.NumberOfNodes();
    std::vector<double> values(size);

    EMPIRE_API_recvDataField(rName.c_str(), size, values.data());

    std::size_t counter=0;
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(rVariable) = values[counter++];
    }
}

void recvDataField_scalar_DefaultName(ModelPart& rModelPart, const Variable<double>& rVariable)
{
    recvDataField_scalar(rModelPart, rVariable.Name(), rVariable);
}


void sendDataField_vector(const ModelPart& rModelPart, const std::string& rName, const Variable< array_1d<double, 3> >& rVariable)
{
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rVariable)) << "Missing nodal solutionstepvariable: " << rVariable.Name() << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    const int size = rModelPart.NumberOfNodes()*3;
    std::vector<double> values(size);

    std::size_t counter=0;
    for (const auto& r_node : rModelPart.Nodes()) {
        const array_1d<double, 3>& r_val = r_node.FastGetSolutionStepValue(rVariable);
        values[counter++] = r_val[0];
        values[counter++] = r_val[1];
        values[counter++] = r_val[2];
    }

    EMPIRE_API_sendDataField(rName.c_str(), size, values.data());
}

void sendDataField_vector_DefaultName(ModelPart& rModelPart, const Variable< array_1d<double, 3> >& rVariable)
{
    sendDataField_vector(rModelPart, rVariable.Name(), rVariable);
}

void recvDataField_vector(ModelPart& rModelPart, const std::string& rName, const Variable< array_1d<double, 3> >& rVariable)
{
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rVariable)) << "Missing nodal solutionstepvariable: " << rVariable.Name() << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    const int size = rModelPart.NumberOfNodes()*3;
    std::vector<double> values(size);

    EMPIRE_API_recvDataField(rName.c_str(), size, values.data());

    std::size_t counter=0;
    for (auto& r_node : rModelPart.Nodes()) {
        array_1d<double, 3>& r_val = r_node.FastGetSolutionStepValue(rVariable);
        r_val[0] = values[counter++];
        r_val[1] = values[counter++];
        r_val[2] = values[counter++];
    }
}

void recvDataField_vector_DefaultName(ModelPart& rModelPart, const Variable< array_1d<double, 3> >& rVariable)
{
    recvDataField_vector(rModelPart, rVariable.Name(), rVariable);
}


void sendDataField_doubleVector(const ModelPart& rModelPart, const std::string& rName, const Variable< array_1d<double, 3> >& rVariable1, const Variable< array_1d<double, 3> >& rVariable2)
{
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rVariable1)) << "Missing nodal solutionstepvariable: " << rVariable1.Name() << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rVariable2)) << "Missing nodal solutionstepvariable: " << rVariable2.Name() << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    const int size = rModelPart.NumberOfNodes()*6;
    std::vector<double> values(size);

    std::size_t counter=0;
    for (const auto& r_node : rModelPart.Nodes()) {
        const array_1d<double, 3>& r_val_1 = r_node.FastGetSolutionStepValue(rVariable1);
        const array_1d<double, 3>& r_val_2 = r_node.FastGetSolutionStepValue(rVariable2);
        values[counter++] = r_val_1[0];
        values[counter++] = r_val_1[1];
        values[counter++] = r_val_1[2];
        values[counter++] = r_val_2[0];
        values[counter++] = r_val_2[1];
        values[counter++] = r_val_2[2];
    }

    EMPIRE_API_sendDataField(rName.c_str(), size, values.data());
}

void sendDataField_doubleVector_DefaultName(ModelPart& rModelPart, const Variable< array_1d<double, 3> >& rVariable1, const Variable< array_1d<double, 3> >& rVariable2)
{
    const std::string name = rVariable1.Name() + rVariable2.Name();
    sendDataField_doubleVector(rModelPart, name, rVariable1, rVariable2);
}

void recvDataField_doubleVector(ModelPart& rModelPart, const std::string& rName, const Variable< array_1d<double, 3> >& rVariable1, const Variable< array_1d<double, 3> >& rVariable2)
{
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rVariable1)) << "Missing nodal solutionstepvariable: " << rVariable1.Name() << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rVariable2)) << "Missing nodal solutionstepvariable: " << rVariable2.Name() << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    const int size = rModelPart.NumberOfNodes()*6;
    std::vector<double> values(size);

    EMPIRE_API_recvDataField(rName.c_str(), size, values.data());

    std::size_t counter=0;
    for (auto& r_node : rModelPart.Nodes()) {
        array_1d<double, 3>& r_val_1 = r_node.FastGetSolutionStepValue(rVariable1);
        array_1d<double, 3>& r_val_2 = r_node.FastGetSolutionStepValue(rVariable2);
        r_val_1[0] = values[counter++];
        r_val_1[1] = values[counter++];
        r_val_1[2] = values[counter++];
        r_val_2[0] = values[counter++];
        r_val_2[1] = values[counter++];
        r_val_2[2] = values[counter++];
    }
}

void recvDataField_doubleVector_DefaultName(ModelPart& rModelPart, const Variable< array_1d<double, 3> >& rVariable1, const Variable< array_1d<double, 3> >& rVariable2)
{
    const std::string name = rVariable1.Name() + rVariable2.Name();
    recvDataField_doubleVector(rModelPart, name, rVariable1, rVariable2);
}

void sendMesh(const ModelPart& rModelPart, const std::string& rName, const bool UseConditions)
{
    // extract information from ModelPart
    const int numNodes = rModelPart.NumberOfNodes();
    const int numElems = (UseConditions) ? (rModelPart.NumberOfConditions()) : (rModelPart.NumberOfElements());

    std::vector<double> nodes(numNodes*3);
    std::vector<int> nodeIDs(numNodes);
    std::size_t node_counter = 0;
    for (const auto& r_node : rModelPart.Nodes()) {
        const auto& r_coords = r_node.GetInitialPosition(); // TODO or current coords?
        nodes[node_counter*3]   = r_coords[0];
        nodes[node_counter*3+1] = r_coords[1];
        nodes[node_counter*3+2] = r_coords[2];
        nodeIDs[node_counter++] = r_node.Id();
    }

    std::vector<int> numNodesPerElem(numElems);
    std::vector<int> elems;
    elems.reserve(numElems*3);
    std::size_t elem_counter = 0;
    if (UseConditions) {
        for (const auto& r_elem : rModelPart.Conditions()) {
            const auto& r_geom = r_elem.GetGeometry();
            numNodesPerElem[elem_counter++] = r_geom.PointsNumber();
            for (const auto& r_node : r_geom) {
                elems.push_back(r_node.Id());
            }
        }
    } else {
        for (const auto& r_elem : rModelPart.Elements()) {
            const auto& r_geom = r_elem.GetGeometry();
            numNodesPerElem[elem_counter++] = r_geom.PointsNumber();
            for (const auto& r_node : r_geom) {
                elems.push_back(r_node.Id());
            }
        }
    }

    EMPIRE_API_sendMesh(rName.c_str(), numNodes, numElems, nodes.data(), nodeIDs.data(), numNodesPerElem.data(), elems.data());
}

void sendMesh_DefaultName(ModelPart& rModelPart, const bool UseConditions)
{
    sendMesh(rModelPart, rModelPart.Name(), UseConditions);
}

template <typename TDouble, typename TInt>
static void createModelPartFromReceivedMesh(const int numNodes, const int numElems, const TDouble* nodes, const TInt* nodeIDs, const TInt* numNodesPerElem, const TInt* elem, ModelPart& rModelPart, const bool UseConditions)
{
    KRATOS_ERROR_IF(rModelPart.NumberOfNodes() > 0) << "ModelPart is not empty, it has nodes!" << std::endl;
    KRATOS_ERROR_IF(rModelPart.NumberOfProperties() > 0) << "ModelPart is not empty, it has properties!" << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    const std::unordered_map<int, std::string> element_name_map = {
        // {1 , "Element3D1N"}, // does not yet exist
        {2 , "Element3D2N"},
        {3 , "Element3D3N"},
        {4 , "Element3D4N"}
    };
    const std::unordered_map<int, std::string> condition_name_map = {
        {1 , "PointCondition3D1N"},
        {2 , "LineCondition3D2N"},
        {3 , "SurfaceCondition3D3N"},
        {4 , "SurfaceCondition3D4N"}
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
        for (int j=0; j<num_nodes_elem; ++j) {
            elem_node_ids[j] = (*elem)[counter++];
        }
        if (UseConditions) {
            rModelPart.CreateNewCondition(condition_name_map.at(num_nodes_elem), i+1, elem_node_ids, p_props);
        } else {
            rModelPart.CreateNewElement(element_name_map.at(num_nodes_elem), i+1, elem_node_ids, p_props);
        }
    }
}

void recvMesh(ModelPart& rModelPart, const std::string& rName, const bool UseConditions, const bool UseRawPointers)
{
    if (UseRawPointers) {
        int numNodes;
        int numElems;
        double** nodes = new double*[1];
        int** nodeIDs = new int*[1];
        int** numNodesPerElem = new int*[1];
        int** elem = new int*[1];

        EMPIRE_API_recvMesh(rName.c_str(), &numNodes, &numElems, nodes, nodeIDs, numNodesPerElem, elem);

        createModelPartFromReceivedMesh(numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elem, rModelPart, UseConditions);

        // deallocating memory
        delete [] nodes[0];
        delete [] nodeIDs[0];
        delete [] numNodesPerElem[0];
        delete [] elem[0];

        delete [] nodes;
        delete [] nodeIDs;
        delete [] numNodesPerElem;
        delete [] elem;

    } else {
        int numNodes;
        int numElems;
        std::vector<double> nodes;
        std::vector<int> nodeIDs;
        std::vector<int> numNodesPerElem;
        std::vector<int> elem;

        EMPIRE_API_recvMesh(rName.c_str(), &numNodes, &numElems, &nodes, &nodeIDs, &numNodesPerElem, &elem);

        createModelPartFromReceivedMesh(numNodes, numElems, &nodes, &nodeIDs, &numNodesPerElem, &elem, rModelPart, UseConditions);
    }
}

void recvMesh_DefaultName(ModelPart& rModelPart, const bool UseConditions, const bool UseRawPointers)
{
    recvMesh(rModelPart, rModelPart.Name().c_str(), UseConditions, UseRawPointers);
}

void SetEchoLevel(const int EchoLevel)
{
    EMPIRE_API_helpers::EchoLevel = EchoLevel;
}

void SetPrintTiming(const bool PrintTiming)
{
    EMPIRE_API_helpers::PrintTiming = PrintTiming;
}

} // helpers namespace

void  AddCustomIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto mEMPIREAPI = m.def_submodule("EMPIRE_API");

    mEMPIREAPI.def("EMPIRE_API_Connect", EMPIRE_API_Connect);
    mEMPIREAPI.def("EMPIRE_API_Disconnect", EMPIRE_API_Disconnect);

    mEMPIREAPI.def("EMPIRE_API_getUserDefinedText", EMPIRE_API_getUserDefinedText);

    mEMPIREAPI.def("EMPIRE_API_sendMesh", EMPIRE_API_Wrappers::sendMesh, py::arg("model_part"), py::arg("name"), py::arg("use_conditions")=false);
    mEMPIREAPI.def("EMPIRE_API_recvMesh", EMPIRE_API_Wrappers::recvMesh, py::arg("model_part"), py::arg("name"), py::arg("use_conditions")=false, py::arg("use_raw_pointers")=false);

    mEMPIREAPI.def("EMPIRE_API_sendMesh", EMPIRE_API_Wrappers::sendMesh_DefaultName, py::arg("model_part"), py::arg("use_conditions")=false);
    mEMPIREAPI.def("EMPIRE_API_recvMesh", EMPIRE_API_Wrappers::recvMesh_DefaultName, py::arg("model_part"), py::arg("use_conditions")=false, py::arg("use_raw_pointers")=false);

    mEMPIREAPI.def("EMPIRE_API_sendDataField", EMPIRE_API_Wrappers::SendArray<true>);
    mEMPIREAPI.def("EMPIRE_API_recvDataField", EMPIRE_API_Wrappers::ReceiveArray<true>);

    mEMPIREAPI.def("EMPIRE_API_sendDataField", EMPIRE_API_Wrappers::sendDataField_scalar);
    mEMPIREAPI.def("EMPIRE_API_recvDataField", EMPIRE_API_Wrappers::recvDataField_scalar);
    mEMPIREAPI.def("EMPIRE_API_sendDataField", EMPIRE_API_Wrappers::sendDataField_scalar_DefaultName);
    mEMPIREAPI.def("EMPIRE_API_recvDataField", EMPIRE_API_Wrappers::recvDataField_scalar_DefaultName);

    mEMPIREAPI.def("EMPIRE_API_sendDataField", EMPIRE_API_Wrappers::sendDataField_vector);
    mEMPIREAPI.def("EMPIRE_API_recvDataField", EMPIRE_API_Wrappers::recvDataField_vector);
    mEMPIREAPI.def("EMPIRE_API_sendDataField", EMPIRE_API_Wrappers::sendDataField_vector_DefaultName);
    mEMPIREAPI.def("EMPIRE_API_recvDataField", EMPIRE_API_Wrappers::recvDataField_vector_DefaultName);

    mEMPIREAPI.def("EMPIRE_API_sendDataField", EMPIRE_API_Wrappers::sendDataField_doubleVector);
    mEMPIREAPI.def("EMPIRE_API_recvDataField", EMPIRE_API_Wrappers::recvDataField_doubleVector);
    mEMPIREAPI.def("EMPIRE_API_sendDataField", EMPIRE_API_Wrappers::sendDataField_doubleVector_DefaultName);
    mEMPIREAPI.def("EMPIRE_API_recvDataField", EMPIRE_API_Wrappers::recvDataField_doubleVector_DefaultName);

    mEMPIREAPI.def("EMPIRE_API_sendSignal_double", EMPIRE_API_Wrappers::SendArray<false>);
    mEMPIREAPI.def("EMPIRE_API_recvSignal_double", EMPIRE_API_Wrappers::ReceiveArray<false>);

    mEMPIREAPI.def("EMPIRE_API_recvConvergenceSignal", EMPIRE_API_recvConvergenceSignal, py::arg("file_name_extension")="default");
    mEMPIREAPI.def("EMPIRE_API_sendConvergenceSignal", EMPIRE_API_sendConvergenceSignal, py::arg("signal"), py::arg("file_name_extension")="default");

    mEMPIREAPI.def("EMPIRE_API_SetEchoLevel", EMPIRE_API_Wrappers::SetEchoLevel);
    mEMPIREAPI.def("EMPIRE_API_PrintTiming",  EMPIRE_API_Wrappers::SetPrintTiming);
}

}  // namespace Python.
} // Namespace Kratos

