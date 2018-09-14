#include "custom_io/hdf5_nodal_data_value_io.h"

#include "utilities/openmp_utils.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"
#include "includes/kratos_parameters.h"
#include "includes/communicator.h"
#include "custom_utilities/registered_variable_lookup.h"
#include "custom_utilities/local_ghost_splitting_utility.h"

namespace Kratos
{
namespace HDF5
{

NodalDataValueIO::NodalDataValueIO(Parameters Settings, File::Pointer pFile)
: mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "prefix": "",
            "list_of_variables": []
        })");
    Settings.ValidateAndAssignDefaults(default_params);
    mPrefix = Settings["prefix"].GetString();
    mVariableNames.resize(Settings["list_of_variables"].size());
    for (std::size_t i = 0; i < mVariableNames.size(); ++i)
        mVariableNames[i] = Settings["list_of_variables"].GetArrayItem(i).GetString();

    KRATOS_CATCH("");
}

namespace {
template <typename TVariable>
class WriteNonHistoricalVariableFunctor;
}

void NodalDataValueIO::WriteNodalResults(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    if (mVariableNames.size() == 0)
        return;

    WriteInfo info;
    std::vector<NodeType*> local_nodes;
    GetLocalNodes(rNodes, local_nodes);

    // Write each variable.
    for (const std::string& r_variable_name : mVariableNames)
        RegisteredVariableLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                                 Variable<Vector<double>>, Variable<Matrix<double>>>(r_variable_name)
            .Execute<WriteNonHistoricalVariableFunctor>(local_nodes, *mpFile, mPrefix, info);

    // Write block partition.
    WritePartitionTable(*mpFile, mPrefix + "/NodalDataValues", info);

    KRATOS_CATCH("");
}

namespace {
template <typename TDataType>
void SetNodalDataValueBuffer(Variable<TDataType> const&,
                             std::vector<NodeType*> const&,
                             Vector<TDataType>&);

void SetNodalDataValueBuffer(Variable<Vector<double>> const&,
                             std::vector<NodeType*> const&,
                             Matrix<double>&);

void SetNodalDataValueBuffer(Variable<Matrix<double>> const&,
                             std::vector<NodeType*> const&,
                             Matrix<double>&);

template <typename TVariable>
class WriteNonHistoricalVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    std::vector<NodeType*>& rNodes,
                    File& rFile,
                    std::string const& rPrefix,
                    WriteInfo& rInfo)
    {
        Vector<typename TVariable::Type> data;
        SetNodalDataValueBuffer(rVariable, rNodes, data);
        rFile.WriteDataSet(rPrefix + "/NodalDataValues/" + rVariable.Name(), data, rInfo);
    }
};

template <>
class WriteNonHistoricalVariableFunctor<Variable<Vector<double>>>
{
public:
    void operator()(Variable<Vector<double>> const& rVariable,
                    std::vector<NodeType*>& rNodes,
                    File& rFile,
                    std::string const& rPrefix,
                    WriteInfo& rInfo)
    {
        Matrix<double> data;
        SetNodalDataValueBuffer(rVariable, rNodes, data);
        rFile.WriteDataSet(rPrefix + "/NodalDataValues/" + rVariable.Name(), data, rInfo);
    }
};

template <>
class WriteNonHistoricalVariableFunctor<Variable<Matrix<double>>>
{
public:
    void operator()(Variable<Matrix<double>> const& rVariable,
                    std::vector<NodeType*>& rNodes,
                    File& rFile,
                    std::string const& rPrefix,
                    WriteInfo& rInfo)
    {
        Matrix<double> data;
        SetNodalDataValueBuffer(rVariable, rNodes, data);
        rFile.WriteDataSet(rPrefix + "/NodalDataValues/" + rVariable.Name(), data, rInfo);
        const int size1 = rNodes.front()->GetValue(rVariable).size1();
        const int size2 = rNodes.front()->GetValue(rVariable).size2();
        rFile.WriteAttribute(
            rPrefix + "/NodalDataValues/" + rVariable.Name(), "Size1", size1);
        rFile.WriteAttribute(
            rPrefix + "/NodalDataValues/" + rVariable.Name(), "Size2", size2);
    }
};

template <typename TDataType>
void SetNodalDataValueBuffer(Variable<TDataType> const& rVariable,
                             std::vector<NodeType*> const& rNodes,
                             Vector<TDataType>& rData)
{
    KRATOS_TRY;

    rData.resize(rNodes.size(), false);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rNodes.size()); ++i)
    {
        const NodeType& r_node = *rNodes[i];
        rData[i] = r_node.GetValue(rVariable);
    }

    KRATOS_CATCH("");
}

void SetNodalDataValueBuffer(Variable<Vector<double>> const& rVariable,
                             std::vector<NodeType*> const& rNodes,
                             Matrix<double>& rData)
{
    KRATOS_TRY;

    rData.resize(rNodes.size(), rNodes.front()->GetValue(rVariable).size(), false);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rNodes.size()); ++i)
    {
        const NodeType& r_node = *rNodes[i];
        Vector<double> const& r_vec = r_node.GetValue(rVariable);
        KRATOS_ERROR_IF(r_vec.size() != rData.size2())
            << "Invalid dimension.\n";
        for (std::size_t j = 0; j < rData.size2(); ++j)
            rData(i, j) = r_vec(j);
    }

    KRATOS_CATCH("");
}

void SetNodalDataValueBuffer(Variable<Matrix<double>> const& rVariable,
                             std::vector<NodeType*> const& rNodes,
                             Matrix<double>& rData)
{
    KRATOS_TRY;

    rData.resize(rNodes.size(), rNodes.front()->GetValue(rVariable).data().size(), false);
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rNodes.size()); ++i)
    {
        const NodeType& r_node = *rNodes[i];
        const auto& r_data = r_node.GetValue(rVariable).data();
        KRATOS_ERROR_IF(r_data.size() != rData.size2())
            << "Invalid dimension.\n";
        for (std::size_t j = 0; j < rData.size2(); ++j)
            rData(i, j) = r_data[j];
    }

    KRATOS_CATCH("");
}

template <typename TVariable>
class ReadNonHistoricalVariableFunctor;
}

void NodalDataValueIO::ReadNodalResults(NodesContainerType& rNodes, Communicator& rComm)
{
    KRATOS_TRY;

    if (mVariableNames.size() == 0)
        return;

    std::vector<NodeType*> local_nodes;
    std::vector<NodeType*> ghost_nodes;
    SplitNodesIntoLocalAndGhost(rNodes, local_nodes, ghost_nodes);
    unsigned start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, mPrefix + "/NodalDataValues");

    // Read local data for each variable.
    for (const std::string& r_variable_name : mVariableNames)
        RegisteredVariableLookup<Variable<array_1d<double, 3>>, Variable<double>, Variable<int>,
                                 Variable<Vector<double>>, Variable<Matrix<double>>>(r_variable_name)
            .Execute<ReadNonHistoricalVariableFunctor>(
                local_nodes, ghost_nodes, rComm, *mpFile, mPrefix, start_index, block_size);

    KRATOS_CATCH("");
}

namespace {
template <typename TDataType>
void SetNodalDataValue(Variable<TDataType> const&,
                       std::vector<NodeType*>&,
                       Vector<TDataType> const&);

void SetNodalDataValue(Variable<Vector<double>> const&,
                       std::vector<NodeType*>&,
                       Matrix<double> const&);

void SetNodalDataValue(Variable<Matrix<double>> const&,
                       std::vector<NodeType*>&,
                       Matrix<double> const&,
                       int,
                       int);

template <typename TDataType>
void ZeroNodalDataValue(Variable<TDataType> const&, std::vector<NodeType*>&);

template <typename TVariable>
class ReadNonHistoricalVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    std::vector<NodeType*>& rLocalNodes,
                    std::vector<NodeType*>& rGhostNodes,
                    Communicator& rComm,
                    File& rFile,
                    std::string const& rPrefix,
                    unsigned StartIndex,
                    unsigned BlockSize)
    {
        Vector<typename TVariable::Type> data;
        rFile.ReadDataSet(rPrefix + "/NodalDataValues/" + rVariable.Name(), data,
                          StartIndex, BlockSize);
        SetNodalDataValue(rVariable, rLocalNodes, data);
        ZeroNodalDataValue(rVariable, rGhostNodes);
        rComm.AssembleNonHistoricalData(rVariable);
    }
};

template <>
class ReadNonHistoricalVariableFunctor<Variable<Vector<double>>>
{
public:
    void operator()(Variable<Vector<double>> const& rVariable,
                    std::vector<NodeType*>& rLocalNodes,
                    std::vector<NodeType*>& rGhostNodes,
                    Communicator& rComm,
                    File& rFile,
                    std::string const& rPrefix,
                    unsigned StartIndex,
                    unsigned BlockSize)
    {
        Matrix<double> data;
        rFile.ReadDataSet(rPrefix + "/NodalDataValues/" + rVariable.Name(), data,
                          StartIndex, BlockSize);
        SetNodalDataValue(rVariable, rLocalNodes, data);
        ZeroNodalDataValue(rVariable, rGhostNodes);
        rComm.AssembleNonHistoricalData(rVariable);
    }
};

template <>
class ReadNonHistoricalVariableFunctor<Variable<Matrix<double>>>
{
public:
    void operator()(Variable<Matrix<double>> const& rVariable,
                    std::vector<NodeType*>& rLocalNodes,
                    std::vector<NodeType*>& rGhostNodes,
                    Communicator& rComm,
                    File& rFile,
                    std::string const& rPrefix,
                    unsigned StartIndex,
                    unsigned BlockSize)
    {
        Matrix<double> data;
        rFile.ReadDataSet(rPrefix + "/NodalDataValues/" + rVariable.Name(), data,
                          StartIndex, BlockSize);
        int size1, size2;
        rFile.ReadAttribute(
            rPrefix + "/NodalDataValues/" + rVariable.Name(), "Size1", size1);
        rFile.ReadAttribute(
            rPrefix + "/NodalDataValues/" + rVariable.Name(), "Size2", size2);
        SetNodalDataValue(rVariable, rLocalNodes, data, size1, size2);
        ZeroNodalDataValue(rVariable, rGhostNodes);
        rComm.AssembleNonHistoricalData(rVariable);
    }
};

template <typename TDataType>
void SetNodalDataValue(Variable<TDataType> const& rVariable,
                       std::vector<NodeType*>& rNodes,
                       Vector<TDataType> const& rData)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rNodes.size() != rData.size())
        << "Number of nodes does not equal data set dimension\n";
    for (std::size_t i = 0; i < rNodes.size(); ++i)
        rNodes[i]->GetValue(rVariable) = rData[i];

    KRATOS_CATCH("");
}

void SetNodalDataValue(Variable<Vector<double>> const& rVariable,
                       std::vector<NodeType*>& rNodes,
                       Matrix<double> const& rData)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rNodes.size() != rData.size1())
        << "Number of nodes does not equal data set dimension\n";
    for (std::size_t i = 0; i < rNodes.size(); ++i)
    {
        Vector<double>& r_vec = rNodes[i]->GetValue(rVariable);
        r_vec.resize(rData.size2(), false);
        for (std::size_t j = 0; j < rData.size2(); ++j)
            r_vec(j) = rData(i, j);
    }

    KRATOS_CATCH("");
}

void SetNodalDataValue(Variable<Matrix<double>> const& rVariable,
                       std::vector<NodeType*>& rNodes,
                       Matrix<double> const& rData,
                       int size1,
                       int size2)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rNodes.size() != rData.size1())
        << "Number of nodes does not equal data set dimension\n";
    for (std::size_t k = 0; k < rNodes.size(); ++k)
    {
        Matrix<double>& r_mat = rNodes[k]->GetValue(rVariable);
        r_mat.resize(size1, size2, false);
        for (int i = 0; i < size1; ++i)
            for (int j = 0; j < size2; ++j)
                r_mat(i, j) = rData(k, i * size2 + j);
    }

    KRATOS_CATCH("");
}

template <typename TDataType>
void ZeroNodalDataValue(Variable<TDataType> const& rVariable,
                        std::vector<NodeType*>& rNodes)
{
    for (auto& p_node : rNodes)
        p_node->GetValue(rVariable) = rVariable.Zero();
}
}

} // namespace HDF5.
} // namespace Kratos.
