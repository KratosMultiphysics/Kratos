#include "custom_io/hdf5_nodal_solution_step_data_io.h"

#include "utilities/openmp_utils.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"
#include "includes/kratos_parameters.h"
#include "includes/communicator.h"
#include "custom_utilities/registered_component_lookup.h"
#include "custom_utilities/local_ghost_splitting_utility.h"

namespace Kratos
{
namespace HDF5
{
namespace
{
template <class TVariableType, class TFileDataType>
void SetDataBuffer(TVariableType const& rVariable,
                   std::vector<NodeType*> const& rNodes,
                   Vector<TFileDataType>& rData,
                   unsigned Step);

template <class TVariableType, class TFileDataType>
void SetNodalSolutionStepData(TVariableType const& rVariable,
                              Vector<TFileDataType> const& rData,
                              std::vector<NodeType*>& rNodes,
                              unsigned Step);

template <typename TVariable>
class WriteNodalVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    std::vector<NodeType*>& rNodes,
                    unsigned Step,
                    File& rFile,
                    std::string const& rPrefix,
                    WriteInfo& rInfo)
    {
        Vector<typename TVariable::Type> data;
        SetDataBuffer(rVariable, rNodes, data, Step);
        rFile.WriteDataSet(rPrefix + "/NodalSolutionStepData/" + rVariable.Name(), data, rInfo);
    }
};

template <typename TVariable>
class ReadNodalVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    std::vector<NodeType*>& rNodes,
                    unsigned Step,
                    File& rFile,
                    std::string const& rPrefix,
                    unsigned StartIndex,
                    unsigned BlockSize)
    {
        Vector<typename TVariable::Type> data;
        rFile.ReadDataSet(rPrefix + "/NodalSolutionStepData/" + rVariable.Name(), data,
                          StartIndex, BlockSize);
        SetNodalSolutionStepData(rVariable, data, rNodes, Step);
    }
};
} // unnamed namespace

NodalSolutionStepDataIO::NodalSolutionStepDataIO(Parameters Settings, File::Pointer pFile)
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
    for (unsigned i = 0; i < mVariableNames.size(); ++i)
        mVariableNames[i] = Settings["list_of_variables"].GetArrayItem(i).GetString();

    KRATOS_CATCH("");
}

void NodalSolutionStepDataIO::WriteNodalResults(NodesContainerType const& rNodes, unsigned Step)
{
    KRATOS_TRY;

    if (mVariableNames.size() == 0)
        return;

    WriteInfo info;
    std::vector<NodeType*> local_nodes;
    GetLocalNodes(rNodes, local_nodes);

    // Write each variable.
    for (const std::string& r_variable_name : mVariableNames)
        RegisteredComponentLookup<Variable<array_1d<double, 3>>,
                                 VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>,
                                 Variable<double>, Variable<int>>(r_variable_name)
            .Execute<WriteNodalVariableFunctor>(local_nodes, Step, *mpFile, mPrefix, info);

    // Write block partition.
    WritePartitionTable(*mpFile, mPrefix + "/NodalSolutionStepData", info);

    KRATOS_CATCH("");
}

void NodalSolutionStepDataIO::ReadNodalResults(NodesContainerType& rNodes, Communicator& rComm, unsigned Step)
{
    KRATOS_TRY;

    if (mVariableNames.size() == 0)
        return;

    std::vector<NodeType*> local_nodes;
    GetLocalNodes(rNodes, local_nodes);
    unsigned start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, mPrefix + "/NodalSolutionStepData");

    // Read local data for each variable.
    for (const std::string& r_variable_name : mVariableNames)
        RegisteredComponentLookup<Variable<array_1d<double, 3>>,
                                 VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>,
                                 Variable<double>, Variable<int>>(r_variable_name)
            .Execute<ReadNodalVariableFunctor>(local_nodes, Step, *mpFile, mPrefix,
                                          start_index, block_size);

    // Synchronize ghost nodes.
    rComm.SynchronizeNodalSolutionStepsData();

    KRATOS_CATCH("");
}

namespace
{
template <class TVariableType, class TFileDataType>
void SetDataBuffer(TVariableType const& rVariable,
                   std::vector<NodeType*> const& rNodes,
                   Vector<TFileDataType>& rData,
                   unsigned Step)
{
    KRATOS_TRY;

    rData.resize(rNodes.size(), false);
    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(rNodes.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
            rData[i] = rNodes[i]->FastGetSolutionStepValue(rVariable, Step);
    }

    KRATOS_CATCH("");
}

template <class TVariableType, class TFileDataType>
void SetNodalSolutionStepData(TVariableType const& rVariable,
                              Vector<TFileDataType> const& rData,
                              std::vector<NodeType*>& rNodes,
                              unsigned Step)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(rData.size() != rNodes.size())
        << "File data block size (" << rData.size()
        << ") is not equal to number of nodes (" << rNodes.size() << ")." << std::endl;

    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(rNodes.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
            rNodes[i]->FastGetSolutionStepValue(rVariable, Step) = rData[i];
    }

    KRATOS_CATCH("");
}
} // unnamed namespace.
} // namespace HDF5.
} // namespace Kratos.