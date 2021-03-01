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

template<typename TVariable>
class CheckVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    const ModelPart& rModelPart,
                    bool& IsVariableFound)
    {
        IsVariableFound = rModelPart.HasNodalSolutionStepVariable(rVariable);
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

void NodalSolutionStepDataIO::WriteNodalResults(ModelPart& rModelPart, unsigned Step)
{
    KRATOS_TRY;

    if (mVariableNames.size() == 0)
        return;

    if (mVariableNames.size() == 1 && mVariableNames[0] == "ALL_VARIABLES_FROM_VARIABLES_LIST") {
        // Here we can change the original list because, solution step variables
        // list is same for all the time steps. So we have to check this only
        // once, and it will be the same for all steps.
        mVariableNames.clear();

        const auto& variable_data_list = rModelPart.GetNodalSolutionStepVariablesList();
        for (const auto& variable_data : variable_data_list) {
            if (variable_data.IsComponent()) {
                const auto& source_variable = variable_data.GetSourceVariable();
                if (std::find(mVariableNames.begin(), mVariableNames.end(), source_variable.Name()) == mVariableNames.end()) {
                    mVariableNames.push_back(source_variable.Name());
                }
            } else {
                mVariableNames.push_back(variable_data.Name());
            }
        }
    }

    WriteInfo info;
    std::vector<NodeType*> local_nodes;
    GetLocalNodes(rModelPart.Nodes(), local_nodes);

    // Write each variable.
    for (const std::string& r_variable_name : mVariableNames)
        RegisteredComponentLookup<Variable<array_1d<double, 3>>,
                                 Variable<double>, Variable<int>>(r_variable_name)
            .Execute<WriteNodalVariableFunctor>(local_nodes, Step, *mpFile, mPrefix, info);

    // Write block partition.
    WritePartitionTable(*mpFile, mPrefix + "/NodalSolutionStepData", info);

    KRATOS_CATCH("");
}

void NodalSolutionStepDataIO::ReadNodalResults(ModelPart& rModelPart, unsigned Step)
{
    KRATOS_TRY;

    if (mVariableNames.size() == 0)
        return;

    const std::string& current_path = mPrefix + "/NodalSolutionStepData";

    if (mVariableNames.size() == 1 && mVariableNames[0] == "ALL_VARIABLES_FROM_FILE") {
        if (mpFile->HasPath(current_path)) {
            // Here we can change the original list because, solution step variables
            // list is same for all the time steps. So we have to check this only
            // once, and it will be the same for all steps.
            mVariableNames.clear();
            const auto& r_variables_list = mpFile->GetDataSetNames(current_path);
            for (const auto& r_variable_name : r_variables_list) {
                bool is_variable_found = false;
                RegisteredComponentLookup<
                    Variable<array_1d<double, 3>>,
                    Variable<double>,
                    Variable<int>
                    >(r_variable_name)
                    .Execute<CheckVariableFunctor>(rModelPart, is_variable_found);

                if (is_variable_found) {
                    mVariableNames.push_back(r_variable_name);
                } else {
                    KRATOS_WARNING("NodalSolutionStepDataIO")
                        << "Skipping reading of " << r_variable_name
                        << " from " << mpFile->GetFileName() << ". It is not found in solution step variables list of "
                        << rModelPart.Name() << ".\n";
                }
            }
        } else {
            KRATOS_WARNING("NodalSolutionStepDataIO")
                << "ALL_VARIABLES_FROM_FILE are specified to be read, but no variable "
                   "data is found at "
                << current_path << " in " << mpFile->GetFileName() <<".\n";
            return;
        }
    }

    std::vector<NodeType*> local_nodes;
    GetLocalNodes(rModelPart.Nodes(), local_nodes);
    unsigned start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, current_path);

    // Read local data for each variable.
    for (const std::string& r_variable_name : mVariableNames)
        RegisteredComponentLookup<Variable<array_1d<double, 3>>,
                                 Variable<double>, Variable<int>>(r_variable_name)
            .Execute<ReadNodalVariableFunctor>(local_nodes, Step, *mpFile, mPrefix,
                                          start_index, block_size);

    // Synchronize ghost nodes.
    rModelPart.GetCommunicator().SynchronizeNodalSolutionStepsData();

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