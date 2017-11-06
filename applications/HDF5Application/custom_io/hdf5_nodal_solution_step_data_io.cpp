#include "custom_io/hdf5_nodal_solution_step_data_io.h"

#include "includes/kratos_components.h"
#include "custom_utilities/hdf5_utils.h"

namespace Kratos
{
namespace HDF5
{
NodalSolutionStepDataIO::NodalSolutionStepDataIO(Parameters& rParams, File::Pointer pFile)
: mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "prefix": "",
            "list_of_variables": []
        })");

    rParams.ValidateAndAssignDefaults(default_params);

    mPrefix = rParams["prefix"].GetString();

    mVariableNames.resize(rParams["list_of_variables"].size());
    for (unsigned i = 0; i < mVariableNames.size(); ++i)
        mVariableNames[i] = rParams["list_of_variables"].GetArrayItem(i).GetString();

    KRATOS_CATCH("");
}

std::string NodalSolutionStepDataIO::GetPrefix() const
{
    return mPrefix;
}

void NodalSolutionStepDataIO::SetPrefix(std::string const& rPrefix)
{
    mPrefix = rPrefix;
}

void NodalSolutionStepDataIO::WriteNodalResults(NodesContainerType const& rNodes, unsigned Step)
{
    KRATOS_TRY;

    std::vector<NodeType*> local_nodes;
    Detail::GetLocalNodes(rNodes, local_nodes, IsPartitionedSimulation());

    // Write each variable.
    for (const std::string& r_variable_name : mVariableNames)
    {
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name))
        {
            const Variable<array_1d<double, 3>>& rVARIABLE =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name);
            Vector<array_1d<double, 3>> data;
            SetDataBuffer(rVARIABLE, local_nodes, data, Step);
            mpFile->WriteDataSet(mPrefix + "/NodalResults/" + r_variable_name, data);
        }
        else if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(
                     r_variable_name))
        {
            const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& rVARIABLE =
                KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                    r_variable_name);
            Vector<double> data;
            SetDataBuffer(rVARIABLE, local_nodes, data, Step);
            mpFile->WriteDataSet(mPrefix + "/NodalResults/" + r_variable_name, data);
        }
        else if (KratosComponents<Variable<double>>::Has(r_variable_name))
        {
            const Variable<double>& rVARIABLE =
                KratosComponents<Variable<double>>::Get(r_variable_name);
            Vector<double> data;
            SetDataBuffer(rVARIABLE, local_nodes, data, Step);
            mpFile->WriteDataSet(mPrefix + "/NodalResults/" + r_variable_name, data);
        }
        else if (KratosComponents<Variable<int>>::Has(r_variable_name))
        {
            const Variable<int>& rVARIABLE =
                KratosComponents<Variable<int>>::Get(r_variable_name);
            Vector<int> data;
            SetDataBuffer(rVARIABLE, local_nodes, data, Step);
            mpFile->WriteDataSet(mPrefix + "/NodalResults/" + r_variable_name, data);
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable type: " << r_variable_name << std::endl;
        }
    }

    // Write block partition.
    Vector<int> dummy(local_nodes.size());
    mpFile->WriteDataPartition(mPrefix + "/NodalResults/Partition", dummy);

    KRATOS_CATCH("");
}

void NodalSolutionStepDataIO::ReadNodalResults(NodesContainerType& rNodes, Communicator& rComm, unsigned Step)
{
    KRATOS_TRY;

    std::vector<NodeType*> local_nodes;
    Detail::GetLocalNodes(rNodes, local_nodes, IsPartitionedSimulation());
    unsigned start_index, block_size;
    std::tie(start_index, block_size) = GetStartIndexAndBlockSize();

    // Read local data for each variable.
    for (const std::string& r_variable_name : mVariableNames)
    {
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name))
        {
            Vector<array_1d<double, 3>> data;
            mpFile->ReadDataSet(mPrefix + "/NodalResults/" + r_variable_name,
                                data, start_index, block_size);
            const Variable<array_1d<double, 3>>& rVARIABLE =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name);
            SetNodalSolutionStepData(rVARIABLE, data, local_nodes, Step);
        }
        else if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(
                     r_variable_name))
        {
            Vector<double> data;
            mpFile->ReadDataSet(mPrefix + "/NodalResults/" + r_variable_name,
                                data, start_index, block_size);
            const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& rVARIABLE =
                KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                    r_variable_name);
            SetNodalSolutionStepData(rVARIABLE, data, local_nodes, Step);
        }
        else if (KratosComponents<Variable<double>>::Has(r_variable_name))
        {
            Vector<double> data;
            mpFile->ReadDataSet(mPrefix + "/NodalResults/" + r_variable_name,
                                data, start_index, block_size);
            const Variable<double>& rVARIABLE =
                KratosComponents<Variable<double>>::Get(r_variable_name);
            SetNodalSolutionStepData(rVARIABLE, data, local_nodes, Step);
        }
        else if (KratosComponents<Variable<int>>::Has(r_variable_name))
        {
            Vector<int> data;
            mpFile->ReadDataSet(mPrefix + "/NodalResults/" + r_variable_name,
                                data, start_index, block_size);
            const Variable<int>& rVARIABLE =
                KratosComponents<Variable<int>>::Get(r_variable_name);
            SetNodalSolutionStepData(rVARIABLE, data, local_nodes, Step);
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable type: " << r_variable_name << std::endl;
        }
    }

    // Synchronize ghost nodes.
    rComm.SynchronizeNodalSolutionStepsData();

    KRATOS_CATCH("");
}

bool NodalSolutionStepDataIO::IsPartitionedSimulation() const
{
    return (mpFile->GetTotalProcesses() > 1);
}

std::tuple<unsigned, unsigned> NodalSolutionStepDataIO::GetStartIndexAndBlockSize() const
{
    KRATOS_TRY;

    std::vector<unsigned> dims =
        mpFile->GetDataDimensions(mPrefix + "/NodalResults/Partition");
    KRATOS_ERROR_IF(dims.size() != 1) << "Invalid partition dimension." << std::endl;
    const unsigned file_partition_size = dims[0] - 1; // Number of procs that wrote the data block.
    unsigned start_index;
    unsigned block_size;
    if (IsPartitionedSimulation() == false)
    {
        start_index = 0;
        // Read the global size of the data block.
        Vector<int> last_partition_index;
        mpFile->ReadDataSet(mPrefix + "/NodalResults/Partition", last_partition_index, file_partition_size, 1);
        block_size = last_partition_index[0];
    }
    else if (mpFile->GetTotalProcesses() == file_partition_size)
    {
        Vector<int> my_partition;
        unsigned my_pid = mpFile->GetPID();
        mpFile->ReadDataSet(mPrefix + "/NodalResults/Partition", my_partition, my_pid, 2);
        start_index = my_partition[0];
        block_size = my_partition[1] - my_partition[0];
    }
    else
    {
        KRATOS_ERROR << "Failed to find a valid data block for reading. Number "
                        "of processors does not match the file partition."
                     << std::endl;
    }

    return std::tuple<unsigned, unsigned>(start_index, block_size);

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.