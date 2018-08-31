#include "custom_io/hdf5_nodal_solution_step_bossak_io.h"

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

NodalSolutionStepBossakIO::NodalSolutionStepBossakIO(Parameters Settings, File::Pointer pFile)
    : NodalSolutionStepDataIO(Settings, pFile)
{
}

namespace {
template <typename TVariable>
class WriteVariableFunctor;
}

void NodalSolutionStepBossakIO::WriteNodalResults(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    if (VariableNames().size() == 0)
        return;

    std::vector<NodeType*> local_nodes;
    GetLocalNodes(rNodes, local_nodes);

    // Write each variable.
    const std::string& prefix = GetPrefix();
    WriteInfo info;
    for (const std::string& r_name : VariableNames())
        RegisteredVariableLookup<Variable<array_1d<double, 3>>,
                                 VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>,
                                 Variable<double>, Variable<int>>(r_name)
            .Execute<WriteVariableFunctor>(local_nodes, GetFile(), prefix,
                                           mAlphaBossak, info);

    // Write block partition.
    WritePartitionTable(GetFile(), prefix + "/NodalSolutionStepData", info);

    KRATOS_CATCH("");
}

namespace
{
template <class TVariableType, class TFileDataType>
void SetDataBuffer(TVariableType const&, std::vector<NodeType*> const&, double, Vector<TFileDataType>&);

template <typename TVariable>
class WriteVariableFunctor
{
public:
    void operator()(TVariable const& rVariable,
                    std::vector<NodeType*>& rNodes,
                    File& rFile,
                    std::string const& rPrefix,
                    double AlphaBossak,
                    WriteInfo& rInfo)
    {
        Vector<typename TVariable::Type> data;
        SetDataBuffer(rVariable, rNodes, AlphaBossak, data);
        rFile.WriteDataSet(rPrefix + "/NodalSolutionStepData/" + rVariable.Name(), data, rInfo);
    }
};

template <class TVariableType, class TFileDataType>
void SetDataBuffer(TVariableType const& rVariable,
                   std::vector<NodeType*> const& rNodes,
                   double AlphaBossak,
                   Vector<TFileDataType>& rData)
{
    KRATOS_TRY;

    rData.resize(rNodes.size(), false);
    if (rVariable == ACCELERATION)
    {
#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rNodes.size()); ++i)
            rData[i] =
                (1.0 - AlphaBossak) * rNodes[i]->FastGetSolutionStepValue(rVariable, 0) +
                AlphaBossak * rNodes[i]->FastGetSolutionStepValue(rVariable, 1);
    }
    else
    {
#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rNodes.size()); ++i)
            rData[i] = rNodes[i]->FastGetSolutionStepValue(rVariable);
    }
    KRATOS_CATCH("");
}
}

void NodalSolutionStepBossakIO::ReadNodalResults(NodesContainerType& rNodes, Communicator& rComm)
{
    KRATOS_TRY;
    NodalSolutionStepDataIO::ReadNodalResults(rNodes, rComm);
    KRATOS_CATCH("");
}

void NodalSolutionStepBossakIO::SetAlphaBossak(double alpha) noexcept
{
    mAlphaBossak = alpha;
}

} // namespace HDF5.
} // namespace Kratos.