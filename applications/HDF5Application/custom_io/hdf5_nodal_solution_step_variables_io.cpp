#include "custom_io/hdf5_nodal_solution_step_variables_io.h"

#include <vector>
#include "includes/kratos_components.h"
#include "utilities/openmp_utils.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{

void WriteVariablesList(File& rFile, std::string const& rPrefix, ModelPart const& rModelPart)
{
    KRATOS_TRY;

    const VariablesList& r_variables_list = rModelPart.GetNodalSolutionStepVariablesList();
    int pos = 0;
    rFile.AddPath(rPrefix + "/NodalSolutionStep/VariablesList");
    for (auto it = r_variables_list.begin(); it != r_variables_list.end(); ++it)
        rFile.WriteAttribute(rPrefix + "/NodalSolutionStep/VariablesList", it->Name(), pos++);

    KRATOS_CATCH("");
}

void ReadAndAssignVariablesList(File& rFile, std::string const& rPrefix, ModelPart& rModelPart)
{
    KRATOS_TRY;

    VariablesList& r_variables_list = rModelPart.GetNodalSolutionStepVariablesList();
    r_variables_list.clear();
    std::vector<std::string> variable_names = rFile.GetAttributeNames(rPrefix + "/NodalSolutionStep/VariablesList");

    // Ensure the variables order is the same as in the original model part.
    std::vector<std::string> ordered_variable_names(variable_names.size());
    for (const auto& r_name : variable_names)
    {
        int pos;
        rFile.ReadAttribute(rPrefix + "/NodalSolutionStep/VariablesList", r_name, pos);
        ordered_variable_names[pos] = r_name;
    }

    // Add the variables to the variables list.
    for (const auto& r_name : ordered_variable_names)
    {
        if (KratosComponents<VariableData>::Has(r_name))
        {
            const VariableData& rVARIABLE = KratosComponents<VariableData>::Get(r_name);
            r_variables_list.Add(rVARIABLE);
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable type: " << r_name << std::endl;
        }
    }

    // Assign variables list to nodes.
#pragma omp parallel
    {
        NodesContainerType::iterator nodes_begin, nodes_end;
        OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), nodes_begin, nodes_end);
        for (auto it = nodes_begin; it != nodes_end; ++it)
            it->SetSolutionStepVariablesList(&r_variables_list);
    }

    KRATOS_CATCH("");
}

void WriteBufferSize(File& rFile, std::string const& rPrefix, int BufferSize)
{
    KRATOS_TRY;

    rFile.AddPath(rPrefix + "/NodalSolutionStep");
    rFile.WriteAttribute(rPrefix + "/NodalSolutionStep", "BufferSize", BufferSize);

    KRATOS_CATCH("");
}

void ReadAndAssignBufferSize(File& rFile, std::string const& rPrefix, ModelPart& rModelPart)
{
    KRATOS_TRY;

    int buffer_size;
    rFile.ReadAttribute(rPrefix + "/NodalSolutionStep", "BufferSize", buffer_size);
    rModelPart.SetBufferSize(buffer_size);

    KRATOS_CATCH("");
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.