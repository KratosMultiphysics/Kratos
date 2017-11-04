#include "custom_io/hdf5_nodal_solution_step_variables_io.h"

#include <vector>
#include "includes/kratos_components.h"

namespace Kratos
{
namespace HDF5
{
namespace Detail
{
NodalSolutionStepVariablesIO::NodalSolutionStepVariablesIO(std::string Prefix, File::Pointer pFile)
: mPrefix(Prefix), mpFile(pFile)
{
}

void NodalSolutionStepVariablesIO::WriteVariablesList(VariablesList const& rVariablesList)
{
    KRATOS_TRY;

    int pos = 0;
    mpFile->AddPath(mPrefix + "/NodalVariablesList");
    for (auto it = rVariablesList.begin(); it != rVariablesList.end(); ++it)
        mpFile->WriteAttribute(mPrefix + "/NodalVariablesList", it->Name(), pos++);

    KRATOS_CATCH("");
}

void NodalSolutionStepVariablesIO::ReadVariablesList(VariablesList& rVariablesList) const
{
    KRATOS_TRY;

    rVariablesList.clear();
    std::vector<std::string> variable_names;
    mpFile->GetAttributeNames(mPrefix + "/NodalVariablesList", variable_names);
    std::vector<std::string> ordered_variable_names(variable_names.size());
    for (const auto& r_name : variable_names)
    {
        int pos;
        mpFile->ReadAttribute(mPrefix + "/NodalVariablesList", r_name, pos);
        ordered_variable_names[pos] = r_name;
    }

    for (const auto& r_name : ordered_variable_names)
    {
        if (KratosComponents<VariableData>::Has(r_name))
        {
            const VariableData& rVARIABLE = KratosComponents<VariableData>::Get(r_name);
            rVariablesList.Add(rVARIABLE);
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable type: " << r_name << std::endl;
        }
    }

    KRATOS_CATCH("");
}

void NodalSolutionStepVariablesIO::WriteBufferSize(int BufferSize)
{
    KRATOS_TRY;

    mpFile->AddPath(mPrefix);
    mpFile->WriteAttribute(mPrefix, "BufferSize", BufferSize);

    KRATOS_CATCH("");
}

int NodalSolutionStepVariablesIO::ReadBufferSize() const
{
    KRATOS_TRY;

    int buffer_size;
    mpFile->ReadAttribute(mPrefix, "BufferSize", buffer_size);
    return buffer_size;

    KRATOS_CATCH("");
}

} // namespace Detail.
} // namespace HDF5.
} // namespace Kratos.