#include "custom_io/hdf5_model_part_io_base.h"

#include "custom_io/hdf5_properties_io.h"
#include "custom_io/hdf5_nodal_solution_step_variables_io.h"
#include "custom_io/hdf5_data_value_container_io.h"

namespace Kratos
{
namespace HDF5
{
ModelPartIOBase::ModelPartIOBase(Parameters Settings, File::Pointer pFile)
: mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "prefix": ""
        })");

    Settings.ValidateAndAssignDefaults(default_params);

    mPrefix = Settings["prefix"].GetString();
    if (mPrefix == "/")
        mPrefix = "";

    KRATOS_CATCH("");
}

std::size_t ModelPartIOBase::ReadNodesNumber()
{
    const std::vector<unsigned> dims = mpFile->GetDataDimensions(mPrefix + "/Nodes/Local/Ids");
    return dims[0];
}

void ModelPartIOBase::ReadProperties(PropertiesContainerType& rProperties)
{
    Internals::PropertiesIO prop_io(mPrefix, mpFile);
    prop_io.ReadProperties(rProperties);
}

void ModelPartIOBase::WriteProperties(Properties const& rProperties)
{
    Internals::PropertiesIO prop_io(mPrefix, mpFile);
    prop_io.WriteProperties(rProperties);
}

void ModelPartIOBase::WriteProperties(PropertiesContainerType const& rProperties)
{
    Internals::PropertiesIO prop_io(mPrefix, mpFile);
    prop_io.WriteProperties(rProperties);
}

void ModelPartIOBase::WriteModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    Internals::NodalSolutionStepVariablesIO nodal_variables_io(mPrefix, mpFile);
    nodal_variables_io.WriteVariablesList(rModelPart);
    nodal_variables_io.WriteBufferSize(rModelPart.GetBufferSize());
    WriteProperties(rModelPart.rProperties());
    Internals::DataValueContainerIO process_info_io(mPrefix + "/ProcessInfo", mpFile);
    process_info_io.WriteDataValueContainer(rModelPart.GetProcessInfo());
    rModelPart.Nodes().Sort(); // Avoid inadvertently reordering partway through
                               // the writing process.
    WriteNodes(rModelPart.Nodes());
    WriteElements(rModelPart.Elements());
    WriteConditions(rModelPart.Conditions());

    KRATOS_CATCH("");
}

void ModelPartIOBase::ReadModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    ReadProperties(rModelPart.rProperties());
    Internals::DataValueContainerIO process_info_io(mPrefix + "/ProcessInfo", mpFile);
    process_info_io.ReadDataValueContainer(rModelPart.GetProcessInfo());
    ReadNodes(rModelPart.Nodes());
    ReadElements(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Elements());
    ReadConditions(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Conditions());
    Internals::NodalSolutionStepVariablesIO nodal_variables_io(mPrefix, mpFile);
    nodal_variables_io.ReadAndAssignVariablesList(rModelPart);
    nodal_variables_io.ReadAndAssignBufferSize(rModelPart);

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
