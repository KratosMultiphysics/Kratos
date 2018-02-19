#include "custom_io/hdf5_model_part_io_base.h"

#include "custom_io/hdf5_properties_io.h"
#include "custom_io/hdf5_nodal_solution_step_variables_io.h"
#include "custom_io/hdf5_data_value_container_io.h"
#include "utilities/builtin_timer.h"

namespace Kratos
{
namespace HDF5
{
ModelPartIOBase::ModelPartIOBase(File::Pointer pFile, std::string const& rPrefix)
: mpFile(pFile), mPrefix(rPrefix)
{
}

std::size_t ModelPartIOBase::ReadNodesNumber()
{
    const std::vector<unsigned> dims = mpFile->GetDataDimensions(mPrefix + "/Nodes/Local/Ids");
    return dims[0];
}

void ModelPartIOBase::ReadProperties(PropertiesContainerType& rProperties)
{
    Internals::ReadProperties(*mpFile, mPrefix, rProperties);
}

void ModelPartIOBase::WriteProperties(Properties const& rProperties)
{
    Internals::WriteProperties(*mpFile, mPrefix, rProperties);
}

void ModelPartIOBase::WriteProperties(PropertiesContainerType const& rProperties)
{
    Internals::WriteProperties(*mpFile, mPrefix, rProperties);
}

void ModelPartIOBase::WriteModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    BuiltinTimer timer;
    Internals::WriteVariablesList(*mpFile, mPrefix, rModelPart);
    Internals::WriteBufferSize(*mpFile, mPrefix, rModelPart.GetBufferSize());
    WriteProperties(rModelPart.rProperties());
    Internals::WriteDataValueContainer(*mpFile, mPrefix + "/ProcessInfo", rModelPart.GetProcessInfo());
    rModelPart.Nodes().Sort(); // Avoid inadvertently reordering partway through
                               // the writing process.
    WriteNodes(rModelPart.Nodes());
    WriteElements(rModelPart.Elements());
    WriteConditions(rModelPart.Conditions());

     if (mpFile->GetEchoLevel() == 1 && mpFile->GetPID() == 0)
        std::cout << "Time to write model part \"" << rModelPart.Name() << "\": " << timer.ElapsedSeconds() << " seconds." << std::endl;

    KRATOS_CATCH("");
}

void ModelPartIOBase::ReadModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    BuiltinTimer timer;
    ReadProperties(rModelPart.rProperties());
    Internals::ReadDataValueContainer(*mpFile, mPrefix + "/ProcessInfo", rModelPart.GetProcessInfo());
    ReadNodes(rModelPart.Nodes());
    ReadElements(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Elements());
    ReadConditions(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Conditions());
    Internals::ReadAndAssignVariablesList(*mpFile, mPrefix, rModelPart);
    Internals::ReadAndAssignBufferSize(*mpFile, mPrefix, rModelPart);

    if (mpFile->GetEchoLevel() == 1 && mpFile->GetPID() == 0)
        std::cout << "Time to read model part \"" << rModelPart.Name() << "\": " << timer.ElapsedSeconds() << " seconds." << std::endl;

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
