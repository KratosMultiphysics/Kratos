#include "custom_io/hdf5_model_part_io.h"

#include "custom_utilities/hdf5_points_data.h"
#include "custom_io/hdf5_nodal_solution_step_variables_io.h"
#include "custom_io/hdf5_data_value_container_io.h"

namespace Kratos
{
namespace HDF5
{
ModelPartIO::ModelPartIO(Parameters Settings, File::Pointer pFile)
: ModelPartIOBase(Settings, pFile)
{
    KRATOS_TRY;

    Check();

    KRATOS_CATCH("");
}

bool ModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    rNodes.clear();

    const std::size_t num_nodes = ReadNodesNumber();
    File& r_file = GetFile();
    Internals::PointsData points;

    points.ReadData(r_file, mPrefix + "/Nodes/Local", 0, num_nodes);
    points.CreateNodes(rNodes);

    return true;
    KRATOS_CATCH("");
}

void ModelPartIO::WriteNodes(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    Internals::PointsData points;
    points.SetData(rNodes);
    File& r_file = GetFile();
    points.WriteData(r_file, mPrefix + "/Nodes/Local");
    
    KRATOS_CATCH("");
}

void ModelPartIO::ReadElements(NodesContainerType& rNodes,
                               PropertiesContainerType& rProperties,
                               ElementsContainerType& rElements)
{
    KRATOS_TRY;

    rElements.clear();

    File& r_file = GetFile();

    Internals::ConnectivitiesInput<ElementType> elem_inputs(mElementIO);
    for (auto& r_item : elem_inputs)
    {
        const unsigned num_elems = r_file.GetDataDimensions(r_item.Path + "/Ids")[0];
        r_item.ReadConnectivities(r_file, rNodes, rProperties, 0, num_elems, rElements);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    Internals::ConnectivitiesOutput<ElementType> elem_outputs(mElementIO, rElements);
    File& r_file = GetFile();
    for (auto& r_item : elem_outputs)
        r_item.WriteConnectivities(r_file);

    KRATOS_CATCH("");
}

void ModelPartIO::ReadConditions(NodesContainerType& rNodes,
                                 PropertiesContainerType& rProperties,
                                 ConditionsContainerType& rConditions)
{
    KRATOS_TRY;

    rConditions.clear();
    File& r_file = GetFile();

    Internals::ConnectivitiesInput<ConditionType> cond_inputs(mConditionIO);
    for (auto& r_item : cond_inputs)
    {
        const unsigned num_conds = r_file.GetDataDimensions(r_item.Path + "/Ids")[0];
        r_item.ReadConnectivities(r_file, rNodes, rProperties, 0, num_conds, rConditions);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    Internals::ConnectivitiesOutput<ConditionType> cond_outputs(mConditionIO, rConditions);
    File& r_file = GetFile();
    for (auto& r_item : cond_outputs)
        r_item.WriteConnectivities(r_file);

    KRATOS_CATCH("");
}

void ModelPartIO::ReadModelPart(ModelPart& rModelPart)
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

void ModelPartIO::Check()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(GetFile().GetTotalProcesses() != 1)
        << "Serial IO expects file access by a single process only." << std::endl;

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
