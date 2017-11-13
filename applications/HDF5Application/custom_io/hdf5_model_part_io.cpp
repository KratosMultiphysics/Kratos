#include "custom_io/hdf5_model_part_io.h"

#include "includes/kratos_components.h"
#include "custom_utilities/hdf5_points_data.h"
#include "custom_utilities/hdf5_connectivities_data.h"
#include "custom_utilities/hdf5_pointer_bins_utility.h"
#include "custom_io/hdf5_properties_io.h"
#include "custom_io/hdf5_nodal_solution_step_variables_io.h"
#include "custom_io/hdf5_data_value_container_io.h"

namespace Kratos
{
namespace HDF5
{
ModelPartIO::ModelPartIO(Parameters& rParams, File::Pointer pFile)
: mpFile(pFile)
{
    KRATOS_TRY;

    Parameters default_params(R"(
        {
            "prefix": "",
            "list_of_elements": [],
            "list_of_conditions": []
        })");

    rParams.ValidateAndAssignDefaults(default_params);

    mPrefix = rParams["prefix"].GetString();
    if (mPrefix == "/")
        mPrefix = "";

    mElementNames.resize(rParams["list_of_elements"].size());
    for (unsigned i = 0; i < mElementNames.size(); ++i)
        mElementNames[i] = rParams["list_of_elements"].GetArrayItem(i).GetString();

    mConditionNames.resize(rParams["list_of_conditions"].size());
    for (unsigned i = 0; i < mConditionNames.size(); ++i)
        mConditionNames[i] = rParams["list_of_conditions"].GetArrayItem(i).GetString();

    Check();

    mElementPointers.resize(mElementNames.size());
    for (unsigned i = 0; i < mElementNames.size(); ++i)
    {
        const ElementType& r_elem = KratosComponents<Element>::Get(mElementNames[i]);
        mElementPointers[i] = &r_elem;
    }

    mConditionPointers.resize(mConditionNames.size());
    for (unsigned i = 0; i < mConditionNames.size(); ++i)
    {
        const ConditionType& r_cond = KratosComponents<Condition>::Get(mConditionNames[i]);
        mConditionPointers[i] = &r_cond;
    }

    KRATOS_CATCH("");
}

bool ModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    rNodes.clear();

    const std::size_t num_nodes = ReadNodesNumber();
    File& r_file = GetFile();
    Detail::PointsData points;

    points.ReadData(r_file, mPrefix + "/Nodes/Local", 0, num_nodes);
    points.CreateNodes(rNodes);

    return true;
    KRATOS_CATCH("");
}

std::size_t ModelPartIO::ReadNodesNumber()
{
    const std::vector<unsigned> dims = GetFile().GetDataDimensions(mPrefix + "/Nodes/Local/Ids");
    return dims[0];
}

void ModelPartIO::WriteNodes(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    Detail::PointsData points;
    points.SetData(rNodes);
    File& r_file = GetFile();
    points.WriteData(r_file, mPrefix + "/Nodes/Local");
    
    KRATOS_CATCH("");
}

void ModelPartIO::ReadProperties(PropertiesContainerType& rProperties)
{
    Detail::PropertiesIO prop_io(mPrefix, mpFile);
    prop_io.ReadProperties(rProperties);
}

void ModelPartIO::WriteProperties(Properties const& rProperties)
{
    Detail::PropertiesIO prop_io(mPrefix, mpFile);
    prop_io.WriteProperties(rProperties);
}

void ModelPartIO::WriteProperties(PropertiesContainerType const& rProperties)
{
    Detail::PropertiesIO prop_io(mPrefix, mpFile);
    prop_io.WriteProperties(rProperties);
}

void ModelPartIO::ReadElements(NodesContainerType& rNodes,
                               PropertiesContainerType& rProperties,
                               ElementsContainerType& rElements)
{
    KRATOS_TRY;

    rElements.clear();

    File& r_file = GetFile();

    for (unsigned i = 0; i < mElementNames.size(); ++i)
    {
        const std::string elem_path = mPrefix + "/Elements/" + mElementNames[i];
        const unsigned num_elems = r_file.GetDataDimensions(elem_path + "/Ids")[0];
        Detail::ConnectivitiesData connectivities;
        connectivities.ReadData(r_file, elem_path, 0, num_elems);
        const ElementType& r_elem = *mElementPointers[i];
        connectivities.CreateElements(r_elem, rNodes, rProperties, rElements);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    const unsigned num_elem_types = mElementNames.size();
    Detail::PointerBinsUtility<ElementType> elem_bins(mElementPointers);
    elem_bins.CreateBins(rElements);
    File& r_file = GetFile();
    for (unsigned i_type = 0; i_type < num_elem_types; ++i_type)
    {
        std::string& r_elem_name = mElementNames[i_type];
        const ElementType* elem_key = mElementPointers[i_type];
        ConstElementsContainerType& r_elems = elem_bins.GetBin(elem_key);
        Detail::ConnectivitiesData connectivities;
        connectivities.SetData(r_elems);
        connectivities.WriteData(r_file, mPrefix + "/Elements/" + r_elem_name);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::ReadConditions(NodesContainerType& rNodes,
                                 PropertiesContainerType& rProperties,
                                 ConditionsContainerType& rConditions)
{
    KRATOS_TRY;

    rConditions.clear();

    File& r_file = GetFile();

    for (unsigned i = 0; i < mConditionNames.size(); ++i)
    {
        const std::string cond_path = mPrefix + "/Conditions/" + mConditionNames[i];
        const unsigned num_conds = r_file.GetDataDimensions(cond_path + "/Ids")[0];
        Detail::ConnectivitiesData connectivities;
        connectivities.ReadData(r_file, cond_path, 0, num_conds);
        const ConditionType& r_cond = *mConditionPointers[i];
        connectivities.CreateConditions(r_cond, rNodes, rProperties, rConditions);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    const unsigned num_cond_types = mConditionNames.size();
    Detail::PointerBinsUtility<ConditionType> cond_bins(mConditionPointers);
    cond_bins.CreateBins(rConditions);
    File& r_file = GetFile();
    for (unsigned i_type = 0; i_type < num_cond_types; ++i_type)
    {
        std::string& r_cond_name = mConditionNames[i_type];
        const ConditionType* cond_key = mConditionPointers[i_type];
        ConstConditionsContainerType& r_conds = cond_bins.GetBin(cond_key);
        Detail::ConnectivitiesData connectivities;
        connectivities.SetData(r_conds);
        connectivities.WriteData(r_file, mPrefix + "/Conditions/" + r_cond_name);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::ReadModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    ReadProperties(rModelPart.rProperties());
    Detail::DataValueContainerIO process_info_io(mPrefix + "/ProcessInfo", mpFile);
    process_info_io.ReadDataValueContainer(rModelPart.GetProcessInfo());
    ReadNodes(rModelPart.Nodes());
    ReadElements(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Elements());
    ReadConditions(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Conditions());
    Detail::NodalSolutionStepVariablesIO nodal_variables_io(mPrefix, mpFile);
    nodal_variables_io.ReadAndAssignVariablesList(rModelPart);
    nodal_variables_io.ReadAndAssignBufferSize(rModelPart);

    KRATOS_CATCH("");
}

void ModelPartIO::WriteModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    Detail::NodalSolutionStepVariablesIO nodal_variables_io(mPrefix, mpFile);
    nodal_variables_io.WriteVariablesList(rModelPart);
    nodal_variables_io.WriteBufferSize(rModelPart.GetBufferSize());
    WriteProperties(rModelPart.rProperties());
    Detail::DataValueContainerIO process_info_io(mPrefix + "/ProcessInfo", mpFile);
    process_info_io.WriteDataValueContainer(rModelPart.GetProcessInfo());
    rModelPart.Nodes().Sort(); // Avoid inadvertently reordering partway through
                               // the writing process.
    WriteNodes(rModelPart.Nodes());
    WriteElements(rModelPart.Elements());
    WriteConditions(rModelPart.Conditions());

    KRATOS_CATCH("");
}

File& ModelPartIO::GetFile() const
{
    return *mpFile;
}

void ModelPartIO::Check()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(GetFile().GetTotalProcesses() != 1)
        << "Serial IO expects file access by a single process only." << std::endl;

    for (const std::string& r_elem_name : mElementNames)
    {
        KRATOS_ERROR_IF_NOT(KratosComponents<Element>::Has(r_elem_name))
            << "Element " << r_elem_name << " not registered." << std::endl;
    }

    for (const std::string& r_cond_name : mConditionNames)
    {
        KRATOS_ERROR_IF_NOT(KratosComponents<Condition>::Has(r_cond_name))
            << "Condition " << r_cond_name << " not registered." << std::endl;
    }

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
