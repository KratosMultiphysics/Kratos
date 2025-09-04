//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "input_output/logger.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_io/hdf5_connectivities_data.h"
#include "custom_io/hdf5_data_value_container_io.h"
#include "custom_io/hdf5_nodal_solution_step_variables_io.h"
#include "custom_io/hdf5_points_data.h"
#include "custom_io/hdf5_properties_io.h"
#include "custom_utilities/container_io_utils.h"
#include "custom_utilities/factor_elements_and_conditions_utility.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"
#include "hdf5_application_variables.h"

// Include base h
#include "custom_io/hdf5_model_part_io.h"

namespace Kratos
{
namespace HDF5
{
ModelPartIO::ModelPartIO(
    Parameters Settings,
    File::Pointer pFile)
    : mpFile(pFile),
      mWriteEntityPropertyIds(true)
{
    Parameters default_params(R"(
        {
            "prefix"           : "",
            "custom_attributes": {}
        })");

    Settings.AddMissingParameters(default_params);
    mPrefix = Settings["prefix"].GetString();

    mCustomAttributes = Settings["custom_attributes"];

    KRATOS_ERROR_IF(mCustomAttributes.Has("__model_part_name"))
        << "The \"__model_part_name\" is a reserved attribute name. Please remove it from custom attributes."
        << " Custom attributes:\n" << mCustomAttributes << "\n";

    KRATOS_ERROR_IF(mPrefix.empty()) << "ModelPartIO requires non-empty prefix.\n";

    KRATOS_ERROR_IF(mPrefix.back() == '/')
        << "The ModelPartIO prefix is always assumed to be representing a group,"
        << " hence the trailing \"/\" is not required [ provided prefix = "
        << mPrefix << " ].\n";
}

ModelPartIO::ModelPartIO(
    const std::string& rPrefix,
    const bool WriteEntityPropertyIds,
    File::Pointer pFile)
    : ModelPartIO(
        [&rPrefix](){Parameters output; output.AddString("prefix", rPrefix); return output; }(),
        pFile)
{
    mWriteEntityPropertyIds = WriteEntityPropertyIds;
}

bool ModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    rNodes.clear();

    Internals::PointsData<Internals::NodesIO> points(mPrefix + "/Nodes/Local", mpFile);
    points.Read(rNodes, Internals::NodesIO{});

    return true;

    KRATOS_CATCH("");
}

std::size_t ModelPartIO::ReadNodesNumber()
{
    const std::vector<unsigned> dims = mpFile->GetDataDimensions(mPrefix + "/Nodes/Local/Ids");
    return dims[0];
}

void ModelPartIO::WriteNodes(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    Internals::PointsData<Internals::NodesIO> points(mPrefix + "/Nodes/Local", mpFile);
    points.Write(rNodes, Internals::NodesIO{}, Parameters(R"({})"));

    KRATOS_CATCH("");
}

void ModelPartIO::ReadProperties(PropertiesContainerType& rProperties)
{
    Internals::ReadProperties(*mpFile, mPrefix, rProperties);
}

void ModelPartIO::WriteProperties(PropertiesContainerType const& rProperties)
{
    Internals::WriteProperties(*mpFile, mPrefix, rProperties);
}

void ModelPartIO::ReadElements(
    NodesContainerType& rNodes,
    PropertiesContainerType& rProperties,
    ElementsContainerType& rElements)
{
    KRATOS_TRY;

    rElements.clear();

    const auto& element_path = mPrefix + "/Elements";

    for (const auto& r_name : mpFile->GetGroupNames(element_path)) {
        Internals::ConnectivitiesData<ElementsContainerType> connectivities(element_path, mpFile);
        connectivities.Read(r_name, rNodes, rProperties, rElements);
    }

    rElements.Unique();

    KRATOS_CATCH("");
}

void ModelPartIO::WriteElements(ElementsContainerType const& rElements)
{
    KRATOS_TRY;

    std::vector<std::string> names;
    std::vector<ElementsContainerType> factored_elements;
    FactorElements(rElements, names, factored_elements);

    const auto& element_path = mPrefix + "/Elements";

    mpFile->AddPath(element_path);
    for (unsigned int i = 0; i < names.size(); ++i) {
        Internals::ConnectivitiesData<ElementsContainerType> connectivities(element_path, mpFile);
        connectivities.Write(factored_elements[i], mWriteEntityPropertyIds);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::ReadConditions(
    NodesContainerType& rNodes,
    PropertiesContainerType& rProperties,
    ConditionsContainerType& rConditions)
{
    KRATOS_TRY;

    rConditions.clear();

    const auto& condition_path = mPrefix + "/Conditions";

    for (const auto& r_name : mpFile->GetGroupNames(condition_path)) {
        Internals::ConnectivitiesData<ConditionsContainerType> connectivities(condition_path, mpFile);
        connectivities.Read(r_name, rNodes, rProperties, rConditions);
    }

    rConditions.Unique();

    KRATOS_CATCH("");
}

void ModelPartIO::WriteConditions(ConditionsContainerType const& rConditions)
{
    KRATOS_TRY;

    std::vector<std::string> names;
    std::vector<ConditionsContainerType> factored_conditions;
    FactorConditions(rConditions, names, factored_conditions);
    const auto& condition_path = mPrefix + "/Conditions";

    mpFile->AddPath(condition_path);
    for (unsigned int i = 0; i < names.size(); ++i) {
        Internals::ConnectivitiesData<ConditionsContainerType> connectivities(condition_path, mpFile);
        connectivities.Write(factored_conditions[i], mWriteEntityPropertyIds);
    }

    KRATOS_CATCH("");
}

void ModelPartIO::WriteModelPart(ModelPart& rModelPart)
{
    KRATOS_TRY;

    BuiltinTimer timer;

    // write model part name as an attribute
    mpFile->AddPath(mPrefix);

    if (!mCustomAttributes.Has("__model_part_name")) {
        mCustomAttributes.AddString("__model_part_name", rModelPart.FullName());
    } else {
        KRATOS_ERROR_IF(!mCustomAttributes["__model_part_name"].IsString() || mCustomAttributes["__model_part_name"].GetString() != rModelPart.FullName())
            << "The provided custom attributes has the reserved \"__model_part_name\" attribute with a different value other than \""
            << rModelPart.FullName() << "\" [ attribute value = " << mCustomAttributes["__model_part_name"] << " ].\n";
    }
    mpFile->WriteAttribute(mPrefix, mCustomAttributes);

    Internals::WriteVariablesList(*mpFile, mPrefix, rModelPart);
    Internals::WriteBufferSize(*mpFile, mPrefix, rModelPart.GetBufferSize());
    WriteProperties(rModelPart.rProperties());
    Internals::WriteDataValueContainer(*mpFile, mPrefix + "/ProcessInfo", rModelPart.GetProcessInfo());

    // Avoid inadvertently reordering partway through
    // the writing process.
    // TODO: Remove after PVS refactoring
    rModelPart.Nodes().Unique();
    rModelPart.Conditions().Unique();
    rModelPart.Elements().Unique();

    WriteNodes(rModelPart.Nodes());
    WriteElements(rModelPart.Elements());
    WriteConditions(rModelPart.Conditions());

    // writing submodel parts
    mpFile->AddPath(mPrefix + "/SubModelParts");
    WriteSubModelParts(rModelPart.SubModelParts(), mPrefix + "/SubModelParts");

    KRATOS_INFO_IF("HDF5Application", mpFile->GetEchoLevel() == 1)
        << "Time to write model part \"" << rModelPart.Name()
        << "\": " << timer.ElapsedSeconds() << " seconds." << std::endl;

    KRATOS_CATCH("");
}

void ModelPartIO::ReadModelPart(ModelPart& rModelPart)
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

    for (const auto& r_sub_model_part_name : mpFile->GetGroupNames(mPrefix + "/SubModelParts")) {
        ReadSubModelParts(rModelPart, mPrefix + "/SubModelParts/" + r_sub_model_part_name);
    }

    this->ReadParitionIndices(rModelPart);

    this->SetCommunicator(rModelPart);

    KRATOS_INFO_IF("HDF5Application", mpFile->GetEchoLevel() == 1)
        << "Time to read model part \"" << rModelPart.Name()
        << "\": " << timer.ElapsedSeconds() << " seconds." << std::endl;

    KRATOS_CATCH("");
}

void ModelPartIO::ReadParitionIndices(ModelPart& rModelPart)
{
    // Do nothing in serial IO.
}

void ModelPartIO::SetCommunicator(ModelPart& rModelPart) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mpFile->GetDataCommunicator().IsDistributed()) << "Using ModelPartIO with an HDF5File which is MPI based.";

    // set the serial communicator
    if (rModelPart.GetCommunicator().IsDistributed()) {
        rModelPart.SetCommunicator(Communicator().Create(mpFile->GetDataCommunicator()));
    }

    // now set local meshes and ghost meshes
    auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    r_local_mesh.Nodes() = rModelPart.Nodes();
    r_local_mesh.Conditions() = rModelPart.Conditions();
    r_local_mesh.Elements() = rModelPart.Elements();

    // set the sub model part communicators recursively.
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        this->SetCommunicator(r_sub_model_part);
    }

    KRATOS_CATCH("");
}

std::vector<std::size_t> ModelPartIO::ReadContainerIds(std::string const& rPath) const
{
    KRATOS_TRY

    const auto [start_index, block_size] = HDF5::StartIndexAndBlockSize(*mpFile, rPath);

    Vector<int> id_buf;
    mpFile->ReadDataSet(rPath + "/Ids", id_buf, start_index, block_size);
    std::vector<std::size_t> ids(id_buf.size());

    IndexPartition<IndexType>(id_buf.size()).for_each([&id_buf, &ids](const auto Index) {
        ids[Index] = id_buf[Index];
    });

    return ids;

    KRATOS_CATCH("");
}

void ModelPartIO::WriteNodeIds(
    const std::string& rPath,
    const NodesContainerType & rNodes)
{
    KRATOS_TRY;

    Vector<int> ids;
    ids.resize(rNodes.size());

    IndexPartition<IndexType>(rNodes.size()).for_each([&rNodes, &ids](const auto Index) {
        ids[Index] = (rNodes.begin() + Index)->Id();
    });

    WriteInfo info;
    mpFile->AddPath(rPath);
    mpFile->WriteDataSet(rPath + "/Ids", ids, info);
    WritePartitionTable(*mpFile, rPath, info);

    KRATOS_CATCH("");
}

void ModelPartIO::WriteSubModelParts(
    const ModelPart::SubModelPartsContainerType& rSubModelPartsContainer,
    const std::string& GroupName)
{
    for (const auto& r_sub_model_part : rSubModelPartsContainer) {
        const std::string sub_model_part_path = GroupName + "/" + r_sub_model_part.Name();

        const auto& r_communicator = r_sub_model_part.GetCommunicator();

        mpFile->AddPath(sub_model_part_path);

        if (r_communicator.GlobalNumberOfNodes() > 0) {
            WriteNodeIds(sub_model_part_path + "/Nodes", r_sub_model_part.Nodes());
        }

        if (r_communicator.GlobalNumberOfElements() > 0) {
            ModelPartIO current_model_part_io(sub_model_part_path, false,  mpFile);
            current_model_part_io.WriteElements(r_sub_model_part.Elements());
        }

        if (r_communicator.GlobalNumberOfConditions() > 0) {
            ModelPartIO current_model_part_io(sub_model_part_path, false, mpFile);
            current_model_part_io.WriteConditions(r_sub_model_part.Conditions());
        }

        WriteSubModelParts(r_sub_model_part.SubModelParts(), sub_model_part_path);
    }
}

void ModelPartIO::ReadSubModelParts(ModelPart& rModelPart, const std::string& rPath)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rPath.rfind("/") == std::string::npos)
        << "The rPath = \"" << rPath << "\" is not a valid hdf5 path. HDF5 path should always start with \"/\"\n.";

    auto& r_sub_model_part = rModelPart.CreateSubModelPart(rPath.substr(rPath.rfind("/") + 1));

    const auto &group_names = mpFile->GetGroupNames(rPath);
    for (const auto &group_name : group_names) {
        const auto &current_path = rPath + "/" + group_name;
        if (group_name == "Nodes") {
            r_sub_model_part.AddNodes(ReadContainerIds(current_path));
        } else if (group_name == "Elements") {
            // iterate over all types of elements
            for (const auto &element_name : mpFile->GetGroupNames(current_path)) {
                r_sub_model_part.AddElements(ReadContainerIds(current_path + "/" + element_name));
            }
        } else if (group_name == "Conditions") {
            // iterate over all types of conditions
            for (const auto &condition_name : mpFile->GetGroupNames(current_path)) {
                r_sub_model_part.AddConditions(ReadContainerIds(current_path + "/" + condition_name));
            }
        } else {
            ReadSubModelParts(r_sub_model_part, current_path);
        }
    }

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.