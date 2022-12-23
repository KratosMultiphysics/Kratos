// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <optional>

// External includes
#include "med.h"

// Project includes
#include "med_model_part_io.h"
#include "includes/model_part_io.h"
#include "utilities/variable_utils.h"

namespace Kratos {

class MedModelPartIO::MedFileHandler
{
public:
    MedFileHandler(
        const std::filesystem::path& rFileName,
        const Kratos::Flags Options) :
            mFileName(rFileName)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(Options.Is(IO::APPEND)) << "Appending to med files is not supported!" << std::endl;
        KRATOS_ERROR_IF(Options.Is(IO::READ) && Options.Is(IO::WRITE)) << "Either reading OR writing is possible, not both!" << std::endl;

        const bool mIsReadMode = !Options.Is(IO::WRITE);

        // Set the mode (consistent with ModelPartIO)
        // read only by default, unless other settings are specified
        med_access_mode open_mode;

        if (mIsReadMode) {
            open_mode = MED_ACC_RDONLY;

            // check if file exists
            KRATOS_ERROR_IF(!std::filesystem::exists(rFileName)) << "File " << rFileName << " does not exist!" << std::endl;

            // basic checks if the file is compatible with the MED library
            med_bool hdf_ok;
            med_bool med_ok;
            const med_err err = MEDfileCompatibility(rFileName.c_str(), &hdf_ok, &med_ok);

            KRATOS_ERROR_IF(err != 0) << "A problem occured while trying to check the compatibility of file " << rFileName << "!" << std::endl;
            KRATOS_ERROR_IF(hdf_ok != MED_TRUE) << "A problem with HDF occured while trying to open file " << rFileName << "!" << std::endl;
            KRATOS_ERROR_IF(med_ok != MED_TRUE) << "A problem with MED occured while trying to open file " << rFileName << "!" << std::endl;
        } else {
            open_mode = MED_ACC_CREAT;
            mMeshName = "Kratos_Mesh";
        }

        mFileHandle = MEDfileOpen(rFileName.c_str(), open_mode);
        KRATOS_ERROR_IF(mFileHandle < 0) << "A problem occured while opening file " << rFileName << "!" << std::endl;

        if (mIsReadMode) {
            // when reading the mesh, it is necessary to querry more information upfront

            const int num_meshes = MEDnMesh(mFileHandle);
            KRATOS_ERROR_IF(num_meshes != 1) << "Expected one mesh, but file " << mFileName << " contains " << num_meshes << " meshes!" << std::endl;

            mMeshName.resize(MED_NAME_SIZE+1);
            med_int space_dim;
            med_int mesh_dim;
            med_mesh_type mesh_type;
            std::string description(MED_COMMENT_SIZE+1, '\0');
            std::string dt_unit(MED_SNAME_SIZE+1, '\0');
            med_sorting_type sorting_type;
            med_int n_step;
            med_axis_type axis_type;
            std::string axis_name(MED_SNAME_SIZE+100, '\0'); // Not sure why the exessive size is required, but segfaults otherwise
            std::string axis_unit(MED_SNAME_SIZE+100, '\0'); // Not sure why the exessive size is required, but segfaults otherwise

            med_err err = MEDmeshInfo(
                mFileHandle,
                1,
                mMeshName.data(),
                &space_dim,
                &mesh_dim,
                &mesh_type,
                description.data(),
                dt_unit.data(),
                &sorting_type,
                &n_step,
                &axis_type,
                axis_name.data(),
                axis_unit.data());

            mMeshName.erase(std::find(mMeshName.begin(), mMeshName.end(), '\0'), mMeshName.end());
            mDimension = space_dim;

            KRATOS_WATCH(space_dim);
            KRATOS_WATCH(mesh_dim);
        }

        KRATOS_CATCH("")
    }

    med_idt GetFileHandle() const
    {
        return mFileHandle;
    }

    const char* GetMeshName() const
    {
        return mMeshName.c_str();
    }

    bool IsReadMode() const
    {
        return mIsReadMode;
    }

    int GetDimension() const
    {
        KRATOS_ERROR_IF_NOT(mDimension.has_value()) << "Dimension can only be querried in read mode!";
        return mDimension.value();
    }

    ~MedFileHandler()
    {
        KRATOS_TRY

        KRATOS_WARNING_IF("MedModelPartIO", MEDfileClose(mFileHandle) < 0) << "Closing of file " << mFileName << " failed!" << std::endl;

        KRATOS_CATCH("")
    }

private:
    std::filesystem::path mFileName;
    med_idt mFileHandle;
    std::string mMeshName;
    bool mIsReadMode;
    std::optional<int> mDimension;
};

MedModelPartIO::MedModelPartIO(const std::filesystem::path& rFileName, const Flags Options)
    : mFileName(rFileName), mOptions(Options)
{
    KRATOS_TRY

    mpFileHandler = Kratos::make_shared<MedFileHandler>(rFileName, Options);

    KRATOS_CATCH("")
}

void MedModelPartIO::ReadModelPart(ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpFileHandler->IsReadMode()) << "MedModelPartIO needs to be created in read mode to read a ModelPart!" << std::endl;

    KRATOS_CATCH("")
}

void MedModelPartIO::WriteModelPart(const ModelPart& rThisModelPart)
{
    KRATOS_TRY

    // TODO what happens if this function is called multiple times?
    // will it overwrite the mesh?
    // or just crash?

    KRATOS_ERROR_IF(mpFileHandler->IsReadMode()) << "MedModelPartIO needs to be created in write mode to write a ModelPart!" << std::endl;

    KRATOS_ERROR_IF_NOT(rThisModelPart.GetProcessInfo().Has(DOMAIN_SIZE)) << "\"DOMAIN_SIZE\" is not defined in ModelPart " << rThisModelPart.FullName() << std::endl;

    const int dimension = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    MEDfileCommentWr(mpFileHandler->GetFileHandle(), "A 2D unstructured mesh : 15 nodes, 12 cells");

    const char axisname[MED_SNAME_SIZE] = "x y";
    const char unitname[MED_SNAME_SIZE] = "cm cm";

    med_err err = MEDmeshCr(
        mpFileHandler->GetFileHandle(),
        mpFileHandler->GetMeshName(),
        dimension , //spacedim
        dimension , //meshdim
        MED_UNSTRUCTURED_MESH,
        "Kratos med", // description
        "",
        MED_SORT_DTIT,
        MED_CARTESIAN,
        axisname,
        unitname);


    KRATOS_CATCH("")
}

void MedModelPartIO::WriteNodes(NodesContainerType const& rThisNodes)
{
    KRATOS_TRY

    const int dimension = mpFileHandler->GetDimension();
    const Vector nodal_coords = VariableUtils().GetCurrentPositionsVector(rThisNodes, dimension);

    // TODO find better solution than to copy
    const std::vector<double> vec_nodal_coords(nodal_coords.begin(), nodal_coords.end());

    med_err err = MEDmeshNodeCoordinateWr(
        mpFileHandler->GetFileHandle(),
        mpFileHandler->GetMeshName(),
        MED_NO_DT,
        MED_NO_IT,
        0.0,
        MED_FULL_INTERLACE,
        rThisNodes.size() / dimension,
        vec_nodal_coords.data());

    KRATOS_CATCH("")
}


void MedModelPartIO::WriteGeometries(GeometryContainerType const& rThisGeometries)
{
    KRATOS_TRY
    KRATOS_ERROR << "MedModelPartIO::WriteGeometries is not yet implemented!" << std::endl;
    KRATOS_CATCH("")
}

void MedModelPartIO::DivideInputToPartitions(SizeType NumberOfPartitions,
                                             const PartitioningInfo& rPartitioningInfo)
{
    // these are not used in ModelPartIO, as the partitioned files are created independently
    std::stringbuf dummy_strbuf;
    auto dummy_stream(Kratos::make_shared<std::iostream>(&dummy_strbuf));

    // ModelPartIO(dummy_stream).DivideInputToPartitions(
    //     NumberOfPartitions,
    //     rPartitioningInfo);
}

void MedModelPartIO::DivideInputToPartitions(Kratos::shared_ptr<std::iostream> * pStreams,
                                             SizeType NumberOfPartitions,
                                             const PartitioningInfo& rPartitioningInfo)
{
    // these are not used in ModelPartIO, streams are passed from outside
    std::stringbuf dummy_strbuf;
    auto dummy_stream(Kratos::make_shared<std::iostream>(&dummy_strbuf));

    // ModelPartIO(dummy_stream).DivideInputToPartitions(
    //     pStreams,
    //     NumberOfPartitions,
    //     rPartitioningInfo);
}

int MedModelPartIO::GetNumberOfMedMeshes() const
{
    return MEDnMesh(mpFileHandler->GetFileHandle());
}

} // namespace Kratos.
