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
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

namespace Kratos {

namespace {

static const std::map<GeometryData::KratosGeometryType, med_geometry_type> KratosToMedGeometryType {
    { GeometryData::KratosGeometryType::Kratos_Point2D,          MED_POINT1 },
    { GeometryData::KratosGeometryType::Kratos_Point3D,          MED_POINT1 },

    { GeometryData::KratosGeometryType::Kratos_Line2D2,          MED_SEG2 },
    { GeometryData::KratosGeometryType::Kratos_Line3D2,          MED_SEG2 },
    { GeometryData::KratosGeometryType::Kratos_Line2D3,          MED_SEG3 },
    { GeometryData::KratosGeometryType::Kratos_Line3D3,          MED_SEG3 },

    { GeometryData::KratosGeometryType::Kratos_Triangle2D3,      MED_TRIA3 },
    { GeometryData::KratosGeometryType::Kratos_Triangle3D3,      MED_TRIA3 },
    { GeometryData::KratosGeometryType::Kratos_Triangle2D6,      MED_TRIA6 },
    { GeometryData::KratosGeometryType::Kratos_Triangle3D6,      MED_TRIA6 },

    { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, MED_QUAD4 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, MED_QUAD4 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, MED_QUAD8 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, MED_QUAD8 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9, MED_QUAD9 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9, MED_QUAD9 },

    { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,    MED_TETRA4 },
    { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10,   MED_TETRA10 },

    { GeometryData::KratosGeometryType::Kratos_Pyramid3D5,       MED_PYRA5 },
    { GeometryData::KratosGeometryType::Kratos_Pyramid3D13,      MED_PYRA13 },

    { GeometryData::KratosGeometryType::Kratos_Prism3D6,         MED_PENTA6 },
    { GeometryData::KratosGeometryType::Kratos_Prism3D15,        MED_PENTA15 },

    { GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,     MED_HEXA8 },
    { GeometryData::KratosGeometryType::Kratos_Hexahedra3D20,    MED_HEXA20 },
    { GeometryData::KratosGeometryType::Kratos_Hexahedra3D27,    MED_HEXA27 }
};

void CheckConnectivitiesSize(
    const std::size_t ExpectedSize,
    const std::vector<med_int>& Conns)
{
    KRATOS_DEBUG_ERROR_IF_NOT(Conns.size() == ExpectedSize) << "Connectivities must have a size of " << ExpectedSize << ", but have " << Conns.size() << "!" << std::endl;
};

std::function<void(std::vector<med_int>&)> GetReorderFunction(const med_geometry_type MedGeomType)
{
    switch (MedGeomType)
    {
    case MED_TRIA3:
        return [](std::vector<med_int>& Connectivities){
            CheckConnectivitiesSize(3, Connectivities);
            std::swap(Connectivities[1], Connectivities[2]);
        };

    case MED_QUAD4:
        return [](std::vector<med_int>& Connectivities){
            CheckConnectivitiesSize(4, Connectivities);
            std::swap(Connectivities[1], Connectivities[3]);
        };

    default:
        return [](std::vector<med_int>& Connectivities){
            // does nothing if no reordering is needed
        };
    }
}

}

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
            med_int space_dim = MEDmeshnAxis(mFileHandle, 1);
            med_int mesh_dim;
            med_mesh_type mesh_type;
            std::string description(MED_COMMENT_SIZE+1, '\0');
            std::string dt_unit(MED_SNAME_SIZE+1, '\0');
            med_sorting_type sorting_type;
            med_int n_step;
            med_axis_type axis_type;
            std::string axis_name(MED_SNAME_SIZE*space_dim+1, '\0');
            std::string axis_unit(MED_SNAME_SIZE*space_dim+1, '\0');

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

    MEDfileCommentWr(mpFileHandler->GetFileHandle(), "A 2D unstructured mesh : 15 nodes, 12 cells");

    constexpr med_int dimension = 3;

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


//     KRATOS_CATCH("")
// }

// void MedModelPartIO::WriteNodes(NodesContainerType const& rThisNodes)
// {
//     KRATOS_TRY

    const Vector nodal_coords = VariableUtils().GetCurrentPositionsVector(rThisModelPart.Nodes(), dimension);

    // TODO find better solution than to copy
    const std::vector<double> vec_nodal_coords(nodal_coords.begin(), nodal_coords.end());

    // add check and warning if ModelPart is empty
    err = MEDmeshNodeCoordinateWr(
        mpFileHandler->GetFileHandle(),
        mpFileHandler->GetMeshName(),
        MED_NO_DT,
        MED_NO_IT,
        0.0,
        MED_FULL_INTERLACE,
        rThisModelPart.NumberOfNodes(),
        vec_nodal_coords.data());

    if (rThisModelPart.NumberOfGeometries() > 0) {
        std::string geometry_name;

        auto it_geom_begin = rThisModelPart.GeometriesBegin();

        auto geom_type_current = it_geom_begin->GetGeometryType();
        auto it_geom_current = it_geom_begin;
        std::size_t n_points_current = it_geom_current->PointsNumber();


        using ConnectivitiesType = std::vector<med_int>;
        using ConnectivitiesVector = std::vector<ConnectivitiesType>;

        ConnectivitiesVector connectivities;
        connectivities.reserve(rThisModelPart.NumberOfGeometries()/3); // assuming that three different types of geometries exist

        auto write_geometries = [this](
            const GeometryData::KratosGeometryType GeomType,
            const std::size_t NumberOfPoints,
            ConnectivitiesVector& Connectivities) {

                const auto med_geom_type = KratosToMedGeometryType.at(GeomType);
                const auto reorder_fct = GetReorderFunction(med_geom_type);

                auto GetMedConnectivities = [&reorder_fct](
                    const std::size_t NumberOfPoints,
                    ConnectivitiesVector& Connectivities) {

                    std::vector<med_int> med_conn(Connectivities.size() * NumberOfPoints);

                    // flatten an reorder the connectivities
                    IndexPartition(Connectivities.size()).for_each([&](const std::size_t i) {
                        reorder_fct(Connectivities[i]);
                        // for (std::size_t p=0; p<NumberOfPoints; ++p) {
                        //     med_conn[i*NumberOfPoints+p] = Connectivities[i][p];
                        // }

                        std::copy(Connectivities[i].begin(), Connectivities[i].end(), med_conn.begin()+(i*NumberOfPoints));
                    });

                    return med_conn;
                };

            const std::vector<med_int> med_conn = GetMedConnectivities(NumberOfPoints, Connectivities);

            if (MEDmeshElementConnectivityWr (
                mpFileHandler->GetFileHandle(),
                mpFileHandler->GetMeshName(),
                MED_NO_DT, MED_NO_IT , 0.0,
                MED_CELL, med_geom_type ,
                MED_NODAL, MED_FULL_INTERLACE,
                Connectivities.size(), med_conn.data()) < 0) {
                    KRATOS_ERROR << "Writing failed!" << std::endl;
            }
        };

        for (std::size_t i=0; i<rThisModelPart.NumberOfGeometries(); ++i) {
            auto this_geom_type = it_geom_current->GetGeometryType();
            if (this_geom_type != geom_type_current) {
                // write the previously collected geometries
                write_geometries(geom_type_current, n_points_current, connectivities);

                connectivities.clear();

                geom_type_current = this_geom_type;
                n_points_current = it_geom_current->PointsNumber();
            }

            ConnectivitiesType conn;
            for (const auto& r_point : it_geom_current->Points()) {
                conn.push_back(r_point.Id());
            }
            connectivities.push_back(conn);

            ++it_geom_current;
        }

        // write the last batch of geometries
        write_geometries(geom_type_current, n_points_current, connectivities);
    }

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

    ModelPartIO(dummy_stream).DivideInputToPartitions(
        NumberOfPartitions,
        rPartitioningInfo);
}

void MedModelPartIO::DivideInputToPartitions(Kratos::shared_ptr<std::iostream> * pStreams,
                                             SizeType NumberOfPartitions,
                                             const PartitioningInfo& rPartitioningInfo)
{
    // these are not used in ModelPartIO, streams are passed from outside
    std::stringbuf dummy_strbuf;
    auto dummy_stream(Kratos::make_shared<std::iostream>(&dummy_strbuf));

    ModelPartIO(dummy_stream).DivideInputToPartitions(
        pStreams,
        NumberOfPartitions,
        rPartitioningInfo);
}

int MedModelPartIO::GetNumberOfMedMeshes() const
{
    return MEDnMesh(mpFileHandler->GetFileHandle());
}

} // namespace Kratos.
