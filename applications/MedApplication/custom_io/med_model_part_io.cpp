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

        if (Options.IsNot(IO::WRITE)) {
            // check if file exists
            KRATOS_ERROR_IF(!std::filesystem::exists(rFileName) && Options.IsNot(IO::WRITE)) << "File " << rFileName << " does not exist!" << std::endl;

            // basic checks if the file is compatible with the MED library
            med_bool hdf_ok;
            med_bool med_ok;
            const med_err err = MEDfileCompatibility(rFileName.c_str(), &hdf_ok, &med_ok);

            KRATOS_ERROR_IF(err != 0) << "A problem occured while trying to check the compatibility of file " << rFileName << "!" << std::endl;
            KRATOS_ERROR_IF(hdf_ok != MED_TRUE) << "A problem with HDF occured while trying to open file " << rFileName << "!" << std::endl;
            KRATOS_ERROR_IF(med_ok != MED_TRUE) << "A problem with MED occured while trying to open file " << rFileName << "!" << std::endl;
        }

        // Set the mode (consistent with ModelPartIO)
        // read only by default, unless other settings are specified
        med_access_mode open_mode = MED_ACC_RDONLY;

        if (Options.Is(IO::APPEND)) {
            open_mode = MED_ACC_RDEXT;
        } else if (Options.Is(IO::WRITE)) {
            open_mode = MED_ACC_RDWR; // or MED_ACC_CREAT?
        }

        mFileHandle = MEDfileOpen(rFileName.c_str(), open_mode);
        KRATOS_ERROR_IF(mFileHandle < 0) << "A problem occured while opening file " << rFileName << "!" << std::endl;

        KRATOS_CATCH("")
    }

    med_idt GetFileHandle() const
    {
        return mFileHandle;
    }

    ~MedFileHandler()
    {
        KRATOS_TRY

        KRATOS_WARNING_IF("MedModelPartIO", MEDfileClose(mFileHandle) < 0) << "Closing of file " << mFileName << " failed!" << std::endl;

        KRATOS_CATCH("")
    }

private:
    med_idt mFileHandle;
    std::filesystem::path mFileName;
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

    // KRATOS_ERROR << "ReadModelPart is not yet implemented ..." << std::endl;

    KRATOS_CATCH("")
}

void MedModelPartIO::WriteModelPart(const ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mOptions.Is(IO::WRITE) || mOptions.Is(IO::APPEND)) << "MedModelPartIO needs to be created in write or append mode to write a ModelPart!" << std::endl;

    KRATOS_ERROR_IF_NOT(rThisModelPart.GetProcessInfo().Has(DOMAIN_SIZE)) << "\"DOMAIN_SIZE\" is not defined in ModelPart " << rThisModelPart.FullName() << std::endl;

    const int dimension = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    MEDfileCommentWr(mpFileHandler->GetFileHandle(),"A 2D unstructured mesh : 15 nodes, 12 cells");

    const char axisname[2*MED_SNAME_SIZE+1] = "x               y               ";
    const char unitname[2*MED_SNAME_SIZE+1] = "cm              cm              ";

    med_err err = MEDmeshCr(
        mpFileHandler->GetFileHandle(),
        rThisModelPart.Name().c_str(), // meshname
        dimension , //spacedim
        dimension , //meshdim
        MED_UNSTRUCTURED_MESH,
        "Kratos med", // description
        "",
        MED_SORT_DTIT,
        MED_CARTESIAN,
        axisname,
        unitname);

    // TODO find better solution than to copy
    const Vector nodal_coords = VariableUtils().GetCurrentPositionsVector(rThisModelPart.Nodes(), dimension);

    std::vector<double> vec_nodal_coords(nodal_coords.begin(), nodal_coords.end());

    err = MEDmeshNodeCoordinateWr(
        mpFileHandler->GetFileHandle(),
        rThisModelPart.Name().c_str(),
        MED_NO_DT,
        MED_NO_IT,
        0.0,
        MED_FULL_INTERLACE,
        vec_nodal_coords.size(),
        vec_nodal_coords.data());

    KRATOS_CATCH("")
}

void MedModelPartIO::DivideInputToPartitions(SizeType NumberOfPartitions,
                                             GraphType const& rDomainsColoredGraph,
                                             PartitionIndicesType const& rNodesPartitions,
                                             PartitionIndicesType const& rElementsPartitions,
                                             PartitionIndicesType const& rConditionsPartitions,
                                             PartitionIndicesContainerType const& rNodesAllPartitions,
                                             PartitionIndicesContainerType const& rElementsAllPartitions,
                                             PartitionIndicesContainerType const& rConditionsAllPartitions)
{
    // these are not used in ModelPartIO, as the partitioned files are created independently
    std::stringbuf dummy_strbuf;
    auto dummy_stream(Kratos::make_shared<std::iostream>(&dummy_strbuf));

    ModelPartIO(dummy_stream).DivideInputToPartitions(
        NumberOfPartitions,
        rDomainsColoredGraph,
        rNodesPartitions,
        rElementsPartitions,
        rConditionsPartitions,
        rNodesAllPartitions,
        rElementsAllPartitions,
        rConditionsAllPartitions);
}

void MedModelPartIO::DivideInputToPartitions(Kratos::shared_ptr<std::iostream> * pStreams,
                                             SizeType NumberOfPartitions,
                                             GraphType const& rDomainsColoredGraph,
                                             PartitionIndicesType const& rNodesPartitions,
                                             PartitionIndicesType const& rElementsPartitions,
                                             PartitionIndicesType const& rConditionsPartitions,
                                             PartitionIndicesContainerType const& rNodesAllPartitions,
                                             PartitionIndicesContainerType const& rElementsAllPartitions,
                                             PartitionIndicesContainerType const& rConditionsAllPartitions)
{
    // these are not used in ModelPartIO, streams are passed from outside
    std::stringbuf dummy_strbuf;
    auto dummy_stream(Kratos::make_shared<std::iostream>(&dummy_strbuf));

    ModelPartIO(dummy_stream).DivideInputToPartitions(
        pStreams,
        NumberOfPartitions,
        rDomainsColoredGraph,
        rNodesPartitions,
        rElementsPartitions,
        rConditionsPartitions,
        rNodesAllPartitions,
        rElementsAllPartitions,
        rConditionsAllPartitions);
}



} // namespace Kratos.
