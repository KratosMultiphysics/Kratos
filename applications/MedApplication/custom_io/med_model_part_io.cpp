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
#include "includes/model_part_io.h"
#include "med_model_part_io.h"

namespace Kratos {

MedModelPartIO::MedModelPartIO(const std::filesystem::path& rFileName, const Flags Options)
    : mFileName(rFileName)
{
    KRATOS_TRY

    mpFileHandler = Kratos::make_shared<MedFileHandler>(rFileName, Options);


    med_int v_hdf_major, v_hdf_minor, v_hdf_release;
    med_int v_med_major, v_med_minor, v_med_release;

    MEDlibraryHdfNumVersion(&v_hdf_major, &v_hdf_minor, &v_hdf_release);
    MEDlibraryNumVersion(&v_med_major, &v_med_minor, &v_med_release);

    KRATOS_WATCH(v_hdf_major)
    KRATOS_WATCH(v_hdf_minor)
    KRATOS_WATCH(v_hdf_release)

    KRATOS_WATCH(v_med_major)
    KRATOS_WATCH(v_med_minor)
    KRATOS_WATCH(v_med_release)

    KRATOS_WATCH(MED_MAJOR_NUM)
    KRATOS_WATCH(MED_NUM_MINEUR)
    KRATOS_WATCH(MED_NUM_RELEASE)

    // TODO check macros vs what comes from the functions (both for HDF abd MED)
    // might indicate that versions are differen, aka different versions of the library loaded at runtime!
    // probably do in the MedApplication class, in register (actually also the filehandle type check could/should be done there ...)

    med_err err;

    med_bool hdf_ok;
    med_bool med_ok;
    err = MEDfileCompatibility(rFileName.c_str(), &hdf_ok, &med_ok);
    KRATOS_WATCH(err)

    KRATOS_ERROR_IF(hdf_ok != MED_TRUE) << "HDF is incompatible" << std::endl;
    KRATOS_ERROR_IF(med_ok != MED_TRUE) << "MED is incompatible" << std::endl;

    KRATOS_WATCH("About to open trhge file")

    // mFileHandle = MEDfileOpen(rFileName.c_str(), MED_ACC_RDWR);
    // if (mFileHandle < 0) {
    //     KRATOS_WATCH("Erreur à la creation du fichier");
    // }

    // med_int v_major;
    // med_int v_minor;
    // med_int v_release;

    // err = MEDfileNumVersionRd(mFileHandle, &v_major, &v_minor, &v_release);
    // KRATOS_WATCH(err)


    // KRATOS_WATCH(mFileHandle)
    // KRATOS_WATCH(v_major)
    // KRATOS_WATCH(v_minor)
    // KRATOS_WATCH(v_release)


    KRATOS_CATCH("")
}

void MedModelPartIO::ReadModelPart(ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR << "ReadModelPart is not yet implemented ..." << std::endl;

    KRATOS_CATCH("")
}

void MedModelPartIO::WriteModelPart(const ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR << "WriteModelPart is not yet implemented ..." << std::endl;

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

} // namespace Kratos.
