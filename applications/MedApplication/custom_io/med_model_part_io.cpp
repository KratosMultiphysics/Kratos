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

// Ensure that the type of the file handle used by Kratos is the same as what is used by MED
// MED will not work properly if this check fails!
static_assert(std::is_same<Kratos::MedModelPartIO::MedFileHandleType, med_idt>::value,
    "Detecting the type of the MED file handle failed, please check the compilation/versions of MED and HDF5!");

namespace Kratos {

MedModelPartIO::MedModelPartIO(const std::filesystem::path& rFileName)
    : mFileName(rFileName)
{
    KRATOS_TRY

    mFileHandle = MEDfileOpen("test1.med", MED_ACC_RDWR);
    if (mFileHandle < 0) {
        KRATOS_WATCH("Erreur à la creation du fichier");
    }


    KRATOS_CATCH("")
}

MedModelPartIO::~MedModelPartIO()
{
    if (MEDfileClose(mFileHandle) < 0)
    std::cout << "Closing failed ...";
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

} // namespace Kratos.
