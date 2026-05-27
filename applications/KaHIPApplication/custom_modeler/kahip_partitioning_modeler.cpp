//     __ __      __  __________  ___                ___            __  _           
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____ 
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __ \ 
//  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /
// /_/ |_\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/ 
//                                  /_/   /_/                                       
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// Project includes
#include "includes/parallel_environment.h"
#include "includes/reorder_consecutive_model_part_io.h"
#include "includes/model_part_io.h"
#include "custom_modeler/kahip_partitioning_modeler.h"
#include "custom_processes/kahip_divide_heterogeneous_input_process.h"
#include "custom_processes/kahip_divide_heterogeneous_input_in_memory_process.h"

namespace Kratos {

// ─── Constructor ──────────────────────────────────────────────────────────────

KaHIPPartitioningModeler::KaHIPPartitioningModeler(Model& rModel, Parameters rParameters)
    : Modeler(rModel, rParameters),
      mpModel(&rModel)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

// ─── GetDefaultParameters ─────────────────────────────────────────────────────

const Parameters KaHIPPartitioningModeler::GetDefaultParameters() const
{
    return Parameters(R"({
        "model_part_name"                            : "Main",
        "input_filename"                             : "",
        "number_of_partitions"                       : 0,
        "dimension"                                  : 3,
        "verbosity"                                  : 0,
        "synchronize_conditions"                     : false,
        "partition_in_memory"                        : false,
        "perform_partitioning"                       : true,
        "skip_timer"                                 : true,
        "ignore_variables_not_in_solution_step_data" : false,
        "data_communicator_name"                     : "",
        "kahip_settings": {
            "preconfiguration": "eco",
            "imbalance":         0.03,
            "seed":              0,
            "suppress_output":   true,
            "num_trials":        1
        }
    })");
}

// ─── SetupModelPart ───────────────────────────────────────────────────────────

void KaHIPPartitioningModeler::SetupModelPart()
{
    KRATOS_TRY

    // ── Retrieve settings ───────────────────────────────────────────────────
    const std::string model_part_name = mParameters["model_part_name"].GetString();
    const std::string input_filename  = mParameters["input_filename"].GetString();
    const int         verbosity       = mParameters["verbosity"].GetInt();
    const bool        sync_conditions = mParameters["synchronize_conditions"].GetBool();
    const bool        in_memory       = mParameters["partition_in_memory"].GetBool();
    const bool        perform_part    = mParameters["perform_partitioning"].GetBool();
    const bool        skip_timer      = mParameters["skip_timer"].GetBool();
    const bool        ignore_vars     = mParameters["ignore_variables_not_in_solution_step_data"].GetBool();
    const std::string comm_name       = mParameters["data_communicator_name"].GetString();

    // Merge verbosity into the process settings so the process can log correctly
    Parameters kahip_settings = mParameters["kahip_settings"].Clone();
    if (!kahip_settings.Has("verbosity")) {
        kahip_settings.AddEmptyValue("verbosity");
    }
    kahip_settings["verbosity"].SetInt(verbosity);

    KRATOS_ERROR_IF(input_filename.empty())
        << "KaHIPPartitioningModeler: 'input_filename' must be specified." << std::endl;

    // ── Get model part and communicator ─────────────────────────────────────
    ModelPart& r_model_part = mpModel->HasModelPart(model_part_name)
                            ? mpModel->GetModelPart(model_part_name)
                            : mpModel->CreateModelPart(model_part_name);

    const DataCommunicator& r_comm = comm_name.empty()
        ? ParallelEnvironment::GetDefaultDataCommunicator()
        : ParallelEnvironment::GetDataCommunicator(comm_name);

    const int mpi_rank = r_comm.Rank();
    const int mpi_size = r_comm.Size();

    // ── IO flags ────────────────────────────────────────────────────────────
    auto io_flags = ModelPartIO::READ;
    if (skip_timer)    io_flags = ModelPartIO::SKIP_TIMER | io_flags;
    if (ignore_vars)   io_flags = ModelPartIO::IGNORE_VARIABLES_ERROR | io_flags;

    // ── Decide number of partitions ─────────────────────────────────────────
    const int user_nparts = mParameters["number_of_partitions"].GetInt();
    const int nparts      = (user_nparts > 0) ? user_nparts : mpi_size;

    // ── Single process or partitioning disabled ─────────────────────────────
    if (!perform_part || mpi_size == 1) {
        KRATOS_INFO_IF("KaHIPPartitioningModeler", verbosity > 0)
            << "Partitioning skipped (single process or perform_partitioning=false)."
            << " Reading '" << input_filename << "' directly." << std::endl;

        ModelPartIO io(input_filename, io_flags);
        io.ReadModelPart(r_model_part);
        return;
    }

    // ── MPI partitioning ────────────────────────────────────────────────────
    KRATOS_INFO_IF("KaHIPPartitioningModeler", verbosity > 0)
        << "Starting KaHIP partitioning of '" << input_filename
        << "' into " << nparts << " partitions"
        << (in_memory ? " (in-memory)" : " (file-based)") << "." << std::endl;

    if (in_memory) {
        // ── In-memory: scatter per-rank streams; no intermediate files ──────
        // ReorderConsecutiveModelPartIO ensures node/element IDs are gapless
        ReorderConsecutiveModelPartIO reorder_io(input_filename, io_flags);
        ModelPartIO serial_io(input_filename, io_flags);

        KaHIPDivideHeterogeneousInputInMemoryProcess process(
            reorder_io, serial_io, r_comm, kahip_settings, sync_conditions);
        process.Execute();

        // After Execute(), serial_io now reads from the local in-memory buffer
        serial_io.ReadModelPart(r_model_part);
    } else {
        // ── File-based: rank 0 partitions, all ranks read their own file ────
        if (mpi_rank == 0) {
            ReorderConsecutiveModelPartIO reorder_io(input_filename, io_flags);

            KaHIPDivideHeterogeneousInputProcess process(
                reorder_io,
                static_cast<std::size_t>(nparts),
                kahip_settings,
                sync_conditions);
            process.Execute();
        }

        // Barrier: wait for rank 0 to finish writing per-rank files
        r_comm.Barrier();

        // Deduce per-rank filename: <base>_partitioned/<base>_<rank>
        // This follows ModelPartIO's DivideInputToPartitions naming convention.
        const std::string base     = input_filename;
        const std::size_t sep      = base.rfind('/');
        const std::string dir_name = (sep != std::string::npos ? base.substr(0, sep + 1) : "")
                                   + (sep != std::string::npos ? base.substr(sep + 1) : base)
                                   + "_partitioned";
        const std::string mesh_name = (sep != std::string::npos ? base.substr(sep + 1) : base);
        const std::string rank_file = dir_name + "/" + mesh_name + "_" + std::to_string(mpi_rank);

        KRATOS_INFO_IF("KaHIPPartitioningModeler", verbosity > 0)
            << "Rank " << mpi_rank << " reading partition file '" << rank_file << "'." << std::endl;

        ModelPartIO rank_io(rank_file, io_flags);
        rank_io.ReadModelPart(r_model_part);
    }

    KRATOS_INFO_IF("KaHIPPartitioningModeler", verbosity > 0)
        << "KaHIPPartitioningModeler::SetupModelPart() done on rank " << mpi_rank << "." << std::endl;

    KRATOS_CATCH("")
}

} // namespace Kratos
