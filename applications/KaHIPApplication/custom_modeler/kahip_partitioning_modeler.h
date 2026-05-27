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

#pragma once

// System includes

// External libraries

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "modeler/modeler.h"

namespace Kratos
{

///@addtogroup KaHIPApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class KaHIPPartitioningModeler
 * @ingroup KaHIPApplication
 * @brief Modeler that partitions a Kratos mesh using KaHIP and loads the local
 *        partition into a @c ModelPart.
 * @details This modeler implements the complete partitioning workflow as a single
 *          @c SetupModelPart() call:
 *
 *          **Single-process or no-partition mode**
 *          The input .mdpa file is read directly into the model part without any
 *          partitioning.
 *
 *          **MPI file-based mode** (@c partition_in_memory = @c false, default)
 *          Rank 0 reads the full mesh, calls
 *          @c KaHIPDivideHeterogeneousInputProcess to partition it, and writes
 *          per-rank .mdpa files. All ranks then read their own file.
 *
 *          **MPI in-memory mode** (@c partition_in_memory = @c true)
 *          Rank 0 partitions, then per-rank data is distributed to all ranks
 *          via @c KaHIPDivideHeterogeneousInputInMemoryProcess (Scatterv). Each
 *          rank reads from its in-memory buffer without touching the filesystem.
 *
 *          **Parameters schema** (all optional, defaults shown):
 *          @code{.json}
 *          {
 *              "model_part_name":                        "Main",
 *              "input_filename":                         "",
 *              "number_of_partitions":                    0,
 *              "dimension":                               3,
 *              "verbosity":                               0,
 *              "synchronize_conditions":                  false,
 *              "partition_in_memory":                     false,
 *              "perform_partitioning":                    true,
 *              "skip_timer":                              true,
 *              "ignore_variables_not_in_solution_step_data": false,
 *              "data_communicator_name":                 "",
 *              "kahip_settings": {
 *                  "preconfiguration":  "eco",
 *                  "imbalance":          0.03,
 *                  "seed":               0,
 *                  "suppress_output":    true,
 *                  "num_trials":         1
 *              }
 *          }
 *          @endcode
 *
 *          When @c number_of_partitions is 0, the MPI communicator size is used.
 *
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KAHIP_APPLICATION) KaHIPPartitioningModeler : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(KaHIPPartitioningModeler);

    using SizeType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor (required by factory registration).
    KaHIPPartitioningModeler() : Modeler() {}

    /**
     * @brief Constructor.
     * @param rModel       The Model containing the model part to populate.
     * @param rParameters  Configuration parameters (see class docs for schema).
     */
    KaHIPPartitioningModeler(Model& rModel, Parameters rParameters);

    /// Destructor.
    ~KaHIPPartitioningModeler() override = default;

    ///@}
    ///@name Modeler interface
    ///@{

    /// Creates a new instance of this modeler (required by the modeler factory).
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<KaHIPPartitioningModeler>(rModel, ModelParameters);
    }

    /**
     * @brief Perform partitioning and load the local partition into the model part.
     * @details Implements the full workflow described in the class documentation.
     *          Called automatically during the analysis-stage setup.
     */
    void SetupModelPart() override;

    /**
     * @brief Return the default parameters schema for this modeler.
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return "KaHIPPartitioningModeler";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    void PrintData(std::ostream& rOStream) const override {}

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr;

    ///@}

}; // class KaHIPPartitioningModeler

///@}

} // namespace Kratos
