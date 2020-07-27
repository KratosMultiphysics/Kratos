//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala

#if !defined(CHIMERA_DISTANCE_CALCULATION_UTILITY )
#define  CHIMERA_DISTANCE_CALCULATION_UTILITY



// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/variational_distance_calculation_process.h"
#include "utilities/parallel_levelset_distance_calculator.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/gather_modelpart_on_all_ranks.h"


#ifdef KRATOS_USING_MPI
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/includes/mpi_communicator.h"
#endif
#include "includes/data_communicator.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Utility for calculating the Distance on a given modelpart
template <int TDim>
class ChimeraDistanceCalculationUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ChimeraDistanceCalculationUtility
    KRATOS_CLASS_POINTER_DEFINITION(ChimeraDistanceCalculationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ChimeraDistanceCalculationUtility() = delete;

    /// Destructor.
    /// Deleted copy constructor
    ChimeraDistanceCalculationUtility(const ChimeraDistanceCalculationUtility& rOther) = delete;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Calculates distance on the whole of rVolumeModelPart from rSkinModelPart
     * @param rVolumeModelPart The background modelpart where distances are calculated.
     * @param rSkinModelPart The skin modelpart from where the distances are calculated
     */
    static inline void CalculateDistance(ModelPart &rVolumeModelPart, ModelPart &rSkinModelPart)
    {

        //typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        //typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        //typedef LinearSolverFactory<SparseSpaceType, LocalSparseSpaceType> LinearSolverFactoryType;
        //typedef LinearSolver<SparseSpaceType, TLocalSpaceType> LinearSolverType;
        //typedef VariationalDistanceCalculationProcess<TDim, TSparseSpaceType, TLocalSpaceType, LinearSolverType> VariationalDistanceCalculationProcessType;
        typedef CalculateDistanceToSkinProcess<TDim> CalculateDistanceToSkinProcessType;
        const int n_vol_nodes = static_cast<int>(rVolumeModelPart.NumberOfNodes());

#pragma omp parallel for
        for (int i_node = 0; i_node < n_vol_nodes; ++i_node)
        {
            auto it_node = rVolumeModelPart.NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(DISTANCE, 0) = 0.0;
            it_node->FastGetSolutionStepValue(DISTANCE, 1) = 0.0;
            it_node->SetValue(DISTANCE, 0.0);
        }
        Model &current_model = rVolumeModelPart.GetModel();
        const DataCommunicator &r_comm =
            rVolumeModelPart.GetCommunicator().GetDataCommunicator();
        ModelPart &r_gathered_skin_mp = r_comm.IsDistributed() ? current_model.CreateModelPart("GatheredSkin") : rSkinModelPart;

#ifdef KRATOS_USING_MPI
        ModelPart& r_root_mp = rVolumeModelPart.GetRootModelPart();
        Communicator::Pointer pnew_comm = Kratos::make_shared< MPICommunicator >(&rVolumeModelPart.GetNodalSolutionStepVariablesList(), r_root_mp.GetCommunicator().GetDataCommunicator());
        r_gathered_skin_mp.SetCommunicator(pnew_comm);
        GatherModelPartOnAllRanksUtility::GatherModelPartOnAllRanks(rSkinModelPart, r_gathered_skin_mp);
        ParallelFillCommunicator(r_gathered_skin_mp).Execute();
#endif
        const int n_skin_nodes = static_cast<int>(r_gathered_skin_mp.NumberOfNodes());
        if(n_skin_nodes != 0 && n_vol_nodes != 0)
            CalculateDistanceToSkinProcessType(rVolumeModelPart, r_gathered_skin_mp).Execute();


        // Parameters amgcl_settings(R"(
        //     {
        //         "solver_type"                   : "amgcl",
        //         "tolerance"                     : 0.001,
        //         "max_iteration"                 : 200,
        //         "krylov_type"                   : "gmres",
        //         "smoother_type"                 : "ilu0",
        //         "verbosity"                     : 0
        //     }
        //     )");
        // const int max_iterations = 1;
        // LinearSolverFactoryType const &linear_solver_factory = KratosComponents<LinearSolverFactoryType>::Get("amgcl");
        // auto amgcl_solver = linear_solver_factory.Create(amgcl_settings);
        // VariationalDistanceCalculationProcessType(rVolumeModelPart, amgcl_solver, max_iterations, VariationalDistanceCalculationProcessType::CALCULATE_EXACT_DISTANCES_TO_PLANE).Execute();

        unsigned int max_level = 100;
		double max_distance = 200;
        auto p_distance_smoother = Kratos::make_shared<ParallelDistanceCalculator<TDim>>();
        p_distance_smoother->CalculateDistances(rVolumeModelPart, DISTANCE, NODAL_AREA, max_level, max_distance);

        VariableUtils().CopyScalarVar(DISTANCE, CHIMERA_DISTANCE, rVolumeModelPart.Nodes());
        current_model.DeleteModelPart("GatheredSkin");
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{


    ///@}

}; // Class ChimeraDistanceCalculationUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // DISTANCE_CALCULATION_UTILITY  defined


