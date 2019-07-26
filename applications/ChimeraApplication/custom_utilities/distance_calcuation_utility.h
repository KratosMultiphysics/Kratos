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

#if !defined(DISTANCE_CALCULATION_UTILITY )
#define  DISTANCE_CALCULATION_UTILITY



// System includes


// External includes
#include "omp.h"


// Project includes
#include "includes/define.h"
#include "processes/variational_distance_calculation_process.h"
#include "utilities/parallel_levelset_distance_calculator.h"
#include "custom_processes/calculate_signed_distance_to_2d_condition_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Utility for calculating the Distance on a given modelpart
template <int TDim, class TSparseSpaceType, class TLocalSpaceType>
class DistanceCalculationUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceCalculationUtility
    KRATOS_CLASS_POINTER_DEFINITION(DistanceCalculationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistanceCalculationUtility() = delete;

    /// Destructor.
    /// Deleted copy constructor
    DistanceCalculationUtility(const DistanceCalculationUtility& rOther) = delete;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Calculates distance on the whole of rBackgroundModelPart from rSkinModelPart
     * @param rBackgroundModelPart The background modelpart where distances are calculated.
     * @param rSkinModelPart The skin modelpart from where the distances are calculated
     */
    static inline void CalculateDistance(ModelPart &rBackgroundModelPart, ModelPart &rSkinModelPart, const double OverlapDistance)
    {
        typedef CalculateDistanceToSkinProcess<TDim> CalculateDistanceToSkinProcessType;
        IndexType nnodes = static_cast<IndexType>(rBackgroundModelPart.NumberOfNodes());
        auto p_distance_smoother = Kratos::make_shared<ParallelDistanceCalculator<TDim>>();
        unsigned int max_level = 100;
		double max_distance = 200;

#pragma omp parallel for
        for (IndexType i_node = 0; i_node < nnodes; ++i_node)
        {
            auto it_node = rBackgroundModelPart.NodesBegin() + i_node;
            double &node_distance = it_node->FastGetSolutionStepValue(DISTANCE);
            node_distance = 0;
        }

        CalculateDistanceToSkinProcessType(rBackgroundModelPart, rSkinModelPart).Execute();

        p_distance_smoother->CalculateDistances(rBackgroundModelPart, DISTANCE, NODAL_AREA, max_level, max_distance);
    }

    /**
     * @brief Calculates distance on the whole of rBackgroundModelPart from rSkinModelPart
     * @param rBackgroundModelPart The background modelpart where distances are calculated.
     * @param rSkinModelPart The skin modelpart from where the distances are calculated
     */
    static inline void CalculateDistanceChimeraApplication(ModelPart &rBackgroundModelPart, ModelPart &rSkinModelPart, const double OverlapDistance)
    {
        IndexType nnodes = static_cast<IndexType>(rBackgroundModelPart.NumberOfNodes());
        auto p_distance_smoother = Kratos::make_shared<ParallelDistanceCalculator<TDim>>();
        unsigned int max_level = 100;
		double max_distance = 200;

#pragma omp parallel for
        for (IndexType i_node = 0; i_node < nnodes; ++i_node)
        {
            auto it_node = rBackgroundModelPart.NodesBegin() + i_node;
            double &node_distance = it_node->FastGetSolutionStepValue(DISTANCE);
            node_distance = 0;
        }

        if(TDim ==2)
            CalculateSignedDistanceTo2DConditionSkinProcess(rSkinModelPart, rBackgroundModelPart).Execute();
        else if(TDim==3)
            CalculateSignedDistanceTo3DConditionSkinProcess(rSkinModelPart, rBackgroundModelPart).Execute();

        p_distance_smoother->CalculateDistances(rBackgroundModelPart, DISTANCE, NODAL_AREA, max_level, max_distance);
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


private:
    ///@name Static Member Variables
    ///@{




    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class DistanceCalculationUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // DISTANCE_CALCULATION_UTILITY  defined


