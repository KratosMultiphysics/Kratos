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


namespace Kratos
{

///@name Kratos Classes
///@{

/// Utility for calculating the Distance on a given modelpart
template <int TDim>
class KRATOS_API(CHIMERA_APPLICATION) ChimeraDistanceCalculationUtility
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
     * @brief Calculates distance on the whole of rBackgroundModelPart from rSkinModelPart
     * @param rBackgroundModelPart The background modelpart where distances are calculated.
     * @param rSkinModelPart The skin modelpart from where the distances are calculated
     */
    static inline void CalculateDistance(ModelPart &rBackgroundModelPart, ModelPart &rSkinModelPart)
    {
        typedef CalculateDistanceToSkinProcess<TDim> CalculateDistanceToSkinProcessType;
        const int nnodes = static_cast<int>(rBackgroundModelPart.NumberOfNodes());

#pragma omp parallel for
        for (int i_node = 0; i_node < nnodes; ++i_node)
        {
            auto it_node = rBackgroundModelPart.NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(DISTANCE, 0) = 0.0;
            it_node->FastGetSolutionStepValue(DISTANCE, 1) = 0.0;
            it_node->SetValue(DISTANCE, 0.0);
        }

        CalculateDistanceToSkinProcessType(rBackgroundModelPart, rSkinModelPart).Execute();

        unsigned int max_level = 100;
		double max_distance = 200;
        auto p_distance_smoother = Kratos::make_shared<ParallelDistanceCalculator<TDim>>();
        p_distance_smoother->CalculateDistances(rBackgroundModelPart, DISTANCE, NODAL_AREA, max_level, max_distance);

        VariableUtils().CopyScalarVar(DISTANCE, CHIMERA_DISTANCE, rBackgroundModelPart.Nodes());
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


