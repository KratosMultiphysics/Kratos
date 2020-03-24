//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_SHALLOW_WATER_UTILITIES_H_INCLUDED
#define KRATOS_SHALLOW_WATER_UTILITIES_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/model_part.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(SHALLOW_WATER_APPLICATION) ShallowWaterUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ShallowWaterUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ShallowWaterUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    /// Destructor.

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ComputeFreeSurfaceElevation(ModelPart& rModelPart);

    void ComputeHeightFromFreeSurface(ModelPart& rModelPart);

    void ComputeVelocity(ModelPart& rModelPart);

    void ComputeMomentum(ModelPart& rModelPart);

    void ComputeAccelerations(ModelPart& rModelPart);

    void FlipScalarVariable(Variable<double>& rOriginVariable, Variable<double>& rDestinationVariable, ModelPart& rModelPart);

    void IdentifySolidBoundary(ModelPart& rModelPart, double SeaWaterLevel, Flags SolidBoundaryFlag);

    void IdentifyWetDomain(ModelPart& rModelPart, Flags WetFlag, double Thickness = 0.0);

    void ResetDryDomain(ModelPart& rModelPart, double Thickness = 0.0);

    template<class TContainerType>
    void DeactivateDryEntities(TContainerType& rContainer, Flags WetFlag)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rContainer.size()); ++i)
        {
            auto it = rContainer.begin() + i;
            it->Set(ACTIVE, it->Is(WetFlag));
        }
    }

    void ComputeVisualizationWaterHeight(ModelPart& rModelPart, Flags WetFlag, double SeaWaterLevel = 0.0);

    void ComputeVisualizationWaterSurface(ModelPart& rModelPart);

    void NormalizeVector(ModelPart& rModelPart, Variable<array_1d<double,3>>& rVariable);

    template<class TVarType>
    void CopyVariableToPreviousTimeStep(ModelPart& rModelPart, TVarType& rVariable)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
        {
            auto const it_node = rModelPart.NodesBegin() + i;
            it_node->FastGetSolutionStepValue(rVariable,1) = it_node->FastGetSolutionStepValue(rVariable);
        }
    }

    void SetMinimumValue(ModelPart& rModelPart, const Variable<double>& rVariable, double MinValue);

    /**
     * @brief Computes the root mean square of a double or component type of non historical variable
     * @param rVariable reference to the variable to compute
     * @param rWeightVariable reference to the weighting variable
     * @param rContainer Reference to the objective container
     */
    template<class TVarType, class TContainerType>
    double RootMeanSquareNonHistorical(
        const TVarType& rVariable,
        const Variable<double>& rWeightVariable,
        TContainerType& rContainer)
    {
        double rms_sum = 0.0;
        double weight_sum = 0.0;

        const auto it_begin = rContainer.begin();

        #pragma omp parallel for reduction(+:rms_sum, weight_sum)
        for (int k = 0; k < static_cast<int>(rContainer.size()); ++k)
        {
            const auto it = it_begin + k;
            rms_sum += it->GetValue(rWeightVariable) * std::pow(it->GetValue(rVariable), 2);
            weight_sum += it->GetValue(rWeightVariable);
        }

        return std::sqrt(rms_sum / weight_sum);
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

}; // Class ShallowWaterUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SHALLOW_WATER_UTILITIES_H_INCLUDED  defined
