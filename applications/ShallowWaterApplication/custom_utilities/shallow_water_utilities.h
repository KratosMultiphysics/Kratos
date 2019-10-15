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

    void UpdatePrimitiveVariables(ModelPart& rModelPart);

    void UpdatePrimitiveVariables(ModelPart& rModelPart, double Epsilon);

    void ComputeAccelerations(ModelPart& rModelPart);

    void FlipScalarVariable(Variable<double>& rOriginVariable, Variable<double>& rDestinationVariable, ModelPart& rModelPart);

    void IdentifySolidBoundary(ModelPart& rModelPart, double SeaWaterLevel, Flags SolidBoundaryFlag);

    void IdentifyWetDomain(ModelPart& rModelPart, Flags WetFlag, double Thickness = 0.0);

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
