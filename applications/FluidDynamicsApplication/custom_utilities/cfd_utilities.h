//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_CFD_UTILITIES_H_INCLUDED)
#define KRATOS_CFD_UTILITIES_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) CFDUtilities
{
    class ReactionBasedYPlus
    {
    public:
        void CalculateData(ModelPart::ConditionType& rCondition);

        void CalculateYPlusAndUTau(double& rYPlus,
                                   array_1d<double, 3>& rUTau,
                                   const array_1d<double, 3>& rNormal);
    };

    class WallFunctionBasedYPlus
    {
    public:
        void CalculateData(ModelPart::ConditionType& rCondition);

        void CalculateYPlusAndUTau(double& rYPlus,
                                   array_1d<double, 3>& rUTau,
                                   const array_1d<double, 3>& rNormal);
    };

public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;
    using ElementType = ModelPart::ElementType;
    using ConditionType = ModelPart::ConditionType;
    using GeometryType = Geometry<NodeType>;

    /// Pointer definition of CFDUtilities
    KRATOS_CLASS_POINTER_DEFINITION(CFDUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    CFDUtilities(){};

    /// Destructor.
    ~CFDUtilities() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void static CalculateReactionBasedYPlus(ModelPart& rModelPart,
                                            const Variable<double>& rKinematicViscosityVariable,
                                            const Variable<array_1d<double, 3>>& rReactionVariable);

    double static CalculateLogarithmicYPlusLimit(const double VonKarman,
                                                 const double Smoothness,
                                                 const int MaxIterations = 20,
                                                 const double Tolerance = 1e-6);

    void static CalculateYPlusAndUtau(double& rYPlus,
                                      double& rUTau,
                                      const double WallVelocity,
                                      const double WallHeight,
                                      const double KinematicViscosity,
                                      const double VonKarman,
                                      const double Smoothness,
                                      const int MaxIterations = 20,
                                      const double Tolerance = 1e-6);

    void static CalculateWallFunctionBasedYPlus(ModelPart& rModelPart,
                                                const double VonKarman = 0.41,
                                                const double Smoothness = 5.2);

    void static CalculateNumberOfNeighbourConditions(ModelPart& rModelPart);

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
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void static CalculateConditionGeometryData(const GeometryType& rGeometry,
                                               const GeometryData::IntegrationMethod& rIntegrationMethod,
                                               Vector& rGaussWeights,
                                               Matrix& rNContainer);

    template <unsigned int TDim>
    void static CalculateNormal(array_1d<double, 3>& rNormal, const ConditionType& rCondition);

    double static CalculateWallHeight(const ConditionType& rCondition,
                                      const array_1d<double, 3>& rNormal);

    template <typename TDataType>
    TDataType static EvaluateInPoint(const GeometryType& rGeometry,
                                     const Variable<TDataType>& rVariable,
                                     const Vector& rShapeFunction,
                                     const int Step = 0);

    void static CalculateReactionBasedUTau(double& YPlus,
                                           array_1d<double, 3>& rUTau,
                                           const GeometryType& rGeometry,
                                           const double Density,
                                           const double KinematicViscosity,
                                           const double WallHeight,
                                           const double Area,
                                           const Variable<array_1d<double, 3>>& rReactionVariable);

    void static CalculateWallFunctionBasedYPlusUTau(double& YPlus,
                                                    array_1d<double, 3>& rUTau,
                                                    const GeometryType& rGeometry,
                                                    const double Density,
                                                    const double KinematicViscosity,
                                                    const double WallHeight,
                                                    const double VonKarman,
                                                    const double Smoothness,
                                                    const Vector& GaussPointShapeFunctions);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    CFDUtilities& operator=(CFDUtilities const& rOther);

    /// Copy constructor.
    CFDUtilities(CFDUtilities const& rOther);

    ///@}

}; // Class CFDUtilities

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const CFDUtilities& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_CFD_UTILITIES_H_INCLUDED  defined
