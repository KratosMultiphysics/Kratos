//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Jordi Cotela
//

#if !defined(KRATOS_VORTICITY_UTILITIES_H )
#define  KRATOS_VORTICITY_UTILITIES_H

// External includes 

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// A set of functions to compute quantities of interest for turbulent flows.
template< std::size_t TDim >
class VorticityUtilities {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VorticityUtilities
    KRATOS_CLASS_POINTER_DEFINITION(VorticityUtilities);

    /// Type for an array of shape function gradient matrices
    typedef Geometry< Node<3> >::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Deleted default constructor
    VorticityUtilities() = delete;
    
    /// Deleted copy constructor.
    VorticityUtilities(VorticityUtilities const& rOther) = delete;

    /// Destructor.
    ~VorticityUtilities();

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    VorticityUtilities& operator=(VorticityUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    static void CalculateQValue(
        const Geometry<Node<3>>& rGeometry,
        const ShapeFunctionDerivativesArrayType& rShapeFunctionsGradients,
        std::vector<double>& rQValues);

    static void CalculateVorticityMagnitude(
        const Geometry<Node<3>>& rGeometry,
        const ShapeFunctionDerivativesArrayType& rShapeFunctionsGradients,
        std::vector<double>& rVorticityMagnitudes);

    static void CalculateVorticityVector(
        const Geometry<Node<3>>& rGeometry,
        const ShapeFunctionDerivativesArrayType& rShapeFunctionsGradients,
        std::vector<array_1d<double,3>>& rVorticities);

    ///@}
private:

    ///@name Private Operations
    ///@{

    static void NodalContributionToVorticityVector(
        const Matrix& rDN_DX,
        const array_1d<double, 3 > & rVelocity,
        const unsigned int iNode,
        array_1d<double,3>& rVorticity);

    ///@}

};  // Class VorticityUtilities

///@}

///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_VORTICITY_UTILITIES_H  defined 


