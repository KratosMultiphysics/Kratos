//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_FLUID_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_FLUID_ELEMENT_DATA_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"

#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

/// Auxiliary class to hold data for elements based on FluidElement
/** TODO: see how this can work for a generic number of stored arguments.
 */
template <unsigned int TDim, unsigned int TNumNodes>
class FluidElementData
{
public:
    ///@name Type Definitions
    ///@{

    typedef array_1d<double, TNumNodes> ScalarDataType;

    typedef boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> VectorDataType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor, accepting the geometry of the element.
    FluidElementData(Geometry<Node<3>>& rGeometry);

    ///@}
    ///@name Public members
    ///@{

    static constexpr unsigned int Dim = TDim;

    static constexpr unsigned int NumNodes = TNumNodes;

    static constexpr unsigned int BlockSize = TDim + 1;

    static constexpr unsigned int LocalSize = TNumNodes * BlockSize;

    ScalarDataType Pressure;

    ScalarDataType Density;

    ScalarDataType Viscosity;

    ScalarDataType MassProjection;

    VectorDataType Velocity;

    VectorDataType MeshVelocity;

    VectorDataType BodyForce;

    VectorDataType MomentumProjection;

    ///@}
    ///@name Public operations
    ///@{

    static int Check(Element& rElement);

    ///@}

private:

    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FluidElementData& operator=(FluidElementData const& rOther);

    /// Copy constructor.
    FluidElementData(FluidElementData const& rOther);

    ///@}

}; // Class FluidElementData

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_FLUID_ELEMENT_DATA_H_INCLUDED  defined
