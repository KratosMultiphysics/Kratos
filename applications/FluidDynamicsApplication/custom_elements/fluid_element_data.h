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
#include "fluid_element_data.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

/// Auxiliary class to hold data for elements based on FluidElement
/** TODO: see how this can work for a generic number of stored arguments.
 */
template <unsigned int TDim, unsigned int TNumNodes>
struct FluidElementData
{
public:
    ///@name Type Definitions
    ///@{

    typedef array_1d<double, TDim> ScalarDataType;

    typedef boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> VectorDataType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor, accepting the geometry of the element.
    FluidElementData(Geometry<Node<3>>& rGeom);

    ///@}
    ///@name Public members
    ///@{

    constexpr static unsigned int Dim = TDim;

    constexpr static unsigned int NumNodes = TNumNodes;

    ScalarDataType Pressure;

    ScalarDataType Density;

    ScalarDataType Viscosity;

    VectorDataType Velocity;

    VectorDataType MeshVelocity;

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
