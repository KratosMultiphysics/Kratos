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

#if !defined(KRATOS_INTEGRATION_POINT_DATA_H_INCLUDED)
#define KRATOS_INTEGRATION_POINT_DATA_H_INCLUDED

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

/// Auxiliary class to evaluate and hold data at integration points for elements
/// based on FluidElement
template <class TElementData>
class IntegrationPointData
{
public:
    ///@name Life Cycle
    ///@{

    IntegrationPointData();

    ///@}
    ///@name Public members
    ///@{

    unsigned int IntegrationPointIndex;

    double Weight;

    array_1d<double, TElementData::NumNodes> N;

    boost::numeric::ublas::bounded_matrix<double, TElementData::NumNodes, TElementData::Dim> DN_DX;

    double Pressure;

    double Density;

    double Viscosity;

    double MassProjection;

    array_1d<double, TElementData::Dim> Velocity;

    array_1d<double, TElementData::Dim> MeshVelocity;

    array_1d<double, TElementData::Dim> ConvectiveVelocity;

    array_1d<double, TElementData::Dim> BodyForce;

    array_1d<double, TElementData::Dim> MomentumProjection;

    ///@}
    ///@name Public Operations
    ///@{

    static void FillIntegrationPointData(IntegrationPointData& rIntegrationPointData,
                                         const TElementData& rElementData,
                                         int IntegrationPointIndex,
                                         const Vector& rWeights,
                                         const Matrix& rNContainer,
                                         const Geometry< Node<3> >::ShapeFunctionsGradientsType& rDN_DX);

    void EvaluateInPoint(double& rResult,
                         const typename TElementData::ScalarDataType& rNodalValues);

    void EvaluateInPoint(array_1d<double, TElementData::Dim>& rResult,
                         const typename TElementData::VectorDataType& rNodalValues);

    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    IntegrationPointData& operator=(IntegrationPointData const& rOther);

    /// Copy constructor.
    IntegrationPointData(IntegrationPointData const& rOther);

    ///@}

}; // struct IntegrationPointData

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_INTEGRATION_POINT_DATA_H_INCLUDED  defined
