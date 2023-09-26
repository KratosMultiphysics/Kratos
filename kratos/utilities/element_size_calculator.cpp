//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//                   Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// Project includes
#include "integration_utilities.h"

// Include base h
#include "element_size_calculator.h"

namespace Kratos {

template< std::size_t TDim, std::size_t TNumNodes >
ElementSizeCalculator<TDim,TNumNodes>::~ElementSizeCalculator() {};

// Triangle2D3 version.
template<>
double ElementSizeCalculator<2,3>::MinimumElementSize(const Geometry<Node >& rGeometry)
{

    /* Calculate node-edge distances */
    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();

    // node 0, edge 12
    double nx = -(y20-y10);
    double ny = x20-x10;
    double Hsq = x10*nx + y10*ny;
    Hsq *= Hsq / (nx*nx + ny*ny);

    // node 1, edge 20
    nx = -y20;
    ny = x20;
    double hsq = x10*nx + y10*ny;
    hsq *= hsq / (nx*nx + ny*ny);
    Hsq = ( hsq < Hsq ) ? hsq : Hsq;

    // node 2, edge 10
    nx = -y10;
    ny = x10;
    hsq = x20*nx + y20*ny;
    hsq *= hsq / (nx*nx + ny*ny);
    Hsq = ( hsq < Hsq ) ? hsq : Hsq;

    return std::sqrt(Hsq);
}
// Triangle2D3 version.
template<>
double ElementSizeCalculator<2,3>::MinimumElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DerivativeNodeIndex > 2)
        << "Invalid DerivativeNodeIndex [ DerivativeNodeIndex = " << DerivativeNodeIndex
        << ", expected DerivativeNodeIndex < 3 ].\n";

    KRATOS_DEBUG_ERROR_IF(DerivativeDirectionIndex > 1)
        << "Invalid DerivativeDirectionIndex [ DerivativeDirectionIndex = " << DerivativeDirectionIndex
        << ", expected DerivativeDirectionIndex < 2 ].\n";

    /* Calculate node-edge distances */
    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double x10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 0);

    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double y10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 1);

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double x20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 0);

    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    const double y20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 1);

    // node 0, edge 12
    double nx = -(y20-y10);
    double nx_derivative = -(y20_derivative - y10_derivative);
    double ny = x20-x10;
    double ny_derivative = x20_derivative - x10_derivative;
    double Hsq = x10*nx + y10*ny;
    Hsq *= Hsq / (nx*nx + ny*ny);
    double Hsq_derivative = HsqDerivative2D(x10, x10_derivative, nx, nx_derivative, y10, y10_derivative, ny, ny_derivative);

    // node 1, edge 20
    nx = -y20;
    nx_derivative = -y20_derivative;
    ny = x20;
    ny_derivative = x20_derivative;
    double hsq = x10*nx + y10*ny;
    hsq *= hsq / (nx*nx + ny*ny);
    double hsq_derivative = HsqDerivative2D(x10, x10_derivative, nx, nx_derivative, y10, y10_derivative, ny, ny_derivative);
    Hsq_derivative = (hsq < Hsq) ? hsq_derivative : Hsq_derivative;
    Hsq = ( hsq < Hsq ) ? hsq : Hsq;

    // node 2, edge 10
    nx = -y10;
    nx_derivative = -y10_derivative;
    ny = x10;
    ny_derivative = x10_derivative;
    hsq = x20*nx + y20*ny;
    hsq *= hsq / (nx*nx + ny*ny);
    hsq_derivative = HsqDerivative2D(x20, x20_derivative, nx, nx_derivative, y20, y20_derivative, ny, ny_derivative);
    Hsq_derivative = (hsq < Hsq) ? hsq_derivative : Hsq_derivative;
    Hsq = ( hsq < Hsq ) ? hsq : Hsq;

    return 0.5 * Hsq_derivative / std::sqrt(Hsq);

    KRATOS_CATCH("");
}

// Triangle2D6 version.
template<>
double ElementSizeCalculator<2,6>::MinimumElementSize(const Geometry<Node >& rGeometry)
{

    const double minimum_element_size = ElementSizeCalculator<2,3>::MinimumElementSize(rGeometry);


    return minimum_element_size;
}

// Triangle2D6 version.
template<>
double ElementSizeCalculator<2,6>::MinimumElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{

    double element_size_derivative = 0.0;

    if (DerivativeNodeIndex < 3)
        element_size_derivative = ElementSizeCalculator<2,3>::MinimumElementSizeDerivative(DerivativeNodeIndex,DerivativeDirectionIndex,rGeometry);

    return element_size_derivative;

}

// Quadrilateral2D4 version.
template<>
double ElementSizeCalculator<2,4>::MinimumElementSize(const Geometry<Node >& rGeometry)
{

    const Node& r_node_0 = rGeometry[0];
    const Node& r_node_1 = rGeometry[1];
    const Node& r_node_2 = rGeometry[2];
    const Node& r_node_3 = rGeometry[3];

    // Calculate face centers
    const double x10 = (r_node_1.X() + r_node_0.X())/2.;
    const double y10 = (r_node_1.Y() + r_node_0.Y())/2.;

    const double x21 = (r_node_2.X() + r_node_1.X())/2.;
    const double y21 = (r_node_2.Y() + r_node_1.Y())/2.;

    const double x32 = (r_node_3.X() + r_node_2.X())/2.;
    const double y32 = (r_node_3.Y() + r_node_2.Y())/2.;

    const double x03 = (r_node_0.X() + r_node_3.X())/2.;
    const double y03 = (r_node_0.Y() + r_node_3.Y())/2.;

    // Distance between face centers (xi direction)
    const double dxi_x = x21 - x03;
    const double dxi_y = y21 - y03;

    const double h2_xi = dxi_x*dxi_x + dxi_y*dxi_y;

    // Distance between face centers (eta direction)
    const double deta_x = x32 - x10;
    const double deta_y = y32 - y10;

    const double h2_eta = deta_x*deta_x + deta_y*deta_y;

    const double h2 = h2_xi < h2_eta ? h2_xi : h2_eta;

    return std::sqrt(h2);
}

template<>
double ElementSizeCalculator<2,4>::MinimumElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DerivativeNodeIndex > 3)
        << "Invalid DerivativeNodeIndex [ DerivativeNodeIndex = " << DerivativeNodeIndex
        << ", expected DerivativeNodeIndex < 4 ].\n";

    KRATOS_DEBUG_ERROR_IF(DerivativeDirectionIndex > 1)
        << "Invalid DerivativeDirectionIndex [ DerivativeDirectionIndex = " << DerivativeDirectionIndex
        << ", expected DerivativeDirectionIndex < 2 ].\n";

    const Node& r_node_0 = rGeometry[0];
    const Node& r_node_1 = rGeometry[1];
    const Node& r_node_2 = rGeometry[2];
    const Node& r_node_3 = rGeometry[3];

    // Calculate face centers
    const double x10 = (r_node_1.X() + r_node_0.X())/2.;
    const double x10_derivative = ((DerivativeNodeIndex == 1) + (DerivativeNodeIndex == 0)) * (DerivativeDirectionIndex == 0) / 2.;

    const double y10 = (r_node_1.Y() + r_node_0.Y())/2.;
    const double y10_derivative = ((DerivativeNodeIndex == 1) + (DerivativeNodeIndex == 0)) * (DerivativeDirectionIndex == 1) / 2.;

    const double x21 = (r_node_2.X() + r_node_1.X())/2.;
    const double x21_derivative = ((DerivativeNodeIndex == 2) + (DerivativeNodeIndex == 1)) * (DerivativeDirectionIndex == 0) / 2.;

    const double y21 = (r_node_2.Y() + r_node_1.Y())/2.;
    const double y21_derivative = ((DerivativeNodeIndex == 2) + (DerivativeNodeIndex == 1)) * (DerivativeDirectionIndex == 1) / 2.;

    const double x32 = (r_node_3.X() + r_node_2.X())/2.;
    const double x32_derivative = ((DerivativeNodeIndex == 3) + (DerivativeNodeIndex == 2)) * (DerivativeDirectionIndex == 0) / 2.;

    const double y32 = (r_node_3.Y() + r_node_2.Y())/2.;
    const double y32_derivative = ((DerivativeNodeIndex == 3) + (DerivativeNodeIndex == 2)) * (DerivativeDirectionIndex == 1) / 2.;

    const double x03 = (r_node_0.X() + r_node_3.X())/2.;
    const double x03_derivative = ((DerivativeNodeIndex == 0) + (DerivativeNodeIndex == 3)) * (DerivativeDirectionIndex == 0) / 2.;

    const double y03 = (r_node_0.Y() + r_node_3.Y())/2.;
    const double y03_derivative = ((DerivativeNodeIndex == 0) + (DerivativeNodeIndex == 3)) * (DerivativeDirectionIndex == 1) / 2.;

    // Distance between face centers (xi direction)
    const double dxi_x = x21 - x03;
    const double dxi_x_derivative = x21_derivative - x03_derivative;

    const double dxi_y = y21 - y03;
    const double dxi_y_derivative = y21_derivative - y03_derivative;

    const double h2_xi = dxi_x*dxi_x + dxi_y*dxi_y;
    const double h2_xi_derivative = 2 * dxi_x * dxi_x_derivative + 2 * dxi_y * dxi_y_derivative;

    // Distance between face centers (eta direction)
    const double deta_x = x32 - x10;
    const double deta_x_derivative = x32_derivative - x10_derivative;

    const double deta_y = y32 - y10;
    const double deta_y_derivative = y32_derivative - y10_derivative;

    const double h2_eta = deta_x*deta_x + deta_y*deta_y;
    const double h2_eta_derivative = 2 * deta_x * deta_x_derivative + 2 * deta_y * deta_y_derivative;

    const double h2 = h2_xi < h2_eta ? h2_xi : h2_eta;
    double h2_derivative = h2_xi < h2_eta ? h2_xi_derivative : h2_eta_derivative;
    if (std::abs(h2_xi - h2_eta) < 1.0e-12)
        h2_derivative = std::min(h2_xi_derivative,h2_eta_derivative);

    return 0.5 * h2_derivative / std::sqrt(h2);

    KRATOS_CATCH("");
}

// Quadrilateral2D9 version.
template<>
double ElementSizeCalculator<2,9>::MinimumElementSize(const Geometry<Node >& rGeometry)
{

    const double minimum_element_size = ElementSizeCalculator<2,4>::MinimumElementSize(rGeometry);

    return minimum_element_size;
}

template<>
double ElementSizeCalculator<2,9>::MinimumElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    double element_size_derivative = 0.0;

    if (DerivativeNodeIndex < 4)
        element_size_derivative = ElementSizeCalculator<2,4>::MinimumElementSizeDerivative(DerivativeNodeIndex,DerivativeDirectionIndex,rGeometry);

    return element_size_derivative;
}

// Tetrahedra3D4 version.
template<>
double ElementSizeCalculator<3,4>::MinimumElementSize(const Geometry<Node >& rGeometry)
{

    /* Calculate distances between each node and the opposite face */
    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double z10 = rGeometry[1].Z() - rGeometry[0].Z();

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    const double z20 = rGeometry[2].Z() - rGeometry[0].Z();

    const double x30 = rGeometry[3].X() - rGeometry[0].X();
    const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    const double z30 = rGeometry[3].Z() - rGeometry[0].Z();

    // face 123
    double nx = (y30-y10)*(z20-z10) - (z30-z10)*(y20-y10);
    double ny = (z30-z10)*(x20-x10) - (x30-x10)*(z20-z10);
    double nz = (x30-x10)*(y20-y10) - (y30-y10)*(x20-x10);
    double Hsq = x10*nx + y10*ny + z10*nz; // scalar product x10*n
    Hsq *= Hsq / (nx*nx + ny*ny + nz*nz); // H^2 = (x10*n)^2 / ||n||^2

    // face 230
    nx = y30*z20 - z30*y20;
    ny = z30*x20 - x30*z20;
    nz = x30*y20 - y30*x20;
    double hsq = x10*nx + y10*ny + z10*nz;
    hsq *= hsq / (nx*nx + ny*ny + nz*nz);
    Hsq = (hsq < Hsq) ? hsq : Hsq;

    // face 301
    nx = y10*z30 - z10*y30;
    ny = z10*x30 - x10*z30;
    nz = x10*y30 - y10*x30;
    hsq = x20*nx + y20*ny + z20*nz;
    hsq *= hsq / (nx*nx + ny*ny + nz*nz);
    Hsq = (hsq < Hsq) ? hsq : Hsq;

    // face 012
    nx = y10*z20 - z10*y20;
    ny = z10*x20 - x10*z20;
    nz = x10*y20 - y10*x20;
    hsq = x30*nx + y30*ny + z30*nz;
    hsq *= hsq / (nx*nx + ny*ny + nz*nz);
    Hsq = (hsq < Hsq) ? hsq : Hsq;
    return std::sqrt(Hsq);
}

template<>
double ElementSizeCalculator<3,4>::MinimumElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DerivativeNodeIndex > 3)
        << "Invalid DerivativeNodeIndex [ DerivativeNodeIndex = " << DerivativeNodeIndex
        << ", expected DerivativeNodeIndex < 4 ].\n";

    KRATOS_DEBUG_ERROR_IF(DerivativeDirectionIndex > 2)
        << "Invalid DerivativeDirectionIndex [ DerivativeDirectionIndex = " << DerivativeDirectionIndex
        << ", expected DerivativeDirectionIndex < 3 ].\n";

    /* Calculate distances between each node and the opposite face */
    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double x10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 0);

    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double y10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 1);

    const double z10 = rGeometry[1].Z() - rGeometry[0].Z();
    const double z10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 2);

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double x20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 0);

    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    const double y20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 1);

    const double z20 = rGeometry[2].Z() - rGeometry[0].Z();
    const double z20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 2);

    const double x30 = rGeometry[3].X() - rGeometry[0].X();
    const double x30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 0);

    const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    const double y30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 1);

    const double z30 = rGeometry[3].Z() - rGeometry[0].Z();
    const double z30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 2);

    // face 123
    double nx = (y30-y10)*(z20-z10) - (z30-z10)*(y20-y10);
    double nx_derivative = (y30_derivative-y10_derivative)*(z20-z10)+(y30-y10)*(z20_derivative-z10_derivative);
    nx_derivative -= ((z30_derivative-z10_derivative)*(y20-y10)+(z30-z10)*(y20_derivative-y10_derivative));
    double ny = (z30-z10)*(x20-x10) - (x30-x10)*(z20-z10);
    double ny_derivative = (z30_derivative-z10_derivative)*(x20-x10) + (z30-z10)*(x20_derivative-x10_derivative);
    ny_derivative -= ((x30_derivative-x10_derivative)*(z20-z10) + (x30-x10)*(z20_derivative-z10_derivative));
    double nz = (x30-x10)*(y20-y10) - (y30-y10)*(x20-x10);
    double nz_derivative = (x30_derivative-x10_derivative)*(y20-y10) + (x30-x10)*(y20_derivative-y10_derivative);
    nz_derivative -= ((y30_derivative-y10_derivative)*(x20-x10) + (y30-y10)*(x20_derivative-x10_derivative));
    double Hsq = x10*nx + y10*ny + z10*nz; // scalar product x10*n
    Hsq *= Hsq / (nx*nx + ny*ny + nz*nz); // H^2 = (x10*n)^2 / ||n||^2
    double Hsq_derivative = HsqDerivative3D(x10, x10_derivative, nx, nx_derivative, y10, y10_derivative, ny, ny_derivative, z10, z10_derivative, nz, nz_derivative);

    // face 230
    nx = y30*z20 - z30*y20;
    nx_derivative = y30_derivative*z20+y30*z20_derivative-z30_derivative*y20-z30*y20_derivative;
    ny = z30*x20 - x30*z20;
    ny_derivative = z30_derivative*x20+z30*x20_derivative-x30_derivative*z20-x30*z20_derivative;
    nz = x30*y20 - y30*x20;
    nz_derivative = x30_derivative*y20+x30*y20_derivative-y30_derivative*x20-y30*x20_derivative;
    double hsq = x10*nx + y10*ny + z10*nz;
    hsq *= hsq / (nx*nx + ny*ny + nz*nz);
    double hsq_derivative = HsqDerivative3D(x10, x10_derivative, nx, nx_derivative, y10, y10_derivative, ny, ny_derivative, z10, z10_derivative, nz, nz_derivative);
    Hsq_derivative = (hsq < Hsq) ? hsq_derivative : Hsq_derivative;
    Hsq = (hsq < Hsq) ? hsq : Hsq;

    // face 301
    nx = y10*z30 - z10*y30;
    nx_derivative = y10_derivative*z30+y10*z30_derivative-z10_derivative*y30-z10*y30_derivative;
    ny = z10*x30 - x10*z30;
    ny_derivative = z10_derivative*x30+z10*x30_derivative-x10_derivative*z30-x10*z30_derivative;
    nz = x10*y30 - y10*x30;
    nz_derivative = x10_derivative*y30+x10*y30_derivative-y10_derivative*x30-y10*x30_derivative;
    hsq = x20*nx + y20*ny + z20*nz;
    hsq *= hsq / (nx*nx + ny*ny + nz*nz);
    hsq_derivative = HsqDerivative3D(x20, x20_derivative, nx, nx_derivative, y20, y20_derivative, ny, ny_derivative, z20, z20_derivative, nz, nz_derivative);
    Hsq_derivative = (hsq < Hsq) ? hsq_derivative : Hsq_derivative;
    Hsq = (hsq < Hsq) ? hsq : Hsq;

    // face 012
    nx = y10*z20 - z10*y20;
    nx_derivative = y10_derivative*z20+y10*z20_derivative-z10_derivative*y20-z10*y20_derivative;
    ny = z10*x20 - x10*z20;
    ny_derivative = z10_derivative*x20+z10*x20_derivative-x10_derivative*z20-x10*z20_derivative;
    nz = x10*y20 - y10*x20;
    nz_derivative = x10_derivative*y20+x10*y20_derivative-y10_derivative*x20-y10*x20_derivative;
    hsq = x30*nx + y30*ny + z30*nz;
    hsq *= hsq / (nx*nx + ny*ny + nz*nz);
    hsq_derivative = HsqDerivative3D(x30, x30_derivative, nx, nx_derivative, y30, y30_derivative, ny, ny_derivative, z30, z30_derivative, nz, nz_derivative);
    Hsq_derivative = (hsq < Hsq) ? hsq_derivative : Hsq_derivative;
    Hsq = (hsq < Hsq) ? hsq : Hsq;

    return 0.5 * Hsq_derivative / std::sqrt(Hsq);

    KRATOS_CATCH("");
}

// Tetrahedra3D4 version.
template<>
double ElementSizeCalculator<3,10>::MinimumElementSize(const Geometry<Node >& rGeometry)
{

    const double minimum_element_size = ElementSizeCalculator<3,4>::MinimumElementSize(rGeometry);

    return minimum_element_size;
}

template<>
double ElementSizeCalculator<3,10>::MinimumElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    double element_size_derivative = 0.0;

    if (DerivativeNodeIndex < 4)
        element_size_derivative = ElementSizeCalculator<3,4>::MinimumElementSizeDerivative(DerivativeNodeIndex,DerivativeDirectionIndex,rGeometry);

    return element_size_derivative;

    KRATOS_CATCH("");
}

// Prism3D6 version.
template<>
double ElementSizeCalculator<3,6>::MinimumElementSize(const Geometry<Node> &rGeometry)
{
    // Get nodes
    const Node& r_node_0 = rGeometry[0];
    const Node& r_node_1 = rGeometry[1];
    const Node& r_node_2 = rGeometry[2];
    const Node& r_node_3 = rGeometry[3];
    const Node& r_node_4 = rGeometry[4];
    const Node& r_node_5 = rGeometry[5];

    // Calculate face centers (top and bottom face)
    const double one_third = 1.0/3.0;
    const array_1d<double,3> low_dseta = one_third * (r_node_0.Coordinates() + r_node_1.Coordinates() + r_node_2.Coordinates() );
    const array_1d<double,3> high_dseta = one_third * (r_node_3.Coordinates() + r_node_4.Coordinates() + r_node_5.Coordinates() );

    // Calculate node-edge distances (top face)
    const double x10 = r_node_1.X() - r_node_0.X();
    const double y10 = r_node_1.Y() - r_node_0.Y();
    const double x20 = r_node_2.X() - r_node_0.X();
    const double y20 = r_node_2.Y() - r_node_0.Y();

    // node 0, edge 12
    double nx = -(y20-y10);
    double ny = x20-x10;
    double Hsq = x10*nx + y10*ny;
    Hsq *= Hsq / (nx*nx + ny*ny);

    // node 1, edge 20
    nx = -y20;
    ny = x20;
    double hsq = x10*nx + y10*ny;
    hsq *= hsq / (nx*nx + ny*ny);
    Hsq = ( hsq < Hsq ) ? hsq : Hsq;

    // node 2, edge 10
    nx = -y10;
    ny = x10;
    hsq = x20*nx + y20*ny;
    hsq *= hsq / (nx*nx + ny*ny);
    Hsq = ( hsq < Hsq ) ? hsq : Hsq;

    // Calculate distance between the face centers (dseta direction)
    const array_1d<double,3> d_dseta = high_dseta - low_dseta;
    const double h2_dseta = d_dseta[0]*d_dseta[0] + d_dseta[1]*d_dseta[1] + d_dseta[2]*d_dseta[2];

    const double h2 = h2_dseta < Hsq?h2_dseta:Hsq;

    return std::sqrt(h2);
}

template<>
double ElementSizeCalculator<3,6>::MinimumElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DerivativeNodeIndex > 5)
        << "Invalid DerivativeNodeIndex [ DerivativeNodeIndex = " << DerivativeNodeIndex
        << ", expected DerivativeNodeIndex < 6 ].\n";

    KRATOS_DEBUG_ERROR_IF(DerivativeDirectionIndex > 2)
        << "Invalid DerivativeDirectionIndex [ DerivativeDirectionIndex = " << DerivativeDirectionIndex
        << ", expected DerivativeDirectionIndex < 3 ].\n";

    // Get nodes
    const Node& r_node_0 = rGeometry[0];
    const Node& r_node_1 = rGeometry[1];
    const Node& r_node_2 = rGeometry[2];
    const Node& r_node_3 = rGeometry[3];
    const Node& r_node_4 = rGeometry[4];
    const Node& r_node_5 = rGeometry[5];

    // Calculate face centers (top and bottom face)
    const double one_third = 1.0/3.0;
    const array_1d<double,3> low_dseta = one_third * (r_node_0.Coordinates() + r_node_1.Coordinates() + r_node_2.Coordinates() );
    array_1d<double,3> low_dseta_derivative = ZeroVector(3);
    low_dseta_derivative[DerivativeDirectionIndex] = one_third * (DerivativeNodeIndex < 3);
    const array_1d<double,3> high_dseta = one_third * (r_node_3.Coordinates() + r_node_4.Coordinates() + r_node_5.Coordinates() );
    array_1d<double,3> high_dseta_derivative = ZeroVector(3);
    high_dseta_derivative[DerivativeDirectionIndex] = one_third * (DerivativeNodeIndex > 2);

    // Calculate node-edge distances (top face)
    const double x10 = r_node_1.X() - r_node_0.X();
    const double x10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 0);
    const double y10 = r_node_1.Y() - r_node_0.Y();
    const double y10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 1);
    const double x20 = r_node_2.X() - r_node_0.X();
    const double x20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 0);
    const double y20 = r_node_2.Y() - r_node_0.Y();
    const double y20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 1);

    // node 0, edge 12
    double nx = -(y20-y10);
    double nx_derivative = -(y20_derivative-y10_derivative);
    double ny = x20-x10;
    double ny_derivative = x20_derivative-x10_derivative;
    double Hsq = x10*nx + y10*ny;
    Hsq *= Hsq / (nx*nx + ny*ny);
    double Hsq_derivative = HsqDerivative2D(x10, x10_derivative, nx, nx_derivative, y10, y10_derivative, ny, ny_derivative);

    // node 1, edge 20
    nx = -y20;
    nx_derivative = -y20_derivative;
    ny = x20;
    ny_derivative = x20_derivative;
    double hsq = x10*nx + y10*ny;
    hsq *= hsq / (nx*nx + ny*ny);
    double hsq_derivative = HsqDerivative2D(x10, x10_derivative, nx, nx_derivative, y10, y10_derivative, ny, ny_derivative);
    Hsq_derivative = (hsq < Hsq) ? hsq_derivative : Hsq_derivative;
    Hsq = ( hsq < Hsq ) ? hsq : Hsq;

    // node 2, edge 10
    nx = -y10;
    nx_derivative = -y10_derivative;
    ny = x10;
    ny_derivative = x10_derivative;
    hsq = x20*nx + y20*ny;
    hsq *= hsq / (nx*nx + ny*ny);
    hsq_derivative = HsqDerivative2D(x20, x20_derivative, nx, nx_derivative, y20, y20_derivative, ny, ny_derivative);
    Hsq_derivative = ( hsq < Hsq ) ? hsq_derivative : Hsq_derivative;
    Hsq = ( hsq < Hsq ) ? hsq : Hsq;

    // Calculate distance between the face centers (dseta direction)
    const array_1d<double,3> d_dseta = high_dseta - low_dseta;
    const array_1d<double,3> d_dseta_derivative = high_dseta_derivative - low_dseta_derivative;
    const double h2_dseta = d_dseta[0]*d_dseta[0] + d_dseta[1]*d_dseta[1] + d_dseta[2]*d_dseta[2];
    const double h2_dseta_derivative = NormSquareDerivative(d_dseta, d_dseta_derivative);

    const double h2 = h2_dseta < Hsq?h2_dseta:Hsq;
    const double h2_derivative = (h2_dseta < Hsq) ? h2_dseta_derivative : Hsq_derivative;

    return 0.5 * h2_derivative / std::sqrt(h2);

    KRATOS_CATCH("");
}

// Hexahedra3D8 version. We use the distance between face centers to compute lengths.
template<>
double ElementSizeCalculator<3,8>::MinimumElementSize(const Geometry<Node >& rGeometry)
{

    const Node& r_node_0 = rGeometry[0];
    const Node& r_node_1 = rGeometry[1];
    const Node& r_node_2 = rGeometry[2];
    const Node& r_node_3 = rGeometry[3];
    const Node& r_node_4 = rGeometry[4];
    const Node& r_node_5 = rGeometry[5];
    const Node& r_node_6 = rGeometry[6];
    const Node& r_node_7 = rGeometry[7];

    // Calculate face centers
    const array_1d<double,3> low_xi  = 0.25 * (r_node_0.Coordinates() + r_node_4.Coordinates() + r_node_3.Coordinates() + r_node_7.Coordinates());
    const array_1d<double,3> high_xi = 0.25 * (r_node_1.Coordinates() + r_node_2.Coordinates() + r_node_6.Coordinates() + r_node_5.Coordinates());
    const array_1d<double,3> low_eta  = 0.25 * (r_node_0.Coordinates() + r_node_1.Coordinates() + r_node_5.Coordinates() + r_node_4.Coordinates());
    const array_1d<double,3> high_eta = 0.25 * (r_node_3.Coordinates() + r_node_7.Coordinates() + r_node_6.Coordinates() + r_node_2.Coordinates());
    const array_1d<double,3> low_dseta  = 0.25 * (r_node_0.Coordinates() + r_node_3.Coordinates() + r_node_2.Coordinates() + r_node_1.Coordinates());
    const array_1d<double,3> high_dseta = 0.25 * (r_node_4.Coordinates() + r_node_5.Coordinates() + r_node_6.Coordinates() + r_node_7.Coordinates());

    // Distance between face centers (xi direction)
    const array_1d<double,3> d_xi = high_xi - low_xi;
    const double h2_xi = d_xi[0]*d_xi[0] + d_xi[1]*d_xi[1] + d_xi[2]*d_xi[2];

    // Distance between face centers (eta direction)
    const array_1d<double,3> d_eta = high_eta - low_eta;
    const double h2_eta = d_eta[0]*d_eta[0] + d_eta[1]*d_eta[1] + d_eta[2]*d_eta[2];

    // Distance between face centers (dseta direction)
    const array_1d<double,3> d_dseta = high_dseta - low_dseta;
    const double h2_dseta = d_dseta[0]*d_dseta[0] + d_dseta[1]*d_dseta[1] + d_dseta[2]*d_dseta[2];

    double h2 = h2_xi < h2_eta ? h2_xi : h2_eta;
    h2 = h2 < h2_dseta ? h2 : h2_dseta;

    return std::sqrt(h2);
}

template<>
double ElementSizeCalculator<3,8>::MinimumElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DerivativeNodeIndex > 7)
        << "Invalid DerivativeNodeIndex [ DerivativeNodeIndex = " << DerivativeNodeIndex
        << ", expected DerivativeNodeIndex < 8 ].\n";

    KRATOS_DEBUG_ERROR_IF(DerivativeDirectionIndex > 2)
        << "Invalid DerivativeDirectionIndex [ DerivativeDirectionIndex = " << DerivativeDirectionIndex
        << ", expected DerivativeDirectionIndex < 3 ].\n";

    const Node& r_node_0 = rGeometry[0];
    const Node& r_node_1 = rGeometry[1];
    const Node& r_node_2 = rGeometry[2];
    const Node& r_node_3 = rGeometry[3];
    const Node& r_node_4 = rGeometry[4];
    const Node& r_node_5 = rGeometry[5];
    const Node& r_node_6 = rGeometry[6];
    const Node& r_node_7 = rGeometry[7];

    const auto derivative_output = [&](const unsigned int A, const unsigned int B, const unsigned int C, const unsigned int D) -> array_1d<double, 3> {
        array_1d<double, 3> result = ZeroVector(3);
        result[DerivativeDirectionIndex] = 0.25 * ((DerivativeNodeIndex == A) || (DerivativeNodeIndex == B) || (DerivativeNodeIndex == C) || (DerivativeNodeIndex == D));
        return result;
    };

    // Calculate face centers
    const array_1d<double,3> low_xi  = 0.25 * (r_node_0.Coordinates() + r_node_4.Coordinates() + r_node_3.Coordinates() + r_node_7.Coordinates());
    const array_1d<double,3> low_xi_derivative = derivative_output(0, 4, 3, 7);

    const array_1d<double,3> high_xi = 0.25 * (r_node_1.Coordinates() + r_node_2.Coordinates() + r_node_6.Coordinates() + r_node_5.Coordinates());
    const array_1d<double,3> high_xi_derivative = derivative_output(1, 2, 6, 5);

    const array_1d<double,3> low_eta  = 0.25 * (r_node_0.Coordinates() + r_node_1.Coordinates() + r_node_5.Coordinates() + r_node_4.Coordinates());
    const array_1d<double,3> low_eta_derivative = derivative_output(0, 1, 5, 4);

    const array_1d<double,3> high_eta = 0.25 * (r_node_3.Coordinates() + r_node_7.Coordinates() + r_node_6.Coordinates() + r_node_2.Coordinates());
    const array_1d<double,3> high_eta_derivative = derivative_output(3, 7, 6, 2);

    const array_1d<double,3> low_dseta  = 0.25 * (r_node_0.Coordinates() + r_node_3.Coordinates() + r_node_2.Coordinates() + r_node_1.Coordinates());
    const array_1d<double,3> low_dseta_derivative = derivative_output(0, 3, 2, 1);

    const array_1d<double,3> high_dseta = 0.25 * (r_node_4.Coordinates() + r_node_5.Coordinates() + r_node_6.Coordinates() + r_node_7.Coordinates());
    const array_1d<double,3> high_dseta_derivative = derivative_output(4, 5, 6, 7);

    // Distance between face centers (xi direction)
    const array_1d<double,3> d_xi = high_xi - low_xi;
    const array_1d<double,3> d_xi_derivative = high_xi_derivative - low_xi_derivative;
    const double h2_xi = d_xi[0]*d_xi[0] + d_xi[1]*d_xi[1] + d_xi[2]*d_xi[2];
    const double h2_xi_derivative = NormSquareDerivative(d_xi, d_xi_derivative);

    // Distance between face centers (eta direction)
    const array_1d<double,3> d_eta = high_eta - low_eta;
    const array_1d<double,3> d_eta_derivative = high_eta_derivative - low_eta_derivative;
    const double h2_eta = d_eta[0]*d_eta[0] + d_eta[1]*d_eta[1] + d_eta[2]*d_eta[2];
    const double h2_eta_derivative = NormSquareDerivative(d_eta, d_eta_derivative);

    // Distance between face centers (dseta direction)
    const array_1d<double,3> d_dseta = high_dseta - low_dseta;
    const array_1d<double,3> d_dseta_derivative = high_dseta_derivative - low_dseta_derivative;
    const double h2_dseta = d_dseta[0]*d_dseta[0] + d_dseta[1]*d_dseta[1] + d_dseta[2]*d_dseta[2];
    const double h2_dseta_derivative = NormSquareDerivative(d_dseta, d_dseta_derivative);

    double h2 = h2_xi < h2_eta ? h2_xi : h2_eta;
    double h2_derivative = h2_xi < h2_eta ? h2_xi_derivative : h2_eta_derivative;
    h2_derivative = h2 < h2_dseta ? h2_derivative : h2_dseta_derivative;
    h2 = h2 < h2_dseta ? h2 : h2_dseta;

    return 0.5 * h2_derivative / std::sqrt(h2);

    KRATOS_CATCH("");
}

// Hexahedra3D27 version. We use the distance between face centers to compute lengths.
template<>
double ElementSizeCalculator<3,27>::MinimumElementSize(const Geometry<Node >& rGeometry)
{

    const double minimum_element_size = ElementSizeCalculator<3,8>::MinimumElementSize(rGeometry);

    return minimum_element_size;
}

template<>
double ElementSizeCalculator<3,27>::MinimumElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    double element_size_derivative = 0.0;

    if (DerivativeNodeIndex < 8)
        element_size_derivative = ElementSizeCalculator<3,8>::MinimumElementSizeDerivative(DerivativeNodeIndex,DerivativeDirectionIndex,rGeometry);

    return element_size_derivative;

    KRATOS_CATCH("");
}

// Triangle2D3 version.
template<>
double ElementSizeCalculator<2,3>::AverageElementSize(const Geometry<Node >& rGeometry)
{

    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();

    return std::sqrt(0.5 * (x10*y20-x20*y10) );
}

template<>
double ElementSizeCalculator<2,3>::AverageElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DerivativeNodeIndex > 2)
        << "Invalid DerivativeNodeIndex [ DerivativeNodeIndex = " << DerivativeNodeIndex
        << ", expected DerivativeNodeIndex < 3 ].\n";

    KRATOS_DEBUG_ERROR_IF(DerivativeDirectionIndex > 1)
        << "Invalid DerivativeDirectionIndex [ DerivativeDirectionIndex = " << DerivativeDirectionIndex
        << ", expected DerivativeDirectionIndex < 2 ].\n";

    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double x10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 0);

    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double y10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 1);

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double x20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 0);

    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    const double y20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 1);

    return 0.5 * 0.5 * (x10_derivative*y20+x10*y20_derivative-x20_derivative*y10-x20*y10_derivative) / std::sqrt(0.5 * (x10*y20-x20*y10) );

    KRATOS_CATCH("");
}

// Triangle2D6 version.
template<>
double ElementSizeCalculator<2,6>::AverageElementSize(const Geometry<Node >& rGeometry)
{

    const double average_element_size = ElementSizeCalculator<2,3>::AverageElementSize(rGeometry);

    return average_element_size;
}

template<>
double ElementSizeCalculator<2,6>::AverageElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    double element_size_derivative = 0.0;

    if (DerivativeNodeIndex < 3)
        element_size_derivative = ElementSizeCalculator<2,3>::AverageElementSizeDerivative(DerivativeNodeIndex,DerivativeDirectionIndex,rGeometry);

    return element_size_derivative;

    KRATOS_CATCH("");
}

// Quadrilateral2D4 version.
template<>
double ElementSizeCalculator<2,4>::AverageElementSize(const Geometry<Node >& rGeometry)
{
    return rGeometry.Length();
}

template<>
double ElementSizeCalculator<2,4>::AverageElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    const double area = rGeometry.Area();
    const double area_derivative = IntegrationUtilities::ComputeArea2DGeometryDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, rGeometry);

    if (area > 0.0) {
        return 0.5 * area_derivative / std::sqrt(area);
    } else if (area < 0.0) {
        return -0.5 * area_derivative / std::sqrt(-area);
    } else {
        return 0.0;
    }

    KRATOS_CATCH("");
}

// Quadrilateral2D9 version.
template<>
double ElementSizeCalculator<2,9>::AverageElementSize(const Geometry<Node >& rGeometry)
{
    return rGeometry.Length();
}

template<>
double ElementSizeCalculator<2,9>::AverageElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    Matrix jacobian(2, 2);
    Matrix shape_function_local_gradients(rGeometry.size(), 2), rdNdXDerivative;
    double detJ_derivative;
    ShapeParameter shape_param;

    shape_param.NodeIndex = DerivativeNodeIndex;
    shape_param.Direction = DerivativeDirectionIndex;

    rGeometry.ShapeFunctionsLocalGradients(shape_function_local_gradients, Point());
    rGeometry.Jacobian(jacobian, Point());

    GeometricalSensitivityUtility geometrical_sensitivity_utility(jacobian, shape_function_local_gradients);
    geometrical_sensitivity_utility.CalculateSensitivity(shape_param, detJ_derivative, rdNdXDerivative);

    const auto detJ = MathUtils<double>::Det2(jacobian);

    if (detJ > 0.0) {
        return 0.5 * detJ_derivative / (std::sqrt(detJ));
    } else if (detJ < 0.0) {
        return -0.5 * detJ_derivative / (std::sqrt(-detJ));
    } else {
        return 0.0;
    }
}

// Tetrahedra3D4 version.
template<>
double ElementSizeCalculator<3,4>::AverageElementSize(const Geometry<Node >& rGeometry)
{

    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double z10 = rGeometry[1].Z() - rGeometry[0].Z();

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    const double z20 = rGeometry[2].Z() - rGeometry[0].Z();

    const double x30 = rGeometry[3].X() - rGeometry[0].X();
    const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    const double z30 = rGeometry[3].Z() - rGeometry[0].Z();

    const double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

    return std::pow(detJ/6.0,1./3.);
}

template<>
double ElementSizeCalculator<3,4>::AverageElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DerivativeNodeIndex > 3)
        << "Invalid DerivativeNodeIndex [ DerivativeNodeIndex = " << DerivativeNodeIndex
        << ", expected DerivativeNodeIndex < 4 ].\n";

    KRATOS_DEBUG_ERROR_IF(DerivativeDirectionIndex > 2)
        << "Invalid DerivativeDirectionIndex [ DerivativeDirectionIndex = " << DerivativeDirectionIndex
        << ", expected DerivativeDirectionIndex < 3 ].\n";

    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double x10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 0);

    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double y10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 1);

    const double z10 = rGeometry[1].Z() - rGeometry[0].Z();
    const double z10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 2);

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double x20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 0);

    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    const double y20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 1);

    const double z20 = rGeometry[2].Z() - rGeometry[0].Z();
    const double z20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 2);

    const double x30 = rGeometry[3].X() - rGeometry[0].X();
    const double x30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 0);

    const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    const double y30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 1);

    const double z30 = rGeometry[3].Z() - rGeometry[0].Z();
    const double z30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 2);


    const double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
    double detJ_derivative = 0.0;
    detJ_derivative += x10_derivative * y20 * z30;
    detJ_derivative += x10 * y20_derivative * z30;
    detJ_derivative += x10 * y20 * z30_derivative;
    detJ_derivative -= x10_derivative * y30 * z20;
    detJ_derivative -= x10 * y30_derivative * z20;
    detJ_derivative -= x10 * y30 * z20_derivative;
    detJ_derivative += y10_derivative * z20 * x30;
    detJ_derivative += y10 * z20_derivative * x30;
    detJ_derivative += y10 * z20 * x30_derivative;
    detJ_derivative -= y10_derivative * x20 * z30;
    detJ_derivative -= y10 * x20_derivative * z30;
    detJ_derivative -= y10 * x20 * z30_derivative;
    detJ_derivative += z10_derivative * x20 * y30;
    detJ_derivative += z10 * x20_derivative * y30;
    detJ_derivative += z10 * x20 * y30_derivative;
    detJ_derivative -= z10_derivative * y20 * x30;
    detJ_derivative -= z10 * y20_derivative * x30;
    detJ_derivative -= z10 * y20 * x30_derivative;

    return (1./3.) * (detJ_derivative/6.0) / std::pow(detJ/6.0, 2./3.);

    KRATOS_CATCH("");
}

// Tetrahedra3D10 version.
template<>
double ElementSizeCalculator<3,10>::AverageElementSize(const Geometry<Node >& rGeometry)
{

    const double average_element_size = ElementSizeCalculator<3,4>::AverageElementSize(rGeometry);

    return average_element_size;
}

template<>
double ElementSizeCalculator<3,10>::AverageElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    double element_size_derivative = 0.0;

    if (DerivativeNodeIndex < 4)
        element_size_derivative = ElementSizeCalculator<3,4>::AverageElementSizeDerivative(DerivativeNodeIndex,DerivativeDirectionIndex,rGeometry);

    return element_size_derivative;


    KRATOS_CATCH("");
}

// Prism3D6 version
template <>
double ElementSizeCalculator<3,6>::AverageElementSize(const Geometry<Node> &rGeometry)
{
    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double z10 = rGeometry[1].Z() - rGeometry[0].Z();

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    const double z20 = rGeometry[2].Z() - rGeometry[0].Z();

    const double x30 = rGeometry[3].X() - rGeometry[0].X();
    const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    const double z30 = rGeometry[3].Z() - rGeometry[0].Z();

    const double detJ = 0.5 * (x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30);
    return std::pow(detJ, 1. / 3.);
}

template<>
double ElementSizeCalculator<3,6>::AverageElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DerivativeNodeIndex > 5)
        << "Invalid DerivativeNodeIndex [ DerivativeNodeIndex = " << DerivativeNodeIndex
        << ", expected DerivativeNodeIndex < 6 ].\n";

    KRATOS_DEBUG_ERROR_IF(DerivativeDirectionIndex > 2)
        << "Invalid DerivativeDirectionIndex [ DerivativeDirectionIndex = " << DerivativeDirectionIndex
        << ", expected DerivativeDirectionIndex < 3 ].\n";

    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double x10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 0);

    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double y10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 1);

    const double z10 = rGeometry[1].Z() - rGeometry[0].Z();
    const double z10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 2);

    const double x20 = rGeometry[2].X() - rGeometry[0].X();
    const double x20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 0);

    const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    const double y20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 1);

    const double z20 = rGeometry[2].Z() - rGeometry[0].Z();
    const double z20_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 2, 0, 2);

    const double x30 = rGeometry[3].X() - rGeometry[0].X();
    const double x30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 0);

    const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    const double y30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 1);

    const double z30 = rGeometry[3].Z() - rGeometry[0].Z();
    const double z30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 2);

    const double detJ = 0.5 * (x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30);
    double detJ_derivative = 0.0;
    detJ_derivative += x10_derivative * y20 * z30;
    detJ_derivative += x10 * y20_derivative * z30;
    detJ_derivative += x10 * y20 * z30_derivative;

    detJ_derivative -= x10_derivative * y30 * z20;
    detJ_derivative -= x10 * y30_derivative * z20;
    detJ_derivative -= x10 * y30 * z20_derivative;

    detJ_derivative += y10_derivative * z20 * x30;
    detJ_derivative += y10 * z20_derivative * x30;
    detJ_derivative += y10 * z20 * x30_derivative;

    detJ_derivative -= y10_derivative * x20 * z30;
    detJ_derivative -= y10 * x20_derivative * z30;
    detJ_derivative -= y10 * x20 * z30_derivative;

    detJ_derivative += z10_derivative * x20 * y30;
    detJ_derivative += z10 * x20_derivative * y30;
    detJ_derivative += z10 * x20 * y30_derivative;

    detJ_derivative -= z10_derivative * y20 * x30;
    detJ_derivative -= z10 * y20_derivative * x30;
    detJ_derivative -= z10 * y20 * x30_derivative;

    detJ_derivative *= 0.5;

    return (1./3.) * detJ_derivative / std::pow(detJ, 2./3.);

    KRATOS_CATCH("");
}

// Hexahedra3D8 version.
template<>
double ElementSizeCalculator<3,8>::AverageElementSize(const Geometry<Node >& rGeometry)
{

    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double z10 = rGeometry[1].Z() - rGeometry[0].Z();

    const double x30 = rGeometry[3].X() - rGeometry[0].X();
    const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    const double z30 = rGeometry[3].Z() - rGeometry[0].Z();

    const double x40 = rGeometry[4].X() - rGeometry[0].X();
    const double y40 = rGeometry[4].Y() - rGeometry[0].Y();
    const double z40 = rGeometry[4].Z() - rGeometry[0].Z();

    const double detJ = x10 * y30 * z40 - x10 * y40 * z30 + y10 * z30 * x40 - y10 * x30 * z40 + z10 * x30 * y40 - z10 * y30 * x40;
    return std::pow(detJ,1./3.);
}

template<>
double ElementSizeCalculator<3,8>::AverageElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(DerivativeNodeIndex > 7)
        << "Invalid DerivativeNodeIndex [ DerivativeNodeIndex = " << DerivativeNodeIndex
        << ", expected DerivativeNodeIndex < 8 ].\n";

    KRATOS_DEBUG_ERROR_IF(DerivativeDirectionIndex > 2)
        << "Invalid DerivativeDirectionIndex [ DerivativeDirectionIndex = " << DerivativeDirectionIndex
        << ", expected DerivativeDirectionIndex < 3 ].\n";

    const double x10 = rGeometry[1].X() - rGeometry[0].X();
    const double x10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 0);
    const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    const double y10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 1);
    const double z10 = rGeometry[1].Z() - rGeometry[0].Z();
    const double z10_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 1, 0, 2);

    const double x30 = rGeometry[3].X() - rGeometry[0].X();
    const double x30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 0);
    const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    const double y30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 1);
    const double z30 = rGeometry[3].Z() - rGeometry[0].Z();
    const double z30_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 3, 0, 2);

    const double x40 = rGeometry[4].X() - rGeometry[0].X();
    const double x40_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 4, 0, 0);
    const double y40 = rGeometry[4].Y() - rGeometry[0].Y();
    const double y40_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 4, 0, 1);
    const double z40 = rGeometry[4].Z() - rGeometry[0].Z();
    const double z40_derivative = EdgeLengthDerivative(DerivativeNodeIndex, DerivativeDirectionIndex, 4, 0, 2);

    const double detJ = x10 * y30 * z40 - x10 * y40 * z30 + y10 * z30 * x40 - y10 * x30 * z40 + z10 * x30 * y40 - z10 * y30 * x40;
    double detJ_derivative = 0.0;

    detJ_derivative += x10_derivative * y30 * z40;
    detJ_derivative += x10 * y30_derivative * z40;
    detJ_derivative += x10 * y30 * z40_derivative;

    detJ_derivative -= x10_derivative * y40 * z30;
    detJ_derivative -= x10 * y40_derivative * z30;
    detJ_derivative -= x10 * y40 * z30_derivative;

    detJ_derivative += y10_derivative * z30 * x40;
    detJ_derivative += y10 * z30_derivative * x40;
    detJ_derivative += y10 * z30 * x40_derivative;

    detJ_derivative -= y10_derivative * x30 * z40;
    detJ_derivative -= y10 * x30_derivative * z40;
    detJ_derivative -= y10 * x30 * z40_derivative;

    detJ_derivative += z10_derivative * x30 * y40;
    detJ_derivative += z10 * x30_derivative * y40;
    detJ_derivative += z10 * x30 * y40_derivative;

    detJ_derivative -= z10_derivative * y30 * x40;
    detJ_derivative -= z10 * y30_derivative * x40;
    detJ_derivative -= z10 * y30 * x40_derivative;

    return (1./3.) * detJ_derivative / std::pow(detJ, 2./3);

    KRATOS_CATCH("");
}

// Hexahedra3D27 version.
template<>
double ElementSizeCalculator<3,27>::AverageElementSize(const Geometry<Node >& rGeometry)
{

    const double average_element_size = ElementSizeCalculator<3,8>::AverageElementSize(rGeometry);

    return average_element_size;
}

template<>
double ElementSizeCalculator<3,27>::AverageElementSizeDerivative(
    const unsigned int DerivativeNodeIndex,
    const unsigned int DerivativeDirectionIndex,
    const Geometry<Node >& rGeometry)
{
    KRATOS_TRY
    double element_size_derivative = 0.0;

    if (DerivativeNodeIndex < 8)
        element_size_derivative = ElementSizeCalculator<3,8>::AverageElementSizeDerivative(DerivativeNodeIndex,DerivativeDirectionIndex,rGeometry);

    return element_size_derivative;

    KRATOS_CATCH("");
}

// Triangle2D3 version.
template<>
double ElementSizeCalculator<2,3>::ProjectedElementSize(
    const Geometry<Node >& rGeometry,
    const array_1d<double,3>& rVelocity)
{

    double Hvel = 0.0;

    const unsigned int NumNodes = 3;

    // Loop over edges looking for maximum 'projected' length
    array_1d<double,3> Edge;
    for(unsigned int i = 0; i < NumNodes; ++i)
    {
        unsigned int j = (i+1) % NumNodes;
        Edge = rGeometry[j] - rGeometry[i];
        double lu = rVelocity[0] * Edge[0];
        for (unsigned int d = 1; d < 2; ++d)
            lu += rVelocity[d] * Edge[d];
        lu = fabs(lu);
        if(Hvel < lu) Hvel = lu;
    }

    if (Hvel > 0.0)
    {
        const double VelNorm = std::sqrt(rVelocity[0]*rVelocity[0] + rVelocity[1]*rVelocity[1] + rVelocity[2]*rVelocity[2]);
        Hvel /= VelNorm;
    }

    return Hvel;
}

// Triangle2D6 version
template<>
double ElementSizeCalculator<2,6>::ProjectedElementSize(
    const Geometry<Node > &rGeometry,
    const array_1d<double,3>& rVelocity)
{
    KRATOS_ERROR << "This function has not been implemented yet." << std::endl;
    return 0.0; // Just to avoid warning during compilations
}

// Quadrilateral2D4 version.
template<>
double ElementSizeCalculator<2,4>::ProjectedElementSize(
    const Geometry<Node >& rGeometry,
    const array_1d<double,3>& rVelocity)
{

    const double Hvel = ElementSizeCalculator<2,3>::ProjectedElementSize(rGeometry,rVelocity);

    return Hvel;
}

// Quadrilateral2D9 version
template<>
double ElementSizeCalculator<2,9>::ProjectedElementSize(
    const Geometry<Node > &rGeometry,
    const array_1d<double,3>& rVelocity)
{
    KRATOS_ERROR << "This function has not been implemented yet." << std::endl;
    return 0.0; // Just to avoid warning during compilations
}

// Tetrahedra3D4 version.
template<>
double ElementSizeCalculator<3,4>::ProjectedElementSize(
    const Geometry<Node >& rGeometry,
    const array_1d<double,3>& rVelocity)
{

    double Hvel = 0.0;

    const unsigned int NumNodes = 4;

    // Loop over edges looking for maximum 'projected' length
    array_1d<double,3> Edge;
    for(unsigned int i = 0; i < NumNodes; ++i)
    {
        for(unsigned int j = i+1; j < NumNodes; ++j)
        {
            Edge = rGeometry[j] - rGeometry[i];
            double lu = rVelocity[0] * Edge[0];
            for (unsigned int d = 1; d < 3; ++d)
                lu += rVelocity[d] * Edge[d];
            lu = fabs(lu);
            if(Hvel < lu) Hvel = lu;
        }
    }

    if (Hvel > 0.0)
    {
        const double VelNorm = std::sqrt(rVelocity[0]*rVelocity[0] + rVelocity[1]*rVelocity[1] + rVelocity[2]*rVelocity[2]);
        Hvel /= VelNorm;
    }

    return Hvel;
}

// Tetrahedra3D10 version.
template<>
double ElementSizeCalculator<3,10>::ProjectedElementSize(
    const Geometry<Node > &rGeometry,
    const array_1d<double,3>& rVelocity)
{
    KRATOS_ERROR << "This function has not been implemented yet." << std::endl;
    return 0.0; // Just to avoid warning during compilations
}

// Prism3D6 version
template<>
double ElementSizeCalculator<3,6>::ProjectedElementSize(
    const Geometry<Node > &rGeometry,
    const array_1d<double,3>& rVelocity)
{
    KRATOS_ERROR << "This function has not been implemented yet." << std::endl;
    return 0.0; // Just to avoid warning during compilations
}

// Hexahedra3D8 version.
template<>
double ElementSizeCalculator<3,8>::ProjectedElementSize(
    const Geometry<Node >& rGeometry,
    const array_1d<double,3>& rVelocity)
{

    double Hvel = 0.0;

    // Logic: define a box given by hexahedra edges 10,30,40 (which I'm assuming to be orthogonal)
    // Transform the given direction to the coordinate system defined by the edges.
    // In this reference frame, place a vector aligned with rVelocity on the origin and see
    // which side of the box intersects the vector (that is: which component of the vector is larger)
    // Then scale the vector to hit the limit of the box.
    // The lenght of the vector is what we are looking for.
    // NOTE: we use absolute values on all checks, this allows us to simplify the problem
    // to a single sector (otherwise we would start by determining in which of the eight sectors
    // given by +-x, +-y, +-z we have to look for the intersection).
    array_1d<double,3> U = rVelocity;
    double u_norm = std::sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);
    if (u_norm > 1e-12)
    {
        //Normalize U
        U /= std::sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);

        array_1d<double,3> v10 = rGeometry[1].Coordinates() - rGeometry[0].Coordinates();
        array_1d<double,3> v30 = rGeometry[3].Coordinates() - rGeometry[0].Coordinates();
        array_1d<double,3> v40 = rGeometry[4].Coordinates() - rGeometry[0].Coordinates();

        // Express U in the coordinate system defined by {v10,v30,v40}
        Matrix Q = ZeroMatrix(3,3);
        for (unsigned int i = 0; i < 3; i++)
        {
            Q(i,0) = v10[i];
            Q(i,1) = v30[i];
            Q(i,2) = v40[i];
        }

        // Invert Matrix Q
        Matrix QInv;
        double det;
        MathUtils<double>::InvertMatrix(Q,QInv,det);

        array_1d<double,3> Uq = ZeroVector(3);
        for (unsigned int i = 0; i < 3; i++)
        {
            for (unsigned int j = 0; j < 3; j++)
            {
                Uq[i] += QInv(i,j)*U[j];
            }

            // Work in absolute values
            Uq[i] = std::fabs(Uq[i]);
        }

        double max_v = Uq[0];
        for (unsigned int d = 1; d < 3; d++)
            if (Uq[d] > max_v)
                max_v = Uq[d];

        double scale = 1.0/max_v;
        Uq *= scale;

        // Undo the transform
        for (unsigned int i = 0; i < 3; i++)
        {
            U[i] = 0.0;
            for (unsigned int j = 0; j < 3; j++)
            {
                U[i] += Q(i,j)*Uq[j];
            }
        }

        // Module
        Hvel = std::sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);

        if (Hvel > 0.0)
        {
            double VelNorm = std::sqrt(rVelocity[0]*rVelocity[0] + rVelocity[1]*rVelocity[1] + rVelocity[2]*rVelocity[2]);
            Hvel /= VelNorm;
        }
    }

    return Hvel;
}

// Hexahedra3D27 version.
template<>
double ElementSizeCalculator<3,27>::ProjectedElementSize(
    const Geometry<Node > &rGeometry,
    const array_1d<double,3>& rVelocity)
{
    KRATOS_ERROR << "This function has not been implemented yet." << std::endl;
    return 0.0; // Just to avoid warning during compilations
}

// Triangle2D3 version.
template<std::size_t TDim, std::size_t TNumNodes>
double ElementSizeCalculator<TDim,TNumNodes>::GradientsElementSize(const BoundedMatrix<double,3,2>& rDN_DX)
{

    double h = 0.0;
    for(unsigned int i = 0; i < 3; ++i){
        double h_inv = 0.0;
        for(unsigned int k = 0; k < 2; ++k){
            h_inv += rDN_DX(i,k)*rDN_DX(i,k);
        }
        h += 1.0/h_inv;
    }
    return sqrt(h)/3.0;
}

// Tetrahedra3D4 version.
template<std::size_t TDim, std::size_t TNumNodes>
double ElementSizeCalculator<TDim,TNumNodes>::GradientsElementSize(const BoundedMatrix<double,4,3>& rDN_DX)
{

    double h = 0.0;
    for(unsigned int i = 0; i < 4; ++i){
        double h_inv = 0.0;
        for(unsigned int k = 0; k < 3; ++k){
            h_inv += rDN_DX(i,k)*rDN_DX(i,k);
        }
        h += 1.0/h_inv;
    }
    return sqrt(h)/4.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class ElementSizeCalculator<2,3>;
template class ElementSizeCalculator<2,6>;
template class ElementSizeCalculator<2,4>;
template class ElementSizeCalculator<2,9>;
template class ElementSizeCalculator<3,4>;
template class ElementSizeCalculator<3,6>;
template class ElementSizeCalculator<3,10>;
template class ElementSizeCalculator<3,8>;
template class ElementSizeCalculator<3,27>;

///////////////////////////////////////////////////////////////////////////////////////////////////

}
