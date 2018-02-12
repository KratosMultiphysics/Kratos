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
//

#include "element_size_calculator.h"

namespace Kratos {

template< std::size_t TDim, std::size_t TNumNodes >
ElementSizeCalculator<TDim,TNumNodes>::~ElementSizeCalculator() {};

// Triangle2D3 version.
template<>
double ElementSizeCalculator<2,3>::MinimumElementSize(const Geometry<Node<3> >& rGeometry) {

    /* Calculate node-edge distances */
    double x10 = rGeometry[1].X() - rGeometry[0].X();
    double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    double x20 = rGeometry[2].X() - rGeometry[0].X();
    double y20 = rGeometry[2].Y() - rGeometry[0].Y();

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

// Quadrilateral2D4 version.
template<>
double ElementSizeCalculator<2,4>::MinimumElementSize(const Geometry<Node<3> >& rGeometry) {

    const Node<3>& r_node_0 = rGeometry[0];
    const Node<3>& r_node_1 = rGeometry[1];
    const Node<3>& r_node_2 = rGeometry[2];
    const Node<3>& r_node_3 = rGeometry[3];

    // Calculate face centers
    double x10 = (r_node_1.X() + r_node_3.X())/2.;
    double y10 = (r_node_1.Y() + r_node_3.Y())/2.;

    double x21 = (r_node_2.X() + r_node_1.X())/2.;
    double y21 = (r_node_2.Y() + r_node_1.Y())/2.;

    double x32 = (r_node_3.X() + r_node_2.X())/2.;
    double y32 = (r_node_3.Y() + r_node_2.Y())/2.;

    double x03 = (r_node_0.X() + r_node_3.X())/2.;
    double y03 = (r_node_0.Y() + r_node_3.Y())/2.;

    // Distance between face centers (xi direction)
    double dxi_x = x21 - x03;
    double dxi_y = y21 - y03;
    
    double h2_xi = dxi_x*dxi_x + dxi_y*dxi_y;

    // Distance between face centers (eta direction)
    double deta_x = x32 - x10;
    double deta_y = y32 - y10;

    double h2_eta = deta_x*deta_x + deta_y*deta_y;

    double h2 = h2_xi < h2_eta ? h2_xi : h2_eta;

    return std::sqrt(h2);
}

// Tetrahedra3D4 version.
template<>
double ElementSizeCalculator<3,4>::MinimumElementSize(const Geometry<Node<3> >& rGeometry) {

    /* Calculate distances between each node and the opposite face */
    double x10 = rGeometry[1].X() - rGeometry[0].X();
    double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    double z10 = rGeometry[1].Z() - rGeometry[0].Z();

    double x20 = rGeometry[2].X() - rGeometry[0].X();
    double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    double z20 = rGeometry[2].Z() - rGeometry[0].Z();

    double x30 = rGeometry[3].X() - rGeometry[0].X();
    double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    double z30 = rGeometry[3].Z() - rGeometry[0].Z();

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

// Hexahedra3D8 version. We use the distance between face centers to compute lengths.
template<>
double ElementSizeCalculator<3,8>::MinimumElementSize(const Geometry<Node<3> >& rGeometry) {

    const Node<3>& r_node_0 = rGeometry[0];
    const Node<3>& r_node_1 = rGeometry[1];
    const Node<3>& r_node_2 = rGeometry[2];
    const Node<3>& r_node_3 = rGeometry[3];
    const Node<3>& r_node_4 = rGeometry[4];
    const Node<3>& r_node_5 = rGeometry[5];
    const Node<3>& r_node_6 = rGeometry[6];
    const Node<3>& r_node_7 = rGeometry[7];

    // Calculate face centers
    array_1d<double,3> low_xi  = 0.25 * (r_node_0.Coordinates() + r_node_4.Coordinates() + r_node_3.Coordinates() + r_node_7.Coordinates());
    array_1d<double,3> high_xi = 0.25 * (r_node_1.Coordinates() + r_node_2.Coordinates() + r_node_6.Coordinates() + r_node_5.Coordinates());
    array_1d<double,3> low_eta  = 0.25 * (r_node_0.Coordinates() + r_node_1.Coordinates() + r_node_5.Coordinates() + r_node_4.Coordinates());
    array_1d<double,3> high_eta = 0.25 * (r_node_3.Coordinates() + r_node_7.Coordinates() + r_node_6.Coordinates() + r_node_2.Coordinates());
    array_1d<double,3> low_dseta  = 0.25 * (r_node_0.Coordinates() + r_node_3.Coordinates() + r_node_2.Coordinates() + r_node_1.Coordinates());
    array_1d<double,3> high_dseta = 0.25 * (r_node_4.Coordinates() + r_node_5.Coordinates() + r_node_6.Coordinates() + r_node_7.Coordinates());

    // Distance between face centers (xi direction)
    array_1d<double,3> d_xi = high_xi - low_xi;
    double h2_xi = d_xi[0]*d_xi[0] + d_xi[1]*d_xi[1] + d_xi[2]*d_xi[2];

    // Distance between face centers (eta direction)
    array_1d<double,3> d_eta = high_eta - low_eta;
    double h2_eta = d_eta[0]*d_eta[0] + d_eta[1]*d_eta[1] + d_eta[2]*d_eta[2];

    // Distance between face centers (dseta direction)
    array_1d<double,3> d_dseta = high_dseta - low_dseta;
    double h2_dseta = d_dseta[0]*d_dseta[0] + d_dseta[1]*d_dseta[1] + d_dseta[2]*d_dseta[2];

    double h2 = h2_xi < h2_eta ? h2_xi : h2_eta;
    h2 = h2 < h2_dseta ? h2 : h2_dseta;

    return std::sqrt(h2);
}

// Triangle2D3 version. 
template<>
double ElementSizeCalculator<2,3>::GradientsElementSize(const bounded_matrix<double,3,2>& rDN_DX) {

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
template<>
double ElementSizeCalculator<3,4>::GradientsElementSize(const bounded_matrix<double,4,3>& rDN_DX) {

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
template class ElementSizeCalculator<2,4>;
template class ElementSizeCalculator<3,4>;
template class ElementSizeCalculator<3,8>;

///////////////////////////////////////////////////////////////////////////////////////////////////

}