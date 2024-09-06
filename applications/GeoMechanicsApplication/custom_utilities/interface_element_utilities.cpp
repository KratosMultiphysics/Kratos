// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "interface_element_utilities.h"
#include "math_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos
{

void InterfaceElementUtilities::CalculateNuMatrix(BoundedMatrix<double, 2, 4>& rNu,
                                                  const Matrix&                Ncontainer,
                                                  const unsigned int&          GPoint)
{
    // Line_interface_2d_2
    rNu(0, 0) = -Ncontainer(GPoint, 0);
    rNu(0, 2) = Ncontainer(GPoint, 1);
    rNu(1, 1) = -Ncontainer(GPoint, 0);
    rNu(1, 3) = Ncontainer(GPoint, 1);
}

void InterfaceElementUtilities::CalculateNuMatrix(BoundedMatrix<double, 2, 8>& rNu,
                                                  const Matrix&                Ncontainer,
                                                  const unsigned int&          GPoint)
{
    // Quadrilateral_interface_2d_4
    rNu(0, 0) = -Ncontainer(GPoint, 0);
    rNu(0, 2) = -Ncontainer(GPoint, 1);
    rNu(1, 1) = -Ncontainer(GPoint, 0);
    rNu(1, 3) = -Ncontainer(GPoint, 1);

    rNu(0, 4) = Ncontainer(GPoint, 2);
    rNu(0, 6) = Ncontainer(GPoint, 3);
    rNu(1, 5) = Ncontainer(GPoint, 2);
    rNu(1, 7) = Ncontainer(GPoint, 3);
}

void InterfaceElementUtilities::CalculateNuMatrix(BoundedMatrix<double, 3, 12>& rNu,
                                                  const Matrix&                 Ncontainer,
                                                  const unsigned int&           GPoint)
{
    // Quadrilateral_interface_3d_4
    rNu(0, 0) = -Ncontainer(GPoint, 0);
    rNu(0, 3) = -Ncontainer(GPoint, 1);
    rNu(1, 1) = -Ncontainer(GPoint, 0);
    rNu(1, 4) = -Ncontainer(GPoint, 1);
    rNu(2, 2) = -Ncontainer(GPoint, 0);
    rNu(2, 5) = -Ncontainer(GPoint, 1);

    rNu(0, 6)  = Ncontainer(GPoint, 2);
    rNu(0, 9)  = Ncontainer(GPoint, 3);
    rNu(1, 7)  = Ncontainer(GPoint, 2);
    rNu(1, 10) = Ncontainer(GPoint, 3);
    rNu(2, 8)  = Ncontainer(GPoint, 2);
    rNu(2, 11) = Ncontainer(GPoint, 3);
}

void InterfaceElementUtilities::CalculateNuMatrix(BoundedMatrix<double, 3, 18>& rNu,
                                                  const Matrix&                 Ncontainer,
                                                  const unsigned int&           GPoint)
{
    // Prism_interface_3d_6
    rNu(0, 0) = -Ncontainer(GPoint, 0);
    rNu(0, 3) = -Ncontainer(GPoint, 1);
    rNu(0, 6) = -Ncontainer(GPoint, 2);
    rNu(1, 1) = -Ncontainer(GPoint, 0);
    rNu(1, 4) = -Ncontainer(GPoint, 1);
    rNu(1, 7) = -Ncontainer(GPoint, 2);
    rNu(2, 2) = -Ncontainer(GPoint, 0);
    rNu(2, 5) = -Ncontainer(GPoint, 1);
    rNu(2, 8) = -Ncontainer(GPoint, 2);

    rNu(0, 9)  = Ncontainer(GPoint, 3);
    rNu(0, 12) = Ncontainer(GPoint, 4);
    rNu(0, 15) = Ncontainer(GPoint, 5);
    rNu(1, 10) = Ncontainer(GPoint, 3);
    rNu(1, 13) = Ncontainer(GPoint, 4);
    rNu(1, 16) = Ncontainer(GPoint, 5);
    rNu(2, 11) = Ncontainer(GPoint, 3);
    rNu(2, 14) = Ncontainer(GPoint, 4);
    rNu(2, 17) = Ncontainer(GPoint, 5);
}

void InterfaceElementUtilities::CalculateNuMatrix(BoundedMatrix<double, 3, 24>& rNu,
                                                  const Matrix&                 Ncontainer,
                                                  const unsigned int&           GPoint)
{
    // Hexahedral_interface_3d_8
    rNu(0, 0)  = -Ncontainer(GPoint, 0);
    rNu(0, 3)  = -Ncontainer(GPoint, 1);
    rNu(0, 6)  = -Ncontainer(GPoint, 2);
    rNu(0, 9)  = -Ncontainer(GPoint, 3);
    rNu(1, 1)  = -Ncontainer(GPoint, 0);
    rNu(1, 4)  = -Ncontainer(GPoint, 1);
    rNu(1, 7)  = -Ncontainer(GPoint, 2);
    rNu(1, 10) = -Ncontainer(GPoint, 3);
    rNu(2, 2)  = -Ncontainer(GPoint, 0);
    rNu(2, 5)  = -Ncontainer(GPoint, 1);
    rNu(2, 8)  = -Ncontainer(GPoint, 2);
    rNu(2, 11) = -Ncontainer(GPoint, 3);

    rNu(0, 12) = Ncontainer(GPoint, 4);
    rNu(0, 15) = Ncontainer(GPoint, 5);
    rNu(0, 18) = Ncontainer(GPoint, 6);
    rNu(0, 21) = Ncontainer(GPoint, 7);
    rNu(1, 13) = Ncontainer(GPoint, 4);
    rNu(1, 16) = Ncontainer(GPoint, 5);
    rNu(1, 19) = Ncontainer(GPoint, 6);
    rNu(1, 22) = Ncontainer(GPoint, 7);
    rNu(2, 14) = Ncontainer(GPoint, 4);
    rNu(2, 17) = Ncontainer(GPoint, 5);
    rNu(2, 20) = Ncontainer(GPoint, 6);
    rNu(2, 23) = Ncontainer(GPoint, 7);
}

void InterfaceElementUtilities::FillPermeabilityMatrix(BoundedMatrix<double, 2, 2>& rPermeabilityMatrix,
                                                       const double& JointWidth,
                                                       const double& Transversal_Permeability)
{
    // Quadrilateral_interface_2d_4
    rPermeabilityMatrix(0, 0) = JointWidth * JointWidth / 12.0;
    rPermeabilityMatrix(1, 1) = Transversal_Permeability;
}

void InterfaceElementUtilities::FillPermeabilityMatrix(BoundedMatrix<double, 3, 3>& rPermeabilityMatrix,
                                                       const double& JointWidth,
                                                       const double& Transversal_Permeability)
{
    // Prism_interface_3d_6 and Hexahedral_interface_3d_8
    rPermeabilityMatrix(0, 0) = JointWidth * JointWidth / 12.0;
    rPermeabilityMatrix(1, 1) = JointWidth * JointWidth / 12.0;
    rPermeabilityMatrix(2, 2) = Transversal_Permeability;
}

void InterfaceElementUtilities::CalculateVoigtVector(array_1d<double, 2>& rVoigtVector)
{
    // Quadrilateral_interface_2d_4
    rVoigtVector[0] = 0.0;
    rVoigtVector[1] = 1.0;
}

void InterfaceElementUtilities::CalculateVoigtVector(array_1d<double, 3>& rVoigtVector)
{
    // Prism_interface_3d_6 and Hexahedral_interface_3d_8
    rVoigtVector[0] = 0.0;
    rVoigtVector[1] = 0.0;
    rVoigtVector[2] = 1.0;
}

void InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(BoundedMatrix<double, 2, 2>& rPermeabilityMatrix,
                                                                const double& JointWidth)
{
    // Quadrilateral_interface_2d_4
    rPermeabilityMatrix(0, 0) = JointWidth * JointWidth / 12.0;
    rPermeabilityMatrix(1, 1) = JointWidth * JointWidth / 12.0;
}

void InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(BoundedMatrix<double, 3, 3>& rPermeabilityMatrix,
                                                                const double& JointWidth)
{
    // Prism_interface_3d_6 and Hexahedral_interface_3d_8
    rPermeabilityMatrix(0, 0) = JointWidth * JointWidth / 12.0;
    rPermeabilityMatrix(1, 1) = JointWidth * JointWidth / 12.0;
    rPermeabilityMatrix(2, 2) = JointWidth * JointWidth / 12.0;
}

Matrix InterfaceElementUtilities::Calculate2DRotationMatrixForLineGeometry(const Geometry<Node>& rGeometry,
                                                                            const array_1d<double, 3>& rLocalCoordinate)
{
    // Since the shape functions depend on one coordinate only
    // for lines, the jacobian only has one column.
    Matrix jacobian;
    rGeometry.Jacobian(jacobian, rLocalCoordinate);
    const auto tangential_vector = GeoMechanicsMathUtilities::Normalized(Vector{column(jacobian, 0)});

    // clang-format off
    Matrix result(2, 2);
    result <<= tangential_vector[0], -tangential_vector[1],
               tangential_vector[1],  tangential_vector[0];
    // clang-format on

    return result;
}

} /* namespace Kratos.*/
