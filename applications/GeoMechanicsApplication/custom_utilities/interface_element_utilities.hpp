// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#if !defined(KRATOS_INTERFACE_ELEMENT_UTILITIES)
#define KRATOS_INTERFACE_ELEMENT_UTILITIES

// Project includes
#include "includes/element.h"
#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos
{

class InterfaceElementUtilities
{
public:
    static inline void CalculateNuMatrix(BoundedMatrix<double, 2, 4>& rNu,
                                         const Matrix&                Ncontainer,
                                         const unsigned int&          GPoint)
    {
        // Line_interface_2d_2
        rNu(0, 0) = -Ncontainer(GPoint, 0);
        rNu(0, 2) = Ncontainer(GPoint, 1);
        rNu(1, 1) = -Ncontainer(GPoint, 0);
        rNu(1, 3) = Ncontainer(GPoint, 1);
    }

    static inline void CalculateNuMatrix(BoundedMatrix<double, 2, 8>& rNu,
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

    static inline void CalculateNuMatrix(BoundedMatrix<double, 3, 12>& rNu,
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

    static inline void CalculateNuMatrix(BoundedMatrix<double, 3, 18>& rNu,
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

    static inline void CalculateNuMatrix(BoundedMatrix<double, 3, 24>& rNu,
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

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void FillPermeabilityMatrix(BoundedMatrix<double, 2, 2>& rPermeabilityMatrix,
                                              const double&                JointWidth,
                                              const double&                Transversal_Permeability)
    {
        // Quadrilateral_interface_2d_4
        rPermeabilityMatrix(0, 0) = JointWidth * JointWidth / 12.0;
        rPermeabilityMatrix(1, 1) = Transversal_Permeability;
    }

    static inline void FillPermeabilityMatrix(BoundedMatrix<double, 3, 3>& rPermeabilityMatrix,
                                              const double&                JointWidth,
                                              const double&                Transversal_Permeability)
    {
        // Prism_interface_3d_6 and Hexahedral_interface_3d_8
        rPermeabilityMatrix(0, 0) = JointWidth * JointWidth / 12.0;
        rPermeabilityMatrix(1, 1) = JointWidth * JointWidth / 12.0;
        rPermeabilityMatrix(2, 2) = Transversal_Permeability;
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateVoigtVector(array_1d<double, 2>& rVoigtVector)
    {
        // Quadrilateral_interface_2d_4
        rVoigtVector[0] = 0.0;
        rVoigtVector[1] = 1.0;
    }

    static inline void CalculateVoigtVector(array_1d<double, 3>& rVoigtVector)
    {
        // Prism_interface_3d_6 and Hexahedral_interface_3d_8
        rVoigtVector[0] = 0.0;
        rVoigtVector[1] = 0.0;
        rVoigtVector[2] = 1.0;
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void CalculateLinkPermeabilityMatrix(BoundedMatrix<double, 2, 2>& rPermeabilityMatrix,
                                                       const double& JointWidth)
    {
        // Quadrilateral_interface_2d_4
        rPermeabilityMatrix(0, 0) = JointWidth * JointWidth / 12.0;
        rPermeabilityMatrix(1, 1) = JointWidth * JointWidth / 12.0;
    }

    static inline void CalculateLinkPermeabilityMatrix(BoundedMatrix<double, 3, 3>& rPermeabilityMatrix,
                                                       const double& JointWidth)
    {
        // Prism_interface_3d_6 and Hexahedral_interface_3d_8
        rPermeabilityMatrix(0, 0) = JointWidth * JointWidth / 12.0;
        rPermeabilityMatrix(1, 1) = JointWidth * JointWidth / 12.0;
        rPermeabilityMatrix(2, 2) = JointWidth * JointWidth / 12.0;
    }

    static Matrix Calculate2DRotationMatrix(const Geometry<Node>& rGeometry)
    {
        array_1d<double, 3> xi{0.0, 0.0, 0.0};
        Matrix              shape_functions_gradients;
        rGeometry.ShapeFunctionsLocalGradients(shape_functions_gradients, xi);
        KRATOS_INFO("ShapeFunctionsGradients") << shape_functions_gradients << std::endl;
        const auto mid_points = CalculateMidPoints(rGeometry);
        KRATOS_INFO("MidPoints") << mid_points << std::endl;
        const auto tangential_vector = CalculateTangentialVector(shape_functions_gradients, mid_points);
        KRATOS_INFO("TangetialVector") << tangential_vector << std::endl;
        // clang-format off
        Matrix rotation_matrix(2, 2);
        rotation_matrix <<= tangential_vector[0], -tangential_vector[1],
                            tangential_vector[1], tangential_vector[0];
        // clang-format on

        KRATOS_INFO("RotationMatrix") << rotation_matrix << std::endl;

        return rotation_matrix;
    }

private:
    static std::vector<array_1d<double, 3>> CalculateMidPoints(const Geometry<Node>& rGeometry)
    {
        std::vector<array_1d<double, 3>> mid_points;

        const auto half_way_point = rGeometry.begin() + rGeometry.PointsNumber() / 2;
        std::transform(rGeometry.begin(), half_way_point, half_way_point, std::back_inserter(mid_points),
                       [](const Node& rNode1, const Node& rNode2) -> array_1d<double, 3> {
            return 0.5 * (rNode1 + rNode2);
        });

        return mid_points;
    }

    static array_1d<double, 3> CalculateTangentialVector(const Matrix& rShapeFunctionGradients,
                                                         const std::vector<array_1d<double, 3>>& rLocations)
    {
        array_1d<double, 3> tangential_vector = ZeroVector(3);
        for (std::size_t i = 0; i < rLocations.size(); ++i) {
            tangential_vector += rLocations[i] * rShapeFunctionGradients(i, 0);
        }
        tangential_vector /= norm_2(tangential_vector);

        return tangential_vector;
    }

}; /* Class InterfaceElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_INTERFACE_ELEMENT_UTILITIES defined */
