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

namespace Kratos
{

class InterfaceElementUtilities
{
public:
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

    //----------------------------------------------------------------------------------------
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

    //----------------------------------------------------------------------------------------
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

    //----------------------------------------------------------------------------------------
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

    //----------------------------------------------------------------------------------------
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

    //----------------------------------------------------------------------------------------
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

    //----------------------------------------------------------------------------------------
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

    //----------------------------------------------------------------------------------------
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
        Matrix rotation_matrix = ZeroMatrix(2, 2);

        array_1d<double,3> xi{0.0, 0.0, 0.0};
        Matrix shape_functions_gradients;
        rGeometry.ShapeFunctionsLocalGradients(shape_functions_gradients, xi);
        KRATOS_INFO("Calculate2DRotationMatrix") << "shape_functions_gradients: " << shape_functions_gradients << std::endl;

        std::vector<array_1d<double, 3>> mid_points;
        for (std::size_t i = 0; i < rGeometry.size()/2; ++i) {
            mid_points.push_back(0.5 * (rGeometry[i] + rGeometry[i + rGeometry.size()/2]));
        }

        KRATOS_INFO("Calculate2DRotationMatrix") << "mid_points: " << mid_points << std::endl;

        array_1d<double, 3> tangential_vector = ZeroVector(3);
        for (std::size_t i = 0; i < mid_points.size(); ++i) {
            tangential_vector += mid_points[i] * shape_functions_gradients(0, i);
        }

        tangential_vector /= norm_2(tangential_vector);
        KRATOS_INFO("Calculate2DRotationMatrix") << "tangential_vector: " << tangential_vector << std::endl;
        auto out_of_plane = array_1d<double, 3> {0.0, 0.0, 1.0};
        array_1d<double, 3> normal_vector = MathUtils<double>::CrossProduct(tangential_vector, out_of_plane);


        rotation_matrix(0,0) = tangential_vector[0];
        rotation_matrix(0,1) = -tangential_vector[1];
        rotation_matrix(1,0) = tangential_vector[1];
        rotation_matrix(1,1) = tangential_vector[0];
        return rotation_matrix;
    }

}; /* Class InterfaceElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_INTERFACE_ELEMENT_UTILITIES defined */
