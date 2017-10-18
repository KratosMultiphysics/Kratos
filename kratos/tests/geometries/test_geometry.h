//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

#if !defined(KRATOS_TEST_GEOMETRY_H_INCLUDED)
#define KRATOS_TEST_GEOMETRY_H_INCLUDED

// System includes
#include <set>

// External includes

// Project includes
#include "includes/node.h"
#include "includes/element.h"
#include "testing/testing.h"
#include "geometries/geometry.h"

namespace Kratos {
namespace Testing {

    constexpr double EPSILON = std::numeric_limits<double>::epsilon();
    constexpr double TOLERANCE = 1e-6;

    /// Factory functions

    /** Generates a point.
     * Generates a point. If no coordinates are provided, random ones between 0 and 1 are selected.
     * Coordinate Z is forced to 0 for 2D element tests.
     * @param  coordx Coordinate X
     * @param  coordy Coordinate Y
     * @return        A point with coordinates coordx, coordy
     */
    template<class TPointType>
    typename TPointType::Pointer GeneratePoint(
        double coordx = (((double) std::rand() / (RAND_MAX)) + 1),
        double coordy = (((double) std::rand() / (RAND_MAX)) + 1),
        double coordz = (((double) std::rand() / (RAND_MAX)) + 1)) {
      static int id = 0;
      return typename TPointType::Pointer(new TPointType(id++, coordx, coordy, coordz));
    }

    /// Auxiliar check functions (from geometry_tester.h)
    /// - All this functions should probably me moved somewhere else.

    /** Gets the corresponding string of the integration method provided.
     * Gets the corresponding string of the integration method provided.
     * @param  geom       Geometry that is used for nothing
     * @param  ThisMethod Input Integration method
     * @return            String with the name of the input integration method
     */
    std::string GetIntegrationName(Geometry<Node<3>>& geom, Geometry<Node<3>>::IntegrationMethod ThisMethod);

    /** Gets the corresponding string of the geometry name.
     * Gets the corresponding string of the geometry name.
     * @param  geom Input Geometry
     * @return      String corresponding to the name of the input geometry
     */
    std::string GetGeometryName(Geometry< Node<3> >& geom);

    /** Computes the linear strain matrix.
     * Computes the linear strain matrix which is useful to verify that
     * a constant strain can be correctly reproduced
     * @param B               [description]
     * @param DN_DX           [description]
     * @param number_of_nodes Nuber of nodes of the geometry
     * @param dimension       Dimension (i.e. 1, 2 or 3)
     */
    void CalculateB(Matrix& B, Matrix& DN_DX, const unsigned int number_of_nodes, const unsigned int dimension);

    /** Verifies the area of the geometry using the integration method.
     * Verifies the area of the geometry using the integration method.
     * @param  geom           Geometry to be tested
     * @param  ThisMethod     Integration method used
     * @param  reference_area Expected area
     * @param  error_msg      Buffer to write the error message
     * @return                Area claculated using the selected integration method.
     */
    double CalculateAreaByIntegration(Geometry<Node<3>>& geom, Geometry<Node<3> >::IntegrationMethod ThisMethod);
    /** Verifies that a displacement field produces the expected strain distribution.
     * Verifies that a displacement field which varies linearly in space, produces the expected strain distribution.
     * This shall be considered a test for shape function derivatives
     * @param geom       Geometry to be tested
     * @param ThisMethod Integration method used
     * @param error_msg  Buffer to write the error message
     */
    void VerifyStrainExactness(Geometry<Node<3>>& geom,  Geometry<Node<3> >::IntegrationMethod ThisMethod);
}
}

#endif // KRATOS_TEST_GEOMETRY_H_INCLUDED defined
