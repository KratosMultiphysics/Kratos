 
// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_CONTACT_UTILITIES)
#define KRATOS_CONTACT_UTILITIES

#include "utilities/math_utils.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "geometries/point.h"

namespace Kratos
{
class ContactUtilities
{
public:
    /**
     * @name Type definitions
     * @{
     */

    typedef Node<3>                                      NodeType;
    typedef Geometry<NodeType>                       GeometryType;
    typedef GeometryData::IntegrationMethod     IntegrationMethod;

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * This
     * @param
     * @return
     */

    static inline void ContactContainerFiller(
            contact_container & contact_container,
            const array_1d<double, 3> & contact_normal,
            const Point<3>& ContactPoint,
            Geometry<Node<3> > & Geom1,
            Geometry<Node<3> > & Geom2,
            const IntegrationMethod IntegrationOrder
            )
    {
        // Define the discrete contact gap
         Point<3> ProjectedPoint;
         const unsigned int number_nodes = Geom1.PointsNumber();
         contact_container.contact_gap.resize(number_nodes);
         contact_container.active_nodes.resize(number_nodes);

         for (unsigned int index = 0; index < number_nodes; index++)
         {
             ContactUtilities::Project(ContactPoint,  Geom1[index], ProjectedPoint, contact_container.contact_gap[index], contact_normal);

             array_1d<double, 3> result;
             bool inside = Geom2.IsInside(ProjectedPoint, result);
             if (inside == true)
             {
                 contact_container.active_nodes[index] = true;
             }
             else
             {
                 contact_container.active_nodes[index] = false;
             }
         }

         // Define the contact area
         contact_container.contact_area = 0.0;

     //     KRATOS_WATCH("-----------------------------------------------------------------------------------------------------------------")
     //     KRATOS_WATCH(Geom1);
     //     KRATOS_WATCH(Geom2);

         double aux_int = 0.0;
         /* Reading integration points */
         const GeometryType::IntegrationPointsArrayType& integration_points = Geom1.IntegrationPoints( IntegrationOrder );
         for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
         {
             Point<3> GaussPoint;
             Point<3> GaussPointLocalCoordinates;
             Point<3> ProjectedGaussPoint;
             GaussPointLocalCoordinates.Coordinate(1) = integration_points[PointNumber].X();
             GaussPointLocalCoordinates.Coordinate(2) = integration_points[PointNumber].Y();
             GaussPointLocalCoordinates.Coordinate(3) = integration_points[PointNumber].Z(); // This is supposed to be 0 always, in 1D and 2D

             array_1d<double, 3> result;
             GaussPoint = Geom1.GlobalCoordinates(result, GaussPointLocalCoordinates);

             double dist_aux;

             ContactUtilities::Project(ContactPoint, GaussPoint,  ProjectedGaussPoint, dist_aux, contact_normal);

             bool inside = Geom2.IsInside(ProjectedGaussPoint, result);

     //         KRATOS_WATCH(inside);
     //         KRATOS_WATCH(result);
     //         KRATOS_WATCH(GaussPoint);
     //         KRATOS_WATCH(ProjectedGaussPoint);

             // Integration weigth
             double IntegrationWeight = integration_points[PointNumber].Weight();
             aux_int += IntegrationWeight;

             if (inside == true)
             {
                 contact_container.contact_area += IntegrationWeight;
             }
         }

         contact_container.contact_area /= aux_int;

//         contact_container.print();
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Project
     * @param
     * @return
     */

    static inline void Project(
            const Point<3>& PointOrigin,
            const Point<3>& PointDestiny,
            Point<3>& PointProjected,
            double& dist,
            const array_1d<double,3>& Normal
            )
    {
         array_1d<double,3> vector_points;
         vector_points[0] = PointDestiny.Coordinate(1) - PointOrigin.Coordinate(1);
         vector_points[1] = PointDestiny.Coordinate(2) - PointOrigin.Coordinate(2);
         vector_points[2] = PointDestiny.Coordinate(3) - PointOrigin.Coordinate(3);

         dist = vector_points[0] * Normal[0]
              + vector_points[1] * Normal[1]
              + vector_points[2] * Normal[2];

         PointProjected.Coordinate(1) = PointDestiny.Coordinate(1) - Normal[0] * dist;
         PointProjected.Coordinate(2) = PointDestiny.Coordinate(2) - Normal[1] * dist;
         PointProjected.Coordinate(3) = PointDestiny.Coordinate(3) - Normal[2] * dist;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * This function calculates the center and radius of the geometry of a condition
     * @param Cond: The pointer to the condition of interest
     * @return Center: The center of the condition
     * @return Radius: The radius of the condition
     * @return Normal: The normal of the condition
     */

    static inline void CenterAndRadius(
            Condition::Pointer pCond,
            Point<3>& Center,
            double& Radius,
            const unsigned int dimension
            )
    {
        Radius = 0.0;
        Center = pCond->GetGeometry().Center();
        array_1d<double,3> & Normal = pCond->GetValue(NORMAL);
        noalias(Normal) = ZeroVector(3);

        // TODO: To calculate the normal I am going to use the Newell's method for quadrilateral, I recommend to find some way to compute in a way that it is possible to have the normal in the nodes and use the Nagata Patch
        if (pCond->GetGeometry().PointsNumber() == 2) // A linear line
        {
            array_1d<double,3> v1,v2;

            // Assuming plane X-Y
            v1[0] = 0.0;
            v1[1] = 0.0;
            v1[2] = 1.0;

            v2[0] = pCond->GetGeometry()[1].X() - pCond->GetGeometry()[0].X();
            v2[1] = pCond->GetGeometry()[1].Y() - pCond->GetGeometry()[0].Y();
            v2[2] = pCond->GetGeometry()[1].Z() - pCond->GetGeometry()[0].Z();

            MathUtils<double>::CrossProduct(Normal,v1,v2);

            double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1]);
            Normal /= NNorm;
        }
        else if (pCond->GetGeometry().PointsNumber() == 3) // A triangle or quadratic line
        {
            if (dimension == 2) // Quadratic line
            {
                boost::numeric::ublas::bounded_matrix<double, 3, 3 >  matrix_coeficients     = ZeroMatrix(3, 3);
                boost::numeric::ublas::bounded_matrix<double, 3, 3 >  inv_matrix_coeficients = ZeroMatrix(3, 3);
                boost::numeric::ublas::bounded_matrix<double, 3, 1 >  vector_coeficients     = ZeroMatrix(3, 1);
                for (unsigned int i = 0; i < 3; i++)
                {
                    matrix_coeficients(i, 0) = pCond->GetGeometry()[i].X() * pCond->GetGeometry()[i].X();
                    matrix_coeficients(i, 1) = pCond->GetGeometry()[i].X();
                    matrix_coeficients(i, 2) = 1.0;

                    vector_coeficients(i, 0) = pCond->GetGeometry()[i].Y();
                }

                StructuralMechanicsMathUtilities::InvMat3x3(matrix_coeficients, inv_matrix_coeficients);

                boost::numeric::ublas::bounded_matrix<double, 3, 1 >  coeficients;
                noalias(coeficients) = prod(inv_matrix_coeficients, vector_coeficients);

                Normal[0] =   2.0 * pCond->GetGeometry()[1].X() * coeficients(0, 0) + coeficients(1, 0);
                Normal[1] = - 1.0;
                Normal[2] =   0.0;

                double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1]);
                Normal /= NNorm;
            }
            else // Triangle
            {
                array_1d<double,3> v1,v2;

                v1[0] = pCond->GetGeometry()[1].X() - pCond->GetGeometry()[0].X();
                v1[1] = pCond->GetGeometry()[1].Y() - pCond->GetGeometry()[0].Y();
                v1[2] = pCond->GetGeometry()[1].Z() - pCond->GetGeometry()[0].Z();

                v2[0] = pCond->GetGeometry()[2].X() - pCond->GetGeometry()[0].X();
                v2[1] = pCond->GetGeometry()[2].Y() - pCond->GetGeometry()[0].Y();
                v2[2] = pCond->GetGeometry()[2].Z() - pCond->GetGeometry()[0].Z();

                MathUtils<double>::CrossProduct(Normal,v1,v2);

                double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1] + Normal[2] * Normal[2]);
                Normal /= NNorm;
            }
        }
        else if (pCond->GetGeometry().PointsNumber() == 4) // A quadrilateral
        {
            // Newell's method
            for(unsigned int i = 0; i < pCond->GetGeometry().PointsNumber(); i++)
            {
                unsigned int index_aux = i + 1;
                if (i == pCond->GetGeometry().PointsNumber() - 1)
                {
                    index_aux = 0;
                }
                Normal[0] += (pCond->GetGeometry()[i].Y() - pCond->GetGeometry()[index_aux].Y()) *
                             (pCond->GetGeometry()[i].Z() - pCond->GetGeometry()[index_aux].Z());

                Normal[1] += (pCond->GetGeometry()[i].Z() - pCond->GetGeometry()[index_aux].Z()) *
                             (pCond->GetGeometry()[i].X() - pCond->GetGeometry()[index_aux].X());

                Normal[2] += (pCond->GetGeometry()[i].X() - pCond->GetGeometry()[index_aux].X()) *
                             (pCond->GetGeometry()[i].Y() - pCond->GetGeometry()[index_aux].Y());
            }

            double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1] + Normal[2] * Normal[2]);
            Normal /= NNorm;
        }
        else // The Newell's method can be used, but nodes must be reordered
        {
            KRATOS_THROW_ERROR( std::logic_error, " There is not any method to calculate the normal for this geometry. Number of nodes: ", pCond->GetGeometry().PointsNumber() );
        }

        // TODO: Add calculation of radius to geometry.h in june (after release)
        // TODO: Add calculation of normal to geometry.h in june (after release)
        for(unsigned int i = 0; i < pCond->GetGeometry().PointsNumber(); i++)
        {
            double dx = Center.Coordinate(1) - pCond->GetGeometry()[i].X();
            double dy = Center.Coordinate(2) - pCond->GetGeometry()[i].Y();
            double dz = Center.Coordinate(3) - pCond->GetGeometry()[i].Z();

            double tmp = dx * dx + dy * dy + dz * dz;

            if(tmp > Radius)
            {
                Radius = tmp;
            }
        }

        Radius = std::sqrt(Radius);
    }
private:
};// class ContactUtilities
}
#endif /* KRATOS_CONTACT_UTILITIES defined */
 
