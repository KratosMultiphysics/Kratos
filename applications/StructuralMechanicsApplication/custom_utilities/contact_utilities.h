 
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
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"

namespace Kratos
{
class ContactUtilities
{
public:
    /**
     * @name Type definitions
     * @{
     */
    
    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef GeometryData::IntegrationMethod         IntegrationMethod;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * This function fills the contact_container for the Mortar condition
     * @param ContactPoint: The destination point
     * @param Geom1 and Geom2: The geometries of the slave and master respectively
     * @param contact_normal1 and contact_normal2: The normals of the slave and master respectively
     * @param IntegrationOrder: The integration order   
     * @return contact_container: Once has been filled
     */

    static inline void ContactContainerFiller(
            contact_container & contact_container,
            const Point<3>& ContactPoint,
            Geometry<Node<3> > & Geom1, // SLAVE
            Geometry<Node<3> > & Geom2, // MASTER
            const array_1d<double, 3> & contact_normal1, // SLAVE
            const array_1d<double, 3> & contact_normal2, // MASTER
            const double ActiveCheckFactor,
            const IntegrationMethod IntegrationOrder
            )
    {
        // Define the basic information
        const unsigned int number_nodes = Geom1.PointsNumber();
        const unsigned int dimension = Geom1.WorkingSpaceDimension();
        
         for (unsigned int index = 0; index < number_nodes; index++)
         {
             if (Geom1[index].Is(ACTIVE) == false)
             {
                Point<3> ProjectedPoint;
                double aux_dist = 0.0;
                if (norm_2(Geom1[index].FastGetSolutionStepValue(NORMAL, 0)) < 1.0e-12)
                {
                    ProjectDirection(Geom2, Geom1[index], ProjectedPoint, aux_dist, contact_normal1);
                }
                else
                {
                    ProjectDirection(Geom2, Geom1[index], ProjectedPoint, aux_dist, Geom1[index].FastGetSolutionStepValue(NORMAL, 0));
                }  
                
                double dist_tol = ActiveCheckFactor * Geom1.Length();
                dist_tol = (dist_tol <= ActiveCheckFactor * Geom2.Length()) ? (ActiveCheckFactor * Geom2.Length()):dist_tol;
                
                array_1d<double, 3> result;
                if (aux_dist < dist_tol)
                {
                    if (Geom2.IsInside(ProjectedPoint, result) == true)
                    {
                        Geom1[index].Set(ACTIVE, true);
                    }
                }
             }
         }

         if (dimension == 2)
         {
             if (number_nodes == 2)
             {
                 contact_container.local_coordinates_slave.clear();
                 contact_container.local_coordinates_slave.resize(2);
                 LocalLine2D2NProcess(contact_container.local_coordinates_slave, Geom1, Geom2, contact_normal1, contact_normal2, IntegrationOrder);  
             }
             else
             {
                 KRATOS_THROW_ERROR( std::logic_error, "NOT IMPLEMENTED. Number of nodes:",  number_nodes);
                 // TODO: IMPLEMENT MORE GEOMETRIES
             }
         }
         else
         {
             KRATOS_THROW_ERROR( std::logic_error, "NOT IMPLEMENTED. Dimension:",  dimension);
             // TODO: IMPLEMENT IN 3D
         }
         
//         contact_container.print();
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Project a point over a line/plane following an arbitrary direction
     * @param Geom: The geometry where to be projected
     * @param PointDestiny: The point to be projected
     * @param Vector: The direction to project
     * @return PointProjected: The point pojected over the plane
     * @return dist: The distance between the point and the plane
     */

    static inline void ProjectDirection(
            const Geometry<Node<3> > & Geom,
            const Point<3>& PointDestiny,
            Point<3>& PointProjected,
            double& dist,
            const array_1d<double,3>& Vector
            )
    {        
        const double tol = 1.0e-15;
        
        array_1d<double,3> Normal;
        
        GeometryNormal(Normal, Geom);
        
        array_1d<double,3> vector_points = Geom.Center() - PointDestiny.Coordinates();

        if (std::abs(inner_prod(Vector, Normal) ) > tol)
        {
            dist = inner_prod(vector_points, Normal)/inner_prod(Vector, Normal); 

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * dist;
        }
        else
        {
            PointProjected.Coordinates() = PointDestiny.Coordinates();
            dist = 0.0;
            std::cout << " The line and the plane are coplanar, something wrong happened " << std::endl;
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Project a point over a plane
     * @param PointOrigin: A point in the plane
     * @param PointDestiny: The point to be projected
     * @param Normal: The normal of the plane
     * @return PointProjected: The point pojected over the plane
     * @return dist: The distance between the point and the plane
     */
    
    static inline void Project(
            const Point<3>& PointOrigin,
            const Point<3>& PointDestiny,
            Point<3>& PointProjected,
            double& dist,
            const array_1d<double,3>& Normal
            )
    {
         array_1d<double,3> vector_points = PointDestiny.Coordinates() - PointOrigin.Coordinates();

         dist = inner_prod(vector_points, Normal); 

         PointProjected.Coordinates() = PointDestiny.Coordinates() - Normal * dist;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * calculates the distance between nodes
     * @param PointOrigin: A point in the plane
     * @param PointDestiny: The point to be projected
     */
    
    static inline double DistancePoints(
            const Point<3>& PointOrigin,
            const Point<3>& PointDestiny
            )
    {
        double dist = std::sqrt((PointOrigin.Coordinate(1) - PointDestiny.Coordinate(1)) * (PointOrigin.Coordinate(1) - PointDestiny.Coordinate(1))
                              + (PointOrigin.Coordinate(2) - PointDestiny.Coordinate(2)) * (PointOrigin.Coordinate(2) - PointDestiny.Coordinate(2))
                              + (PointOrigin.Coordinate(3) - PointDestiny.Coordinate(3)) * (PointOrigin.Coordinate(3) - PointDestiny.Coordinate(3)));
        
        return dist;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This function calculates the local coordinates of the projected line for the mortar condition
     * @param Geom1 and Geom2: The geometries of the slave and master respectively
     * @param contact_normal1 and contact_normal2: The normals of the slave and master respectively
     * @param IntegrationOrder: The integration order   
     * @return local_coordinates_slave and local_coordinates_master: The local coordinates of the pojected points of the line
     */
    
    static inline void LocalLine2D2NProcess(
        std::vector<double> & local_coordinates_slave,
        Geometry<Node<3> > & Geom1, // SLAVE
        Geometry<Node<3> > & Geom2, // MASTER
        const array_1d<double, 3> & contact_normal1, // SLAVE
        const array_1d<double, 3> & contact_normal2, // MASTER
        const IntegrationMethod & IntegrationOrder
    )
    {   
        // Define auxiliar values
        Point<3> ProjectedPoint;
        
        // Domain 1
        for (unsigned int index_1 = 0; index_1 < Geom1.PointsNumber(); index_1++)
        {
             double aux_dist = 0.0;
             if (norm_2(Geom1[index_1].FastGetSolutionStepValue(NORMAL, 0)) < 1.0e-12)
             {
                 ProjectDirection(Geom2, Geom1[index_1], ProjectedPoint, aux_dist, contact_normal1);
             }
             else
             {
                 ProjectDirection(Geom2, Geom1[index_1], ProjectedPoint, aux_dist, Geom1[index_1].FastGetSolutionStepValue(NORMAL, 0));
             }  

             array_1d<double, 3> projected_local_coor_1;
             const bool in_out_1 = Geom2.IsInside(ProjectedPoint, projected_local_coor_1);
             array_1d<double, 3> local_coor_1 = Geom1.PointLocalCoordinates(local_coor_1, Geom1[index_1]);
                              
             if (in_out_1 == true)
             {
                 local_coordinates_slave[index_1]  = local_coor_1[0];
             }
             else
             {
                 // Domain 2
                 double proj_dist = 2.0;
                 for (unsigned int index_2 = 0; index_2 < Geom2.PointsNumber(); index_2++)
                 {
                     ProjectDirection(Geom1,  Geom2[index_2], ProjectedPoint, aux_dist, - contact_normal1);

                     array_1d<double, 3> projected_local_coor_2;
                     const bool in_out_2 = Geom1.IsInside(ProjectedPoint, projected_local_coor_2);

                     if (in_out_2 == true)
                     {
                         if( std::abs( projected_local_coor_2[0] - local_coor_1[0] ) <= proj_dist )
                         {
                             proj_dist = std::abs(projected_local_coor_2[0] - local_coor_1[0] );
                             local_coordinates_slave[index_1]  = projected_local_coor_2[0];
                         }
                     }
                }
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * This function calculates the center and radius of the geometry of a condition
     * @param Cond: The pointer to the condition of interest
     * @return Center: The center of the condition
     * @return Radius: The radius of the condition
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
        
        // TODO: Add calculation of radius to geometry.h 
        array_1d<double, 3> aux_vector;
        for(unsigned int i = 0; i < pCond->GetGeometry().PointsNumber(); i++)
        {
            noalias(aux_vector) = Center.Coordinates() - pCond->GetGeometry()[i].Coordinates();;
            
            double tmp = inner_prod(aux_vector, aux_vector);

            if(tmp > Radius)
            {
                Radius = tmp;
            }
        }

        Radius = std::sqrt(Radius);
        
        ConditionNormal(pCond);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * This function calculates the normal of a condition
     * @param Cond: The pointer to the condition of interest
     */

    static inline void ConditionNormal(Condition::Pointer pCond)
    {
        // TODO: Add calculation of normal to geometry.h 
        array_1d<double,3> & Normal = pCond->GetValue(NORMAL);
        
        GeometryNormal(Normal, pCond->GetGeometry());
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This function calculates the normal of a geometry
     * @param Cond: The pointer to the condition of interest
     */

    static inline void GeometryNormal(
            array_1d<double,3> & Normal,
            const Geometry<Node<3> > & Geom
            )
    {
        const unsigned int dimension = Geom.WorkingSpaceDimension();
        
        noalias(Normal) = ZeroVector(3);
        
        // TODO: To calculate the normal I am going to use the Newell's method for quadrilateral, I recommend to find some way to compute in a way that it is possible to have the normal in the nodes and use the Nagata Patch
        if (Geom.PointsNumber() == 2) // A linear line
        {
            array_1d<double,3> v1,v2;

            // Assuming plane X-Y
            noalias(v1) = Geom[1].Coordinates() - Geom[0].Coordinates();
            
            noalias(v2) = ZeroVector(3);
            v2[2] = 1.0;
            
            MathUtils<double>::CrossProduct(Normal,v1,v2);

            double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1]);
            Normal /= NNorm;
        }
        else if (Geom.PointsNumber() == 3) // A triangle or quadratic line
        {
            if (dimension == 2) // Quadratic line
            {
                boost::numeric::ublas::bounded_matrix<double, 3, 3 >  matrix_coeficients     = ZeroMatrix(3, 3);
                boost::numeric::ublas::bounded_matrix<double, 3, 3 >  inv_matrix_coeficients = ZeroMatrix(3, 3);
                boost::numeric::ublas::bounded_matrix<double, 3, 1 >  vector_coeficients     = ZeroMatrix(3, 1);
                for (unsigned int i = 0; i < 3; i++)
                {
                    matrix_coeficients(i, 0) = Geom[i].X() * Geom[i].X();
                    matrix_coeficients(i, 1) = Geom[i].X();
                    matrix_coeficients(i, 2) = 1.0;

                    vector_coeficients(i, 0) = Geom[i].Y();
                }

                StructuralMechanicsMathUtilities::InvMat3x3(matrix_coeficients, inv_matrix_coeficients);

                boost::numeric::ublas::bounded_matrix<double, 3, 1 >  coeficients;
                noalias(coeficients) = prod(inv_matrix_coeficients, vector_coeficients);

                Normal[0] =   2.0 * Geom[1].X() * coeficients(0, 0) + coeficients(1, 0);
                Normal[1] = - 1.0;
                Normal[2] =   0.0;

                double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1]);
                Normal /= NNorm;
            }
            else // Triangle
            {
                array_1d<double,3> v1,v2;

                noalias(v1) = Geom[1].Coordinates() - Geom[0].Coordinates();

                noalias(v2) = Geom[2].Coordinates() - Geom[0].Coordinates();

                MathUtils<double>::CrossProduct(Normal,v1,v2);

                double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1] + Normal[2] * Normal[2]);
                Normal /= NNorm;
            }
        }
        else if (Geom.PointsNumber() == 4) // A quadrilateral
        {
            // Newell's method
            for(unsigned int i = 0; i < Geom.PointsNumber(); i++)
            {
                unsigned int index_aux = i + 1;
                if (i == Geom.PointsNumber() - 1)
                {
                    index_aux = 0;
                }
                
                Normal[0] += (Geom[i].Y() - Geom[index_aux].Y()) *
                             (Geom[i].Z() - Geom[index_aux].Z());

                Normal[1] += (Geom[i].Z() - Geom[index_aux].Z()) *
                             (Geom[i].X() - Geom[index_aux].X());

                Normal[2] += (Geom[i].X() - Geom[index_aux].X()) *
                             (Geom[i].Y() - Geom[index_aux].Y());
            }

            double NNorm = std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1] + Normal[2] * Normal[2]);
            Normal /= NNorm;
        }
        else // The Newell's method can be used, but nodes must be reordered
        {
            KRATOS_THROW_ERROR( std::logic_error, " There is not any method to calculate the normal for this geometry. Number of nodes: ", Geom.PointsNumber() );
        }
        
//         KRATOS_WATCH(Normal);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * It computes the mean of the normal in the condition in all the nodes
     * @param ModelPart: The model part to compute
     * @return The modelparts with the normal computed
     */
    
    static inline void ComputeNodesMeanNormalModelPart(ModelPart & ModelPart)
    {
        // Tolerance
        const double tol = 1.0e-14;
        
        // Initialize normal vectors
        const array_1d<double,3> ZeroNormal = ZeroVector(3);
        
        NodesArrayType& pNode  = ModelPart.Nodes();
        NodesArrayType::iterator it_node_begin = pNode.ptr_begin();
        NodesArrayType::iterator it_node_end   = pNode.ptr_end();
        
        for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
        {
            if (node_it->Is(INTERFACE))
            {
                noalias(node_it->GetValue(NORMAL)) = ZeroNormal;
            }
        }
        
        // Sum all the nodes normals
        ConditionsArrayType& pCond  = ModelPart.Conditions();
        ConditionsArrayType::iterator it_cond_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_cond_end   = pCond.ptr_end();
        
        for(ConditionsArrayType::iterator cond_it = it_cond_begin; cond_it!=it_cond_end; cond_it++)
        {
            if (cond_it->Is(ACTIVE) || cond_it->Is(MASTER))
            {
                ConditionNormal(*(cond_it.base()));
                
                const array_1d<double, 3> & rNormal = cond_it->GetValue(NORMAL);
                
                for (unsigned int i = 0; i < cond_it->GetGeometry().PointsNumber(); i++)
                {
                    noalias( cond_it->GetGeometry()[i].GetValue(NORMAL) ) += rNormal;
                }
            }
        }
        
        // Normalize normal vectors
        it_node_begin = pNode.ptr_begin();
        it_node_end   = pNode.ptr_end();
        
        for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
        {
            if (node_it->Is(INTERFACE))
            {
                const double norm = norm_2(node_it->GetValue(NORMAL));
                
                if (norm > tol)
                {
                    node_it->FastGetSolutionStepValue(NORMAL)  /= norm;
                }
            }
        }
      
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * It changes from active to inactive and viceversa the nodes 
     * @param ModelPart: The model part to compute
     * @param cn: Kind of penalty, not necessarily the YOUNG_MODULUS
     * @return The modelparts with the conditions changed
     */
    
    static inline void ReComputeActiveInactive(
        ModelPart & ModelPart, 
        const double cn
        )
    {
        NodesArrayType& pNode  = ModelPart.Nodes();
        NodesArrayType::iterator it_node_begin = pNode.ptr_begin();
        NodesArrayType::iterator it_node_end   = pNode.ptr_end();
        
        for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
        {
            const array_1d<double,3> lagrange_multiplier = node_it->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER, 0);
            const array_1d<double,3>        nodal_normal = node_it->GetValue(NORMAL); 
            const double lambda_n = inner_prod(lagrange_multiplier, nodal_normal);
            
            const double check = lambda_n - cn * node_it->GetValue(WEIGHTED_GAP); 
            
            if (check >= 0.0)
            {
                node_it->Set(ACTIVE, false);
            }
            else
            {
                node_it->Set(ACTIVE, true);
            }
        }
        
    }

private:
};// class ContactUtilities
}
#endif /* KRATOS_CONTACT_UTILITIES defined */
 
