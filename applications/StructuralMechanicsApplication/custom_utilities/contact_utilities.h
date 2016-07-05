 
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
     * This
     * @param
     * @return
     */

    static inline void ContactContainerFiller(
            contact_container & contact_container,
            const Point<3>& ContactPoint,
            Geometry<Node<3> > & Geom1,
            Geometry<Node<3> > & Geom2,
            const array_1d<double, 3> & contact_normal1,
            const array_1d<double, 3> & contact_normal2,
            const IntegrationMethod IntegrationOrder
            )
    {        
        // Define the discrete contact gap
         Point<3> ProjectedPoint;
         const unsigned int number_nodes = Geom1.PointsNumber();
         const unsigned int dimension = Geom1.WorkingSpaceDimension();
         contact_container.contact_gap.resize(number_nodes);
         contact_container.active_nodes_slave.resize(number_nodes);

         for (unsigned int index = 0; index < number_nodes; index++)
         {
             Project(ContactPoint,  Geom1[index], ProjectedPoint, contact_container.contact_gap[index], contact_normal1);
             
             array_1d<double, 3> result;
             contact_container.active_nodes_slave[index] =  Geom2.IsInside(ProjectedPoint, result);
         }
         
         std::vector<bool> active_nodes_master;
         active_nodes_master.resize(Geom2.PointsNumber());

         double aux_value;
         for (unsigned int index = 0; index < Geom2.PointsNumber(); index++)
         {
             Project(Geom1[0],  Geom2[index], ProjectedPoint, aux_value, contact_normal2);
             
             array_1d<double, 3> result;
             active_nodes_master[index] = Geom1.IsInside(ProjectedPoint, result);
         }

//          KRATOS_WATCH("-----------------------------------------------------------------------------------------------------------------")
//          KRATOS_WATCH(Geom1);
//          KRATOS_WATCH(Geom2);
         
         /* Reading integration points slave condition */
         const GeometryType::IntegrationPointsArrayType& integration_points1 = Geom1.IntegrationPoints( IntegrationOrder );
         std::vector<bool> active_gauss_slave;
         active_gauss_slave.resize(integration_points1.size());
         
         for (unsigned int PointNumber = 0; PointNumber < integration_points1.size(); PointNumber++)
         {
             Point<3> GaussPoint;
             Point<3> GaussPointLocalCoordinates;
             Point<3> ProjectedGaussPoint;
             GaussPointLocalCoordinates.Coordinates() = integration_points1[PointNumber].Coordinates();

             array_1d<double, 3> result;
             GaussPoint = Geom1.GlobalCoordinates(result, GaussPointLocalCoordinates);

             double dist_aux;

             Project(ContactPoint, GaussPoint,  ProjectedGaussPoint, dist_aux, contact_normal1);

             active_gauss_slave[PointNumber] = Geom2.IsInside(ProjectedGaussPoint, result);
             
//              KRATOS_WATCH(inside);
//              KRATOS_WATCH(result);
//              KRATOS_WATCH(GaussPoint);
//              KRATOS_WATCH(ProjectedGaussPoint);
         }
         
         /* Reading integration points master condition */
         const GeometryType::IntegrationPointsArrayType& integration_points2 = Geom2.IntegrationPoints( IntegrationOrder );
         std::vector<bool> active_gauss_master;
         active_gauss_master.resize(integration_points2.size());
         
         for (unsigned int PointNumber = 0; PointNumber < integration_points2.size(); PointNumber++)
         {
             Point<3> GaussPoint;
             Point<3> GaussPointLocalCoordinates;
             Point<3> ProjectedGaussPoint;
             GaussPointLocalCoordinates.Coordinates() = integration_points2[PointNumber].Coordinates();

             array_1d<double, 3> result;
             GaussPoint = Geom2.GlobalCoordinates(result, GaussPointLocalCoordinates);

             double dist_aux;

             Project(ContactPoint, GaussPoint,  ProjectedGaussPoint, dist_aux, contact_normal1);

             active_gauss_master[PointNumber] =  Geom1.IsInside(ProjectedGaussPoint, result);
         }

         if (dimension == 2)
         {
             if (number_nodes == 2)
             {
                 contact_container.local_coordinates_slave.resize(2, false);
                 LocalLine2D2NProcess(contact_container.active_nodes_slave, active_gauss_slave, Geom1, contact_container.local_coordinates_slave[0], contact_container.local_coordinates_slave[1], IntegrationOrder);
                
                 contact_container.local_coordinates_master.resize(2, false);
                 LocalLine2D2NProcess(active_nodes_master, active_gauss_master, Geom2, contact_container.local_coordinates_master[0], contact_container.local_coordinates_master[1], IntegrationOrder);
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
         noalias(vector_points) = PointDestiny.Coordinates() - PointOrigin.Coordinates();

         dist = inner_prod(vector_points, Normal); 

         PointProjected.Coordinates() = PointDestiny.Coordinates() - Normal * dist;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * 
     * @param 
     * @return 
     */

    static inline void LocalLine2D2NProcess(
        const std::vector<bool> & active_nodes,
        const std::vector<bool> & active_gauss,
        Geometry<Node<3> > & Geom,
        Geometry<Point<3> > & GeomOut,
        double & coord1,
        double & coord2,
        const IntegrationMethod & IntegrationOrder
    )
    {   
        Point<3>::Pointer pPoint1, pPoint2;
        
        bool point1_assigned = false;
        bool point2_assigned = false;
        
        if (active_nodes[0] == true)
        {
            pPoint1->Coordinates() = Geom[0].Coordinates();
            coord1 = - 1.0;
            point1_assigned = true;
        }
        
        if (active_nodes[1] == true)
        {
            pPoint2->Coordinates() = Geom[1].Coordinates();
            coord2 = 1.0;
            point2_assigned = true;
        }
        
        if ((point1_assigned && point2_assigned) == false)
        {
            const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( IntegrationOrder );
            const unsigned int number_integration_points = integration_points.size();
            for (unsigned int PointNumber = 0; PointNumber < number_integration_points; PointNumber++)
            {
                Point<3> GaussPointLocalCoordinates;
                
                if (point1_assigned == false)
                {
                    if (active_gauss[PointNumber] == true)
                    {
                        GaussPointLocalCoordinates.Coordinates() = integration_points[PointNumber].Coordinates();

                        array_1d<double, 3> result;
                        *(pPoint1) = Geom.GlobalCoordinates(result, GaussPointLocalCoordinates);
                        
                        coord1 = GaussPointLocalCoordinates.Coordinate(1);
                        
                        point1_assigned = true;
                    }
                }
                
                if (point2_assigned == false)
                {
                    if (active_gauss[number_integration_points - PointNumber - 1] == true)
                    {
                        GaussPointLocalCoordinates.Coordinates() = integration_points[number_integration_points - PointNumber - 1].Coordinates();

                        array_1d<double, 3> result;
                        *(pPoint2) = Geom.GlobalCoordinates(result, GaussPointLocalCoordinates);
                         
                        coord2 = GaussPointLocalCoordinates.Coordinate(1);
                        
                        point2_assigned = true;
                    }
                }
                
                if ((point1_assigned && point2_assigned) == true)
                {
                    break;
                }
            }
        }

        GeomOut = Line2D2< Point<3> >(pPoint1, pPoint2);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * 
     * @param 
     * @return 
     */

    static inline void LocalLine2D2NProcess(
        const std::vector<bool> & active_nodes,
        const std::vector<bool> & active_gauss,
        Geometry<Node<3> > & Geom,
        double & coord1,
        double & coord2,
        const IntegrationMethod & IntegrationOrder
    )
    {   
        bool point1_assigned = false;
        bool point2_assigned = false;
        
        if (active_nodes[0] == true)
        {
            coord1 = - 1.0;
            point1_assigned = true;
        }
        
        if (active_nodes[1] == true)
        {
            coord2 = 1.0;
            point2_assigned = true;
        }
        
        if ((point1_assigned && point2_assigned) == false)
        {
            const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( IntegrationOrder );
            const unsigned int number_integration_points = integration_points.size();
            for (unsigned int PointNumber = 0; PointNumber < number_integration_points; PointNumber++)
            {
                Point<3> GaussPointLocalCoordinates;
                
                if (point1_assigned == false)
                {
                    if (active_gauss[PointNumber] == true)
                    {
                        coord1 = integration_points[PointNumber].Coordinate(1);
                        point1_assigned = true;
                    }
                }
                
                if (point2_assigned == false)
                {
                    if (active_gauss[PointNumber] == true)
                    {
                        coord2 = integration_points[number_integration_points - PointNumber - 1].Coordinate(1);
                        point2_assigned = true;
                    }
                }
                
                if ((point1_assigned && point2_assigned) == true)
                {
                    break;
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
        
        ConditionNormal(pCond, dimension);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * This function calculates the normal of a condition
     * @param Cond: The pointer to the condition of interest
     */

    static inline void ConditionNormal(
            Condition::Pointer pCond,
            const unsigned int dimension
            )
    {
        // TODO: Add calculation of normal to geometry.h 
        array_1d<double,3> & Normal = pCond->GetValue(NORMAL);
        noalias(Normal) = ZeroVector(3);

        // TODO: To calculate the normal I am going to use the Newell's method for quadrilateral, I recommend to find some way to compute in a way that it is possible to have the normal in the nodes and use the Nagata Patch
        if (pCond->GetGeometry().PointsNumber() == 2) // A linear line
        {
            array_1d<double,3> v1,v2;

            // Assuming plane X-Y
            noalias(v1) = ZeroVector(3);
            v1[2] = 1.0;

            noalias(v2) = pCond->GetGeometry()[1].Coordinates() - pCond->GetGeometry()[0].Coordinates();
            
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

                noalias(v1) = pCond->GetGeometry()[1].Coordinates() - pCond->GetGeometry()[0].Coordinates();

                noalias(v2) = pCond->GetGeometry()[2].Coordinates() - pCond->GetGeometry()[0].Coordinates();

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
        const unsigned int dimension = ModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
               
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
                noalias(node_it->FastGetSolutionStepValue(NORMAL)) = ZeroNormal;
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
                ConditionNormal(*(cond_it.base()), dimension);
                
                const array_1d<double, 3> & rNormal = cond_it->GetValue(NORMAL);
//                 KRATOS_WATCH(rNormal);
                for (unsigned int i = 0; i < cond_it->GetGeometry().PointsNumber(); i++)
                {
                    noalias( cond_it->GetGeometry()[i].FastGetSolutionStepValue(NORMAL) ) += rNormal;
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
                const double norm = norm_2(node_it->FastGetSolutionStepValue(NORMAL));
                if (norm > tol)
                {
                    node_it->FastGetSolutionStepValue(NORMAL)  /= norm;
                }
//                 KRATOS_WATCH(node_it->FastGetSolutionStepValue(NORMAL) );
            }
        }
      
    }

private:
};// class ContactUtilities
}
#endif /* KRATOS_CONTACT_UTILITIES defined */
 
