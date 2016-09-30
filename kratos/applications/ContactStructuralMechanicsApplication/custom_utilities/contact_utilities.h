 
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
#include "contact_structural_mechanics_application_variables.h"
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
            const double ActiveCheckFactor
            )
    {
        // The original code was split into two methods to be more handy.. nothing actually is changed
        InitializeActiveInactiveSets( Geom1, Geom2, contact_normal1, contact_normal2, ActiveCheckFactor );
        GenerateMortarSegmentsProcess( contact_container, Geom1, Geom2, contact_normal1, contact_normal2 );
    }
    
    static inline void ContactContainerFiller(
            contact_container & contact_container,
            const Point<3>& ContactPoint,
            Condition::Pointer & pCond_1,       // SLAVE
            const Condition::Pointer & pCond_2, // MASTER
            const array_1d<double, 3> & contact_normal1, // SLAVE
            const array_1d<double, 3> & contact_normal2, // MASTER
            const double ActiveCheckFactor
            )
    {

        ContactContainerFiller(contact_container, ContactPoint, pCond_1->GetGeometry(), pCond_2->GetGeometry(), contact_normal1, contact_normal2, ActiveCheckFactor);
        
        Geometry<Node<3> > & Geom1 = pCond_1->GetGeometry();
        const unsigned int number_nodes = Geom1.PointsNumber();
        
        for (unsigned int index = 0; index < number_nodes; index++)
        {
            if (Geom1[index].Is(ACTIVE) == true)
            {
                pCond_1->Set(ACTIVE, true);
                break;
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline void InitializeActiveInactiveSets(
        Geometry<Node<3> > & Geom1, // SLAVE
        Geometry<Node<3> > & Geom2, // MASTER
        const array_1d<double, 3> & contact_normal1, // SLAVE
        const array_1d<double, 3> & contact_normal2, // MASTER
        const double ActiveCheckFactor
        )
    {
        // Define the basic information
        const unsigned int number_nodes = Geom1.PointsNumber();
        
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
                
                double dist_tol = ActiveCheckFactor;    // the actual gap tolerance is user-define instead of being a factor of the length
//                double dist_tol = ActiveCheckFactor * Geom1.Length();
//                dist_tol = (dist_tol <= ActiveCheckFactor * Geom2.Length()) ? (ActiveCheckFactor * Geom2.Length()):dist_tol;
                
                array_1d<double, 3> result;
                // NOTE: We don't use std::abs() because if the aux_dist is negative is penetrating, in fact we just consider dist_tol > 0 to have some tolerance and for the static schemes
                if (aux_dist <= dist_tol && Geom2.IsInside(ProjectedPoint, result))
                {
//                    // For debug purpose // NOTE: Look for using echo_level
//                    if (aux_dist < 0.0)
//                    {
//                        std::cout << "Penetration in node: " << Geom1[index].Id() << " of " << aux_dist << " m" << std::endl;
//                    }    
                    Geom1[index].Set(ACTIVE, true);
                    Geom1[index].GetSolutionStepValue( IS_ACTIVE_SET ) = true;
                }
                else  
                {
                }
             }
         }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    static inline void GenerateMortarSegmentsProcess(
        contact_container & contact_container,
        Geometry<Node<3> > & Geom1, // SLAVE
        Geometry<Node<3> > & Geom2, // MASTER
        const array_1d<double, 3> & contact_normal1, // SLAVE
        const array_1d<double, 3> & contact_normal2 // MASTER
        )
    {
        // Define the basic information
        const unsigned int number_nodes = Geom1.PointsNumber();
        const unsigned int dimension = Geom1.WorkingSpaceDimension();
        
        if (dimension == 2)
        {
            if (number_nodes == 2)
            {
                contact_container.local_coordinates_slave.clear();
                contact_container.local_coordinates_slave.resize(2);
                LocalLine2D2NProcess(contact_container.local_coordinates_slave, Geom1, Geom2, contact_normal1, contact_normal2);  
            }
            else
            {
                KRATOS_THROW_ERROR( std::logic_error, "NOT IMPLEMENTED. Number of nodes:",  number_nodes);
                // TODO: IMPLEMENT MORE GEOMETRIES
            }
        }
        else   // In 3D, we won't use mortar segments. We use colocation integration instead
        {
            contact_container.local_coordinates_slave.clear();
            contact_container.local_coordinates_slave.resize(0);
        }
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

        if( norm_2( Vector ) <= tol && norm_2( Normal ) >= tol )
        {
            dist = inner_prod(vector_points, Normal)/norm_2(Normal);

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * dist;
            std::cout << " :: Warning: Zero projection vector. Projection using the condition vector instead." << std::endl;
        }
        else if (std::abs(inner_prod(Vector, Normal) ) >= tol)
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
    
    static inline void ProjectCoordDirection(
            const Geometry<Node<3> > & Geom,
            const GeometryType::CoordinatesArrayType& CoordDestiny,
            GeometryType::CoordinatesArrayType& CoordProjected,
            double& dist,
            const array_1d<double,3>& Vector
            )
    {        
        const double tol = 1.0e-15;
        
        array_1d<double,3> Normal;
        
        GeometryNormal(Normal, Geom);
        
        array_1d<double,3> vector_points = Geom.Center() - CoordDestiny;

        if( norm_2( Vector ) <= tol && norm_2( Normal ) >= tol )
        {
            dist = inner_prod(vector_points, Normal)/norm_2(Normal);

            CoordProjected = CoordDestiny + Vector * dist;
            std::cout << " :: Warning: Zero projection vector. Projection using the condition normal vector instead." << std::endl;
        }
        else if (std::abs(inner_prod(Vector, Normal) ) >= tol)
        {
            dist = inner_prod(vector_points, Normal)/inner_prod(Vector, Normal); 

            CoordProjected = CoordDestiny + Vector * dist;
        }
        else
        {
            CoordProjected = CoordDestiny;
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
        const array_1d<double, 3> & contact_normal2 // MASTER
    )
    {   
        // Define auxiliar values
        Point<3> ProjectedPoint;
        
        // Domain 1
        for (unsigned int index_1 = 0; index_1 < Geom1.PointsNumber(); index_1++)
        {
             double aux_dist = 0.0;
             if (norm_2(Geom1[index_1].GetValue(NORMAL)) < 1.0e-12)
             {
                 ProjectDirection(Geom2, Geom1[index_1], ProjectedPoint, aux_dist, contact_normal1);
             }
             else
             {
                 ProjectDirection(Geom2, Geom1[index_1], ProjectedPoint, aux_dist, Geom1[index_1].GetValue(NORMAL));
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
     * This function calculates the normal of a condition
     * @param Cond: The pointer to the condition of interest
     */

    static inline array_1d<double,3> GaussPointNormal(
        const Vector N,
        const Geometry<Node<3> > & Geom
        )
    {
        array_1d<double,3> normal = ZeroVector(3);
        for( unsigned int iNode = 0; iNode < Geom.PointsNumber(); ++iNode )
        {
            normal += N[iNode] * Geom[iNode].GetValue(NORMAL); // The opposite direction
        }
        
        normal = normal/norm_2(normal); // It is suppossed to be already unitary (just in case)
        
        return normal;
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
        array_1d<double, 3> v1, v2;
        noalias(Normal) = ZeroVector(3);
        
        // Geom normal is the sum of all nodal normals
        // nodal normal = tangent_eta (or e3 in 2D) x tangent_xi
        for ( unsigned int i = 0; i < Geom.PointsNumber( ); ++i )
        {
            NodalTangents( v1, v2, Geom, i );
            Normal += MathUtils<double>::CrossProduct( v2, v1 );
        }
        
        Normal /= norm_2( Normal );
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This function calculates the tangents at a given node
     */

    static inline void NodalTangents(
            array_1d<double,3> & t1,    // tangent in xi direction
            array_1d<double,3> & t2,    // tangent in eta direction in 3D or simply e3 cartesian vector in 2D
            const Geometry<Node<3> > & Geom,
            const unsigned int i_node
            )
    {
        const unsigned int dimension = Geom.WorkingSpaceDimension( );
        const unsigned int local_dim = Geom.LocalSpaceDimension( );
        Matrix J_i = ZeroMatrix( dimension, local_dim ); 
        Geom.Jacobian( J_i, Geom[i_node].Coordinates( ) );

        t1 = column(J_i, 0);
        
        if( dimension == 2 )
        {
            t2 = ZeroVector(3);
            t2[2] = 1.0;
            t1[2] = 0.0;
        }
        else if( dimension == 3 )
        {
            t2 = column(J_i, 1);
        }
        else
        {
            std::cout << "\033[31m" << "Error in: " << __PRETTY_FUNCTION__ << "\033[0m" << std::endl;
            KRATOS_THROW_ERROR( std::logic_error, "Can't calculate nodal tangents. Bad dimension provided. Dimension = ", dimension );
        }
        
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
//         const unsigned int dimension = ModelPart.ConditionsBegin()->WorkingSpaceDimension();
        
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
        
        // Applied laziness - MUST be calculated BEFORE normalizing the normals
        ComputeDeltaNodesMeanNormalModelPart( ModelPart );
        
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
                    node_it->GetValue(NORMAL)  /= norm;
                }
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * It computes the directional derivative of the normal in the condition in all the nodes
     * @param ModelPart: The model part to compute
     * @return The modelparts with the normal computed
     */

    static inline void ComputeDeltaNodesMeanNormalModelPart(ModelPart & ModelPart)
    {
        // Conditions
        ConditionsArrayType& pCond  = ModelPart.Conditions();
        ConditionsArrayType::iterator it_cond_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_cond_end   = pCond.ptr_end();

        // Nodes
        NodesArrayType& pNode  = ModelPart.Nodes();
        NodesArrayType::iterator it_node_begin = pNode.ptr_begin();
        NodesArrayType::iterator it_node_end   = pNode.ptr_end();

        // Tolerance
        const double tol = 1.0e-14;

        // Initialize directional derivative
        const unsigned int dimension = it_cond_begin->WorkingSpaceDimension( ); 
        const Matrix ZeroDeltaNormal = ZeroMatrix( dimension, dimension );

        Matrix Delta_ne_adj  = Matrix(dimension, dimension);
        Matrix Ce = Matrix(dimension, dimension);
        
        const Matrix I = IdentityMatrix(dimension, dimension);

        // Initialize the normal and tangential directional derivatives
        for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
        {
            if (node_it->Is(INTERFACE))
            {
                node_it->GetValue(DELTA_NORMAL) = ZeroDeltaNormal;
            }
        }

        // Sum the directional derivatives of the adjacent segments
        for(ConditionsArrayType::iterator cond_it = it_cond_begin; cond_it!=it_cond_end; cond_it++)
        {
            if (cond_it->Is(ACTIVE) || cond_it->Is(MASTER))
            {
                const array_1d<double, 3> ne = cond_it->GetValue(NORMAL);   // normalized condition normal
                Matrix ne_o_ne = subrange( outer_prod( ne, ne ), 0, dimension, 0, dimension );

                const unsigned int num_nodes = cond_it->GetGeometry( ).PointsNumber();
                if( dimension == 2 )
                {
                    const double ne_norm = cond_it->GetGeometry( ).Length( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                    
                    Delta_ne_adj( 0, 0 ) =  0.0;
                    Delta_ne_adj( 0, 1 ) = -1.0;
                    Delta_ne_adj( 1, 0 ) =  1.0;
                    Delta_ne_adj( 1, 1 ) =  0.0;
                    
                    Ce = prod( I - ne_o_ne, Delta_ne_adj ) / ne_norm;     // In 2D, Delta_ne_adj is node-independent => evaluated outside the nodes loop
                    
                    for (unsigned int i = 0; i < num_nodes; i++)
                    {
                        NodeType& node_j = cond_it->GetGeometry( )[i];
                        
                        // -/+ 0.5 are the values of DN_Dxi for linear line elements at nodes 1 and 2 - no need to call the function
                        double DN_De_j = 0.0;
                        if( i == 0 )
                        {
                            DN_De_j = -0.5;
                        }
                        else if( i == 1 )
                        {
                            DN_De_j =  0.5;
                        }
                        else
                        {
                            KRATOS_THROW_ERROR( std::logic_error, "DELTA_NORMAL is not yet defined for higher order 1D elements. Number of nodes: ", num_nodes );
                        }
                        
                        node_j.GetValue(DELTA_NORMAL) += Ce * DN_De_j;
                    }
                }
                else if ( dimension == 3 )
                {
                    const double ne_norm = cond_it->GetGeometry( ).Area( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                    
                    for (unsigned int i = 0; i < num_nodes; i++)
                    {
                        Matrix J = ZeroMatrix( 3, 2 ); // Jacobian [ 3D global x 2D local ]
                        array_1d<double, 2> DN_De_j        = ZeroVector( 2 );  // Shape functions derivatives for node j [ DN_Dxi_j, DN_Deta_j ]
                        array_1d<double, 3> local_coords_j = ZeroVector( 3 );
                        
                        NodeType& node_j = cond_it->GetGeometry( )[i];
                        
                        if( num_nodes == 3 )    // linear triangle element
                        {
                            if( i == 0 )
                            {
                                local_coords_j[0] = 0.0;
                                local_coords_j[1] = 0.0;
                                DN_De_j( 0 ) = -1.0;
                                DN_De_j( 1 ) = -1.0;
                            }
                            else if( i == 1 )
                            {
                                local_coords_j[0] = 1.0;
                                local_coords_j[1] = 0.0;
                                DN_De_j( 0 ) = 1.0;
                                DN_De_j( 1 ) = 0.0;
                            }
                            else // i == 2
                            {
                                local_coords_j[0] = 0.0;
                                local_coords_j[1] = 1.0;
                                DN_De_j( 0 ) = 0.0;
                                DN_De_j( 1 ) = 1.0;
                            }
                        }
                        else if( num_nodes == 4 )    // linear quad element - FIXME: this will screw things up if the user defines a 4-node tri elem
                        {
                            if( i == 0 )
                            {
                                local_coords_j[0] = -1.0;
                                local_coords_j[1] = -1.0;
                                DN_De_j( 0 ) = -0.5;
                                DN_De_j( 1 ) = -0.5;
                            }
                            else if( i == 1 )
                            {
                                local_coords_j[0] =  1.0;
                                local_coords_j[1] = -1.0;
                                DN_De_j( 0 ) =  0.5;
                                DN_De_j( 1 ) = -0.5;
                            }
                            else if( i == 2 )
                            {
                                local_coords_j[0] =  1.0;
                                local_coords_j[1] =  1.0;
                                DN_De_j( 0 ) =  0.5;
                                DN_De_j( 1 ) =  0.5;
                            }
                            else // i == 3
                            {
                                local_coords_j[0] = -1.0;
                                local_coords_j[1] =  1.0;
                                DN_De_j( 0 ) = -0.5;
                                DN_De_j( 1 ) =  0.5;
                            }
                        }
                        else
                        {
                            KRATOS_THROW_ERROR( std::logic_error, "DELTA_NORMAL is not yet defined for higher order 2D elements. Number of nodes: ", cond_it->GetGeometry( ).PointsNumber() );
                        }
                        
                        cond_it->GetGeometry( ).Jacobian( J, local_coords_j );
                        
                        /*
                         * NOTES:
                         * Delta_ne_adj here is Delta_n_hat_e in Popp's thesis equation A.3 in appendix A.
                         * In 2D, it is also like this, but I split it into a node-dependent DN_De_j and a node-independent Delta_ne_adj
                         * In 3D, Delta_ne_adj is node-dependent => evaluated completely inside the nodes loop
                         * In both cases, Ce is the elemental contribution in Popp's thesis equation A.2 in appendix A
                         * 
                         * ::::: DERIVATION OF Delta_ne_adj Matrix in 3D ::::: 
                         * 
                         *   [ SUM( N,xi * Delta_x ) x SUM( N,eta * x ) ] + [ SUM( N,xi * x ) x SUM( N,eta * Delta_x ) ]
                         * = [ SUM( N,xi * Delta_x ) x J_eta            ] + [            J_xi x SUM( N,eta * Delta_x ) ]
                         * = [ SUM( N,xi * Delta_x ) x J_eta            ] - [ SUM( N,eta * Delta_x ) x J_xi ]
                         * SUM( N,xi * Delta_x ) is the consistent assembly of N,xi in blocks. Similarily for N,eta
                         * therefore, for node j we only care about N,xi(j) and N,eta(j) nodal blocks
                         * = [ N,xi(j) * Delta_x(j) x J_eta ] - [ N,eta(j) * Delta_x(j) x J_xi ]
                         * = [ N,xi(j) * ( ones(3) x J_eta ) - N,eta(j) *( ones(3) x J_xi ) ] * Delta_x(j)
                         * = Delta_ne_adj * Delta_x
                         */
                        
                        Delta_ne_adj(0,0) = 0.0;
                        Delta_ne_adj(0,1) = +J(2,1) * DN_De_j(0) - J(2,0) * DN_De_j(1); 
                        Delta_ne_adj(0,2) = -J(1,1) * DN_De_j(0) + J(1,0) * DN_De_j(1); 
                        Delta_ne_adj(1,0) = -J(2,1) * DN_De_j(0) + J(2,0) * DN_De_j(1); 
                        Delta_ne_adj(1,1) = 0.0;                   
                        Delta_ne_adj(1,2) = +J(0,1) * DN_De_j(0) - J(0,0) * DN_De_j(1); 
                        Delta_ne_adj(2,0) = +J(1,1) * DN_De_j(0) - J(1,0) * DN_De_j(1); 
                        Delta_ne_adj(2,1) = -J(0,1) * DN_De_j(0) + J(0,0) * DN_De_j(1); 
                        Delta_ne_adj(2,2) = 0.0;
                        
                        Ce = prod( I - ne_o_ne, Delta_ne_adj ) / ne_norm;
                        node_j.GetValue(DELTA_NORMAL) += Ce;
                    }
                }
                else
                {
                    KRATOS_THROW_ERROR( std::logic_error, "Bad dimension provided to calculate DELTA_NORMAL. Dimension = ", dimension )
                }
            }
        }

        // Normalize normal directional derivatives
        it_node_begin = pNode.ptr_begin();
        it_node_end   = pNode.ptr_end();

        for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
        {
            if (node_it->Is(INTERFACE))
            {
                const array_1d<double, 3> & nj = node_it->GetValue(NORMAL); // nodal non-normalized normal (this function is called before normalization)
                
                Matrix nj_o_nj = subrange( outer_prod( nj, nj ), 0, dimension, 0, dimension );
                const double nj_norm = norm_2( nj );
                const double nj_norm_3 = nj_norm * nj_norm * nj_norm;
                
                if ( nj_norm_3 > tol )
                {
                    const Matrix Cj = I / nj_norm - nj_o_nj / nj_norm_3;
                    node_it->GetValue(DELTA_NORMAL) = prod( Cj, node_it->GetValue(DELTA_NORMAL) );
                }
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Calculates the determinant of the jacobian of the contact element  
     * @param Jacobian: The element's jacobian
     * @return The determinant of the provided jacobian
     */
    static inline const double ContactElementDetJacobian( const Matrix& J )
    {
        Matrix JTJ = prod( trans(J), J );
        if( J.size2( ) == 1 )
        {
            return std::sqrt( JTJ(0,0) );
        }
        else if( J.size2( ) == 2 )
        {
            return std::sqrt( JTJ(0,0) * JTJ(1,1) - JTJ(1,0) * JTJ(0,1) );
        }
        else
        {
            KRATOS_THROW_ERROR( std::logic_error, "Illegal local dimension for contact element. Dimension = ", J.size2( ) );
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
        ModelPart & rModelPart, 
        const double cn
        )
    {
        NodesArrayType& pNode  = rModelPart.Nodes();
        NodesArrayType::iterator it_node_begin = pNode.ptr_begin();
        NodesArrayType::iterator it_node_end   = pNode.ptr_end();
        
        for(NodesArrayType::iterator node_it = it_node_begin; node_it!=it_node_end; node_it++)
        {
            if (node_it->Is(INTERFACE))
            {
                const array_1d<double,3> lagrange_multiplier = node_it->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0);
                const array_1d<double,3>        nodal_normal = node_it->GetValue(NORMAL); 
                const double lambda_n = inner_prod(lagrange_multiplier, nodal_normal);

                const double check = lambda_n - cn * node_it->GetValue(WEIGHTED_GAP);

                if (check <= 0.0)
                {
                    node_it->Set(ACTIVE, false);
                    node_it->GetSolutionStepValue( IS_ACTIVE_SET ) = false;
                }
                else
                {
                    node_it->Set(ACTIVE, true);
                    node_it->GetSolutionStepValue( IS_ACTIVE_SET ) = true;
                }
                

//                /// DEBUG ///
//                if( node_it->Id() == 616 || node_it->Id() == 629 )
//                {
//                    DEBUG_MSG( "Recomputing the active/inactive sets using PDASS" )
//                    KRATOS_WATCH( node_it->Id( ) )
//                    LOG_VECTOR3( lagrange_multiplier )
//                    LOG_VECTOR3( node_it->GetValue(NORMAL) )
//                    LOG_SCALAR( node_it->GetValue(WEIGHTED_GAP) )
//                    LOG_SCALAR( check )
//                    LOG_SCALAR( lambda_n )
//                    LOG_SCALAR( node_it->Is( ACTIVE ) )
//                }
//                /// DEBUG ///
            }
        }
    }
    
private:
};// class ContactUtilities
}
#endif /* KRATOS_CONTACT_UTILITIES defined */
 
