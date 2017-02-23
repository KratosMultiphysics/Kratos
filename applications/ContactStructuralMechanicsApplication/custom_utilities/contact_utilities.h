 
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
#include "utilities/openmp_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
class ContactUtilities
{
public:
    ///@name Type Definitions
    ///@{
    
    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
    ///@}
    ///@name Operations
    ///@{

    
    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef GeometryData::IntegrationMethod         IntegrationMethod;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;



    /**
     * This function fills the contact_container for the Mortar condition
     * @param ContactPoint: The destination point
     * @param Geom1 and Geom2: The geometries of the slave and master respectively
     * @param contact_normal1 and contact_normal2: The normals of the slave and master respectively
     * @param IntegrationOrder: The integration order   
     * @return contact_container: Once has been filled
     */

    static inline bool ContactContainerFiller(
            GeometryType& Geom1, // SLAVE
            GeometryType& Geom2, // MASTER
            const array_1d<double, 3> & contact_normal1, // SLAVE
            const array_1d<double, 3> & contact_normal2, // MASTER
            const double ActiveCheckFactor
            )
    {
        // Define the basic information
        const unsigned int number_nodes = Geom1.PointsNumber();
        
        bool cond_active = false;
        
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
              
                // TODO: Think about this
                double dist_tol = ActiveCheckFactor;    // The actual gap tolerance is user-define instead of being a factor of the length
//                 double dist_tol = ActiveCheckFactor * Geom1.Length();
//                 dist_tol = (dist_tol <= ActiveCheckFactor * Geom2.Length()) ? (ActiveCheckFactor * Geom2.Length()):dist_tol;
                
                array_1d<double, 3> result;
                // NOTE: We don't use std::abs() because if the aux_dist is negative is penetrating, in fact we just consider dist_tol > 0 to have some tolerance and for the static schemes
                if (aux_dist <= dist_tol && Geom2.IsInside(ProjectedPoint, result))
                { 
                    Geom1[index].Set(ACTIVE, true);
                    cond_active = true;
                }
             }
             else
             {
                 cond_active = true;
             }
         }
         
         return cond_active;
    }
    
    static inline void ContactContainerFiller(
            std::vector<contact_container> *& ConditionPointers,
            Condition::Pointer & pCond_1,       // SLAVE
            const Condition::Pointer & pCond_2, // MASTER
            const array_1d<double, 3> & contact_normal1, // SLAVE
            const array_1d<double, 3> & contact_normal2, // MASTER
            const double ActiveCheckFactor
            )
    {
        const bool cond_active = ContactContainerFiller(pCond_1->GetGeometry(), pCond_2->GetGeometry(), contact_normal1, contact_normal2, ActiveCheckFactor);
        
        if (cond_active == true)
        {
            pCond_1->Set(ACTIVE, true);
            contact_container aux_contact_container;
            aux_contact_container.condition   = pCond_2;
            aux_contact_container.active_pair = true;
            ConditionPointers->push_back(aux_contact_container);
        }
    }
    
    /**
     * Project a point over a line/plane following an arbitrary direction
     * @param Geom: The geometry where to be projected
     * @param PointDestiny: The point to be projected
     * @param Vector: The direction to project
     * @return PointProjected: The point pojected over the plane
     * @return dist: The distance between the point and the plane
     */

    static inline void ProjectDirection(
            const GeometryType& Geom,
            const Point<3>& PointDestiny,
            Point<3>& PointProjected,
            double& dist,
            const array_1d<double,3>& Vector
            )
    {        
        const double tol = 1.0e-15;
        
        array_1d<double,3> Normal;
        
        GeometryNormal(Normal, Geom);
        
//         const array_1d<double,3> vector_points = Geom.Center() - PointDestiny.Coordinates();
        const array_1d<double,3> vector_points = Geom[0].Coordinates() - PointDestiny.Coordinates();

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
            const GeometryType& Geom,
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
    
    /**
     * Projects iteratively to get the coordinate
     * @param GeomOrigin: The origin geometry
     * @param PointDestiny: The destination point
     * @return ResultingPoint: The distance between the point and the plane
     * @return Inside: True is inside, false not
     */
    
    static inline bool ProjectIterative(
        GeometryType& GeomOrigin,
        const GeometryType::CoordinatesArrayType& PointDestiny,
        GeometryType::CoordinatesArrayType& ResultingPoint
        )
    {
        Matrix J = ZeroMatrix( GeomOrigin.LocalSpaceDimension(), GeomOrigin.LocalSpaceDimension() );

        ResultingPoint.clear();

        Vector DeltaXi = ZeroVector( GeomOrigin.LocalSpaceDimension() );

        GeometryType::CoordinatesArrayType CurrentGlobalCoords( ZeroVector( 3 ) );

        Vector NOrigin;
        
        //Newton iteration:
        double tol = 1.0e-8;

        unsigned int maxiter = 30;

        for ( unsigned int k = 0; k < maxiter; k++ )
        {
            GeomOrigin.ShapeFunctionsValues( NOrigin, ResultingPoint );
            
            const array_1d<double,3> normal_xi = GaussPointNormal(NOrigin, GeomOrigin);
            
            CurrentGlobalCoords = ZeroVector( 3 );
            GeomOrigin.GlobalCoordinates( CurrentGlobalCoords, ResultingPoint );
            
            double dist;
            ProjectCoordDirection(GeomOrigin, PointDestiny, CurrentGlobalCoords, dist, normal_xi);
            
            noalias( CurrentGlobalCoords ) = PointDestiny - CurrentGlobalCoords;
            GeomOrigin.InverseOfJacobian( J, ResultingPoint );
            noalias( DeltaXi ) = prod( J, CurrentGlobalCoords );
            noalias( ResultingPoint ) += DeltaXi;

            if ( norm_2( DeltaXi ) < tol )
            {
                return true;
            }
        }

        return false;
    }
    
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
    
    /**
     * Calculates the distance between nodes
     * @param PointOrigin: A point in the plane
     * @param PointDestiny: The point to be projected
     */
    
    static inline double DistancePoints(
            const Point<3>& PointOrigin,
            const Point<3>& PointDestiny
            )
    {
        const double dist = std::sqrt((PointOrigin.Coordinate(1) - PointDestiny.Coordinate(1)) * (PointOrigin.Coordinate(1) - PointDestiny.Coordinate(1))
                                    + (PointOrigin.Coordinate(2) - PointDestiny.Coordinate(2)) * (PointOrigin.Coordinate(2) - PointDestiny.Coordinate(2))
                                    + (PointOrigin.Coordinate(3) - PointDestiny.Coordinate(3)) * (PointOrigin.Coordinate(3) - PointDestiny.Coordinate(3)));
        
        return dist;
    }

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
            
            const double tmp = inner_prod(aux_vector, aux_vector);

            if(tmp > Radius)
            {
                Radius = tmp;
            }
        }

        Radius = std::sqrt(Radius);
        
        ConditionNormal(pCond);
    }

    /**
     * This function calculates the normal of a condition
     * @param Cond: The pointer to the condition of interest
     */

    static inline void ConditionNormal(Condition::Pointer pCond)
    {
        // TODO: Add calculation of normal to geometry.h 
        array_1d<double,3> & Normal     = pCond->GetValue(NORMAL);
        array_1d<double,3> & TangentXi  = pCond->GetValue(TANGENT_XI);
        array_1d<double,3> & TangentEta = pCond->GetValue(TANGENT_ETA);
        
        GeometryNormal(Normal, TangentXi, TangentEta, pCond->GetGeometry());
    }

    /**
     * This function calculates the normal of a condition
     * @param Cond: The pointer to the condition of interest
     */

    static inline array_1d<double,3> GaussPointNormal(
        const Vector N,
        const GeometryType & Geom
        )
    {
        array_1d<double,3> normal = ZeroVector(3);
        for( unsigned int iNode = 0; iNode < Geom.PointsNumber(); ++iNode )
        {
            normal += N[iNode] * Geom[iNode].GetValue(NORMAL); 
        }
        
        normal = normal/norm_2(normal); // It is suppossed to be already unitary (just in case)
        
        return normal;
    }
    
    /**
     * This function calculates the normal of a geometry
     * @param Cond: The pointer to the condition of interest
     */

    static inline void GeometryNormal(
            array_1d<double,3> & Normal,
            array_1d<double,3> & TangentXi,
            array_1d<double,3> & TangentEta,
            const GeometryType & Geom
            )
    {
        noalias(Normal) = ZeroVector(3);
        
        // Geom normal is the sum of all nodal normals
        // nodal normal = tangent_eta (or e3 in 2D) x tangent_xi
        for ( unsigned int i = 0; i < Geom.PointsNumber( ); ++i )
        {
            NodalTangents( TangentXi, TangentEta, Geom, i );
            Normal += MathUtils<double>::CrossProduct( TangentEta, TangentXi );
        }
        
        Normal     /= norm_2( Normal );
        TangentXi  /= norm_2( TangentXi );
        TangentEta /= norm_2( TangentEta );
    }
        
    static inline void GeometryNormal(
            array_1d<double,3> & Normal,
            const GeometryType & Geom
            )
    {
        array_1d<double,3> TangentXi, TangentEta;
        
        GeometryNormal(Normal, TangentXi, TangentEta, Geom);
    }
    
    /**
     * This function calculates the tangents at a given node
     */

    static inline void NodalTangents(
            array_1d<double,3> & t1,    // tangent in xi direction
            array_1d<double,3> & t2,    // tangent in eta direction in 3D or simply e3 cartesian vector in 2D
            const GeometryType & Geom,
            const unsigned int i_node
            )
    {
        const unsigned int dimension = Geom.WorkingSpaceDimension( );
        const unsigned int local_dim = Geom.LocalSpaceDimension( );
        Matrix J_i = ZeroMatrix( dimension, local_dim ); 
        Geom.Jacobian( J_i, Geom[i_node].Coordinates( ) );

        if( dimension == 2 )
        {
            // t1 (axis-direction)
            t1[0] = J_i(0, 0);
            t1[1] = J_i(1, 0);
            t1[2] = 0.0;
            // t2 (z-direction)
            t2[0] = 0.0;
            t2[1] = 0.0;
            t2[2] = 1.0;
        }
        else // We can use just else, if it is not 2D it is 3D 
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                // Using the Jacobian tangent directions 
                t1[i] = J_i(i, 0);
                t2[i] = J_i(i, 1);
            }
        }
        
    }
    
    /**
     * It computes the mean of the normal in the condition in all the nodes
     * @param ModelPart: The model part to compute
     * @return The modelparts with the normal computed
     */
    
    static inline void ComputeNodesMeanNormalModelPart(ModelPart & rModelPart) 
    {
        // Tolerance
        const double tol = 1.0e-14;

        // Initialize normal vectors
        const array_1d<double,3> ZeroVect = ZeroVector(3);
        
        NodesArrayType& pNode = rModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
//         #pragma omp parallel for // NOTE: Giving problems!!
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->GetValue(NODAL_AREA)           = 0.0;
            noalias(itNode->GetValue(NORMAL))      = ZeroVect;
            noalias(itNode->GetValue(TANGENT_XI))  = ZeroVect; 
            noalias(itNode->GetValue(TANGENT_ETA)) = ZeroVect; 
        }
        
        // Sum all the nodes normals
        ConditionsArrayType& pCond = rModelPart.Conditions();
        auto numConditions = pCond.end() - pCond.begin();
        
//         #pragma omp parallel for // NOTE: Don't parallelize, you are accesing to the nodes (try with atomic)
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            if (itCond->Is(ACTIVE) || itCond->Is(MASTER))
            {
                ConditionNormal(*(itCond.base()));
                
                const unsigned int number_nodes = itCond->GetGeometry().PointsNumber();
                const double & rArea = itCond->GetGeometry().Area()/number_nodes;
                const array_1d<double, 3> & rNormal     = itCond->GetValue(NORMAL);
                const array_1d<double, 3> & rTangentXi  = itCond->GetValue(TANGENT_XI);
                const array_1d<double, 3> & rTangentEta = itCond->GetValue(TANGENT_ETA);
                
                for (unsigned int i = 0; i < number_nodes; i++)
                {
                    itCond->GetGeometry()[i].GetValue(NODAL_AREA)             += rArea;
                    noalias( itCond->GetGeometry()[i].GetValue(NORMAL) )      += rArea * rNormal;
                    noalias( itCond->GetGeometry()[i].GetValue(TANGENT_XI) )  += rArea * rTangentXi;
                    noalias( itCond->GetGeometry()[i].GetValue(TANGENT_ETA) ) += rArea * rTangentEta;
                }
            }
        }
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;

            const double total_area        = itNode->GetValue(NODAL_AREA);
            itNode->GetValue(NORMAL)      /= total_area;
            itNode->GetValue(TANGENT_XI)  /= total_area;
            itNode->GetValue(TANGENT_ETA) /= total_area;
        }
        
//         // Applied laziness - MUST be calculated BEFORE normalizing the normals
//         ComputeDeltaNodesMeanNormalModelPart( rModelPart ); // NOTE: The area pondering?, after or before

        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;

//             const double total_area        = itNode->GetValue(NODAL_AREA);
//             itNode->GetValue(NORMAL)      /= total_area;
//             itNode->GetValue(TANGENT_XI)  /= total_area;
//             itNode->GetValue(TANGENT_ETA) /= total_area;

            const double norm_normal     = norm_2(itNode->GetValue(NORMAL));
            const double norm_tangentxi  = norm_2(itNode->GetValue(TANGENT_XI));
            const double norm_tangenteta = norm_2(itNode->GetValue(TANGENT_ETA));
            
            if (norm_normal > tol)
            {
                itNode->GetValue(NORMAL)      /= norm_normal;
                itNode->GetValue(TANGENT_XI)  /= norm_tangentxi;
                itNode->GetValue(TANGENT_ETA) /= norm_tangenteta;
            }
        }
    }

    /**
     * It computes the directional derivative of the normal in the condition in all the nodes
     * @param ModelPart: The model part to compute
     * @return The modelparts with the normal computed
     */

    static inline void ComputeDeltaNodesMeanNormalModelPart(ModelPart & rModelPart)
    {
        // TODO: Add parallelization

        // Conditions
        ConditionsArrayType& pCond  = rModelPart.Conditions();
        ConditionsArrayType::iterator it_cond_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_cond_end   = pCond.ptr_end();

        // Tolerance
        const double tol = 1.0e-14;

        // Initialize directional derivative
        const unsigned int dimension = it_cond_begin->WorkingSpaceDimension( ); 
        const Matrix ZeroDeltaNormal = ZeroMatrix( dimension, dimension );

        Matrix Delta_ne_adj  = Matrix(dimension, dimension);
        Matrix Ce = Matrix(dimension, dimension);
        
        const Matrix I = IdentityMatrix(dimension, dimension);

        NodesArrayType& pNode = rModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->GetValue(DELTA_NORMAL) = ZeroDeltaNormal;
        }
        
        // Sum the directional derivatives of the adjacent segments
        for(ConditionsArrayType::iterator itCond = it_cond_begin; itCond!=it_cond_end; itCond++)
        {
            if (itCond->Is(ACTIVE) || itCond->Is(MASTER)) 
            {
                const array_1d<double, 3> ne = itCond->GetValue(NORMAL);   // normalized condition normal
                Matrix ne_o_ne = subrange( outer_prod( ne, ne ), 0, dimension, 0, dimension );

                const unsigned int num_nodes = itCond->GetGeometry( ).PointsNumber();
                if( dimension == 2 )
                {
                    if (num_nodes == 2)
                    {
                        const double ne_norm = itCond->GetGeometry( ).Length( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                        
                        Delta_ne_adj( 0, 0 ) =  0.0;
                        Delta_ne_adj( 0, 1 ) = -1.0;
                        Delta_ne_adj( 1, 0 ) =  1.0;
                        Delta_ne_adj( 1, 1 ) =  0.0;
                        
                        Ce = prod( I - ne_o_ne, Delta_ne_adj ) / ne_norm;     // In 2D, Delta_ne_adj is node-independent => evaluated outside the nodes loop
                        
                        for (unsigned int i = 0; i < num_nodes; i++)
                        {
                            NodeType& node_j = itCond->GetGeometry( )[i];
                            
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
                            node_j.GetValue(DELTA_NORMAL) += Ce * DN_De_j;
                        }
                    }
                    else
                    {
                        KRATOS_ERROR << "DELTA_NORMAL is not yet defined for higher order 1D elements. Number of nodes: " << num_nodes << std::endl;
                    }
                }
                else if ( dimension == 3 )
                {
                    const double ne_norm = itCond->GetGeometry( ).Area( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                    
                    for (unsigned int i = 0; i < num_nodes; i++)
                    {
                        Matrix J = ZeroMatrix( 3, 2 ); // Jacobian [ 3D global x 2D local ]
                        array_1d<double, 2> DN_De_j        = ZeroVector( 2 );  // Shape functions derivatives for node j [ DN_Dxi_j, DN_Deta_j ]
                        array_1d<double, 3> local_coords_j = ZeroVector( 3 );
                        
                        NodeType& node_j = itCond->GetGeometry( )[i];
                        
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
                            KRATOS_ERROR << "DELTA_NORMAL is not yet defined for higher order 2D elements. Number of nodes: " << itCond->GetGeometry( ).PointsNumber() << std::endl;
                        }
                        
                        itCond->GetGeometry( ).Jacobian( J, local_coords_j );
                        
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
                    KRATOS_ERROR << "Bad dimension provided to calculate DELTA_NORMAL. Dimension = " << dimension << std::endl;
                }
            }
        }
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            const array_1d<double, 3> & nj = itNode->GetValue(NORMAL); // nodal non-normalized normal (this function is called before normalization)
            
            Matrix nj_o_nj = subrange( outer_prod( nj, nj ), 0, dimension, 0, dimension );
            const double nj_norm = norm_2( nj );
            const double nj_norm_3 = nj_norm * nj_norm * nj_norm;
            
            if ( nj_norm_3 > tol )
            {
                const Matrix Cj = I / nj_norm - nj_o_nj / nj_norm_3;
                itNode->GetValue(DELTA_NORMAL) = prod( Cj, itNode->GetValue(DELTA_NORMAL) );
            }
        }
    }
    
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
            KRATOS_ERROR << "Illegal local dimension for contact element. Dimension = " << J.size2( ) << std::endl;
        }
    }
    
    /**
     * It changes from active to inactive and viceversa the nodes 
     * @param ModelPart: The model part to compute
     * @return The modelparts with the conditions changed
     */
    
    static inline void ReComputeActiveInactive(ModelPart & rModelPart)
    {
        // TODO: If works make it parallell (it is difficult, be careful with the repeated nodes) 
       
        ConditionsArrayType& pCond = rModelPart.GetSubModelPart("Contact").Conditions();
        ConditionsArrayType::iterator it_cond_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_cond_end   = pCond.ptr_end();
        
        for(ConditionsArrayType::iterator itCond = it_cond_begin; itCond!=it_cond_end; itCond++)
        {
            bool condition_is_active = false; // It is supposed to be always defined, and with this only the slave conditions will be taken into account
            if( (itCond)->IsDefined(ACTIVE) == true)
            {
                condition_is_active = (itCond)->Is(ACTIVE);
            }
            
            if ( condition_is_active == true )
            {
                // Recompute Active/Inactive nodes
                double cn = 0.0;
                if (itCond->GetProperties().Has(NORMAL_AUGMENTATION_FACTOR) == true)
                {
                    cn = itCond->GetProperties().GetValue(NORMAL_AUGMENTATION_FACTOR); 
                }
                
                double ct = 0.0;
                if (itCond->GetProperties().Has(TANGENT_AUGMENTATION_FACTOR) == true)
                {
                    ct = itCond->GetProperties().GetValue(TANGENT_AUGMENTATION_FACTOR); 
                }

                GeometryType & CondGeometry = itCond->GetGeometry();
                
                for(unsigned int itNode = 0; itNode!=CondGeometry.PointsNumber(); itNode++)
                {
                    if (CondGeometry[itNode].Is(VISITED) == false)
                    {
                        const double mu = CondGeometry[itNode].GetValue(WEIGHTED_FRICTION);
                        const double gn = CondGeometry[itNode].GetValue(WEIGHTED_GAP);
                        const array_1d<double,3> lagrange_multiplier = CondGeometry[itNode].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0);
                        const array_1d<double,3>        nodal_normal = CondGeometry[itNode].GetValue(NORMAL); 
                        
//                         const double lambda_n = inner_prod(lagrange_multiplier, nodal_normal);
//                         const double augmented_normal_presssure = lambda_n + cn * gn;           
                        
                        double augmented_normal_presssure = inner_prod(lagrange_multiplier, nodal_normal);
                        if (gn < 0.0) // NOTE: Penetration
                        {
                            augmented_normal_presssure += cn * gn;     
                        }
                        
                        if (augmented_normal_presssure < 0.0) // NOTE: This could be conflictive (< or <=)
                        {
                            CondGeometry[itNode].Set(ACTIVE, true);
                        }
                        else
                        {
                            CondGeometry[itNode].Set(ACTIVE, false);
//                             CondGeometry[itNode].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0) = ZeroVector(3);
                        }
                        
                        const array_1d<double, 3> nodal_tangent_xi  = CondGeometry[itNode].GetValue(TANGENT_XI); 
                        const array_1d<double, 3> nodal_tangent_eta = CondGeometry[itNode].GetValue(TANGENT_ETA);
                        const double tangent_xi_lm  = inner_prod(nodal_tangent_xi,  lagrange_multiplier);
                        const double tangent_eta_lm = inner_prod(nodal_tangent_eta, lagrange_multiplier);
                        const double lambda_t = std::sqrt(tangent_xi_lm * tangent_xi_lm + tangent_eta_lm * tangent_eta_lm); 
                        
                        const double augmented_tangent_presssure = std::abs(lambda_t + ct * CondGeometry[itNode].GetValue(WEIGHTED_SLIP)) + mu * augmented_normal_presssure;
                        
                        if (augmented_tangent_presssure < 0.0) // TODO: Check if it is minor equal or just minor
                        {
                            CondGeometry[itNode].Set(SLIP, false);
                        }
                        else
                        {
                            CondGeometry[itNode].Set(SLIP, true);
                        }
                        
                        CondGeometry[itNode].Set(VISITED, true);
                    }
                }
            }
        }
    }
    
    /**
     * It changes from active to inactive and viceversa the nodes 
     * @param ModelPart: The model part to compute
     * @return The modelparts with the conditions changed
     */
    
    static inline void ReComputeActiveInactiveALMFrictionless(ModelPart & rModelPart)
    {
        // TODO: If works make it parallell (it is difficult, be careful with the repeated nodes) 
       
        ConditionsArrayType& pCond = rModelPart.GetSubModelPart("Contact").Conditions();
        ConditionsArrayType::iterator it_cond_begin = pCond.ptr_begin();
        ConditionsArrayType::iterator it_cond_end   = pCond.ptr_end();
        
        for(ConditionsArrayType::iterator itCond = it_cond_begin; itCond!=it_cond_end; itCond++)
        {
            bool condition_is_active = false; // It is supposed to be always defined, and with this only the slave conditions will be taken into account
            if( (itCond)->IsDefined(ACTIVE) == true)
            {
                condition_is_active = (itCond)->Is(ACTIVE);
            }
            
            if ( condition_is_active == true )
            {
                // Recompute Active/Inactive nodes
                double epsilon = 0.0;
                if (itCond->GetProperties().Has(PENALTY_FACTOR) == true)
                {
                    epsilon = itCond->GetProperties().GetValue(PENALTY_FACTOR); 
                }
                
                double k = 0.0;
                if (itCond->GetProperties().Has(SCALE_FACTOR) == true)
                {
                    k = itCond->GetProperties().GetValue(SCALE_FACTOR); 
                }

                GeometryType & CondGeometry = itCond->GetGeometry();
                
                for(unsigned int itNode = 0; itNode!=CondGeometry.PointsNumber(); itNode++)
                {
                    if (CondGeometry[itNode].Is(VISITED) == false)
                    {
                        const double gn = CondGeometry[itNode].GetValue(WEIGHTED_GAP);
                        const double augmented_normal_presssure = k * CondGeometry[itNode].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS, 0) + epsilon * gn;     
                        
                        if (augmented_normal_presssure < 0.0) // NOTE: This could be conflictive (< or <=)
                        {
                            CondGeometry[itNode].Set(ACTIVE, true);
                        }
                        else
                        {
                            CondGeometry[itNode].Set(ACTIVE, false);
//                             CondGeometry[itNode].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS, 0) = 0.0;
                        }
                        
//                         // Debug 
//                         std::cout << CondGeometry[itNode].Id() << " Gap: " << gn  << " Pressure: " << augmented_normal_presssure << " Active: " << CondGeometry[itNode].Is(ACTIVE) << std::endl;
                        
                        CondGeometry[itNode].Set(VISITED, true);
                    }
                }
            }
        }
    }
    
    /**
     * It resets the visited status in all the nodes
     * @param ModelPart: The model part to compute
     * @return The modelparts with the nodes changed
     */
    
    static inline void ResetVisited(ModelPart & rModelPart)
    {
        NodesArrayType& pNode = rModelPart.GetSubModelPart("Contact").Nodes();
        
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;

            if (itNode->Is(SLAVE))
            {   
                itNode->Set(VISITED, false);
            }
        }
    }
    
    /**
     * It resets the value of the weighted gap and slip 
     * @param ModelPart: The model part to compute
     * @return The modelparts with the nodes changed
     */
    
    static inline void ResetWeightedValues(ModelPart & rModelPart)
    {
        NodesArrayType& pNode = rModelPart.GetSubModelPart("Contact").Nodes();
        
        const bool DLM = rModelPart.GetSubModelPart("Contact").Is(MODIFIED);
        
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;

            if (itNode->Is(SLAVE))
            {   
                if (itNode->Is(ACTIVE) == false)
                {
                    itNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, 0) = ZeroVector(3);
                    if (DLM == true)
                    {
                        itNode->FastGetSolutionStepValue(DOUBLE_LM, 0) = ZeroVector(3);
                    }
                }
                itNode->GetValue(WEIGHTED_GAP)      = 0.0;
                itNode->GetValue(WEIGHTED_SLIP)     = 0.0;
                itNode->GetValue(WEIGHTED_FRICTION) = 0.0;
            }
        }
    }
    
    /**
     * It resets the value of the weighted gap 
     * @param ModelPart: The model part to compute
     * @return The modelparts with the nodes changed
     */
    
    static inline void ResetWeightedALMFrictionlessValues(ModelPart & rModelPart)
    {
        NodesArrayType& pNode = rModelPart.GetSubModelPart("Contact").Nodes();
        
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;

            if (itNode->Is(SLAVE))
            {   
                if (itNode->Is(ACTIVE) == false)
                {
                    itNode->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS, 0) = 0.0;
                }
                itNode->GetValue(WEIGHTED_GAP) = 0.0;
            }
        }
    }
    
    /**
     * It calculates the matrix of coordinates of a geometry
     * @param nodes: The geometry to calculate
     * @param current: If we calculate the current coordinates or the initial ones
     * @return Coordinates: The matrix containing the coordinates of the geometry
     */
    
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline bounded_matrix<double, TNumNodes, TDim> GetCoordinates(
        const GeometryType& nodes,
        const bool current = true
        )
    {
        /* DEFINITIONS */            
        bounded_matrix<double, TNumNodes, TDim> Coordinates;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            array_1d<double, 3> coord = nodes[iNode].Coordinates();
            
            if (current == false)
            {
                coord -= nodes[iNode].FastGetSolutionStepValue(DISPLACEMENT, 0);
            }

            for (unsigned int iDof = 0; iDof < TDim; iDof++)
            {
                Coordinates(iNode, iDof) = coord[iDof];
            }
        }
        
        return Coordinates;
    }

    /**
     * It calculates the vector of a variable of a geometry
     * @param nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return VarVector: The vector containing the variables of the geometry
     */
    
    template< unsigned int TNumNodes >
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& nodes,
        const Variable<double>& rVarName,
        unsigned int step
        )
    {
        /* DEFINITIONS */        
        array_1d<double, TNumNodes> VarVector;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            VarVector[iNode] = nodes[iNode].FastGetSolutionStepValue(rVarName, step);
        }
        
        return VarVector;
    }

    template< unsigned int TNumNodes >
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& nodes,
        const Variable<double>& rVarName
        )
    {
        /* DEFINITIONS */        
        array_1d<double, TNumNodes> VarVector;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            VarVector[iNode] = nodes[iNode].GetValue(rVarName);
        }
        
        return VarVector;
    }
    
    /**
     * It calculates the matrix of a variable of a geometry
     * @param nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return VarMatrix: The matrix containing the variables of the geometry
     */
    
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline Matrix GetVariableMatrix(
        const GeometryType& nodes,
        const Variable<array_1d<double,3> >& rVarName,
        unsigned int step
        )
    {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, TDim> VarMatrix;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            const array_1d<double, 3> Value = nodes[iNode].FastGetSolutionStepValue(rVarName, step);
            for (unsigned int iDof = 0; iDof < TDim; iDof++)
            {
                VarMatrix(iNode, iDof) = Value[iDof];
            }
        }
        
        return VarMatrix;
    }

    template< unsigned int TDim, unsigned int TNumNodes>
    static inline Matrix GetVariableMatrix(
        const GeometryType& nodes,
        const Variable<array_1d<double,3> >& rVarName
        )
    {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, TDim> VarMatrix;
        
        for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
        {
            const array_1d<double, 3> Value = nodes[iNode].GetValue(rVarName);
            for (unsigned int iDof = 0; iDof < TDim; iDof++)
            {
                VarMatrix(iNode, iDof) = Value[iDof];
            }
        }
        
        return VarMatrix;
    }
    
private:
};// class ContactUtilities
}
#endif /* KRATOS_CONTACT_UTILITIES defined */
 
