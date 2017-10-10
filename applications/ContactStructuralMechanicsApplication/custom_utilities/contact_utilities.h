// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CONTACT_UTILITIES)
#define KRATOS_CONTACT_UTILITIES

#include "utilities/math_utils.h"
#include "contact_structural_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "geometries/point.h"
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
    
    // General type definitions
    typedef Node<3>                                              NodeType;
    typedef Point<3>                                            PointType;
    typedef PointType::CoordinatesArrayType          CoordinatesArrayType;
    typedef Geometry<NodeType>                               GeometryType;
    typedef Geometry<PointType>                         GeometryPointType;
    typedef GeometryData::IntegrationMethod             IntegrationMethod;
    typedef ModelPart::NodesContainerType                  NodesArrayType;
    typedef ModelPart::ConditionsContainerType        ConditionsArrayType;
    
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

    /**
     * Project a point over a line/plane following an arbitrary direction
     * @param Geom: The geometry where to be projected
     * @param PointDestiny: The point to be projected
     * @param Normal: The normal of the geometry
     * @param Vector: The direction to project
     * @return PointProjected: The point pojected over the plane
     * @return Distance: The distnace between surfaces
     */

    static inline double FastProjectDirection(
        const GeometryType& Geom,
        const PointType& PointDestiny,
        PointType& PointProjected,
        const array_1d<double,3>& Normal,
        const array_1d<double,3>& Vector
        )
    {    
        // We define the tolerance
        const double tolerance = std::numeric_limits<double>::epsilon();
        
        // We define the distance
        double distance = 0.0;
        
        const array_1d<double,3> vector_points = Geom[0].Coordinates() - PointDestiny.Coordinates();

        if( norm_2( Vector ) < tolerance && norm_2( Normal ) > tolerance )
        {
            distance = inner_prod(vector_points, Normal)/norm_2(Normal);

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * distance;
            std::cout << " :: Warning: Zero projection vector. Projection using the condition vector instead." << std::endl;
        }
        else if (std::abs(inner_prod(Vector, Normal) ) > tolerance)
        {
            distance = inner_prod(vector_points, Normal)/inner_prod(Vector, Normal); 

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * distance;
        }
        else
        {
            PointProjected.Coordinates() = PointDestiny.Coordinates();
            std::cout << " The line and the plane are coplanar, something wrong happened " << std::endl;
        }
        
        return distance;
    }
    
    /**
     * Projects iteratively to get the coordinate
     * @param GeomOrigin: The origin geometry
     * @param PointDestiny: The destination point
     * @return ResultingPoint: The distance between the point and the plane
     * @return Inside: True is inside, false not
     */
    
    static inline bool ProjectIterativeLine2D(
        GeometryType& GeomOrigin,
        const GeometryType::CoordinatesArrayType& PointDestiny,
        GeometryType::CoordinatesArrayType& ResultingPoint,
        const array_1d<double, 3>& Normal,
        const double Tolerance = 1.0e-8,
        double DeltaXi = 0.5
        )
    {
//         ResultingPoint.clear();
        
        double old_delta_xi = 0.0;

        array_1d<double, 3> current_global_coords;

        array_1d<array_1d<double, 3>, 2> normals;
        normals[0] = GeomOrigin[0].GetValue(NORMAL);
        normals[1] = GeomOrigin[1].GetValue(NORMAL);
        
        bounded_matrix<double,2,2> X;
        bounded_matrix<double,2,1> DN;
        for(unsigned int i=0; i<2;i++)
        {
            X(0,i) = GeomOrigin[i].X();
            X(1,i) = GeomOrigin[i].Y();
        }

        Matrix J = ZeroMatrix( 1, 1 );
        
        //Newton iteration:

        const unsigned int max_iter = 20;

        for ( unsigned int k = 0; k < max_iter; k++ )
        {
            array_1d<double, 2> N_origin;
            N_origin[0] = 0.5 * ( 1.0 - ResultingPoint[0]);
            N_origin[1] = 0.5 * ( 1.0 + ResultingPoint[0]);
            
            array_1d<double,3> normal_xi = ZeroVector(3);
            for( unsigned int i_node = 0; i_node < 2; ++i_node )
            {
                normal_xi += N_origin[i_node] * normals[i_node]; 
            }
            
            normal_xi = normal_xi/norm_2(normal_xi); 
            
            current_global_coords = ZeroVector( 3 );
            for( unsigned int i_node = 0; i_node < 2; ++i_node )
            {
                current_global_coords += N_origin[i_node] * GeomOrigin[i_node].Coordinates(); 
            }
            
            const array_1d<double,3> VectorPoints = GeomOrigin.Center() - PointDestiny;
            const double distance = inner_prod(VectorPoints, Normal)/inner_prod(-normal_xi, Normal); 
            const array_1d<double, 3> current_destiny_global_coords = PointDestiny - normal_xi * distance;
            
            // Derivatives of shape functions
            Matrix ShapeFunctionsGradients;
            ShapeFunctionsGradients = GeomOrigin.ShapeFunctionsLocalGradients(ShapeFunctionsGradients, ResultingPoint );
            noalias(DN) = prod(X,ShapeFunctionsGradients);

            noalias(J) = prod(trans(DN),DN); // TODO: Add the non linearity concerning the normal
            Vector RHS = prod(trans(DN),subrange(current_destiny_global_coords - current_global_coords,0,2));
            
            old_delta_xi = DeltaXi;
            DeltaXi = RHS[0]/J(0, 0);
            
            ResultingPoint[0] += DeltaXi;
            
            if (ResultingPoint[0] <= -1.0)
            {
                ResultingPoint[0] = -1.0;
            }
            else if (ResultingPoint[0] >= 1.0)
            {
                ResultingPoint[0] = 1.0;
            }
            
            if ( std::abs(DeltaXi - old_delta_xi) < Tolerance )
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
     * @return Distance: The distance between the point and the plane
     */
    
    static inline void Project(
        const PointType& PointOrigin,
        const PointType& PointDestiny,
        PointType& PointProjected,
        double& Distance,
        const array_1d<double,3>& Normal
        )
    {
        array_1d<double,3> vector_points = PointDestiny.Coordinates() - PointOrigin.Coordinates();

        Distance = inner_prod(vector_points, Normal); 

        PointProjected.Coordinates() = PointDestiny.Coordinates() - Normal * Distance;
    }
    
        
    /**
     * Project a point over a plane (avoiding some steps)
     * @param PointOrigin: A point in the plane
     * @param PointDestiny: The point to be projected
     * @param Normal: The normal of the plane
     * @return PointProjected: The point pojected over the plane
     */
    
    static inline PointType FastProject(
        const PointType& PointOrigin,
        const PointType& PointDestiny,
        const array_1d<double,3>& Normal
        )
    {
        array_1d<double,3> vector_points = PointDestiny.Coordinates() - PointOrigin.Coordinates();

        const double distance = inner_prod(vector_points, Normal); 
        
        PointType point_projected;
        point_projected.Coordinates() = PointDestiny.Coordinates() - Normal * distance;
        
        return point_projected;
    }
    
    /**
     * This functions checks if the length of the line is to short, with the potential of provoque ill condition in the dual LM formulation
     * @param GeometryLine: The line to be checked
     * @param Tolerance: The threshold length
     * @return True if the line is too short, false otherwise
     */
    
    static inline bool LengthCheck(
        const GeometryPointType& GeometryLine,
        const double Tolerance = 1.0e-6
        )
    {
        const double lx = GeometryLine[0].X() - GeometryLine[1].X();
        const double ly = GeometryLine[0].Y() - GeometryLine[1].Y();

        const double length = std::sqrt(lx * lx + ly * ly);
        
        return (length < Tolerance) ? true : false;
    }
    
    /**
     * This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param GeometryTriangle: The triangle to be checked
     * @return True if the triangle is in bad shape, false otherwise
     */
    
    static inline bool HeronCheck(const GeometryPointType GeometryTriangle)
    {
        return HeronCheck(GeometryTriangle[0], GeometryTriangle[1], GeometryTriangle[2]);
    }
    
    /**
     * This functions checks if the semiperimeter is smaller than any of the sides of the triangle
     * @param PointOrig1: The triangle first point
     * @param PointOrig2: The triangle second point
     * @param PointOrig3: The triangle third point
     * @return True if the triangle is in bad shape, false otherwise
     */
    
    static inline bool HeronCheck(        
        const PointType& PointOrig1,
        const PointType& PointOrig2,
        const PointType& PointOrig3
        )
    {
        const double a = MathUtils<double>::Norm3(PointOrig1.Coordinates()-PointOrig2.Coordinates());
        const double b = MathUtils<double>::Norm3(PointOrig2.Coordinates()-PointOrig3.Coordinates());
        const double c = MathUtils<double>::Norm3(PointOrig3.Coordinates()-PointOrig1.Coordinates());
      
        const double s = 0.5 * (a + b + c);
        const double A2 = s * (s - a) * (s - b) * (s - c);
        
        const bool Check = A2 <= 0.0 ? true : false;  // We consider as bad shaped the ones with no area or negative A2 (semiperimeter smaller than any side)
        
//         // Debug
//         std::cout << Check << " A2: " << A2 << std::endl;
//         if (Check == true)
//         {
//             std::cout << "Warning:: The triangle is in bad shape" << std::endl;
//             std::cout << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << PointOrig1.X() << "," << PointOrig1.Y() << "," << PointOrig1.Z()  << "},{" << PointOrig2.X() << "," << PointOrig2.Y() << "," << PointOrig2.Z()  << "},{" << PointOrig3.X() << "," << PointOrig3.Y() << "," << PointOrig3.Z()  << "}}]}]" << std::endl;
//         }
        
        return Check;
    }
    
    /**
     * This function scales the points according to a factor (to increase the bounding box)
     * @param PointToScale: The point to scale
     * @param Center: The reference point
     * @param ScaleFactor: The factor considered to "grow" the node
     */
    
    template<class TPointType>
    static inline void ScaleNode(
        TPointType& PointToScale,
        const array_1d<double, 3>& Normal,
        const double& LengthSearch
        )
    {        
        PointToScale.Coordinates() = PointToScale.Coordinates() + Normal * LengthSearch;
    }
    
    /**
     * Calculates the distance between nodes
     * @param PointOrigin: The first node
     * @param PointDestiny: The second node
     */
    
    static inline double DistancePoints(
        const GeometryType::CoordinatesArrayType& PointOrigin,
        const GeometryType::CoordinatesArrayType& PointDestiny
        )
    {
        return std::sqrt((PointOrigin[0] - PointDestiny[0]) * (PointOrigin[0] - PointDestiny[0])
                       + (PointOrigin[1] - PointDestiny[1]) * (PointOrigin[1] - PointDestiny[1])
                       + (PointOrigin[2] - PointDestiny[2]) * (PointOrigin[2] - PointDestiny[2]));
    }

    /**
     * This function calculates the normal in a specific GP with a given shape function
     * @param N: The shape function considered
     * @param Geom: The geometry of condition of interest
     */

    static inline array_1d<double,3> GaussPointNormal(
        const Vector& N,
        const GeometryType& Geom
        )
    {
        array_1d<double,3> normal = ZeroVector(3);
        for( unsigned int i_node = 0; i_node < Geom.PointsNumber(); ++i_node )
        {
            normal += N[i_node] * Geom[i_node].GetValue(NORMAL); 
        }
        
        if (norm_2(normal) > std::numeric_limits<double>::epsilon())
        {
            normal = normal/norm_2(normal); // It is suppossed to be already unitary (just in case)
        }
        
        return normal;
    }
    
    /**
     * It computes the mean of the normal in the condition in all the nodes
     * @param ModelPart: The model part to compute
     * @return The modelparts with the normal computed
     */
    
    static inline void ComputeNodesMeanNormalModelPart(ModelPart& rModelPart) 
    {
        // Tolerance
        const double& tolerance = std::numeric_limits<double>::epsilon();

        // Initialize normal vectors
        const array_1d<double,3> zero_vect = ZeroVector(3);
        
        NodesArrayType& nodes_array = rModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size()); 
        
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            it_node->SetValue(NORMAL, zero_vect);
        }
        
        // Sum all the nodes normals
        ConditionsArrayType& conditions_array = rModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            if (it_cond->Is(SLAVE) || it_cond->Is(MASTER) || it_cond->Is(ACTIVE))
            {
                // Aux coordinates
                CoordinatesArrayType aux_coords;
                aux_coords = it_cond->GetGeometry().PointLocalCoordinates(aux_coords, it_cond->GetGeometry().Center());
                array_1d<double, 3>& rNormal = it_cond->GetValue(NORMAL);
                rNormal = it_cond->GetGeometry().Normal(aux_coords);
                
                const unsigned int number_nodes = it_cond->GetGeometry().PointsNumber();
                
                for (unsigned int i = 0; i < number_nodes; i++)
                {
                    #pragma omp critical
                    noalias( it_cond->GetGeometry()[i].GetValue(NORMAL) ) += rNormal;
                }
            }
        }
        
        if (rModelPart.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] == true)
        {
            // Applied laziness - MUST be calculated BEFORE normalizing the normals
            ComputeDeltaNodesMeanNormalModelPart( rModelPart ); 
        }

        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;

            const double norm_normal = norm_2(it_node->GetValue(NORMAL));
            
            if (norm_normal > tolerance)
            {
                it_node->GetValue(NORMAL) /= norm_normal;
            }
        }
    }
    
    /**
     * It computes the mean of the normal in the condition in all the nodes using the area to weight it
     * @param ModelPart: The model part to compute
     * @return The modelparts with the normal computed
     */
    
    static inline void ComputeNodesMeanNormalAreaWeightedModelPart(ModelPart& rModelPart) 
    {
        // Tolerance
        const double& tolerance = std::numeric_limits<double>::epsilon();

        // Initialize normal vectors
        const array_1d<double,3> zero_vect = ZeroVector(3);
        
        NodesArrayType& nodes_array = rModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size()); 
        
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            it_node->SetValue(NODAL_AREA, 0.0);
            it_node->SetValue(NORMAL, zero_vect);
        }
        
        // Aux coordinates
        CoordinatesArrayType aux_coords;
        aux_coords.clear();
        
        // Sum all the nodes normals
        ConditionsArrayType& conditions_array = rModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            if (it_cond->Is(SLAVE) || it_cond->Is(MASTER) || it_cond->Is(ACTIVE))
            {
                aux_coords = it_cond->GetGeometry().PointLocalCoordinates(aux_coords, it_cond->GetGeometry().Center());
                array_1d<double, 3>& rNormal = it_cond->GetValue(NORMAL);
                rNormal = it_cond->GetGeometry().Normal(aux_coords);
                
                const unsigned int number_nodes = it_cond->GetGeometry().PointsNumber();
                const double & rArea = it_cond->GetGeometry().Area()/number_nodes;
                
                for (unsigned int i = 0; i < number_nodes; i++)
                {
                    #pragma omp atomic
                    it_cond->GetGeometry()[i].GetValue(NODAL_AREA)        += rArea;
                    #pragma omp critical
                    noalias( it_cond->GetGeometry()[i].GetValue(NORMAL) ) += rArea * rNormal;
                }
            }
        }
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;

            const double& total_area = it_node->GetValue(NODAL_AREA);
            if (total_area > tolerance)
            {
                it_node->GetValue(NORMAL) /= total_area;
            }
        }
        
        if (rModelPart.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] == true)
        {
            // Applied laziness - MUST be calculated BEFORE normalizing the normals
            ComputeDeltaNodesMeanNormalModelPart( rModelPart ); 
        }

        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;

            const double norm_normal = norm_2(it_node->GetValue(NORMAL));
            
            if (norm_normal > tolerance)
            {
                it_node->GetValue(NORMAL) /= norm_normal;
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
        ConditionsArrayType& p_cond  = rModelPart.Conditions();
        ConditionsArrayType::iterator it_cond_begin = p_cond.ptr_begin();
        ConditionsArrayType::iterator it_cond_end   = p_cond.ptr_end();

        // Tolerance
        const double& tolerance = std::numeric_limits<double>::epsilon();

        // Initialize directional derivative
        const unsigned int dimension = it_cond_begin->WorkingSpaceDimension( ); 
        const Matrix zero_delta_normal = ZeroMatrix( dimension, dimension );

        Matrix Delta_ne_adj  = Matrix(dimension, dimension);
        Matrix Ce = Matrix(dimension, dimension);
        
        const Matrix I = IdentityMatrix(dimension, dimension);

        NodesArrayType& nodes_array = rModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size()); 
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            it_node->SetValue(DELTA_NORMAL, zero_delta_normal);
        }
        
        // Sum the directional derivatives of the adjacent segments
        for(ConditionsArrayType::iterator it_cond = it_cond_begin; it_cond!=it_cond_end; it_cond++)
        {
            if (it_cond->Is(ACTIVE) || it_cond->Is(MASTER)) 
            {
                const array_1d<double, 3>& ne = it_cond->GetValue(NORMAL);   // normalized condition normal
                Matrix ne_o_ne = subrange( outer_prod( ne, ne ), 0, dimension, 0, dimension );

                const unsigned int num_nodes = it_cond->GetGeometry( ).PointsNumber();
                if( dimension == 2 )
                {
                    if (num_nodes == 2)
                    {
                        const double ne_norm = it_cond->GetGeometry( ).Length( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                        
                        Delta_ne_adj( 0, 0 ) =  0.0;
                        Delta_ne_adj( 0, 1 ) = -1.0;
                        Delta_ne_adj( 1, 0 ) =  1.0;
                        Delta_ne_adj( 1, 1 ) =  0.0;
                        
                        Ce = prod( I - ne_o_ne, Delta_ne_adj ) / ne_norm;     // In 2D, Delta_ne_adj is node-independent => evaluated outside the nodes loop
                        
                        for (unsigned int i = 0; i < num_nodes; i++)
                        {
                            NodeType& node_j = it_cond->GetGeometry( )[i];
                            
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
                    const double ne_norm = it_cond->GetGeometry( ).Area( ); // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
                    
                    for (unsigned int i = 0; i < num_nodes; i++)
                    {
                        Matrix J = ZeroMatrix( 3, 2 ); // Jacobian [ 3D global x 2D local ]
                        array_1d<double, 2> DN_De_j        = ZeroVector( 2 );  // Shape functions derivatives for node j [ DN_Dxi_j, DN_Deta_j ]
                        array_1d<double, 3> local_coords_j = ZeroVector( 3 );
                        
                        NodeType& node_j = it_cond->GetGeometry( )[i];
                        
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
                            KRATOS_ERROR << "DELTA_NORMAL is not yet defined for higher order 2D elements. Number of nodes: " << it_cond->GetGeometry( ).PointsNumber() << std::endl;
                        }
                        
                        it_cond->GetGeometry( ).Jacobian( J, local_coords_j );
                        
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
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            const array_1d<double, 3>& nj = it_node->GetValue(NORMAL); // nodal non-normalized normal (this function is called before normalization)
            
            Matrix nj_o_nj = subrange( outer_prod( nj, nj ), 0, dimension, 0, dimension );
            const double nj_norm = norm_2( nj );
            const double nj_norm_3 = nj_norm * nj_norm * nj_norm;
            
            if ( nj_norm_3 > tolerance )
            {
                const Matrix Cj = I / nj_norm - nj_o_nj / nj_norm_3;
                it_node->SetValue(DELTA_NORMAL, prod( Cj, it_node->GetValue(DELTA_NORMAL) ));
            }
        }
    }
    
    /**
     * It calculates the center updated in u_n+1/2
     * @param ThisGeometry: The geometry to calculate
     * @return point: The center in u_n+1/2 (Newmark)
     */
    
    static inline Point<3> GetHalfJumpCenter(
        GeometryType& ThisGeometry,
        const double& DeltaTime
        )
    {
        Point<3> center = ThisGeometry.Center();
        
        // Initialize variables
        Vector N;
        GeometryType::CoordinatesArrayType local_point;
        
        // Get shape functions
        ThisGeometry.PointLocalCoordinates( local_point, center.Coordinates() );
        ThisGeometry.ShapeFunctionsValues( N, local_point );
        
        const Matrix new_delta_disp = 0.25 * DeltaTime * (GetVariableMatrix(ThisGeometry, VELOCITY, 0) + GetVariableMatrix(ThisGeometry, VELOCITY, 1)) + 0.125 * DeltaTime * DeltaTime * GetVariableMatrix(ThisGeometry, ACCELERATION, 1);
        
        const Vector new_delta_disp_center = prod(trans(new_delta_disp), N);
        
        center.Coordinates() += new_delta_disp_center;
        
        return center;
    }
    
    /**
     * It calculates the matrix of coordinates of a geometry
     * @param nodes: The geometry to calculate
     * @param current: If we calculate the current coordinates or the initial ones
     * @return coordinates: The matrix containing the coordinates of the geometry
     */
    
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline bounded_matrix<double, TNumNodes, TDim> GetCoordinates(
        const GeometryType& Nodes,
        const bool Current = true,
        const unsigned int Step = 0
        )
    {
        /* DEFINITIONS */            
        bounded_matrix<double, TNumNodes, TDim> coordinates;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            array_1d<double, 3> coord;
            
            if (Current == true)
            {
                coord = Nodes[i_node].Coordinates();
            }
            else
            {
                coord = Nodes[i_node].GetInitialPosition();
                
                if (Step > 0)
                {
                    coord += Nodes[i_node].FastGetSolutionStepValue(DISPLACEMENT, Step);
                }
            }

            for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
            {
                coordinates(i_node, i_dof) = coord[i_dof];
            }
        }
        
        return coordinates;
    }

    /**
     * It calculates the vector of an historical variable of a geometry
     * @param Nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return var_vector: The vector containing the variables of the geometry
     */
    
    template< unsigned int TNumNodes >
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& Nodes,
        const Variable<double>& rVarName,
        const unsigned int Step = 0
        )
    {
        /* DEFINITIONS */        
        array_1d<double, TNumNodes> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            var_vector[i_node] = Nodes[i_node].FastGetSolutionStepValue(rVarName, Step);
        }
        
        return var_vector;
    }
    
    /**
     * It calculates the vector of an historical variable of a geometry
     * @param Nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return var_vector: The vector containing the variables of the geometry
     */
        
    template< unsigned int TNumNodes >
    static inline bounded_matrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& Nodes,
        const Variable<double>& rVarName,
        const unsigned int Step = 0
        )
    {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, 1> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            var_vector(i_node, 0) = Nodes[i_node].FastGetSolutionStepValue(rVarName, Step);
        }
        
        return var_vector;
    }

    /**
     * It calculates the vector of a non-historical variable of a geometry
     * @param Nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @return var_vector: The vector containing the variables of the geometry
     */
        
    template< unsigned int TNumNodes >
    static inline array_1d<double, TNumNodes> GetVariableVector(
        const GeometryType& Nodes,
        const Variable<double>& rVarName
        )
    {
        /* DEFINITIONS */        
        array_1d<double, TNumNodes> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            var_vector[i_node] = Nodes[i_node].GetValue(rVarName);
        }
        
        return var_vector;
    }
    
    /**
     * It calculates the vector of a non-historical variable of a geometry
     * @param Nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @return var_vector: The vector containing the variables of the geometry
     */
    
    template< unsigned int TNumNodes >
    static inline bounded_matrix<double, TNumNodes, 1> GetVariableVectorMatrix(
        const GeometryType& Nodes,
        const Variable<double>& rVarName
        )
    {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, 1> var_vector;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            var_vector(i_node, 0) = Nodes[i_node].GetValue(rVarName);
        }
        
        return var_vector;
    }
    
    /**
     * It calculates the matrix of a variable of a geometry
     * @param Nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return var_matrix: The matrix containing the variables of the geometry
     */
    
    static inline Matrix GetVariableMatrix(
        const GeometryType& Nodes,
        const Variable<array_1d<double,3> >& rVarName,
        const unsigned int& Step
        )
    {
        /* DEFINITIONS */        
        const std::size_t num_nodes = Nodes.size();
        const std::size_t dim = Nodes.WorkingSpaceDimension();
        Matrix var_matrix(num_nodes, dim);
        
        for (unsigned int i_node = 0; i_node < num_nodes; i_node++)
        {
            const array_1d<double, 3> value = Nodes[i_node].FastGetSolutionStepValue(rVarName, Step);
            for (unsigned int i_dof = 0; i_dof < dim; i_dof++)
            {
                var_matrix(i_node, i_dof) = value[i_dof];
            }
        }
        
        return var_matrix;
    }
    
    /**
     * It calculates the matrix of a variable of a geometry
     * @param Nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @param step: The step where calculate
     * @return var_matrix: The matrix containing the variables of the geometry
     */
    
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline Matrix GetVariableMatrix(
        const GeometryType& Nodes,
        const Variable<array_1d<double,3> >& rVarName,
        const unsigned int Step
        )
    {
        /* DEFINITIONS */        
        Matrix var_matrix(TNumNodes, TDim);
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            const array_1d<double, 3> value = Nodes[i_node].FastGetSolutionStepValue(rVarName, Step);
            for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
            {
                var_matrix(i_node, i_dof) = value[i_dof];
            }
        }
        
        return var_matrix;
    }

    /**
     * It calculates the matrix of a non-historical variable of a geometry
     * @param Nodes: The geometry to calculate
     * @param rVarName: The name of the variable to calculate
     * @return var_matrix: The matrix containing the variables of the geometry
     */
        
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline Matrix GetVariableMatrix(
        const GeometryType& Nodes,
        const Variable<array_1d<double,3> >& rVarName
        )
    {
        /* DEFINITIONS */        
        Matrix var_matrix(TNumNodes, TDim);
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            const array_1d<double, 3>& value = Nodes[i_node].GetValue(rVarName);
            for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
            {
                var_matrix(i_node, i_dof) = value[i_dof];
            }
        }
        
        return var_matrix;
    }
    
    /**
     * It calculates the matrix containing the absolute value of another matrix
     * @param InputMatrix: The original matrix
     * @return AbsMatrix: The matrix containing the absolute value of another matrix
     */
        
    template< unsigned int TDim, unsigned int TNumNodes>
    static inline bounded_matrix<double, TNumNodes, TDim> GetAbsMatrix(const bounded_matrix<double, TNumNodes, TDim> InputMatrix)
    {
        /* DEFINITIONS */        
        bounded_matrix<double, TNumNodes, TDim> AbsMatrix;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
        {
            for (unsigned int i_dof = 0; i_dof < TDim; i_dof++)
            {
                AbsMatrix(i_node, i_dof) = std::abs(InputMatrix(i_node, i_dof));
            }
        }
        
        return AbsMatrix;
    }
    
private:
};// class ContactUtilities

///@name Explicit Specializations
///@{

}
#endif /* KRATOS_CONTACT_UTILITIES defined */
