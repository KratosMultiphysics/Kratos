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
#include "utilities/mortar_utilities.h"

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
    typedef Point                                               PointType;
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
     * Project a point over a plane
     * @param PointOrigin A point in the plane
     * @param PointDestiny The point to be projected
     * @param Normal The normal of the plane
     * @param PointProjected The point pojected over the plane
     * @param Distance The distance between the point and the plane
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
     * This function scales the points according to a factor (to increase the bounding box)
     * @param PointToScale The point to scale
     * @param Normal The normal of the point
     * @param LengthSearch The factor considered to "grow" the node
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
     * @param PointOrigin The first node
     * @param PointDestiny The second node
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
     * It computes the mean of the normal in the condition in all the nodes
     * @param rModelPart The model part to compute
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
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;
            it_node->SetValue(NORMAL, zero_vect);
        }
        
        // Sum all the nodes normals
        ConditionsArrayType& conditions_array = rModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for
        for(int i = 0; i < num_conditions; ++i) 
        {
            auto it_cond = conditions_array.begin() + i;
            GeometryType& this_geometry = it_cond->GetGeometry();
            
            // Aux coordinates
            CoordinatesArrayType aux_coords;
            aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_geometry.Center());
            it_cond->SetValue(NORMAL, this_geometry.UnitNormal(aux_coords));
//             const array_1d<double, 3>& normal = it_cond->GetValue(NORMAL);
            
            const unsigned int number_nodes = this_geometry.PointsNumber();
            
            for (unsigned int i = 0; i < number_nodes; ++i)
            {
                auto& this_node = this_geometry[i];
                aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_node.Coordinates());
                const array_1d<double, 3>& normal = this_geometry.UnitNormal(aux_coords);
                auto& aux_normal = this_node.GetValue(NORMAL);
                for (unsigned int index = 0; index < 3; ++index)
                {
                    #pragma omp atomic
                    aux_normal[index] += normal[index];
                }
            }
        }
        
        if (rModelPart.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] == true)
        {
            // Applied laziness - MUST be calculated BEFORE normalizing the normals
            ComputeDeltaNodesMeanNormalModelPart( rModelPart ); 
        }

        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;

            const double norm_normal = norm_2(it_node->GetValue(NORMAL));
            if (norm_normal > tolerance) it_node->GetValue(NORMAL) /= norm_normal;
            else KRATOS_ERROR << "WARNING:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
        }
    }
    
    /**
     * It computes the mean of the normal in the condition in all the nodes using the area to weight it
     * @param rModelPart The model part to compute
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
        for(int i = 0; i < num_nodes; ++i) 
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
        for(int i = 0; i < num_conditions; ++i) 
        {
            auto it_cond = conditions_array.begin() + i;
            GeometryType& this_geometry = it_cond->GetGeometry();
            
            aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_geometry.Center());
            it_cond->SetValue(NORMAL, this_geometry.UnitNormal(aux_coords));
            const array_1d<double, 3>& normal = it_cond->GetValue(NORMAL);
            
            const unsigned int number_nodes = this_geometry.PointsNumber();
            const double & rArea = this_geometry.Area()/number_nodes;
            
            for (unsigned int i = 0; i < number_nodes; ++i)
            {
                auto& this_node = this_geometry[i];
                double& nodal_area = this_node.GetValue(NODAL_AREA);
                #pragma omp atomic
                nodal_area += rArea;
                auto& aux_normal = this_node.GetValue(NORMAL);
                for (unsigned int index = 0; index < 3; ++index)
                {
                    #pragma omp atomic
                    aux_normal[index] += normal[index];
                }
            }
        }
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;
            const double& total_area = it_node->GetValue(NODAL_AREA);
            if (total_area > tolerance) it_node->GetValue(NORMAL) /= total_area;
        }
        
        if (rModelPart.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] == true)
        {
            // Applied laziness - MUST be calculated BEFORE normalizing the normals
            ComputeDeltaNodesMeanNormalModelPart( rModelPart ); 
        }

        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;

            const double norm_normal = norm_2(it_node->GetValue(NORMAL));
            
            if (norm_normal > tolerance) it_node->GetValue(NORMAL) /= norm_normal;
            else KRATOS_ERROR << "WARNING:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
        }
    }

    /**
     * It computes the directional derivative of the normal in the condition in all the nodes
     * @param rModelPart The model part to compute
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
        for(int i = 0; i < num_nodes; ++i) 
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
                        
                        for (unsigned int i = 0; i < num_nodes; ++i)
                        {
                            NodeType& node_j = it_cond->GetGeometry( )[i];
                            
                            // -/+ 0.5 are the values of DN_Dxi for linear line elements at nodes 1 and 2 - no need to call the function
                            const double DN_De_j = ( i == 0 ) ? -0.5 : 0.5;
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
                    
                    for (unsigned int i = 0; i < num_nodes; ++i)
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
        for(int i = 0; i < num_nodes; ++i) 
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
     * @param ThisGeometry The geometry to calculate
     * @param DeltaTime The increment of time considered
     * @return point: The center in u_n+1/2 (Newmark)
     */
    
    static inline array_1d<double, 3> GetHalfJumpCenter(
        GeometryType& ThisGeometry,
        const double& DeltaTime
        )
    {
        PointType center = ThisGeometry.Center();
        
        // Initialize variables
        Vector N;
        GeometryType::CoordinatesArrayType local_point;
        
        // Get shape functions
        ThisGeometry.PointLocalCoordinates( local_point, center.Coordinates() );
        ThisGeometry.ShapeFunctionsValues( N, local_point );
        
        const Matrix new_delta_disp = 0.25 * DeltaTime * (GetVariableMatrix(ThisGeometry, VELOCITY, 0) + GetVariableMatrix(ThisGeometry, VELOCITY, 1)) + 0.125 * DeltaTime * DeltaTime * GetVariableMatrix(ThisGeometry, ACCELERATION, 1);
        
        const Vector new_delta_disp_center = prod(trans(new_delta_disp), N);
        
        for (unsigned int i = 0; i < new_delta_disp_center.size(); ++i)
        {
            center.Coordinates()[i] += new_delta_disp_center[i];
        }
        
        return center.Coordinates();
    }
    
         
    /** 
     * It calculates the matrix of a variable of a geometry 
     * @param Nodes The geometry to calculate 
     * @param rVarName The name of the variable to calculate 
     * @param Step The step where calculate 
     * @return var_matrix: The matrix containing the variables of the geometry 
     */ 
     
    static inline Matrix GetVariableMatrix( 
        const GeometryType& Nodes, 
        const Variable<array_1d<double,3> >& rVarName, 
        const unsigned int Step 
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
    
private:
};// class ContactUtilities

///@name Explicit Specializations
///@{

}
#endif /* KRATOS_CONTACT_UTILITIES defined */
