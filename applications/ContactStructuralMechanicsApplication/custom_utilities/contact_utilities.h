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
     * This function scales the points according to a factor (to increase the bounding box)
     * @param PointToScale The point to scale
     * @param Normal The normal of the point
     * @param LengthSearch The factor considered to "grow" the node
     */
    
    template<class TPointType>
    static inline void ScaleNode(
        TPointType& PointToScale,
        const array_1d<double, 3>& Normal,
        const double LengthSearch
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
     * It calculates the center updated in u_n+1 or u_n+1/2
     * @param rThisModelPart The modelpart to update
     * @param DeltaTime The increment of time considered
     */
    
    static inline void ComputeStepJump(
        ModelPart& rThisModelPart,
        const double DeltaTime,
        const bool HalfJump = true
        )
    {
        // Time constants 
        const double velocity_constant = HalfJump ? 0.25 : 0.5;     
        const double acceleration_constant = HalfJump ? 0.125 : 0.5;
        
        // Iterate over the nodes
        NodesArrayType& nodes_array = rThisModelPart.Nodes();
    
        // We compute the half jump
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) 
        {
            auto it_node = nodes_array.begin() + i;
            array_1d<double, 3> new_delta_disp = velocity_constant * DeltaTime * (it_node->FastGetSolutionStepValue(VELOCITY) + it_node->FastGetSolutionStepValue(VELOCITY, 1)) + acceleration_constant * std::pow(DeltaTime, 2) * it_node->FastGetSolutionStepValue(ACCELERATION);
            if (it_node->IsFixed(DISPLACEMENT_X)) new_delta_disp[0] = 0.0;
            if (it_node->IsFixed(DISPLACEMENT_Y)) new_delta_disp[1] = 0.0;
            if (it_node->IsFixed(DISPLACEMENT_Z)) new_delta_disp[2] = 0.0;
            it_node->SetValue(DELTA_COORDINATES, new_delta_disp);
        }
    }
    
    /**
     * It calculates the center updated in u_n+1/2
     * @param ThisGeometry The geometry to calculate
     * @return point: The center in u_n+1/2 (Newmark)
     */
    
    static inline array_1d<double, 3> GetHalfJumpCenter(GeometryType& ThisGeometry)
    {
        array_1d<double, 3> center = (ThisGeometry.Center()).Coordinates();
        
        // Initialize variables
        Vector N;
        GeometryType::CoordinatesArrayType local_point;
        
        // Get shape functions
        ThisGeometry.PointLocalCoordinates( local_point, center );
        ThisGeometry.ShapeFunctionsValues( N, local_point );
        
    #ifdef KRATOS_DEBUG
        KRATOS_ERROR_IF(ThisGeometry[0].Has(DELTA_COORDINATES) == false) << "WARNING:: Please call ComputeStepJump() first" << std::endl;
    #endif

        const Vector new_delta_disp_center = prod(trans(GetVariableMatrix(ThisGeometry, DELTA_COORDINATES)), N);
        
        for (unsigned int i = 0; i < new_delta_disp_center.size(); ++i)
            center[i] += new_delta_disp_center[i];
        
        return center;
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
        const Variable<array_1d<double,3> >& rVarName
        ) 
    { 
        /* DEFINITIONS */         
        const std::size_t num_nodes = Nodes.size(); 
        const std::size_t dim = Nodes.WorkingSpaceDimension(); 
        Matrix var_matrix(num_nodes, dim); 
         
        for (unsigned int i_node = 0; i_node < num_nodes; i_node++) 
        { 
            const array_1d<double, 3> value = Nodes[i_node].GetValue(rVarName); 
            for (unsigned int i_dof = 0; i_dof < dim; i_dof++) 
                var_matrix(i_node, i_dof) = value[i_dof]; 
        } 
         
        return var_matrix; 
    } 
    
private:
};// class ContactUtilities

///@name Explicit Specializations
///@{

}
#endif /* KRATOS_CONTACT_UTILITIES defined */
