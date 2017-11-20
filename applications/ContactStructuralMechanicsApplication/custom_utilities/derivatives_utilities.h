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

#if !defined(KRATOS_DERIVATIVES_UTILITIES)
#define KRATOS_DERIVATIVES_UTILITIES

#include "contact_structural_mechanics_application_variables.h"

/* Includes */
#include "includes/model_part.h"
#include "includes/mortar_classes.h"

/* Utilities */
#include "utilities/mortar_utilities.h"
#include "utilities/math_utils.h"

/* Geometries */
#include "geometries/point.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"

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
    
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
class DerivativesUtilities
{
public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Properties PropertiesType;
    
    typedef typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N, typename std::conditional<TNumNodes == 3, PointBelongsTriangle3D3N, PointBelongsQuadrilateral3D4N>::type>::type BelongType;
    
    typedef PointBelong<TNumNodes>                                                 PointBelongType;
    
    typedef Geometry<PointBelongType>                                      GeometryPointBelongType;
    
    typedef array_1d<PointBelongType,TDim>                                      ConditionArrayType;
    
    typedef typename std::vector<ConditionArrayType>                        ConditionArrayListType;
    
    typedef Line2D2<PointType>                                                            LineType;
    
    typedef Triangle3D3<PointType>                                                    TriangleType;
    
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type  DecompositionType;
    
    typedef typename std::conditional<TFrictional == true, DerivativeDataFrictional<TDim, TNumNodes>, DerivativeData<TDim, TNumNodes> >::type DerivativeDataType;
    
    typedef MortarKinematicVariablesWithDerivatives<TDim, TNumNodes>              GeneralVariables;
    
    typedef DualLagrangeMultiplierOperatorsWithDerivatives<TDim, TNumNodes, TFrictional>    AeData;
    
    typedef MortarOperatorWithDerivatives<TDim, TNumNodes, TFrictional>    MortarConditionMatrices;
    
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
     * This method is used to compute the directional derivatives of the jacobian determinant
     * @param DecompGeom The triangle used to decompose the geometry
     * @param rVariables The kinematic variables
     * @param rDerivativeData The derivatives container
     */
    static inline void CalculateDeltaDetjSlave(
        const DecompositionType& DecompGeom,
        GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData
        )
    {
        if (TDim == 2)
        {
            // Fill up the elements corresponding to the slave DOFs - the rest remains zero
            for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
            {
                rDerivativeData.DeltaDetjSlave[i    ] = rVariables.jSlave( 0, 0 ) * rVariables.DNDeSlave( i_slave, 0) / rVariables.DetjSlave;
                rDerivativeData.DeltaDetjSlave[i + 1] = rVariables.jSlave( 1, 0 ) * rVariables.DNDeSlave( i_slave, 0) / rVariables.DetjSlave;
            }
        }
        else
        {
            const array_1d<double, 3>& x1cell = DecompGeom[0].Coordinates();
            const array_1d<double, 3>& x2cell = DecompGeom[1].Coordinates();
            const array_1d<double, 3>& x3cell = DecompGeom[2].Coordinates();
            
            const array_1d<double, 3>& x21cell = x2cell - x1cell;
            const array_1d<double, 3>& x31cell = x3cell - x1cell;
            
            array_1d<double, 3> aux_cross_product;
            MathUtils<double>::CrossProduct(aux_cross_product, x21cell, x31cell);
            aux_cross_product /= rVariables.DetjSlave;
//             aux_cross_product /= norm_2(aux_cross_product);
            
            for ( unsigned int i_node = 0; i_node < 2 * TNumNodes; ++i_node )
            {
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof) 
                {
                    const bounded_matrix<double, 3, 3>& local_delta_vertex = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof];
                    array_1d<double, 3> aux_delta_cross_product1, aux_delta_cross_product2;
                    
                    MathUtils<double>::CrossProduct(aux_delta_cross_product1, row(local_delta_vertex, 1) - row(local_delta_vertex, 0), x31cell);
                    MathUtils<double>::CrossProduct(aux_delta_cross_product2, x21cell, row(local_delta_vertex, 2) - row(local_delta_vertex, 0));
                    
                    rDerivativeData.DeltaDetjSlave[i_node * TDim + i_dof] = inner_prod(aux_cross_product, aux_delta_cross_product1) + inner_prod(aux_cross_product, aux_delta_cross_product2);
                }
            }
        }
    }
    
    /**
     * This method is used to compute the local increment of the normal
     * @param CondGeometry The geometry where the delta normal is computed
     * @param NodeIndex The index of the node of the geometry considered
     * @return The matrix containing the delta normals
     * NOTE: Not the mean, look in the contact utilities 
     */
    static inline bounded_matrix<double, TDim, TDim> LocalDeltaNormal(
        const GeometryType& CondGeometry,
        const unsigned int NodeIndex
        )
    {
        // Tolerance
        const double tolerance = std::numeric_limits<double>::epsilon();
            
        bounded_matrix<double, TDim, TDim> DeltaNeAdj;
        bounded_matrix<double, TDim, TDim> Ce;
        
        const bounded_matrix<double, TDim, TDim> I = IdentityMatrix(TDim, TDim);
        
        bounded_matrix<double, TDim, TDim> delta_normal = ZeroMatrix(TDim,TDim);
        
        // Normalized condition normal
        PointType auxiliar_center;
        auxiliar_center.Coordinates() = ZeroVector(3);
        const array_1d<double, 3>& Ne = CondGeometry.UnitNormal(auxiliar_center);
        bounded_matrix<double, TDim, TDim> NeoNe = subrange( outer_prod( Ne, Ne ), 0, TDim, 0, TDim );
        
        // Auxiliar value
        // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
        double NeNorm = (TDim == 2) ? CondGeometry.Length( ) : CondGeometry.Area( );
        
        if (TDim == 2)
        {                
            DeltaNeAdj( 0, 0 ) =  0.0;
            DeltaNeAdj( 0, 1 ) = -1.0;
            DeltaNeAdj( 1, 0 ) =  1.0;
            DeltaNeAdj( 1, 1 ) =  0.0;
            
            const double DNDej = (NodeIndex == 0) ? - 0.5 : 0.5;
            
            Ce = prod( I - NeoNe, DeltaNeAdj ) / NeNorm; // In 2D, DeltaNeAdj is node-independent => evaluated outside the nodes loop
            
            delta_normal = - 2.0 * Ce * DNDej; // NOTE: Check why - 2???!!!, it was the only wayto ensure the same value as the symbolic. You will need to repeat this in 3D            
        //         delta_normal = Ce * DNDej;     
        }
        else
        {
            MatrixType J = ZeroMatrix( 3, 2 ); // Jacobian [ 3D global x 2D local ]
            array_1d<double, 2> DNDej;
            array_1d<double, 3> LocalCoordsj;
            
            if( TNumNodes == 3 )    // linear triangle element
            {
                if( NodeIndex == 0 )
                {
                    LocalCoordsj[0] = 0.0;
                    LocalCoordsj[1] = 0.0;
                    DNDej[0] = - 1.0;
                    DNDej[1] = - 1.0;
                }
                else if( NodeIndex == 1 )
                {
                    LocalCoordsj[0] = 1.0;
                    LocalCoordsj[1] = 0.0;
                    DNDej[0] = 1.0;
                    DNDej[1] = 0.0;
                }
                else // NodeIndex == 2
                {
                    LocalCoordsj[0] = 0.0;
                    LocalCoordsj[1] = 1.0;
                    DNDej[0] = 0.0;
                    DNDej[1] = 1.0;
                }
            }
            else if( TNumNodes == 4 )    // linear quad element 
            {
                if( NodeIndex == 0 )
                {
                    LocalCoordsj[0] = - 1.0;
                    LocalCoordsj[1] = - 1.0;
                    DNDej[0] = - 0.5;
                    DNDej[1] = - 0.5;
                }
                else if( NodeIndex == 1 )
                {
                    LocalCoordsj[0] =   1.0;
                    LocalCoordsj[1] = - 1.0;
                    DNDej[0] =   0.5;
                    DNDej[1] = - 0.5;
                }
                else if( NodeIndex == 2 )
                {
                    LocalCoordsj[0] =  1.0;
                    LocalCoordsj[1] =  1.0;
                    DNDej[0] =  0.5;
                    DNDej[1] =  0.5;
                }
                else // NodeIndex == 3
                {
                    LocalCoordsj[0] = - 1.0;
                    LocalCoordsj[1] =   1.0;
                    DNDej[0] = - 0.5;
                    DNDej[1] =   0.5;
                }
            }
            
            CondGeometry.Jacobian( J, LocalCoordsj );
            
            DeltaNeAdj(0,0) = 0.0;
            DeltaNeAdj(0,1) = +J(2,1) * DNDej[0] - J(2,0) * DNDej[1]; 
            DeltaNeAdj(0,2) = -J(1,1) * DNDej[0] + J(1,0) * DNDej[1]; 
            DeltaNeAdj(1,0) = -J(2,1) * DNDej[0] + J(2,0) * DNDej[1]; 
            DeltaNeAdj(1,1) = 0.0;                   
            DeltaNeAdj(1,2) = +J(0,1) * DNDej[0] - J(0,0) * DNDej[1]; 
            DeltaNeAdj(2,0) = +J(1,1) * DNDej[0] - J(1,0) * DNDej[1]; 
            DeltaNeAdj(2,1) = -J(0,1) * DNDej[0] + J(0,0) * DNDej[1]; 
            DeltaNeAdj(2,2) = 0.0;
            
            Ce = prod( I - NeoNe, DeltaNeAdj ) / NeNorm;
            delta_normal = Ce;
        }
        
        NeNorm = norm_2( Ne );
        const double NeNorm3 = NeNorm * NeNorm * NeNorm;
        
        if ( NeNorm3 > tolerance )
        {
            const bounded_matrix<double, TDim, TDim> Cj = I / NeNorm - NeoNe / NeNorm3;
            delta_normal = prod( Cj, delta_normal );
        }
            
        return delta_normal; 
    }
    
    /**
     * Calculates the increment of the normal in the slave condition
     * @param rDerivativeData The class containing the derivatives
     * @param SlaveGeometry The geometry of the slave side
     */
    
    static inline void CalculateDeltaNormalSlave(
        DerivativeDataType& rDerivativeData,
        const GeometryType& SlaveGeometry
        )
    {
        if (TDim == 2)
        {
            for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave)
            {
                const bounded_matrix<double, TDim, TDim>& delta_normal = SlaveGeometry[i_slave].GetValue(DELTA_NORMAL);
//                 const bounded_matrix<double, TDim, TDim> delta_normal = this->LocalDeltaNormal(SlaveGeometry, i_slave);
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof) 
                {
                    row(rDerivativeData.DeltaNormalSlave[i_slave * TDim + i_dof], i_slave) = trans(column(delta_normal, i_dof)); 
                }
            }
        }
    }
    
    /**
     * Calculates the increment of the normal and in the master condition
     * @param rDerivativeData The class containing the derivatives
     * @param MasterGeometry The geometry of the master side
     */
    
    static inline void CalculateDeltaNormalMaster(
        DerivativeDataType& rDerivativeData,
        const GeometryType& MasterGeometry
        )
    {
        if (TDim == 2)
        {
            for ( unsigned int i_master = 0; i_master < 2; ++i_master )
            {
//                const bounded_matrix<double, 2, 2>& delta_normal = GetGeometry[i_master].GetValue(DELTA_NORMAL);
                const bounded_matrix<double, 2, 2> delta_normal = LocalDeltaNormal(MasterGeometry, i_master);
                for (unsigned i_dof = 0; i_dof < 2; ++i_dof) 
                {
                    row(rDerivativeData.DeltaNormalMaster[i_master * 2 + i_dof], i_master) = trans(column(delta_normal, i_dof));
                }
            }
        }
    }
    
    /**
     * This method is used to compute the directional derivatives of the cell vertex
     * @param rVariables The kinematic variables
     * @param rDerivativeData The derivatives container
     * @param TheseBelongs The belongs list used in the derivatives
     * @param ConsiderNormalVariation If consider the normal derivative
     * @param SlaveGeometry The slave geometry
     * @param MasterGeometry The master geometry
     * @param Normal The normal vector of the slave geometry
     * The  procedure will be the following in order to compute the derivative of the clipping 
     * The expression of the clipping is the following:
     *      xclipp = xs1 - num/denom * diff3
     * Being:
     *      diff1 = xs1 - xs2;
     *      diff2 = xe2 - xs2;
     *      diff3 = xe1 - xs1;
     *      num = (diff1 x diff2) · n0
     *      denom = (diff3 x diff2) · n0
     * The derivative can be defined then as:
     *     delta_xclipp = delta_xs1 - (delta_num * denom - num * delta_denom)/delta_denom**2 * diff3 - num/ denom * delta_diff3
     * And here:
     *     delta_num = num · delta_n0 + n0 · (delta_diff1 x diff2 + diff1 x delta_diff2)
     *     delta_denom = denom · delta_n0 + n0 · (delta_diff3 x diff2 + diff3 x delta_diff2)
     */
    static inline void CalculateDeltaCellVertex(
        GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        const array_1d<BelongType, TDim>& TheseBelongs,
        const bool ConsiderNormalVariation,
        GeometryType& SlaveGeometry,
        GeometryType& MasterGeometry,
        const array_1d<double, 3>& Normal
        )
    {
        // The Normal and delta Normal in the center of the element
        bounded_matrix<double, TDim, TDim> delta_normal;
        
        const double aux_nodes_coeff = static_cast<double>(TNumNodes);
        
        const Point& slave_center = SlaveGeometry.Center().Coordinates();
        
//     #ifdef KRATOS_DEBUG
//         for (unsigned i_triangle = 0; i_triangle < 3; ++i_triangle) 
//         {
//             KRATOS_WATCH(static_cast<unsigned int>(TheseBelongs[i_triangle]));
//         }
//     #endif
        
        for (unsigned i_triangle = 0; i_triangle < 3; ++i_triangle) 
        {
            if (TheseBelongs[i_triangle] >= 2 * TNumNodes) // It belongs to an intersection
            {    
                // We compute the indexes
                unsigned int belong_index_slave_start, belong_index_slave_end, belong_index_master_start, belong_index_master_end;
                ConvertAuxHashIndex(static_cast<unsigned int>(TheseBelongs[i_triangle]), belong_index_slave_start, belong_index_slave_end, belong_index_master_start, belong_index_master_end);
                
                // The coordinates should be in the projected plane
                const array_1d<double, 3>& xs1 = MortarUtilities::FastProject(slave_center, SlaveGeometry[belong_index_slave_start], Normal).Coordinates(); // Start coordinates of the first segment
                const array_1d<double, 3>& xe1 = MortarUtilities::FastProject(slave_center, SlaveGeometry[belong_index_slave_end], Normal).Coordinates(); // End coordinates of the first segment
                const array_1d<double, 3>& xs2 = MortarUtilities::FastProject(slave_center, MasterGeometry[belong_index_master_start], Normal).Coordinates(); // Start coordinates of the second segment
                const array_1d<double, 3>& xe2 = MortarUtilities::FastProject(slave_center, MasterGeometry[belong_index_master_end], Normal).Coordinates(); // End coordinates of the second segment
                
//                 const array_1d<double, 3>& xs1 = SlaveGeometry[belong_index_slave_start].Coordinates();   // Start coordinates of the first segment
//                 const array_1d<double, 3>& xe1 = SlaveGeometry[belong_index_slave_end].Coordinates();     // End coordinates of the first segment
//                 const array_1d<double, 3>& xs2 = MasterGeometry[belong_index_master_start].Coordinates(); // Start coordinates of the second segment
//                 const array_1d<double, 3>& xe2 = MasterGeometry[belong_index_master_end].Coordinates();   // End coordinates of the second segment
                
                // We define the array containing the indexes of the vertexes
                array_1d<unsigned int, 4> belong_indexes;
                belong_indexes[0] = belong_index_slave_start;
                belong_indexes[1] = belong_index_slave_end;
                belong_indexes[2] = belong_index_master_start + TNumNodes;
                belong_indexes[3] = belong_index_master_end   + TNumNodes;
                
                // We define the diffs between the extremes of segements
                const array_1d<double, 3> diff1 = xs1 - xs2;
                const array_1d<double, 3> diff2 = xe2 - xs2;
                const array_1d<double, 3> diff3 = xe1 - xs1;
                
                // We compute the denominator and numerator of the clipping
                array_1d<double, 3> aux_num, aux_denom;
                MathUtils<double>::CrossProduct(aux_num,   diff1, diff2);
                MathUtils<double>::CrossProduct(aux_denom, diff3, diff2);
                const double num   = inner_prod(aux_num,   Normal);
                const double denom = inner_prod(aux_denom, Normal);
                
                for (unsigned int i_belong = 0; i_belong < 4; ++i_belong)
                {
                    // The index of the node
                    const unsigned int belong_index = belong_indexes[i_belong];
                    
                    // We compute the delta normal
                    if (ConsiderNormalVariation == true && belong_index < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index) * (1.0/aux_nodes_coeff);
                    else delta_normal = ZeroMatrix(3, 3);
                    
                    for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                    {                    
                        bounded_matrix<double, 3, 3>& local_delta_vertex = rDerivativeData.DeltaCellVertex[belong_index * TDim + i_dof];
                        
                        // Special cases (slave nodes)
                        if (i_belong == 0) // First node of the slave
                        {
                            const double coeff = 1.0 + num/denom;
                            LocalDeltaVertex(local_delta_vertex, Normal, delta_normal, i_dof, i_triangle, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, coeff);
                        }
                        else if (i_belong == 1) // Second node of the slave
                        {
                            const double coeff = - num/denom;
                            LocalDeltaVertex(local_delta_vertex, Normal, delta_normal, i_dof, i_triangle, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, coeff);
                        }
                        
                        // We define some auxiliar coefficients
                        const double coeff1 = - 1.0/denom;
                        const double coeff2 = num/std::pow(denom, 2);
                        
                        // We add the part corresponding purely to delta normal
                        if (ConsiderNormalVariation == true)
                        {
                            noalias(row(local_delta_vertex, i_triangle)) += diff3 * coeff1 * inner_prod(aux_num,  trans(column(delta_normal, i_dof))); 
                            noalias(row(local_delta_vertex, i_triangle)) += diff3 * coeff2 * inner_prod(aux_denom, trans(column(delta_normal, i_dof))); 
                        }
                        
                        // We compute the delta diffs
                        const array_1d<double, 3> delta_diff1 = (i_belong == 0) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0) : (i_belong == 2) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - 1.0) : ZeroVector(3);
                        const array_1d<double, 3> delta_diff2 = (i_belong == 3) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0) : (i_belong == 2) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - 1.0) : ZeroVector(3);
                        const array_1d<double, 3> delta_diff3 = (i_belong == 1) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0) : (i_belong == 0) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - 1.0) : ZeroVector(3);
                        
                        // We compute now the delta num and denom
                        array_1d<double, 3> aux_cross_product;
                        MathUtils<double>::CrossProduct(aux_cross_product, delta_diff1, diff2);
                        double delta_num = inner_prod(aux_cross_product, Normal);
                        MathUtils<double>::CrossProduct(aux_cross_product, diff1, delta_diff2);
                        delta_num += inner_prod(aux_cross_product, Normal);
                        
                        MathUtils<double>::CrossProduct(aux_cross_product, delta_diff3, diff2);
                        double delta_denom = inner_prod(aux_cross_product, Normal);
                        MathUtils<double>::CrossProduct(aux_cross_product, diff3, delta_diff2);
                        delta_denom += inner_prod(aux_cross_product, Normal);
                        
                        // Finally we add the contributions of delta num and denom
                        noalias(row(local_delta_vertex, i_triangle)) += coeff1 * diff3 * delta_num; 
                        noalias(row(local_delta_vertex, i_triangle)) += coeff2 * diff3 * delta_denom; 
                    }
                }
            }
            else // It belongs to a master/slave node
            {
                const unsigned int belong_index = static_cast<unsigned int>(TheseBelongs[i_triangle]);
                
                if (ConsiderNormalVariation == true && belong_index < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index) * (1.0/aux_nodes_coeff);
                else delta_normal = ZeroMatrix(3, 3);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    bounded_matrix<double, 3, 3>& local_delta_vertex = rDerivativeData.DeltaCellVertex[belong_index * TDim + i_dof];
                    
                    LocalDeltaVertex(local_delta_vertex, Normal, delta_normal, i_dof, i_triangle, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry);
                }
            }
        }
    }
    
    /**
     * Calculates the increment of the shape functions
     * @param rVariables The kinematic variables
     * @param rDerivativeData The derivatives container
     * @param SlaveGeometry The geometry of the slave side
     * @param MasterGeometry The geometry of the master side
     * @param SlaveNormal The normal of the slave side
     * @param DecompGeom The triangle used to decompose the geometry
     * @param LocalPointDecomp The local coordinates in the decomposed geometry
     * @param LocalPointParent The local coordinates in the slave geometry
     * @param ConsiderNormalVariation If consider the normal derivative
     */
    
    static inline void CalculateDeltaN1(
        const GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        GeometryType& SlaveGeometry,
        GeometryType& MasterGeometry,
        const array_1d<double, 3> SlaveNormal,
        const DecompositionType& DecompGeom,
        const PointType& LocalPointDecomp,
        const PointType& LocalPointParent,
        const bool ConsiderNormalVariation = false
        )
    {
        /* Shape functions */
        const VectorType& N1 = rVariables.NSlave;

        /* Local gradients */
        const MatrixType& DNDe1 = rVariables.DNDeSlave;
        
        // The Normal and delta Normal in the center of the element
        bounded_matrix<double, TDim, TDim> delta_normal = ZeroMatrix(TDim, TDim);
        
        /* Shape function decomposition */
        VectorType N_decomp;
        DecompGeom.ShapeFunctionsValues( N_decomp, LocalPointDecomp.Coordinates() );
        
        if (TDim == 2)
        {
            // TODO: Finish this!!!!
        }
        else
        {
            for ( unsigned int i_node = 0; i_node < 2 * TNumNodes; ++i_node)
            {
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof) 
                {
                    array_1d<double, 3> aux_RHS1 = ZeroVector(3);
                    
                    // The vertex cell contribution
                    const auto& local_delta_cell = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof]; 
                    for(std::size_t i_belong = 0; i_belong < 3; ++i_belong)
                    {
                        noalias(aux_RHS1) += N_decomp[i_belong] * row(local_delta_cell, i_belong);
                    }
                    
                    // Local contribution
                    const array_1d<double, 3>& aux_delta_node = LocalDeltaVertex( SlaveNormal,  delta_normal, i_dof, i_node, ConsiderNormalVariation, SlaveGeometry, MasterGeometry );
                    if (i_node < TNumNodes) noalias(aux_RHS1) -= N1[i_node] * aux_delta_node;
                    
                    // We compute the delta coordinates 
                    array_1d<double, 2> aux_delta_coords1;
                    DeltaPointLocalCoordinates(aux_delta_coords1, aux_RHS1, rVariables.DNDeSlave, SlaveGeometry, SlaveNormal);
                    
                    // Now we can compute the delta shape functions
                    auto& delta_n1 = rDerivativeData.DeltaN1[i_node * TDim + i_dof];
                    noalias(delta_n1) = (aux_delta_coords1[0] * column(DNDe1, 0) + aux_delta_coords1[1] * column(DNDe1, 1));
                }
            }
        }
    }
    
    /**
     * Calculates the increment of the shape functions
     * @param rVariables The kinematic variables
     * @param rDerivativeData The derivatives container
     * @param SlaveGeometry The geometry of the slave side
     * @param MasterGeometry The geometry of the master side
     * @param SlaveNormal The normal of the slave side
     * @param MasterNormal The normal of the master side
     * @param DecompGeom The triangle used to decompose the geometry
     * @param LocalPointDecomp The local coordinates in the decomposed geometry
     * @param LocalPointParent The local coordinates in the slave geometry
     * @param ConsiderNormalVariation If consider the normal derivative
     * @param DualLM If the dual Lm formulation is considered
     */
    
    static inline void CalculateDeltaN(
        const GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        GeometryType& SlaveGeometry,
        GeometryType& MasterGeometry,
        const array_1d<double, 3> SlaveNormal,
        const array_1d<double, 3> MasterNormal,
        const DecompositionType& DecompGeom,
        const PointType& LocalPointDecomp,
        const PointType& LocalPointParent,
        const bool ConsiderNormalVariation = false,
        const bool DualLM = false
        )
    {
        /* Shape functions */
        const VectorType& N1 = rVariables.NSlave;
        const VectorType& N2 = rVariables.NMaster;
        
        /* Local gradients */
        const MatrixType& DNDe1 = rVariables.DNDeSlave;
        const MatrixType& DNDe2 = rVariables.DNDeMaster;
        
        // The Normal and delta Normal in the center of the element
        bounded_matrix<double, TDim, TDim> delta_normal = ZeroMatrix(TDim, TDim);
        
        /* Shape function decomposition */
        VectorType N_decomp;
        DecompGeom.ShapeFunctionsValues( N_decomp, LocalPointDecomp.Coordinates() );
        
        if (TDim == 2)
        {
            // TODO: Finish this!!!!
        }
        else
        {
            for ( unsigned int i_node = 0; i_node < 2 * TNumNodes; ++i_node)
            {
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof) 
                {
                    array_1d<double, 3> aux_RHS1 = ZeroVector(3);
                    
                    // The vertex cell contribution
                    const auto& local_delta_cell = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof]; 
                    for(std::size_t i_belong = 0; i_belong < 3; ++i_belong)
                    {
                        noalias(aux_RHS1) += N_decomp[i_belong] * row(local_delta_cell, i_belong);
                    }
                    array_1d<double, 3> aux_RHS2 = aux_RHS1;
                    
                    // Local contribution
                    const array_1d<double, 3>& aux_delta_node = LocalDeltaVertex( SlaveNormal,  delta_normal, i_dof, i_node, ConsiderNormalVariation, SlaveGeometry, MasterGeometry );
                    if (i_node < TNumNodes) noalias(aux_RHS1) -= N1[i_node] * aux_delta_node;
                    else noalias(aux_RHS2) -= N2[i_node - TNumNodes] * aux_delta_node;
                    
                    // We compute the delta coordinates 
                    array_1d<double, 2> aux_delta_coords1, aux_delta_coords2;
                    DeltaPointLocalCoordinates(aux_delta_coords1, aux_RHS1, rVariables.DNDeSlave, SlaveGeometry, SlaveNormal);
                    DeltaPointLocalCoordinates(aux_delta_coords2, aux_RHS2, rVariables.DNDeMaster, MasterGeometry, SlaveNormal);
                    
                    // Now we can compute the delta shape functions
                    auto& delta_n1 = rDerivativeData.DeltaN1[i_node * TDim + i_dof];
                    noalias(delta_n1) = (aux_delta_coords1[0] * column(DNDe1, 0) + aux_delta_coords1[1] * column(DNDe1, 1));
                    
                    auto& delta_n2 = rDerivativeData.DeltaN2[i_node * TDim + i_dof];
                    noalias(delta_n2) = (aux_delta_coords2[0] * column(DNDe2, 0) + aux_delta_coords2[1] * column(DNDe2, 1));
                    
                    // The derivatives of the dual shape function 
                    auto& delta_phi = rDerivativeData.DeltaPhi[i_node * TDim + i_dof];
                    noalias(delta_phi) = prod(rDerivativeData.DeltaAe[i_node * TDim + i_dof], N1);
                    if (DualLM == true) noalias(delta_phi) += prod(rDerivativeData.Ae, delta_n1);
                    else noalias(delta_phi) += delta_n1;
                }
            }
        }
    }
    
    /**
     * Returns a matrix with the increment of displacements, that can be used for compute the Jacobian reference (current) configuration
     * @param DeltaPosition The matrix with the increment of displacements 
     * @param LocalCoordinates The array containing the local coordinates of the exact integration segment
     */
    
    Matrix& CalculateDeltaPosition(
        Matrix& DeltaPosition,
        const GeometryType& ThisGeometry,
        const ConditionArrayType& LocalCoordinates
        )
    {
        KRATOS_TRY;

        DeltaPosition = ZeroMatrix(TDim, TDim);

        for ( unsigned int i_node = 0; i_node < TNumNodes; ++i_node )
        {
            const array_1d<double, 3 > delta_displacement = ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);
            
            for ( unsigned int j_node = 0; j_node < TDim; ++j_node )
            {
                Vector N;
                ThisGeometry.ShapeFunctionsValues( N, LocalCoordinates[j_node].Coordinates() );

                for ( unsigned int j_dim = 0; j_dim < TDim; ++j_dim )
                {
                    DeltaPosition(j_node, j_dim) += N[i_node] * delta_displacement[j_dim];
                }
            }
        }

        return DeltaPosition;

        KRATOS_CATCH( "" );
    }
    
    /**
     * Returns a matrix with the increment of displacements
     * @param DeltaPosition The matrix with the increment of displacements 
     * @param ThisGeometry The geometry considered 
     */
    
    static inline Matrix& CalculateDeltaPosition(
        Matrix& DeltaPosition,
        const GeometryType& ThisGeometry
        )
    {
        KRATOS_TRY;

        DeltaPosition = ZeroMatrix(TNumNodes, TDim);

        for ( unsigned int i_node = 0; i_node < TNumNodes; ++i_node )
        {
            const array_1d<double, 3 > delta_displacement = ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);
            
            for ( unsigned int i_dim = 0; i_dim < TDim; ++i_dim )
            {
                DeltaPosition(i_node, i_dim) += delta_displacement[i_dim];
            }
        }

        return DeltaPosition;

        KRATOS_CATCH( "" );
    }
    
    /**
     * Returns a vector with the increment of displacements
     */
    
    static inline void CalculateDeltaPosition(
        VectorType& DeltaPosition,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const unsigned int IndexNode
        )
    {
        KRATOS_TRY;

        if (IndexNode < TNumNodes)
        {
            DeltaPosition = SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1);
        }
        else
        {
            const unsigned int index_master = IndexNode - TNumNodes;
            DeltaPosition = MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1);
        }

        KRATOS_CATCH( "" );
    }
    
    /**
     * Returns a vector with the increment of displacements
     */
    
    static inline void CalculateDeltaPosition(
        VectorType& DeltaPosition,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const unsigned int IndexNode,
        const unsigned int iDoF
        )
    {
        KRATOS_TRY;

        DeltaPosition = ZeroVector(3);
        
        if (IndexNode < TNumNodes)
        {
            DeltaPosition[iDoF] = (SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        }
        else
        {
            const unsigned int index_master = IndexNode - TNumNodes;
            DeltaPosition[iDoF] = (MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        }

        KRATOS_CATCH( "" );
    }
    
    /**
     * Returns a double with the increment of displacements
     */
    
    static inline void CalculateDeltaPosition(
        double& DeltaPosition,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const unsigned int IndexNode,
        const unsigned int iDoF
        )
    {
        KRATOS_TRY;
    
        if (IndexNode < TNumNodes)
        {
            DeltaPosition = (SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        }
        else
        {
            const unsigned int index_master = IndexNode - TNumNodes;
            DeltaPosition = (MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        }

        KRATOS_CATCH( "" );
    }
    
    /**
     * Calculate Ae and DeltaAe matrices
     * @param SlaveGeometry The geometry of the slave side
     * @param SlaveNormal The normal of the slave side
     * @param pMasterCondition The pointer to the master side
     * @param rDerivativeData The derivatives container
     * @param rVariables The kinematic variables
     * @param ConsiderNormalVariation If consider the normal derivative
     * @param ConditionsPointsSlave The points that configure the exact decomposition of the geometry
     * @param ThisIntegrationMethod The integration method considered
     * @param AxiSymCoeff The axisymmetric coefficient
     */
    static inline bool CalculateAeAndDeltaAe(
        GeometryType& SlaveGeometry,
        const array_1d<double, 3>& SlaveNormal,
        Condition::Pointer pMasterCondition,
        DerivativeDataType& rDerivativeData,
        GeneralVariables& rVariables,
        const bool ConsiderNormalVariation,
        ConditionArrayListType& ConditionsPointsSlave,
        IntegrationMethod ThisIntegrationMethod,
        const double AxiSymCoeff = 1.0
        )
    {
        // We initilize the Ae components
        AeData rAeData;
        rAeData.Initialize();
        
        rDerivativeData.InitializeDeltaAeComponents();

        // Initialize general variables for the current master element
        rVariables.Initialize();
        
        // Master geometry
        GeometryType& master_geometry = pMasterCondition->GetGeometry();
        
        for (unsigned int i_geom = 0; i_geom < ConditionsPointsSlave.size(); ++i_geom)
        {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (unsigned int i_node = 0; i_node < TDim; ++i_node)
            {
                PointType global_point;
                SlaveGeometry.GlobalCoordinates(global_point, ConditionsPointsSlave[i_geom][i_node]);
                points_array[i_node] = boost::make_shared<PointType>(global_point);
                belong_array[i_node] = ConditionsPointsSlave[i_geom][i_node].GetBelong();
            }
            
            DecompositionType decomp_geom( points_array );
            
            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, SlaveGeometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
            
            if (bad_shape == false)
            {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );
                
                // Integrating the mortar operators
                for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); ++point_number )
                {
                    // We reset the derivatives
                    rDerivativeData.ResetDerivatives();
                                
                    // We compute the local coordinates 
                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                    PointType local_point_parent;
                    PointType gp_global;
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    SlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);
                    
                    // Calculate the kinematic variables
                    // We compute the current configuration
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    SlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);
                    
                    // SLAVE KINEMATIC COMPUTATIONS
                    SlaveGeometry.ShapeFunctionsValues( rVariables.NSlave, local_point_parent.Coordinates() );
                    SlaveGeometry.ShapeFunctionsLocalGradients( rVariables.DNDeSlave, local_point_parent );
                    rVariables.PhiLagrangeMultipliers = rVariables.NSlave;
                    
                    rVariables.jSlave = decomp_geom.Jacobian( rVariables.jSlave, local_point_decomp.Coordinates());
                    rVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );

                    // MASTER KINEMATIC COMPUTATIONS
                    PointType projected_gp_global;
                    array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, SlaveGeometry);
                    
                    GeometryType::CoordinatesArrayType slave_gp_global;
                    SlaveGeometry.GlobalCoordinates( slave_gp_global, local_point_parent );
                    MortarUtilities::FastProjectDirection( master_geometry, slave_gp_global, projected_gp_global, SlaveNormal, -gp_normal ); // The opposite direction
                    
                    GeometryType::CoordinatesArrayType projected_gp_local;
                    master_geometry.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                    master_geometry.ShapeFunctionsValues( rVariables.NMaster,    projected_gp_local );         
                    master_geometry.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );
                    rVariables.jMaster = master_geometry.Jacobian( rVariables.jMaster, projected_gp_local);
                    
                    // Update the derivative of the integration vertex (just in 3D)
                    if (TDim == 3) CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, ConsiderNormalVariation, SlaveGeometry, master_geometry, SlaveNormal);
                                    
                    // Update the derivative of DetJ
                    CalculateDeltaDetjSlave(decomp_geom, rVariables, rDerivativeData); 
                    
                    // Update the derivatives of the shape functions and the gap
                    CalculateDeltaN1(rVariables, rDerivativeData, SlaveGeometry, master_geometry, SlaveNormal, decomp_geom, local_point_decomp, local_point_parent, ConsiderNormalVariation);
                    
                    // Integrate
                    const double integration_weight = AxiSymCoeff * integration_points_slave[point_number].Weight();
                    
                    rAeData.CalculateDeltaAeComponents(rVariables, rDerivativeData, integration_weight);
                }
            }
        }
        
        return rAeData.CalculateDeltaAe(rDerivativeData); 
    }
    
private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * This method is used to compute the directional derivatives of the cell vertex (locally)
     * @param Normal The normal of the slave surface
     * @param DeltaNormal The derivative of the normal vector
     * @param iDoF The DoF computed index
     * @param iBelong The belong (intersection, node, etc..) index
     * @param ConsiderNormalVariation If the normal variation is considered
     * @param SlaveGeometry The geometry of the slave side
     * @param MasterGeometry The geometry of the master side
     * @param Coeff The coefficient considered in proportion
     * @return The local vertex derivative
     */
    static inline array_1d<double, 3> LocalDeltaVertex(
        const array_1d<double, 3>& Normal,
        const bounded_matrix<double, TDim, TDim>& DeltaNormal,
        const unsigned int iDoF,
        const unsigned int iBelong,
        const bool ConsiderNormalVariation,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        double Coeff = 1.0
        )
    {
        // We create the auxiliar array
        array_1d<double, 3> aux_delta_vertex = ZeroVector(3);
        
        // This is the coefficient of the center contribution
        const double auxiliar_coeff = 1.0/static_cast<double>(TNumNodes);
        
        //  We initialize some values
        const array_1d<double, 3>& coords_center = SlaveGeometry.Center().Coordinates();
        const array_1d<double, 3>& coords_node = (iBelong < TNumNodes) ? SlaveGeometry[iBelong].Coordinates() : MasterGeometry[iBelong - TNumNodes].Coordinates();
//         double delta_position;
//         CalculateDeltaPosition(delta_position, SlaveGeometry, MasterGeometry, iBelong, iDoF);
//         Coeff = (std::abs(delta_position) > 0.0) ? Coeff : 0.0;
    
        // The corresponding part to the nodal coordinates
        array_1d<double, 3> aux_der = ZeroVector(3);
        aux_der[iDoF] = 1.0;
        aux_delta_vertex += aux_der;
        
        // The corresponding part to the normal
        const double coordsxdeltanormal = (ConsiderNormalVariation == true) ? inner_prod(coords_node - coords_center, column(DeltaNormal, iDoF)) : 0.0;
        
        const double factor_belong = (iBelong < TNumNodes) ? (1.0 - auxiliar_coeff) : 1.0; 
        const double deltacoordsxnormal =  factor_belong * Normal[iDoF];
        aux_delta_vertex += - Normal * (deltacoordsxnormal + coordsxdeltanormal);
        
        // The corresponding part to delta normal
        const double coordsxnormal = - inner_prod(coords_node - coords_center, Normal);
        if (ConsiderNormalVariation == true) aux_delta_vertex += coordsxnormal * trans(column(DeltaNormal, iDoF));
        
        return Coeff * aux_delta_vertex;
    }
    
    /**
     * This method is used to compute the directional derivatives of the cell vertex (locally)
     * @param DeltaVertexMatrix The whole delta vertex matrix
     * @param Normal The normal of the slave surface
     * @param DeltaNormal The derivative of the normal vector
     * @param iDoF The DoF computed index
     * @param iTriangle The triangle point index
     * @param iBelong The belong (intersection, node, etc..) index
     * @param ConsiderNormalVariation If the normal variation is considered
     * @param SlaveGeometry The geometry of the slave side
     * @param MasterGeometry The geometry of the master side
     * @param Coeff The coefficient considered in proportion
     */
    static inline void LocalDeltaVertex(
        bounded_matrix<double, 3, 3>& DeltaVertexMatrix,
        const array_1d<double, 3>& Normal,
        const bounded_matrix<double, TDim, TDim>& DeltaNormal,
        const unsigned int iDoF,
        const unsigned int iTriangle,
        const unsigned int iBelong,
        const bool ConsiderNormalVariation,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        double Coeff = 1.0
        )
    {
        noalias(row(DeltaVertexMatrix, iTriangle)) += LocalDeltaVertex( Normal,  DeltaNormal, iDoF, iBelong, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, Coeff);
    }
    
    /**
     * This method computes the equivalent indexes to the auxiliar hash
     * @param AuxIndex The auxiliar index to decompose
     * @param iBelongSlaveStart The index of the first/slave segment and first node
     * @param iBelongSlaveEnd The index of the first/slave segment and end node
     * @param iBelongMasterStart The index of the second/master segment and first node
     * @param iBelongMasterEnd The index of the second/master segment and end node
     */
    static inline void ConvertAuxHashIndex(
        const unsigned int AuxIndex,
        unsigned int& iBelongSlaveStart, 
        unsigned int& iBelongSlaveEnd, 
        unsigned int& iBelongMasterStart, 
        unsigned int& iBelongMasterEnd
        )
    {
        unsigned int index_to_decompose = AuxIndex - 2 * TNumNodes;
    
        iBelongMasterEnd = index_to_decompose/10000;
        index_to_decompose = std::fmod(index_to_decompose, 10000);
        iBelongMasterStart = index_to_decompose/1000;
        index_to_decompose = std::fmod(index_to_decompose, 1000);
        iBelongSlaveEnd = index_to_decompose/100;
        index_to_decompose = std::fmod(index_to_decompose, 100);
        iBelongSlaveStart = index_to_decompose/10;
    }
    
    /**
     * This method computes the increment of local coordinates
     * @param rResult The solution obtained
     * @param DeltaPoint The increment of position in the points
     * @param ThisGeometry The geometry considered
     * @param ThisNormal The normal of the geometry
     */
    static inline void DeltaPointLocalCoordinates(
        array_1d<double, 2>& rResult,
        const array_1d<double, 3>& DeltaPoint,
        const MatrixType& rDNDe,
        const GeometryType& ThisGeometry,
        const array_1d<double, 3>& ThisNormal
        )
    {
        bounded_matrix<double, 3, TNumNodes> X;
        for(unsigned int i = 0; i < TNumNodes; ++i)
        {
            X(0, i) = ThisGeometry[i].X();
            X(1, i) = ThisGeometry[i].Y();
            X(2, i) = ThisGeometry[i].Z();
        }
        
        const bounded_matrix<double, 3, 2> DN = prod(X, rDNDe);
        
        double det_j;
        const bounded_matrix<double, 2, 2> J = prod(trans(DN),DN);
        const bounded_matrix<double, 2, 2> invJ = MathUtils<double>::InvertMatrix<2>(J, det_j);
        
        const array_1d<double, 2> res = prod(trans(DN), DeltaPoint);
        noalias(rResult) = prod(invJ, res);

//         bounded_matrix<double, 3, 3> L;
//         for(unsigned int i = 0; i < 3; ++i) 
//         {
//             for(unsigned int j = 0; j < 2; ++j) L(i, j) = DN(i, j);
//             L(i, 2) = ThisNormal[i];
//         }
// 
//         double det_L;
//         const bounded_matrix<double, 3, 3> invL = MathUtils<double>::InvertMatrix<3>(L, det_L);
//         array_1d<double, 3> aux = prod(invL, DeltaPoint);
//         rResult[0] = aux[0];
//         rResult[1] = aux[1];
    }
    
    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@} 
};// class DerivativesUtilities

}
#endif /* KRATOS_DERIVATIVES_UTILITIES defined */
 
