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
    
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
class DerivativesUtilities
{
public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef Vector                                                                                     VectorType;

    typedef Matrix                                                                                     MatrixType;

    typedef std::size_t                                                                                 IndexType;

    typedef Geometry<NodeType>                                                                       GeometryType;

    typedef Geometry<NodeType>::PointsArrayType                                                    NodesArrayType;

    typedef Properties                                                                             PropertiesType;
    
    typedef typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N, typename std::conditional<TNumNodes == 3, PointBelongsTriangle3D3N, PointBelongsQuadrilateral3D4N>::type>::type BelongType;
    
    typedef PointBelong<TNumNodes>                                                                PointBelongType;
    
    typedef Geometry<PointBelongType>                                                     GeometryPointBelongType;
    
    typedef array_1d<PointBelongType,TDim>                                                     ConditionArrayType;
    
    typedef typename std::vector<ConditionArrayType>                                       ConditionArrayListType;
    
    typedef Line2D2<PointType>                                                                           LineType;
    
    typedef Triangle3D3<PointType>                                                                   TriangleType;
    
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type                 DecompositionType;
    
    typedef typename std::conditional<TFrictional == true, DerivativeDataFrictional<TDim, TNumNodes, TNormalVariation>, DerivativeData<TDim, TNumNodes, TNormalVariation> >::type DerivativeDataType;
    
    typedef MortarKinematicVariablesWithDerivatives<TDim, TNumNodes>                             GeneralVariables;
    
    typedef DualLagrangeMultiplierOperatorsWithDerivatives<TDim, TNumNodes, TFrictional, TNormalVariation> AeData;
    
    typedef MortarOperatorWithDerivatives<TDim, TNumNodes, TFrictional, TNormalVariation> MortarConditionMatrices;
    
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
        const GeneralVariables& rVariables,
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
                    const auto& local_delta_vertex = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof];
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
     * @param Jacobian The jacobian on the GP
     * @param DNDe The local gradient
     * @return The matrix containing the delta normals
     * NOTE: Not the mean, look in the contact utilities 
     */
    static inline array_1d<array_1d<double, 3>, TDim * TNumNodes> GPDeltaNormal(
        const Matrix& Jacobian,
        const Matrix& DNDe
        )
    {                
        // Tangent directions
        array_1d<double,3> j0(3, 0.0), j1(3, 0.0);
        
        // Using the Jacobian tangent directions
        if (TDim == 2)
        {
            j1[2] = 1.0;
            for (unsigned int i_dim = 0; i_dim < 2; ++i_dim)
            {
                j0[i_dim]  = Jacobian(i_dim, 0);
            } 
        }
        else
        {
            for (unsigned int i_dim = 0; i_dim < 3; ++i_dim)
            {
                j0[i_dim] = Jacobian(i_dim, 0);
                j1[i_dim] = Jacobian(i_dim, 1);
            } 
        }
        
        array_1d<double, 3> normal;;
        MathUtils<double>::CrossProduct(normal, j0, j1);
        const double area_normal_norm = norm_2(normal);
    #ifdef KRATOS_DEBUG
        KRATOS_ERROR_IF(area_normal_norm < std::numeric_limits<double>::epsilon()) << "ZERO NORMAL: " << area_normal_norm << std::endl;
    #endif
        const array_1d<double, 3> unit_normal = normal/area_normal_norm;
        
        array_1d<array_1d<double, 3>, TDim * TNumNodes> delta_normal;
        for ( unsigned int i_node = 0; i_node < TNumNodes; ++i_node )
        {
            for(unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
            {
                const unsigned int i_dof = i_node * TDim + i_dim;

                array_1d<double,3> delta_j0(3, 0.0), delta_j1(3, 0.0);

                delta_j0[i_dim] += DNDe(i_node, 0);
                if (TDim == 3) delta_j1[i_dim] += DNDe(i_node, 1);
                
                array_1d<double, 3> aux_delta_normal0, aux_delta_normal1;
                MathUtils<double>::CrossProduct(aux_delta_normal0, j0, delta_j1);
                MathUtils<double>::CrossProduct(aux_delta_normal1, delta_j0, j1);
                const array_1d<double, 3> aux_delta_normal = aux_delta_normal0 + aux_delta_normal1;
                const double delta_norm = inner_prod(aux_delta_normal, normal);
                
                delta_normal[i_dof] = (aux_delta_normal + unit_normal * delta_norm)/area_normal_norm;
            }
        }

        return delta_normal; 
    }
    
    /**
     * It computes the delta normal of the center of the geoemtry
     * @param rThisGeometry The geometry where the delta normal is computed
     */
    static inline array_1d<array_1d<double, 3>, TDim * TNumNodes> DeltaNormalCenter(GeometryType& rThisGeometry)
    {
        // We compute the gradient and jacobian
        GeometryType::CoordinatesArrayType point_local;  
        rThisGeometry.PointLocalCoordinates( point_local, (rThisGeometry.Center()).Coordinates( ) ) ;
        Matrix jacobian;
        jacobian = rThisGeometry.Jacobian( jacobian, point_local);
        Matrix gradient;
        rThisGeometry.ShapeFunctionsLocalGradients( gradient, point_local );
        
        // We compute the previous normal (TODO: Think about to save it instead)
        const array_1d<double, 3>& previous_normal = PreviousNormalGeometry(rThisGeometry, point_local);
        
        // Now we compute the normal + DeltaNormal
        array_1d<double, 3> aux_delta_normal0(3, 0.0);
        const auto& aux_delta_normal_0 = GPDeltaNormal(jacobian, gradient);
        
        for ( unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
        {
            const array_1d<double, 3> delta_disp = rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT, 1);
            for ( unsigned int i_dof = 0; i_dof < TDim; ++i_dof)
            {
                const auto& aux_delta_normal = aux_delta_normal_0[i_node * TDim + i_dof];
                aux_delta_normal0 += delta_disp[i_dof] * aux_delta_normal;
            }
        }
        
        array_1d<double, 3> calculated_normal_geometry = aux_delta_normal0 + previous_normal;
        calculated_normal_geometry /= norm_2(calculated_normal_geometry);
        
        // We compute the diff matrix to compute the auxiliar matrix later
        const array_1d<double, 3> diff_vector = calculated_normal_geometry - previous_normal;  
        
        // Computing auxiliar matrix
        const bounded_matrix<double, 3, 3>& renormalizer_matrix = ComputeRenormalizerMatrix(diff_vector, aux_delta_normal0);
        array_1d<array_1d<double, 3>, TDim * TNumNodes> normalized_delta_normal_0;
        for ( unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
        {
            for ( unsigned int i_dof = 0; i_dof < TDim; ++i_dof)
            {
                const auto& aux_delta_normal = aux_delta_normal_0[i_node * TDim + i_dof];
                normalized_delta_normal_0[i_node * TDim + i_dof] = prod(renormalizer_matrix, aux_delta_normal);
            }
        }
        
        return normalized_delta_normal_0;
    }

    /**
     * Calculates the increment of the normal and in the master condition
     * @param rDeltaNormal The derivative of the normal
     * @param rThisGeometry The geometry of the master side
     */
    
    static inline void CalculateDeltaNormal(
        array_1d<bounded_matrix<double, TNumNodes, TDim>, TNumNodes * TDim>& rDeltaNormal,
        GeometryType& rThisGeometry
        )
    {
        bounded_matrix<double, TNumNodes, TDim> aux_normal_geometry = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(rThisGeometry,  NORMAL, 1);
        bounded_matrix<double, TNumNodes, TDim> aux_delta_normal_geometry = ZeroMatrix(TNumNodes, TDim);
        
        for ( unsigned int i_geometry = 0; i_geometry < TNumNodes; ++i_geometry )
        {
            // We compute the gradient and jacobian
            GeometryType::CoordinatesArrayType point_local;  
            rThisGeometry.PointLocalCoordinates( point_local, rThisGeometry[i_geometry].Coordinates( ) ) ;
            Matrix jacobian;
            jacobian = rThisGeometry.Jacobian( jacobian, point_local);
            Matrix gradient;
            rThisGeometry.ShapeFunctionsLocalGradients( gradient, point_local );
            
            // We compute the delta normal of the node
            const auto& delta_normal_node = GPDeltaNormal(jacobian, gradient);
            for ( unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            {
                const array_1d<double, 3> delta_disp = rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT, 1);
                for ( unsigned int i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    const auto& aux_delta_normal = subrange(delta_normal_node[i_node * TDim + i_dof], 0, TDim);
                    row(aux_delta_normal_geometry, i_geometry) += delta_disp[i_dof] * aux_delta_normal;
                }
            }
        }
        
        bounded_matrix<double, TNumNodes, TDim> calculated_normal_geometry = aux_delta_normal_geometry + aux_normal_geometry;
        for ( unsigned int i_geometry = 0; i_geometry < TNumNodes; ++i_geometry ) row(calculated_normal_geometry, i_geometry) /= norm_2(row(calculated_normal_geometry, i_geometry));
        
        // We compute the diff matrix to compute the auxiliar matrix later
        const bounded_matrix<double, TNumNodes, TDim> diff_matrix = calculated_normal_geometry - aux_normal_geometry;  
        
        // We iterate over the nodes of the geometry
        for ( unsigned int i_geometry = 0; i_geometry < TNumNodes; ++i_geometry )
        {
            // Computing auxiliar matrix
            const bounded_matrix<double, TDim, TDim> renormalizer_matrix = ComputeRenormalizerMatrix(diff_matrix, aux_delta_normal_geometry, i_geometry);
            
            // We compute the gradient and jacobian
            GeometryType::CoordinatesArrayType point_local;  
            rThisGeometry.PointLocalCoordinates( point_local, rThisGeometry[i_geometry].Coordinates( ) ) ;
            Matrix jacobian;
            jacobian = rThisGeometry.Jacobian( jacobian, point_local);
            Matrix gradient;
            rThisGeometry.ShapeFunctionsLocalGradients( gradient, point_local );
            
            // We compute the delta normal of the node
            const auto& delta_normal_node = GPDeltaNormal(jacobian, gradient);
            for ( unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            {
                for ( unsigned int i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    const auto& aux_delta_normal = subrange(delta_normal_node[i_node * TDim + i_dof], 0, TDim);
                    row(rDeltaNormal[i_node * TDim + i_dof], i_geometry) = prod(renormalizer_matrix, aux_delta_normal);
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
        const GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        const array_1d<BelongType, TDim>& TheseBelongs,
        const NormalDerivativesComputation ConsiderNormalVariation,
        GeometryType& SlaveGeometry,
        GeometryType& MasterGeometry,
        const array_1d<double, 3>& Normal
        )
    {
        // The Normal and delta Normal in the center of the element
        const array_1d<array_1d<double, 3>, TDim * TNumNodes>& all_delta_normal = DeltaNormalCenter(SlaveGeometry);
        array_1d<double, 3> delta_normal;
        
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
                double distance;
                const array_1d<double, 3>& xs1 = MortarUtilities::FastProject(slave_center, SlaveGeometry[belong_index_slave_start], Normal, distance).Coordinates(); // Start coordinates of the first segment
                const array_1d<double, 3>& xe1 = MortarUtilities::FastProject(slave_center, SlaveGeometry[belong_index_slave_end], Normal, distance).Coordinates(); // End coordinates of the first segment
                const array_1d<double, 3>& xs2 = MortarUtilities::FastProject(slave_center, MasterGeometry[belong_index_master_start], Normal, distance).Coordinates(); // Start coordinates of the second segment
                const array_1d<double, 3>& xe2 = MortarUtilities::FastProject(slave_center, MasterGeometry[belong_index_master_end], Normal, distance).Coordinates(); // End coordinates of the second segment
                
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
                    
                    for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                    {                    
                        // We get the delta normal
                        if (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION && belong_index < TNumNodes) delta_normal = all_delta_normal[belong_index * TDim + i_dof] * (1.0/aux_nodes_coeff);
                        else delta_normal = ZeroVector(3);
                    
                        auto& local_delta_vertex = rDerivativeData.DeltaCellVertex[belong_index * TDim + i_dof];
                        
                        // Special cases (slave nodes)
                        if (i_belong == 0) // First node of the slave
                        {
                            const double coeff = 1.0 + num/denom;
                            noalias(row(local_delta_vertex, i_triangle)) += LocalDeltaVertex( Normal,  delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, coeff);
                        }
                        else if (i_belong == 1) // Second node of the slave
                        {
                            const double coeff = - num/denom;
                            noalias(row(local_delta_vertex, i_triangle)) += LocalDeltaVertex( Normal,  delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, coeff);
                        }
                        
                        // We define some auxiliar coefficients
                        const double coeff1 = - 1.0/denom;
                        const double coeff2 = num/std::pow(denom, 2);
                        
                        // We add the part corresponding purely to delta normal
                        if (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION)
                        {
                            noalias(row(local_delta_vertex, i_triangle)) += diff3 * coeff1 * inner_prod(aux_num,  delta_normal); 
                            noalias(row(local_delta_vertex, i_triangle)) += diff3 * coeff2 * inner_prod(aux_denom, delta_normal); 
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
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    // We get the delta normal
                    if (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION && belong_index < TNumNodes) delta_normal = all_delta_normal[belong_index * TDim + i_dof] * (1.0/aux_nodes_coeff);
                    else delta_normal = ZeroVector(3);
                    
                    auto& local_delta_vertex = rDerivativeData.DeltaCellVertex[belong_index * TDim + i_dof];
                    noalias(row(local_delta_vertex, i_triangle)) += LocalDeltaVertex( Normal,  delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry);
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
        const NormalDerivativesComputation ConsiderNormalVariation = NO_DERIVATIVES_COMPUTATION
        )
    {
        /* Shape functions */
        const VectorType& N1 = rVariables.NSlave;

        /* Local gradients */
        const MatrixType& DNDe1 = rVariables.DNDeSlave;
        
        // The Normal and delta Normal in the center of the element
        const array_1d<array_1d<double, 3>, TDim * TNumNodes>& all_delta_normal = DeltaNormalCenter(SlaveGeometry);
        
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
                    // We get the delta normal
                    const array_1d<double, 3>& delta_normal = (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION && i_node < TNumNodes) ? all_delta_normal[i_node * TDim + i_dof] : ZeroVector(3);
                    
                    // We compute the residuals
                    array_1d<double, 3> aux_RHS1 = ZeroVector(3);
                    
                    // The vertex cell contribution
                    const auto& local_delta_cell = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof]; 
                    for(std::size_t i_belong = 0; i_belong < 3; ++i_belong)
                    {
                        noalias(aux_RHS1) += N_decomp[i_belong] * row(local_delta_cell, i_belong);
                    }
                    
                    // Local contribution
                    const array_1d<double, 3>& aux_delta_node = LocalDeltaVertex( SlaveNormal, delta_normal, i_dof, i_node, ConsiderNormalVariation, SlaveGeometry, MasterGeometry );
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
        const NormalDerivativesComputation ConsiderNormalVariation = NO_DERIVATIVES_COMPUTATION,
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
        const array_1d<array_1d<double, 3>, TDim * TNumNodes>& all_delta_normal = DeltaNormalCenter(SlaveGeometry);
        
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
                    // We get the delta normal
                    const array_1d<double, 3>& delta_normal = (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION && i_node < TNumNodes) ? all_delta_normal[i_node * TDim + i_dof] : ZeroVector(3);
                    
                    // We compute the residuals
                    array_1d<double, 3> aux_RHS1 = ZeroVector(3);
                    
                    // The vertex cell contribution
                    const auto& local_delta_cell = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof]; 
                    for(std::size_t i_belong = 0; i_belong < 3; ++i_belong)
                    {
                        noalias(aux_RHS1) += N_decomp[i_belong] * row(local_delta_cell, i_belong);
                    }
                    array_1d<double, 3> aux_RHS2 = aux_RHS1;
                    
                    // Local contribution
                    const array_1d<double, 3>& aux_delta_node = LocalDeltaVertex( SlaveNormal, delta_normal, i_dof, i_node, ConsiderNormalVariation, SlaveGeometry, MasterGeometry );
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
                    if (DualLM == true) 
                    {
                        noalias(delta_phi) = prod(rDerivativeData.DeltaAe[i_node * TDim + i_dof], N1);
                        noalias(delta_phi) += prod(rDerivativeData.Ae, delta_n1);
                    }
                    else 
                    {
                        noalias(delta_phi) = delta_n1;
                    }
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
     * @param MasterGeometry The master side geometry
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
        GeometryType& MasterGeometry,
        DerivativeDataType& rDerivativeData,
        GeneralVariables& rVariables,
        const NormalDerivativesComputation ConsiderNormalVariation,
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
        
        for (unsigned int i_geom = 0; i_geom < ConditionsPointsSlave.size(); ++i_geom)
        {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (unsigned int i_node = 0; i_node < TDim; ++i_node)
            {
                PointType global_point;
                SlaveGeometry.GlobalCoordinates(global_point, ConditionsPointsSlave[i_geom][i_node]);
                points_array[i_node] = PointType::Pointer( new PointType(global_point) );
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
                    MortarUtilities::FastProjectDirection( MasterGeometry, slave_gp_global, projected_gp_global, SlaveNormal, -gp_normal ); // The opposite direction
                    
                    GeometryType::CoordinatesArrayType projected_gp_local;
                    MasterGeometry.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                    MasterGeometry.ShapeFunctionsValues( rVariables.NMaster,    projected_gp_local );         
                    MasterGeometry.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );
                    rVariables.jMaster = MasterGeometry.Jacobian( rVariables.jMaster, projected_gp_local);
                    
                    // Update the derivative of the integration vertex (just in 3D)
                    if (TDim == 3) CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, SlaveNormal);
                                    
                    // Update the derivative of DetJ
                    CalculateDeltaDetjSlave(decomp_geom, rVariables, rDerivativeData); 
                    
                    // Update the derivatives of the shape functions and the gap
                    CalculateDeltaN1(rVariables, rDerivativeData, SlaveGeometry, MasterGeometry, SlaveNormal, decomp_geom, local_point_decomp, local_point_parent, ConsiderNormalVariation);
                    
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
        const array_1d<double, 3>& DeltaNormal,
        const unsigned int iDoF,
        const unsigned int iBelong,
        const NormalDerivativesComputation ConsiderNormalVariation,
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
    
        // The corresponding part to the nodal coordinates
        array_1d<double, 3> aux_der = ZeroVector(3);
        aux_der[iDoF] = 1.0;
        aux_delta_vertex += aux_der;
        
        // The corresponding part to the normal
        const double coordsxdeltanormal = (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION) ? inner_prod(coords_node - coords_center, DeltaNormal) : 0.0;
        
        const double factor_belong = (iBelong < TNumNodes) ? (1.0 - auxiliar_coeff) : 1.0; 
        const double deltacoordsxnormal =  factor_belong * Normal[iDoF];
        aux_delta_vertex += - Normal * (deltacoordsxnormal + coordsxdeltanormal);
        
        // The corresponding part to delta normal
        const double coordsxnormal = - inner_prod(coords_node - coords_center, Normal);
        if (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION) aux_delta_vertex += coordsxnormal * DeltaNormal;
        
        return Coeff * aux_delta_vertex;
    }
    
    /**
     * This method computes the auxiliar matrix used to keep unitary the normal
     * @param DiffVector The auxiliar vector of difference of two normal vectors
     * @param DeltaNormal The vector containing the delta normal
     * @return The auxiliar matrix computed
     */
    static inline bounded_matrix<double, 3, 3> ComputeRenormalizerMatrix(
        const array_1d<double, 3>& DiffVector, 
        const array_1d<double, 3>& DeltaNormal
        ) 
    {
        for (unsigned int itry = 0; itry < 3; ++itry)
        {
            if (DeltaNormal[itry] > std::numeric_limits<double>::epsilon())
            {
                bounded_matrix<double, 3, 3> aux_matrix;
                
                const unsigned int aux_index_1 = itry == 2 ? 0 : itry + 1;
                const unsigned int aux_index_2 = itry == 2 ? 1 : (itry == 1 ? 0 : 2);
                
                const double diff = DeltaNormal[aux_index_1] + DeltaNormal[aux_index_2];
                const double coeff = DeltaNormal[itry];
                
                aux_matrix(0, aux_index_1) = 1.0; 
                aux_matrix(0, aux_index_2) = 1.0;
                aux_matrix(1, aux_index_1) = 1.0; 
                aux_matrix(1, aux_index_2) = 1.0;
                aux_matrix(2, aux_index_1) = 1.0; 
                aux_matrix(2, aux_index_2) = 1.0;

                aux_matrix(0, itry) = (DiffVector[0] - diff)/coeff;
                aux_matrix(1, itry) = (DiffVector[1] - diff)/coeff;
                aux_matrix(2, itry) = (DiffVector[2] - diff)/coeff;
                
                return aux_matrix;
            }
        }
        
        return IdentityMatrix(3, 3);
    }
    
    /**
     * This method computes the auxiliar matrix used to keep unitary the normal
     * @param DiffMatrix The auxiliar matrix of difference of two normal matrices
     * @param DeltaNormal The matrix containing the delta normal
     * @param iGeometry The index of the node of the geoemtry computed
     * @return The auxiliar matrix computed
     */
    static inline bounded_matrix<double, 3, 3> ComputeRenormalizerMatrix(
        const bounded_matrix<double, TNumNodes, TDim>& DiffMatrix, 
        const bounded_matrix<double, TNumNodes, TDim>& DeltaNormal,
        const unsigned int iGeometry
        )
    {        
        for (unsigned int itry = 0; itry < 3; ++itry)
        {            
            if (DeltaNormal(iGeometry, itry) > std::numeric_limits<double>::epsilon())
            {
                bounded_matrix<double, 3, 3> aux_matrix;
                
                const unsigned int aux_index_1 = itry == 2 ? 0 : itry + 1;
                const unsigned int aux_index_2 = itry == 2 ? 1 : (itry == 1 ? 0 : 2);
                
                const double diff = DeltaNormal(iGeometry, aux_index_1) + DeltaNormal(iGeometry, aux_index_2);
                const double coeff = DeltaNormal(iGeometry, itry);
                
                aux_matrix(0, aux_index_1) = 1.0; 
                aux_matrix(0, aux_index_2) = 1.0;
                aux_matrix(1, aux_index_1) = 1.0; 
                aux_matrix(1, aux_index_2) = 1.0;
                aux_matrix(2, aux_index_1) = 1.0; 
                aux_matrix(2, aux_index_2) = 1.0;

                aux_matrix(0, itry) = (DiffMatrix(iGeometry, 0) - diff)/coeff;
                aux_matrix(1, itry) = (DiffMatrix(iGeometry, 1) - diff)/coeff;
                aux_matrix(2, itry) = (DiffMatrix(iGeometry, 2) - diff)/coeff;
                
                return aux_matrix;
            }
        }
        
        return IdentityMatrix(3, 3);
    }
    
    /**
     * This method computes the normal in the previous configuration
     * @param rThisGeometry The geometry where compute
     * @param PointLocal The local coordinates of the point
     * @return The normal in the previous configuration
     */
    static inline array_1d<double, 3> PreviousNormalGeometry(
        const GeometryType& rThisGeometry,
        const GeometryType::CoordinatesArrayType& PointLocal
        )
    {
        // We compute the previous normal in the geometry 
        Matrix previous_jacobian, delta_position;
        delta_position = CalculateDeltaPosition(delta_position, rThisGeometry);
        previous_jacobian = rThisGeometry.Jacobian( previous_jacobian, PointLocal, delta_position);
        
        // We define the normal and tangents
        array_1d<double,3> tangent_xi(3, 0.0), tangent_eta(3, 0.0);
        
        // Using the Jacobian tangent directions
        if (TDim == 2)
        {
            tangent_eta[2] = 1.0;
            for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
            {
                tangent_xi[i_dim]  = previous_jacobian(i_dim, 0);
            } 
        }
        else
        {
            for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
            {
                tangent_xi[i_dim]  = previous_jacobian(i_dim, 0);
                tangent_eta[i_dim] = previous_jacobian(i_dim, 1);
            } 
        }

        array_1d<double, 3> previous_normal;
        MathUtils<double>::CrossProduct(previous_normal, tangent_xi, tangent_eta);
        const double norm_normal = norm_2(previous_normal);
        previous_normal /= norm_normal;
        KRATOS_ERROR_IF(norm_normal < std::numeric_limits<double>::epsilon()) << "ERROR: The normal norm is zero or almost zero. Norm. normal: " << norm_normal << std::endl;
        
        return previous_normal;
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
        // Tolerance
        const double tolerance = std::numeric_limits<double>::epsilon();
        
        bounded_matrix<double, 3, TNumNodes> X;
        for(unsigned int i = 0; i < TNumNodes; ++i)
        {
            X(0, i) = ThisGeometry[i].X();
            X(1, i) = ThisGeometry[i].Y();
            X(2, i) = ThisGeometry[i].Z();
        }
        
        const bounded_matrix<double, 3, 2> DN = prod(X, rDNDe);
        
        const bounded_matrix<double, 2, 2> J = prod(trans(DN),DN);
        double det_j = MathUtils<double>::DetMat<bounded_matrix<double, 2, 2>>(J);
        const bounded_matrix<double, 2, 2> invJ = (std::abs(det_j) < tolerance) ? ZeroMatrix(2,2) : MathUtils<double>::InvertMatrix<2>(J, det_j);
        
    #ifdef KRATOS_DEBUG
        if (std::abs(det_j) < tolerance) std::cout << "WARNING: CANNOT INVERT JACOBIAN TO COMPUTE DELTA COORDINATES" << std::endl;
    #endif
        
        const array_1d<double, 2> res = prod(trans(DN), DeltaPoint);
        noalias(rResult) = prod(invJ, res);

//         bounded_matrix<double, 3, 3> L;
//         for(unsigned int i = 0; i < 3; ++i) 
//         {
//             for(unsigned int j = 0; j < 2; ++j) L(i, j) = DN(i, j);
//             L(i, 2) = ThisNormal[i];
//         }
// 
//         double det_L = MathUtils<double>::DetMat<bounded_matrix<double, 3, 3>>(L);
//         const bounded_matrix<double, 3, 3> invL = (std::abs(det_L) < tolerance) ? ZeroMatrix(3,3) : MathUtils<double>::InvertMatrix<3>(L, det_L);
//         array_1d<double, 3> aux = prod(invL, DeltaPoint);
//         rResult[0] = aux[0];
//         rResult[1] = aux[1];
//         #ifdef KRATOS_DEBUG
//             if (std::abs(det_L) < tolerance) std:::cout << "WARNING: CANNOT INVERT JACOBIAN TO COMPUTE DELTA COORDINATES" << std::endl;
//         #endif
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
 
