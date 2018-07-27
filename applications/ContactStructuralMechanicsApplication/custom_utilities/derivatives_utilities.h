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

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"

/* Includes */
#include "includes/model_part.h"
#include "includes/mortar_classes.h"

/* Utilities */
#include "utilities/geometrical_projection_utilities.h"
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
/**
 * @class DerivativesUtilities 
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This utilities are used in order to compute the directional derivatives during mortar contact
 * @details The derivatives take the same argument templates than the contact conditions
 * @author Vicente Mataix Ferrandiz
 * @author Gabriel Valdes Alonzo 
 */
template< std::size_t TDim, std::size_t TNumNodes, bool TFrictional, bool TNormalVariation>
class DerivativesUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// The vector considered
    typedef Vector                                                                                     VectorType;

    /// The dense matrix considered
    typedef Matrix                                                                                     MatrixType;

    /// The index type
    typedef std::size_t                                                                                 IndexType;

    /// The geometry of nodes
    typedef Geometry<NodeType>                                                                       GeometryType;

    /// The array of nodes contained in a geometry
    typedef Geometry<NodeType>::PointsArrayType                                                    NodesArrayType;

    /// The Properties type
    typedef Properties                                                                             PropertiesType;

    /// The belong type (for derivatives definition)
    typedef typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N, typename std::conditional<TNumNodes == 3, PointBelongsTriangle3D3N, PointBelongsQuadrilateral3D4N>::type>::type BelongType;

    /// The points used for derivatives definition
    typedef PointBelong<TNumNodes>                                                                PointBelongType;

    /// A geometry defined by the point belongs (the points used for derivatives definition)
    typedef Geometry<PointBelongType>                                                     GeometryPointBelongType;

    /// An array of belong point to define the geometries of belong points
    typedef array_1d<PointBelongType,TDim>                                                     ConditionArrayType;

    /// The definition of an array pf conditions
    typedef typename std::vector<ConditionArrayType>                                       ConditionArrayListType;

    /// The line definition
    typedef Line2D2<PointType>                                                                           LineType;

    /// The triangle definition
    typedef Triangle3D3<PointType>                                                                   TriangleType;

    /// The geometry for decomposition (line in 2D and triangle for 3D)
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type                 DecompositionType;

    /// The derivative data type
    typedef typename std::conditional<TFrictional == true, DerivativeDataFrictional<TDim, TNumNodes, TNormalVariation>, DerivativeData<TDim, TNumNodes, TNormalVariation> >::type DerivativeDataType;

    /// The kinematic variables
    typedef MortarKinematicVariablesWithDerivatives<TDim, TNumNodes>                             GeneralVariables;

    /// The dual LM operators
    typedef DualLagrangeMultiplierOperatorsWithDerivatives<TDim, TNumNodes, TFrictional, TNormalVariation> AeData;

    /// The mortar operators
    typedef MortarOperatorWithDerivatives<TDim, TNumNodes, TFrictional, TNormalVariation> MortarConditionMatrices;

    /// Pointer definition of DerivativesUtilities
    KRATOS_CLASS_POINTER_DEFINITION( DerivativesUtilities );
    
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
     * @brief This method is used to compute the directional derivatives of the jacobian determinant
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
        if (TDim == 2) {
            // Fill up the elements corresponding to the slave DOFs - the rest remains zero
            for ( IndexType i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim ) {
                rDerivativeData.DeltaDetjSlave[i    ] = rVariables.jSlave( 0, 0 ) * rVariables.DNDeSlave( i_slave, 0) / rVariables.DetjSlave;
                rDerivativeData.DeltaDetjSlave[i + 1] = rVariables.jSlave( 1, 0 ) * rVariables.DNDeSlave( i_slave, 0) / rVariables.DetjSlave;
            }
        } else {
            const array_1d<double, 3>& x1cell = DecompGeom[0].Coordinates();
            const array_1d<double, 3>& x2cell = DecompGeom[1].Coordinates();
            const array_1d<double, 3>& x3cell = DecompGeom[2].Coordinates();

            const array_1d<double, 3> x21cell = x2cell - x1cell;
            const array_1d<double, 3> x31cell = x3cell - x1cell;

            array_1d<double, 3> aux_cross_product;
            MathUtils<double>::CrossProduct(aux_cross_product, x21cell, x31cell);
            aux_cross_product /= rVariables.DetjSlave;
//             aux_cross_product /= norm_2(aux_cross_product);

            for ( IndexType i_node = 0; i_node < 2 * TNumNodes; ++i_node ) {
                for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
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
     * @brief This method is used to compute the local increment of the normal
     * @param Jacobian The jacobian on the GP
     * @param DNDe The local gradient
     * @return The matrix containing the delta normals
     * @note Not the mean, look in the contact utilities
     */
    static inline array_1d<array_1d<double, 3>, TDim * TNumNodes> GPDeltaNormal(
        const Matrix& Jacobian,
        const Matrix& DNDe
        )
    {
        // Tangent directions
        array_1d<double,3> j0 = ZeroVector(3), j1 = ZeroVector(3);

        // Using the Jacobian tangent directions
        if (TDim == 2) {
            j1[2] = 1.0;
            for (IndexType i_dim = 0; i_dim < 2; ++i_dim)
                j0[i_dim]  = Jacobian(i_dim, 0);
        } else {
            for (IndexType i_dim = 0; i_dim < 3; ++i_dim) {
                j0[i_dim] = Jacobian(i_dim, 0);
                j1[i_dim] = Jacobian(i_dim, 1);
            }
        }

        array_1d<double, 3> normal;;
        MathUtils<double>::CrossProduct(normal, j0, j1);
        const double area_normal_norm = norm_2(normal);

        KRATOS_DEBUG_ERROR_IF(area_normal_norm < std::numeric_limits<double>::epsilon()) << "ZERO NORMAL: " << area_normal_norm << std::endl;

        const array_1d<double, 3> unit_normal = normal/area_normal_norm;

        array_1d<array_1d<double, 3>, TDim * TNumNodes> delta_normal;
        array_1d<double,3> delta_j0, delta_j1;
        const array_1d<double, 3> zero_array(3, 0.0);
        for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node ) {
            for(IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                const IndexType i_dof = i_node * TDim + i_dim;

                delta_j0 = zero_array;
                delta_j1 = zero_array;

                delta_j0[i_dim] += DNDe(i_node, 0);
                if (TDim == 3) 
                    delta_j1[i_dim] += DNDe(i_node, 1);

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
     * @brief It computes the delta normal of the center of the geometry
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
        const array_1d<double, 3> previous_normal = PreviousNormalGeometry(rThisGeometry, point_local);

        // Now we compute the normal + DeltaNormal
        array_1d<double, 3> aux_delta_normal0 = ZeroVector(3);
        const array_1d<array_1d<double, 3>, TDim * TNumNodes> aux_delta_normal_0 = GPDeltaNormal(jacobian, gradient);

        for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3> delta_disp = rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT, 1);
            for ( IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                const array_1d<double, 3>& aux_delta_normal = aux_delta_normal_0[i_node * TDim + i_dof];
                aux_delta_normal0 += delta_disp[i_dof] * aux_delta_normal;
            }
        }

        array_1d<double, 3> calculated_normal_geometry = aux_delta_normal0 + previous_normal;
        calculated_normal_geometry /= norm_2(calculated_normal_geometry);

        // We compute the diff matrix to compute the auxiliar matrix later
        const array_1d<double, 3> diff_vector = calculated_normal_geometry - previous_normal;

        // Computing auxiliar matrix
        const BoundedMatrix<double, 3, 3> renormalizer_matrix = ComputeRenormalizerMatrix(diff_vector, aux_delta_normal0);
        array_1d<array_1d<double, 3>, TDim * TNumNodes> normalized_delta_normal_0;
        for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            for ( IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                const array_1d<double, 3>& aux_delta_normal = aux_delta_normal_0[i_node * TDim + i_dof];
                normalized_delta_normal_0[i_node * TDim + i_dof] = prod(renormalizer_matrix, aux_delta_normal);
            }
        }

        return normalized_delta_normal_0;
    }

    /**
     * @brief Calculates the increment of the normal and in the master condition
     * @param rDeltaNormal The derivative of the normal
     * @param rThisGeometry The geometry of the master side
     */

    static inline void CalculateDeltaNormal(
        array_1d<BoundedMatrix<double, TNumNodes, TDim>, TNumNodes * TDim>& rDeltaNormal,
        GeometryType& rThisGeometry
        )
    {
        BoundedMatrix<double, TNumNodes, TDim> aux_normal_geometry = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(rThisGeometry,  NORMAL, 1);
        BoundedMatrix<double, TNumNodes, TDim> aux_delta_normal_geometry = ZeroMatrix(TNumNodes, TDim);

        for ( IndexType i_geometry = 0; i_geometry < TNumNodes; ++i_geometry ) {
            // We compute the gradient and jacobian
            GeometryType::CoordinatesArrayType point_local;
            rThisGeometry.PointLocalCoordinates( point_local, rThisGeometry[i_geometry].Coordinates( ) ) ;
            Matrix jacobian;
            jacobian = rThisGeometry.Jacobian( jacobian, point_local);
            Matrix gradient;
            rThisGeometry.ShapeFunctionsLocalGradients( gradient, point_local );

            // We compute the delta normal of the node
            const array_1d<array_1d<double, 3>, TDim * TNumNodes> delta_normal_node = GPDeltaNormal(jacobian, gradient);
            for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                const array_1d<double, 3> delta_disp = rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT, 1);
                for ( IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                    const array_1d<double, TDim>& aux_delta_normal = subrange(delta_normal_node[i_node * TDim + i_dof], 0, TDim);
                    row(aux_delta_normal_geometry, i_geometry) += delta_disp[i_dof] * aux_delta_normal;
                }
            }
        }

        BoundedMatrix<double, TNumNodes, TDim> calculated_normal_geometry = aux_delta_normal_geometry + aux_normal_geometry;
        for ( IndexType i_geometry = 0; i_geometry < TNumNodes; ++i_geometry )
            row(calculated_normal_geometry, i_geometry) /= norm_2(row(calculated_normal_geometry, i_geometry));

        // We compute the diff matrix to compute the auxiliar matrix later
        const BoundedMatrix<double, TNumNodes, TDim> diff_matrix = calculated_normal_geometry - aux_normal_geometry;

        // We iterate over the nodes of the geometry
        for ( IndexType i_geometry = 0; i_geometry < TNumNodes; ++i_geometry ) {
            // Computing auxiliar matrix
            const BoundedMatrix<double, TDim, TDim> renormalizer_matrix = (TDim == 3) ? ComputeRenormalizerMatrix(diff_matrix, aux_delta_normal_geometry, i_geometry) : IdentityMatrix(2, 2);

            // We compute the gradient and jacobian
            GeometryType::CoordinatesArrayType point_local;
            rThisGeometry.PointLocalCoordinates( point_local, rThisGeometry[i_geometry].Coordinates( ) ) ;
            Matrix jacobian;
            jacobian = rThisGeometry.Jacobian( jacobian, point_local);
            Matrix gradient;
            rThisGeometry.ShapeFunctionsLocalGradients( gradient, point_local );

            // Auxilary terms
            array_1d<double, TDim> aux_delta_normal;
            
            // We compute the delta normal of the node
            const array_1d<array_1d<double, 3>, TDim * TNumNodes> delta_normal_node = GPDeltaNormal(jacobian, gradient);
            for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                for ( IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                    array_1d<double, TDim> aux_delta_normal = subrange(delta_normal_node[i_node * TDim + i_dof], 0, TDim);
                    row(rDeltaNormal[i_node * TDim + i_dof], i_geometry) = prod(renormalizer_matrix, aux_delta_normal);
                }
            }
        }
    }

    /**
     * @brief This method is used to compute the directional derivatives of the cell vertex
     * @param rVariables The kinematic variables
     * @param rDerivativeData The derivatives container
     * @param TheseBelongs The belongs list used in the derivatives
     * @param ConsiderNormalVariation If consider the normal derivative
     * @param SlaveGeometry The slave geometry
     * @param MasterGeometry The master geometry
     * @param Normal The normal vector of the slave geometry
     * @details The  procedure will be the following in order to compute the derivative of the clipping
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
        const array_1d<array_1d<double, 3>, TDim * TNumNodes> all_delta_normal = DeltaNormalCenter(SlaveGeometry);
	array_1d<double, 3> zero_array(3, 0.0);
        array_1d<double, 3> delta_normal;

        const double aux_nodes_coeff = static_cast<double>(TNumNodes);

        const PointType slave_center = SlaveGeometry.Center();

//     #ifdef KRATOS_DEBUG
//         for (unsigned i_triangle = 0; i_triangle < 3; ++i_triangle)
//             KRATOS_WATCH(static_cast<IndexType>(TheseBelongs[i_triangle]));
//     #endif

        for (IndexType i_triangle = 0; i_triangle < 3; ++i_triangle) {
            if (TheseBelongs[i_triangle] >= 2 * TNumNodes) { // It belongs to an intersection
                // We compute the indexes
                IndexType belong_index_slave_start, belong_index_slave_end, belong_index_master_start, belong_index_master_end;
                ConvertAuxHashIndex(static_cast<IndexType>(TheseBelongs[i_triangle]), belong_index_slave_start, belong_index_slave_end, belong_index_master_start, belong_index_master_end);

                // The coordinates should be in the projected plane
                double distance;
                const array_1d<double, 3> xs1 = GeometricalProjectionUtilities::FastProject(slave_center, SlaveGeometry[belong_index_slave_start], Normal, distance).Coordinates(); // Start coordinates of the first segment
                const array_1d<double, 3> xe1 = GeometricalProjectionUtilities::FastProject(slave_center, SlaveGeometry[belong_index_slave_end], Normal, distance).Coordinates(); // End coordinates of the first segment
                const array_1d<double, 3> xs2 = GeometricalProjectionUtilities::FastProject(slave_center, MasterGeometry[belong_index_master_start], Normal, distance).Coordinates(); // Start coordinates of the second segment
                const array_1d<double, 3> xe2 = GeometricalProjectionUtilities::FastProject(slave_center, MasterGeometry[belong_index_master_end], Normal, distance).Coordinates(); // End coordinates of the second segment

                // We define the array containing the indexes of the vertexes
                array_1d<IndexType, 4> belong_indexes;
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

                for (IndexType i_belong = 0; i_belong < 4; ++i_belong) {
                    // The index of the node
                    const IndexType belong_index = belong_indexes[i_belong];

                    for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                        // We get the delta normal
                        if (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION && belong_index < TNumNodes) delta_normal = all_delta_normal[belong_index * TDim + i_dof] * (1.0/aux_nodes_coeff);
                        else delta_normal = zero_array;

                        auto& local_delta_vertex = rDerivativeData.DeltaCellVertex[belong_index * TDim + i_dof];

                        // Special cases (slave nodes)
                        if (i_belong == 0) { // First node of the slave
                            const double coeff = 1.0 + num/denom;
                            noalias(row(local_delta_vertex, i_triangle)) += LocalDeltaVertex( Normal,  delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, coeff);
                        } else if (i_belong == 1) { // Second node of the slave
                            const double coeff = - num/denom;
                            noalias(row(local_delta_vertex, i_triangle)) += LocalDeltaVertex( Normal,  delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, coeff);
                        }

                        // We define some auxiliar coefficients
                        const double coeff1 = - 1.0/denom;
                        const double coeff2 = num/std::pow(denom, 2);

                        // We add the part corresponding purely to delta normal
                        if (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION) {
                            noalias(row(local_delta_vertex, i_triangle)) += diff3 * coeff1 * inner_prod(aux_num,  delta_normal);
                            noalias(row(local_delta_vertex, i_triangle)) += diff3 * coeff2 * inner_prod(aux_denom, delta_normal);
                        }

                        // We compute the delta diffs
                        const array_1d<double, 3> delta_diff1 = (i_belong == 0) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0) : (i_belong == 2) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - 1.0) : zero_array;
                        const array_1d<double, 3> delta_diff2 = (i_belong == 3) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0) : (i_belong == 2) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - 1.0) : zero_array;
                        const array_1d<double, 3> delta_diff3 = (i_belong == 1) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0) : (i_belong == 0) ? LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - 1.0) : zero_array;

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
            } else { // It belongs to a master/slave node
                const IndexType belong_index = static_cast<IndexType>(TheseBelongs[i_triangle]);

                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof) {
                    // We get the delta normal
                    if (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION && belong_index < TNumNodes) 
                        delta_normal = all_delta_normal[belong_index * TDim + i_dof] * (1.0/aux_nodes_coeff);
                    else 
                        delta_normal = zero_array;

                    auto& local_delta_vertex = rDerivativeData.DeltaCellVertex[belong_index * TDim + i_dof];
                    noalias(row(local_delta_vertex, i_triangle)) += LocalDeltaVertex( Normal,  delta_normal, i_dof, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry);
                }
            }
        }
    }

    /**
     * @brief Calculates the increment of the shape functions
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
        // Auxiliar zero array
        const array_1d<double, 3> zero_array(3, 0.0);

        /* Shape functions */
        const VectorType& N1 = rVariables.NSlave;

        /* Local gradients */
        const MatrixType& DNDe1 = rVariables.DNDeSlave;

        // The Normal and delta Normal in the center of the element
        const array_1d<array_1d<double, 3>, TDim * TNumNodes> all_delta_normal = DeltaNormalCenter(SlaveGeometry);

        /* Shape function decomposition */
        VectorType N_decomp;
        DecompGeom.ShapeFunctionsValues( N_decomp, LocalPointDecomp.Coordinates() );

        if (TDim == 3) { // NOTE: This is not used in 2D
            for ( IndexType i_node = 0; i_node < 2 * TNumNodes; ++i_node) {
                for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                    // We get the delta normal
                    const array_1d<double, 3> delta_normal = (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION && i_node < TNumNodes) ? all_delta_normal[i_node * TDim + i_dof] : zero_array;

                    // We compute the residuals
                    array_1d<double, 3> aux_RHS1(3, 0.0);

                    // The vertex cell contribution
                    const auto& local_delta_cell = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof];
                    for(std::size_t i_belong = 0; i_belong < 3; ++i_belong) {
                        noalias(aux_RHS1) += N_decomp[i_belong] * row(local_delta_cell, i_belong);
                    }

                    // Local contribution
                    const array_1d<double, 3> aux_delta_node = LocalDeltaVertex( SlaveNormal, delta_normal, i_dof, i_node, ConsiderNormalVariation, SlaveGeometry, MasterGeometry );
                    if (i_node < TNumNodes) 
                        noalias(aux_RHS1) -= N1[i_node] * aux_delta_node;

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
     * @brief Calculates the increment of the shape functions
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
        // Auxiliar zero array
        const array_1d<double, 3> zero_array(3, 0.0);

        /* Shape functions */
        const VectorType& N1 = rVariables.NSlave;
        const VectorType& N2 = rVariables.NMaster;

        /* Local gradients */
        const MatrixType& DNDe1 = rVariables.DNDeSlave;
        const MatrixType& DNDe2 = rVariables.DNDeMaster;

        // The Normal and delta Normal in the center of the element
        const array_1d<array_1d<double, 3>, TDim * TNumNodes> all_delta_normal = DeltaNormalCenter(SlaveGeometry);

        /* Shape function decomposition */
        VectorType N_decomp;
        DecompGeom.ShapeFunctionsValues( N_decomp, LocalPointDecomp.Coordinates() );

        /* Tolerance */
        const double tolerance = std::numeric_limits<double>::epsilon();

        if (TDim == 2) {
            array_1d<PointType, TNumNodes> projected_in_slave, projected_in_master;

            for (IndexType i_mortar_node = 0; i_mortar_node < TNumNodes; ++i_mortar_node) {
                // Projecting points in opposite geometry, defining mortar nodes
                GeometricalProjectionUtilities::FastProjectDirection( SlaveGeometry,  MasterGeometry[i_mortar_node], projected_in_slave[i_mortar_node],  SlaveNormal, MasterNormal );
                GeometricalProjectionUtilities::FastProjectDirection( MasterGeometry, SlaveGeometry[i_mortar_node],  projected_in_master[i_mortar_node], MasterNormal, SlaveNormal );
            }

            for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {

                    array_1d<double, TNumNodes> DeltaXi_slave  = ZeroVector(TNumNodes);
                    array_1d<double, TNumNodes> DeltaXi_master = ZeroVector(TNumNodes);
                    double DeltaXi1 = 0.0, DeltaXi2 = 0.0;

                    for (IndexType i_mortar_node = 0; i_mortar_node < TNumNodes; ++i_mortar_node) {
                        // Auxiliary point for local coordinates
                        PointType aux_point_slave, aux_point_master;

                        // Compute DeltaXi on the slave side if point is inside geometry (nodes of geometry excluded)
                        if (SlaveGeometry.IsInside( projected_in_slave[i_mortar_node].Coordinates(), aux_point_slave.Coordinates() ) && (norm_2(projected_in_slave[i_mortar_node].Coordinates() - SlaveGeometry[0].Coordinates()) > tolerance) && (norm_2(projected_in_slave[i_mortar_node].Coordinates() - SlaveGeometry[1].Coordinates()) > tolerance))
                            DeltaXi_slave[i_mortar_node] = LocalDeltaSegmentN1( all_delta_normal, SlaveNormal, SlaveGeometry, MasterGeometry, N1, DNDe1, i_mortar_node, i_node, i_dof, ConsiderNormalVariation );

                        // Compute DeltaXi on the master side if point is inside geometry (nodes of geometry excluded)
                        if (MasterGeometry.IsInside( projected_in_master[i_mortar_node].Coordinates(), aux_point_master.Coordinates() ) && (norm_2(projected_in_master[i_mortar_node].Coordinates() - MasterGeometry[0].Coordinates()) > tolerance) && (norm_2(projected_in_master[i_mortar_node].Coordinates() - MasterGeometry[1].Coordinates()) > tolerance))
                            DeltaXi_master[i_mortar_node] = LocalDeltaSegmentN2( all_delta_normal, SlaveNormal, SlaveGeometry, MasterGeometry, N2, DNDe2, i_mortar_node, i_node, i_dof, ConsiderNormalVariation );
                    }

                    // Evaluate DeltaXi1 expression
                    DeltaXi1 = inner_prod(N_decomp, DeltaXi_slave);

                    // Evaluate DeltaXi2 expression
                    DeltaXi2 = inner_prod(N_decomp, DeltaXi_master);

                    // Multiply for DNDe for obtaining DeltaN

                    auto& delta_n1 = rDerivativeData.DeltaN1[i_node * TDim + i_dof];
                    noalias(delta_n1) = DeltaXi1 * column(DNDe1, 0);

                    auto& delta_n2 = rDerivativeData.DeltaN2[i_node * TDim + i_dof];
                    noalias(delta_n2) = DeltaXi2 * column(DNDe2, 0);

                    // The derivatives of the dual shape function
                    auto& delta_phi = rDerivativeData.DeltaPhi[i_node * TDim + i_dof];
                    if (DualLM == true) {
                        noalias(delta_phi) = prod(rDerivativeData.Ae, delta_n1);
                        if (i_node >= TNumNodes)
                            noalias(delta_phi) += prod(rDerivativeData.DeltaAe[i_node * TDim + i_dof], N1);
                    } else {
                        noalias(delta_phi) = delta_n1;
                    }
                }
            }
        } else {
            for ( IndexType i_node = 0; i_node < 2 * TNumNodes; ++i_node) {
                for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                    // We get the delta normal
                    const array_1d<double, 3> delta_normal = (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION && i_node < TNumNodes) ? all_delta_normal[i_node * TDim + i_dof] : zero_array;

                    // We compute the residuals
                    array_1d<double, 3> aux_RHS1(3, 0.0);

                    // The vertex cell contribution
                    const auto& local_delta_cell = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof];
                    for(std::size_t i_belong = 0; i_belong < 3; ++i_belong) {
                        noalias(aux_RHS1) += N_decomp[i_belong] * row(local_delta_cell, i_belong);
                    }
                    array_1d<double, 3> aux_RHS2 = aux_RHS1;

                    // Local contribution
                    const array_1d<double, 3> aux_delta_node = LocalDeltaVertex( SlaveNormal, delta_normal, i_dof, i_node, ConsiderNormalVariation, SlaveGeometry, MasterGeometry );
                    if (i_node < TNumNodes) 
                        noalias(aux_RHS1) -= N1[i_node] * aux_delta_node;
                    else 
                        noalias(aux_RHS2) -= N2[i_node - TNumNodes] * aux_delta_node;

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
                    if (DualLM == true) {
                        noalias(delta_phi) = prod(rDerivativeData.DeltaAe[i_node * TDim + i_dof], N1);
                        noalias(delta_phi) += prod(rDerivativeData.Ae, delta_n1);
                    } else {
                        noalias(delta_phi) = delta_n1;
                    }
                }
            }
        }
    }

    /**
     * @brief Returns a matrix with the increment of displacements, that can be used for compute the Jacobian reference (current) configuration
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

        for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node ) {
            const array_1d<double, 3 > delta_displacement = ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);

            for ( IndexType j_node = 0; j_node < TDim; ++j_node ) {
                Vector N;
                ThisGeometry.ShapeFunctionsValues( N, LocalCoordinates[j_node].Coordinates() );

                for ( IndexType j_dim = 0; j_dim < TDim; ++j_dim )
                    DeltaPosition(j_node, j_dim) += N[i_node] * delta_displacement[j_dim];
            }
        }

        return DeltaPosition;

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Returns a matrix with the increment of displacements
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

        for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node ) {
            const array_1d<double, 3 > delta_displacement = ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);

            for ( IndexType i_dim = 0; i_dim < TDim; ++i_dim )
                DeltaPosition(i_node, i_dim) += delta_displacement[i_dim];
        }

        return DeltaPosition;

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Returns a vector with the increment of displacements
     * @param DeltaPosition The resulting vector with the increment of position
     * @param SlaveGeometry The slave geometry
     * @param MasterGeometry The master geometry
     * @param IndexNode The node index
     */

    static inline void CalculateDeltaPosition(
        VectorType& DeltaPosition,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const IndexType IndexNode
        )
    {
        KRATOS_TRY;

        if (IndexNode < TNumNodes) {
            DeltaPosition = SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1);
        } else {
            const IndexType index_master = IndexNode - TNumNodes;
            DeltaPosition = MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1);
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Returns a vector with the increment of displacements
     * @param DeltaPosition The resulting vector with the increment of position
     * @param SlaveGeometry The slave geometry
     * @param MasterGeometry The master geometry
     * @param IndexNode The node index
     * @param iDoF The degree of freedom index
     */

    static inline void CalculateDeltaPosition(
        VectorType& DeltaPosition,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const IndexType IndexNode,
        const IndexType iDoF
        )
    {
        KRATOS_TRY;

        DeltaPosition = ZeroVector(3);

        if (IndexNode < TNumNodes) {
            DeltaPosition[iDoF] = (SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        } else {
            const IndexType index_master = IndexNode - TNumNodes;
            DeltaPosition[iDoF] = (MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Returns a double with the increment of displacements
     * @param DeltaPosition The resulting double with the increment of position
     * @param SlaveGeometry The slave geometry
     * @param MasterGeometry The master geometry
     * @param IndexNode The node index
     * @param iDoF The degree of freedom index
     */

    static inline void CalculateDeltaPosition(
        double& DeltaPosition,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const IndexType IndexNode,
        const IndexType iDoF
        )
    {
        KRATOS_TRY;

        if (IndexNode < TNumNodes) {
            DeltaPosition = (SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        } else {
            const IndexType index_master = IndexNode - TNumNodes;
            DeltaPosition = (MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Calculate Ae and DeltaAe matrices
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

        for (IndexType i_geom = 0; i_geom < ConditionsPointsSlave.size(); ++i_geom) {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                PointType global_point;
                SlaveGeometry.GlobalCoordinates(global_point, ConditionsPointsSlave[i_geom][i_node]);
                points_array[i_node] = PointType::Pointer( new PointType(global_point) );
                belong_array[i_node] = ConditionsPointsSlave[i_geom][i_node].GetBelong();
            }

            DecompositionType decomp_geom( points_array );

            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, SlaveGeometry.Length() * 1.0e-12) : MortarUtilities::HeronCheck(decomp_geom);

            if (bad_shape == false) {
                const GeometryType::IntegrationPointsArrayType integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );

                // Integrating the mortar operators
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                    // We reset the derivatives
                    rDerivativeData.ResetDerivatives();

                    // We compute the local coordinates
                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                    PointType local_point_parent;
                    PointType gp_global;

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
                    GeometricalProjectionUtilities::FastProjectDirection( MasterGeometry, slave_gp_global, projected_gp_global, SlaveNormal, -gp_normal ); // The opposite direction

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
     * @brief This method is used to compute the directional derivatives of the cell vertex (locally)
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
        const IndexType iDoF,
        const IndexType iBelong,
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
        const array_1d<double, 3> coords_center = SlaveGeometry.Center().Coordinates();
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
        if (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION) 
            aux_delta_vertex += coordsxnormal * DeltaNormal;

        return Coeff * aux_delta_vertex;
    }

    /**
     * @brief This method computes the auxiliar matrix used to keep unitary the normal
     * @param DiffVector The auxiliar vector of difference of two normal vectors
     * @param DeltaNormal The vector containing the delta normal
     * @return The auxiliar matrix computed
     */
    static inline BoundedMatrix<double, 3, 3> ComputeRenormalizerMatrix(
        const array_1d<double, 3>& DiffVector,
        const array_1d<double, 3>& DeltaNormal
        )
    {
        for (IndexType itry = 0; itry < 3; ++itry) {
            if (DeltaNormal[itry] > std::numeric_limits<double>::epsilon()) {
                BoundedMatrix<double, 3, 3> aux_matrix;

                const IndexType aux_index_1 = itry == 2 ? 0 : itry + 1;
                const IndexType aux_index_2 = itry == 2 ? 1 : (itry == 1 ? 0 : 2);

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
     * @brief This method computes the auxiliar matrix used to keep unitary the normal
     * @param DiffMatrix The auxiliar matrix of difference of two normal matrices
     * @param DeltaNormal The matrix containing the delta normal
     * @param iGeometry The index of the node of the geometry computed
     * @return The auxiliar matrix computed
     */
    static inline BoundedMatrix<double, 3, 3> ComputeRenormalizerMatrix(
        const BoundedMatrix<double, TNumNodes, TDim>& DiffMatrix,
        const BoundedMatrix<double, TNumNodes, TDim>& DeltaNormal,
        const IndexType iGeometry
        )
    {
        for (IndexType itry = 0; itry < 3; ++itry) {
            if (DeltaNormal(iGeometry, itry) > std::numeric_limits<double>::epsilon()) {
                BoundedMatrix<double, 3, 3> aux_matrix;

                const IndexType aux_index_1 = itry == 2 ? 0 : itry + 1;
                const IndexType aux_index_2 = itry == 2 ? 1 : (itry == 1 ? 0 : 2);

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
     * @brief This method computes the normal in the previous configuration
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
        array_1d<double,3> tangent_xi = ZeroVector(3), tangent_eta = ZeroVector(3);

        // Using the Jacobian tangent directions
        if (TDim == 2) {
            tangent_eta[2] = 1.0;
            for (IndexType i_dim = 0; i_dim < TDim; i_dim++)
                tangent_xi[i_dim]  = previous_jacobian(i_dim, 0);
        } else {
            for (IndexType i_dim = 0; i_dim < TDim; i_dim++) {
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
     * @brief This method computes the equivalent indexes to the auxiliar hash
     * @param AuxIndex The auxiliar index to decompose
     * @param iBelongSlaveStart The index of the first/slave segment and first node
     * @param iBelongSlaveEnd The index of the first/slave segment and end node
     * @param iBelongMasterStart The index of the second/master segment and first node
     * @param iBelongMasterEnd The index of the second/master segment and end node
     */
    static inline void ConvertAuxHashIndex(
        const IndexType AuxIndex,
        IndexType& iBelongSlaveStart,
        IndexType& iBelongSlaveEnd,
        IndexType& iBelongMasterStart,
        IndexType& iBelongMasterEnd
        )
    {
        IndexType index_to_decompose = AuxIndex - 2 * TNumNodes;

        iBelongMasterEnd = index_to_decompose/10000;
        index_to_decompose = std::fmod(index_to_decompose, 10000);
        iBelongMasterStart = index_to_decompose/1000;
        index_to_decompose = std::fmod(index_to_decompose, 1000);
        iBelongSlaveEnd = index_to_decompose/100;
        index_to_decompose = std::fmod(index_to_decompose, 100);
        iBelongSlaveStart = index_to_decompose/10;
    }

    /**
     * @brief This method computes the increment of local coordinates
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

        BoundedMatrix<double, 3, TNumNodes> X;
        for(IndexType i = 0; i < TNumNodes; ++i) {
            X(0, i) = ThisGeometry[i].X();
            X(1, i) = ThisGeometry[i].Y();
            X(2, i) = ThisGeometry[i].Z();
        }

        const BoundedMatrix<double, 3, 2> DN = prod(X, rDNDe);

        const BoundedMatrix<double, 2, 2> J = prod(trans(DN),DN);
        double det_j = MathUtils<double>::DetMat<BoundedMatrix<double, 2, 2>>(J);
        const BoundedMatrix<double, 2, 2> invJ = (std::abs(det_j) < tolerance) ? ZeroMatrix(2,2) : MathUtils<double>::InvertMatrix<2>(J, det_j);

    #ifdef KRATOS_DEBUG
        if (std::abs(det_j) < tolerance) 
            KRATOS_WARNING("Jacobian invert") << "WARNING: CANNOT INVERT JACOBIAN TO COMPUTE DELTA COORDINATES" << std::endl;
    #endif

        const array_1d<double, 2> res = prod(trans(DN), DeltaPoint);
        noalias(rResult) = prod(invJ, res);

//         BoundedMatrix<double, 3, 3> L;
//         for(IndexType i = 0; i < 3; ++i)
//         {
//             for(IndexType j = 0; j < 2; ++j) L(i, j) = DN(i, j);
//             L(i, 2) = ThisNormal[i];
//         }
//
//         double det_L = MathUtils<double>::DetMat<BoundedMatrix<double, 3, 3>>(L);
//         const BoundedMatrix<double, 3, 3> invL = (std::abs(det_L) < tolerance) ? ZeroMatrix(3,3) : MathUtils<double>::InvertMatrix<3>(L, det_L);
//         array_1d<double, 3> aux = prod(invL, DeltaPoint);
//         rResult[0] = aux[0];
//         rResult[1] = aux[1];
//         #ifdef KRATOS_DEBUG
//             if (std::abs(det_L) < tolerance) 
//                 KRATOS_WARNING("Jacobian invert") << "WARNING: CANNOT INVERT JACOBIAN TO COMPUTE DELTA COORDINATES" << std::endl;
//         #endif
    }

    /**
     * @brief This method is used to compute the directional derivatives of the mortar segment on the slave side
     * @param DeltaNormal All derivatives of normals in points of both geometries
     * @param SlaveNormal Normal in the center of slave geometry
     * @param SlaveGeometry Slave geometry where to compute the value
     * @param MasterGeometry Master node projected to the Slave Geometry
     * @param N1 Shape function at SlaveGeometry
     * @param DNDe1 Gradient of shape function N1
     * @param MortarNode Index of mortar node where computation is being carried
     * @param iNode Actual node where computation is being carried
     * @param iDoF Direction in which computation is carried
     * @param ConsiderNormalVariation Flag to determine if consider delta_normal
     * @return The mortar node derivative
     */
    static inline double LocalDeltaSegmentN1(
        const array_1d<array_1d<double, 3>, TDim * TNumNodes>& DeltaNormal,
        const array_1d<double, 3> SlaveNormal,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const VectorType& N1,
        const MatrixType& DNDe1,
        const IndexType MortarNode,
        const IndexType iNode,
        const IndexType iDoF,
        const NormalDerivativesComputation ConsiderNormalVariation
        )
    {
        array_1d<double, TDim> Xa, DXa;
        BoundedMatrix<double, TDim, TNumNodes> X1, DX1;
        BoundedMatrix<double, 3, TNumNodes> n1, Dn1;

        // Projected node coordinates
        Xa[0] = MasterGeometry[MortarNode].X();
        Xa[1] = MasterGeometry[MortarNode].Y();

        // Projected node derivatives
        DXa = ZeroVector(TDim);
        DXa[iDoF] = 1.0;

        // Slave element nodes coordinates and normals with respective derivatives
        for(IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            X1(0, i_node) = SlaveGeometry[i_node].X();
            X1(1, i_node) = SlaveGeometry[i_node].Y();

            for(IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
               if(i_dof == iDoF)
                   DX1(i_dof, i_node) = 1.0;
               else
                   DX1(i_dof, i_node) = 0.0;
            }

            column(n1, i_node)  = SlaveNormal;
            column(Dn1, i_node) = (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION) ? DeltaNormal[i_node * TDim + iDoF] : ZeroVector(3);
        }

        // Computation of DeltaXi_a
        const double denom = -inner_prod(row(X1, 0), column(DNDe1, 0))*(N1[iNode]*n1(1,iNode)) + inner_prod(row(X1, 1), column(DNDe1,0))*(N1[iNode]*n1(0,iNode)) - (N1[iNode]*X1(0,iNode) - Xa[0])*inner_prod(row(n1, 1), column(DNDe1, 0)) + (N1[iNode]*X1(1,iNode) - Xa[1])*inner_prod(row(n1, 0), column(DNDe1, 0));
        const double num = (N1[iNode]*DX1(0,iNode) - DXa[0])*(N1[iNode]*n1(1,iNode)) - (N1[iNode]*DX1(1,iNode) - DXa[1])*(N1[iNode]*n1(0,iNode)) + (N1[iNode]*X1(0,iNode) - Xa[0])*(N1[iNode]*Dn1(1,iNode)) - (N1[iNode]*X1(1,iNode) - Xa[1])*(N1[iNode]*Dn1(0,iNode));

        return num/denom;
    }

    /**
     * @brief This method is used to compute the directional derivatives of the mortar segment on the master side
     * @param DeltaNormal All derivatives of normals in points of both geometries
     * @param SlaveNormal Normal in the center of slave geometry
     * @param SlaveGeometry Slave geometry where to compute the value
     * @param MasterGeometry Master node projected to the Slave Geometry
     * @param N2 Shape function at MasterGeometry
     * @param DNDe2 Gradient of shape function N2
     * @param MortarNode Index of mortar node where computation is being carried
     * @param iNode Actual node where computation is being carried
     * @param iDoF Direction in which computation is carried
     * @param ConsiderNormalVariation flag to determine if consider delta_normal
     * @return The mortar node derivative
     */
     static inline double LocalDeltaSegmentN2(
         const array_1d<array_1d<double, 3>, TDim * TNumNodes>& DeltaNormal,
         const array_1d<double, 3> SlaveNormal,
         const GeometryType& SlaveGeometry,
         const GeometryType& MasterGeometry,
         const VectorType& N2,
         const MatrixType& DNDe2,
         const IndexType MortarNode,
         const IndexType iNode,
         const IndexType iDoF,
         const NormalDerivativesComputation ConsiderNormalVariation
         )
     {
         array_1d<double, TDim> Xa, DXa;
         array_1d<double, 3>    na, Dna;
         BoundedMatrix<double, TDim, TNumNodes> X2, DX2;

         // Projected node coordinates
         Xa[0] = SlaveGeometry[MortarNode].X();
         Xa[1] = SlaveGeometry[MortarNode].Y();

         // Projected node derivatives
         DXa = ZeroVector(TDim);
         DXa[iDoF] = 1.0;

         // Projected normal and derivative
         na = SlaveNormal;
         Dna = (ConsiderNormalVariation != NO_DERIVATIVES_COMPUTATION) ? DeltaNormal[MortarNode * TDim + iDoF]: ZeroVector(3);

         // Slave element nodes coordinates and derivatives
         for(IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
             X2(0, i_node) = MasterGeometry[i_node].X();
             X2(1, i_node) = MasterGeometry[i_node].Y();

             for(IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                if(i_dof == iDoF)
                    DX2(i_dof, i_node) = 1.0;
                else
                    DX2(i_dof, i_node) = 0.0;
             }
         }

         // Computation of DeltaXi_a
         const double lhs = -1.0/(inner_prod(row(X2, 0), column(DNDe2,0))*na[1] - inner_prod(row(X2, 1), column(DNDe2,0))*na[0]);
         const double rhs = (N2[iNode]*DX2(0, iNode)-DXa[0])*na[1] - (N2[iNode]*DX2(1, iNode)-DXa[1])*na[0] + (N2[iNode]*X2(0, iNode)-Xa[0])*Dna[1] - (N2[iNode]*X2(1, iNode)-Xa[1])*Dna[0];

         return lhs*rhs;
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
