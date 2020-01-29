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

// System includes

// External includes

// Project includes
#include "custom_utilities/derivatives_utilities.h"
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/mortar_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{
template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaDetjSlave(
    const DecompositionType& DecompGeom,
    const GeneralVariables& rVariables,
    DerivativeDataType& rDerivativeData
    )
{
    KRATOS_TRY

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

        array_1d<double, 3> aux_delta_cross_product1, aux_delta_cross_product2;
        for ( IndexType i_node = 0; i_node < (TNumNodesMaster + TNumNodes); ++i_node ) {
            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                const auto& r_local_delta_vertex = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof];

                MathUtils<double>::CrossProduct(aux_delta_cross_product1, row(r_local_delta_vertex, 1) - row(r_local_delta_vertex, 0), x31cell);
                MathUtils<double>::CrossProduct(aux_delta_cross_product2, x21cell, row(r_local_delta_vertex, 2) - row(r_local_delta_vertex, 0));

                rDerivativeData.DeltaDetjSlave[i_node * TDim + i_dof] = inner_prod(aux_cross_product, aux_delta_cross_product1) + inner_prod(aux_cross_product, aux_delta_cross_product2);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline array_1d<array_1d<double, 3>, TDim * TNumNodes> DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::GPDeltaNormalSlave(
    const Matrix& rJacobian,
    const Matrix& rDNDe
    )
{
    KRATOS_TRY

    // Tangent directions
    array_1d<double,3> j0 = ZeroVector(3), j1 = ZeroVector(3);

    // Using the Jacobian tangent directions
    if (TDim == 2) {
        j1[2] = 1.0;
        for (IndexType i_dim = 0; i_dim < 2; ++i_dim)
            j0[i_dim]  = rJacobian(i_dim, 0);
    } else {
        for (IndexType i_dim = 0; i_dim < 3; ++i_dim) {
            j0[i_dim] = rJacobian(i_dim, 0);
            j1[i_dim] = rJacobian(i_dim, 1);
        }
    }

    array_1d<double, 3> normal;;
    MathUtils<double>::CrossProduct(normal, j0, j1);
    const double area_normal_norm = norm_2(normal);

    KRATOS_DEBUG_ERROR_IF(area_normal_norm < ZeroTolerance) << "ZERO NORMAL: " << area_normal_norm << std::endl;

    const array_1d<double, 3> unit_normal = normal/area_normal_norm;

    array_1d<array_1d<double, 3>, TDim * TNumNodes> delta_normal;
    array_1d<double,3> delta_j0, delta_j1;
    const array_1d<double, 3> zero_array = ZeroVector(3);
    array_1d<double, 3> aux_delta_normal0, aux_delta_normal1;
    for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node ) {
        for(IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
            const IndexType i_dof = i_node * TDim + i_dim;

            delta_j0 = zero_array;
            delta_j1 = zero_array;

            delta_j0[i_dim] += rDNDe(i_node, 0);
            if (TDim == 3)
                delta_j1[i_dim] += rDNDe(i_node, 1);

            MathUtils<double>::CrossProduct(aux_delta_normal0, j0, delta_j1);
            MathUtils<double>::CrossProduct(aux_delta_normal1, delta_j0, j1);
            const array_1d<double, 3> r_aux_delta_normal = aux_delta_normal0 + aux_delta_normal1;
            const double delta_norm = inner_prod(r_aux_delta_normal, normal);

            delta_normal[i_dof] = (r_aux_delta_normal + unit_normal * delta_norm)/area_normal_norm;
        }
    }

    return delta_normal;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline array_1d<array_1d<double, 3>, TDim * TNumNodesMaster> DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::GPDeltaNormalMaster(
    const Matrix& rJacobian,
    const Matrix& rDNDe
    )
{
    KRATOS_TRY

    // Tangent directions
    array_1d<double,3> j0 = ZeroVector(3), j1 = ZeroVector(3);

    // Using the Jacobian tangent directions
    if (TDim == 2) {
        j1[2] = 1.0;
        for (IndexType i_dim = 0; i_dim < 2; ++i_dim)
            j0[i_dim]  = rJacobian(i_dim, 0);
    } else {
        for (IndexType i_dim = 0; i_dim < 3; ++i_dim) {
            j0[i_dim] = rJacobian(i_dim, 0);
            j1[i_dim] = rJacobian(i_dim, 1);
        }
    }

    array_1d<double, 3> normal;;
    MathUtils<double>::CrossProduct(normal, j0, j1);
    const double area_normal_norm = norm_2(normal);

    KRATOS_DEBUG_ERROR_IF(area_normal_norm < ZeroTolerance) << "ZERO NORMAL: " << area_normal_norm << std::endl;

    const array_1d<double, 3> unit_normal = normal/area_normal_norm;

    array_1d<array_1d<double, 3>, TDim * TNumNodesMaster> delta_normal;
    array_1d<double,3> delta_j0, delta_j1, r_aux_delta_normal;
    const array_1d<double, 3> zero_array = ZeroVector(3);
    array_1d<double, 3> aux_delta_normal0, aux_delta_normal1;
    for ( IndexType i_node = 0; i_node < TNumNodesMaster; ++i_node ) {
        for(IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
            const IndexType i_dof = i_node * TDim + i_dim;

            noalias(delta_j0) = zero_array;
            noalias(delta_j1) = zero_array;

            delta_j0[i_dim] += rDNDe(i_node, 0);
            if (TDim == 3)
                delta_j1[i_dim] += rDNDe(i_node, 1);

            MathUtils<double>::CrossProduct(aux_delta_normal0, j0, delta_j1);
            MathUtils<double>::CrossProduct(aux_delta_normal1, delta_j0, j1);
            noalias(r_aux_delta_normal) = aux_delta_normal0 + aux_delta_normal1;
            const double delta_norm = inner_prod(r_aux_delta_normal, normal);

            delta_normal[i_dof] = (r_aux_delta_normal + unit_normal * delta_norm)/area_normal_norm;
        }
    }

    return delta_normal;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
array_1d<array_1d<double, 3>, TDim * TNumNodes> DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::DeltaNormalCenter(const GeometryType& rThisGeometry)
{
    KRATOS_TRY

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
    const array_1d<array_1d<double, 3>, TDim * TNumNodes> aux_delta_normal_0 = GPDeltaNormalSlave(jacobian, gradient);

    for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const array_1d<double, 3> delta_disp = rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT, 1);
        for ( IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
            const array_1d<double, 3>& r_aux_delta_normal = aux_delta_normal_0[i_node * TDim + i_dof];
            aux_delta_normal0 += delta_disp[i_dof] * r_aux_delta_normal;
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
            const array_1d<double, 3>& r_aux_delta_normal = aux_delta_normal_0[i_node * TDim + i_dof];
            normalized_delta_normal_0[i_node * TDim + i_dof] = prod(renormalizer_matrix, r_aux_delta_normal);
        }
    }

    return normalized_delta_normal_0;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaNormalSlave(
    array_1d<BoundedMatrix<double, TNumNodes, TDim>, TNumNodes * TDim>& rDeltaNormal,
    GeometryType& rThisGeometry
    )
{
    KRATOS_TRY

    BoundedMatrix<double, TNumNodes, TDim> aux_normal_geometry = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(rThisGeometry,  NORMAL, 1);
    BoundedMatrix<double, TNumNodes, TDim> aux_delta_normal_geometry = ZeroMatrix(TNumNodes, TDim);

    // Auxiliar variables
    GeometryType::CoordinatesArrayType point_local;
    Matrix jacobian, gradient;
    array_1d<array_1d<double, 3>, TDim * TNumNodes> delta_normal_node;
    array_1d<double, TDim> aux_delta_normal;

    // Iterate over geometry nodes
    for ( IndexType i_geometry = 0; i_geometry < TNumNodes; ++i_geometry ) {
        // We compute the gradient and jacobian
        rThisGeometry.PointLocalCoordinates( point_local, rThisGeometry[i_geometry].Coordinates( ));
        jacobian = rThisGeometry.Jacobian( jacobian, point_local);
        rThisGeometry.ShapeFunctionsLocalGradients( gradient, point_local );

        // We compute the delta normal of the node
        delta_normal_node = GPDeltaNormalSlave(jacobian, gradient);
        for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3> delta_disp = rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT, 1);
            for ( IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                const array_1d<double, 3>& r_delta_normal = delta_normal_node[i_node * TDim + i_dof];
                for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                    aux_delta_normal[i_dim] = r_delta_normal[i_dim];
                }
                for (IndexType j_dof = 0; j_dof < TDim; ++j_dof) {
                    aux_delta_normal_geometry(i_geometry, j_dof) += delta_disp[i_dof] * aux_delta_normal[j_dof];
                }
            }
        }
    }

    BoundedMatrix<double, TNumNodes, TDim> calculated_normal_geometry = aux_delta_normal_geometry + aux_normal_geometry;
    for ( IndexType i_geometry = 0; i_geometry < TNumNodes; ++i_geometry )
        row(calculated_normal_geometry, i_geometry) /= norm_2(row(calculated_normal_geometry, i_geometry));

    // We compute the diff matrix to compute the auxiliar matrix later
    const BoundedMatrix<double, TNumNodes, TDim> diff_matrix = calculated_normal_geometry - aux_normal_geometry;

    // Auxilary terms
    BoundedMatrix<double, TDim, TDim> renormalizer_matrix;

    // We iterate over the nodes of the geometry
    for ( IndexType i_geometry = 0; i_geometry < TNumNodes; ++i_geometry ) {
        // Computing auxiliar matrix
        noalias(renormalizer_matrix) = (TDim == 3) ? ComputeRenormalizerMatrix(diff_matrix, aux_delta_normal_geometry, i_geometry) : IdentityMatrix(2, 2);

        // We compute the gradient and jacobian
        rThisGeometry.PointLocalCoordinates( point_local, rThisGeometry[i_geometry].Coordinates( ) ) ;
        jacobian = rThisGeometry.Jacobian( jacobian, point_local);
        rThisGeometry.ShapeFunctionsLocalGradients( gradient, point_local );

        // We compute the delta normal of the node
        delta_normal_node = GPDeltaNormalSlave(jacobian, gradient);
        for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            for ( IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                const array_1d<double, 3>& r_delta_normal = delta_normal_node[i_node * TDim + i_dof];
                for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                    aux_delta_normal[i_dim] = r_delta_normal[i_dim];
                }
                row(rDeltaNormal[i_node * TDim + i_dof], i_geometry) = prod(renormalizer_matrix, aux_delta_normal);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaNormalMaster(
    array_1d<BoundedMatrix<double, TNumNodesMaster, TDim>, TNumNodesMaster * TDim>& rDeltaNormal,
    const GeometryType& rThisGeometry
    )
{
    KRATOS_TRY

    BoundedMatrix<double, TNumNodesMaster, TDim> aux_normal_geometry = MortarUtilities::GetVariableMatrix<TDim,TNumNodesMaster>(rThisGeometry,  NORMAL, 1);
    BoundedMatrix<double, TNumNodesMaster, TDim> aux_delta_normal_geometry = ZeroMatrix(TNumNodesMaster, TDim);

    // We iterate over the nodes of the geometry
    Matrix jacobian, gradient;
    GeometryType::CoordinatesArrayType point_local;
    array_1d<double, TDim> aux_delta_normal;
    for ( IndexType i_geometry = 0; i_geometry < TNumNodesMaster; ++i_geometry ) {
        // We compute the gradient and jacobian
        rThisGeometry.PointLocalCoordinates( point_local, rThisGeometry[i_geometry].Coordinates( ) ) ;
        jacobian = rThisGeometry.Jacobian( jacobian, point_local);
        rThisGeometry.ShapeFunctionsLocalGradients( gradient, point_local );

        // We compute the delta normal of the node
        const array_1d<array_1d<double, 3>, TDim * TNumNodesMaster> delta_normal_node = GPDeltaNormalMaster(jacobian, gradient);
        for ( IndexType i_node = 0; i_node < TNumNodesMaster; ++i_node) {
            const array_1d<double, 3> delta_disp = rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT, 1);
            for ( IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                const array_1d<double, 3>& delta_normal = delta_normal_node[i_node * TDim + i_dof];
                for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                    aux_delta_normal[i_dim] = delta_normal[i_dim];
                }
                for (IndexType j_dof = 0; j_dof < TDim; ++j_dof) {
                    aux_delta_normal_geometry(i_geometry, j_dof) += delta_disp[i_dof] * aux_delta_normal[j_dof];
                }
            }
        }
    }

    BoundedMatrix<double, TNumNodesMaster, TDim> calculated_normal_geometry = aux_delta_normal_geometry + aux_normal_geometry;
    for ( IndexType i_geometry = 0; i_geometry < TNumNodesMaster; ++i_geometry )
        row(calculated_normal_geometry, i_geometry) /= norm_2(row(calculated_normal_geometry, i_geometry));

    // We compute the diff matrix to compute the auxiliar matrix later
    const BoundedMatrix<double, TNumNodesMaster, TDim> diff_matrix = calculated_normal_geometry - aux_normal_geometry;

    // We iterate over the nodes of the geometry
    BoundedMatrix<double, TDim, TDim> renormalizer_matrix;
    array_1d<array_1d<double, 3>, TDim * TNumNodesMaster> delta_normal_node;
    for ( IndexType i_geometry = 0; i_geometry < TNumNodesMaster; ++i_geometry ) {
        // Computing auxiliar matrix
        noalias(renormalizer_matrix) = (TDim == 3) ? ComputeRenormalizerMatrix(diff_matrix, aux_delta_normal_geometry, i_geometry) : IdentityMatrix(2, 2);

        // We compute the gradient and jacobian
        rThisGeometry.PointLocalCoordinates( point_local, rThisGeometry[i_geometry].Coordinates( ) ) ;

        jacobian = rThisGeometry.Jacobian( jacobian, point_local);
        rThisGeometry.ShapeFunctionsLocalGradients( gradient, point_local );

        // We compute the delta normal of the node
        delta_normal_node = GPDeltaNormalMaster(jacobian, gradient);
        for ( IndexType i_node = 0; i_node < TNumNodesMaster; ++i_node) {
            for ( IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                const array_1d<double, 3>& delta_normal = delta_normal_node[i_node * TDim + i_dof];
                for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                    aux_delta_normal[i_dim] = delta_normal[i_dim];
                }
                row(rDeltaNormal[i_node * TDim + i_dof], i_geometry) = prod(renormalizer_matrix, aux_delta_normal);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaCellVertex(
    const GeneralVariables& rVariables,
    DerivativeDataType& rDerivativeData,
    const array_1d<BelongType, TDim>& rTheseBelongs,
    const NormalDerivativesComputation ConsiderNormalVariation,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    const array_1d<double, 3>& rNormal
    )
{
    KRATOS_TRY

    // The Normal and delta Normal in the center of the element
    const array_1d<array_1d<double, 3>, TDim * TNumNodes> all_delta_normal = DeltaNormalCenter(rSlaveGeometry);
    array_1d<double, 3> zero_array = ZeroVector(3);
    array_1d<double, 3> delta_normal;

    const double aux_nodes_coeff = static_cast<double>(TNumNodes);

    const PointType slave_center = rSlaveGeometry.Center();

//     #ifdef KRATOS_DEBUG
//         for (unsigned i_triangle = 0; i_triangle < 3; ++i_triangle)
//             KRATOS_WATCH(static_cast<IndexType>(TheseBelongs[i_triangle]));
//     #endif

    double distance;
    array_1d<double, 3> xs1, xe1, xs2, xe2;
    array_1d<double, 3> diff1, diff2, diff3;
    array_1d<double, 3> delta_diff1, delta_diff2, delta_diff3;
    array_1d<double, 3> aux_num, aux_denom, aux_cross_product;
    array_1d<IndexType, 4> belong_indexes;
    IndexType belong_index_slave_start, belong_index_slave_end, belong_index_master_start, belong_index_master_end;
    for (IndexType i_triangle = 0; i_triangle < 3; ++i_triangle) {
        if (static_cast<IndexType>(rTheseBelongs[i_triangle]) >= (TNumNodesMaster + TNumNodes)) { // It belongs to an intersection
            // We compute the indexes
            ConvertAuxHashIndex(static_cast<IndexType>(rTheseBelongs[i_triangle]), belong_index_slave_start, belong_index_slave_end, belong_index_master_start, belong_index_master_end);

            // The coordinates should be in the projected plane
            noalias(xs1) = GeometricalProjectionUtilities::FastProject(slave_center, rSlaveGeometry[belong_index_slave_start], rNormal, distance).Coordinates(); // Start coordinates of the first segment
            noalias(xe1) = GeometricalProjectionUtilities::FastProject(slave_center, rSlaveGeometry[belong_index_slave_end], rNormal, distance).Coordinates(); // End coordinates of the first segment
            noalias(xs2) = GeometricalProjectionUtilities::FastProject(slave_center, rMasterGeometry[belong_index_master_start], rNormal, distance).Coordinates(); // Start coordinates of the second segment
            noalias(xe2) = GeometricalProjectionUtilities::FastProject(slave_center, rMasterGeometry[belong_index_master_end], rNormal, distance).Coordinates(); // End coordinates of the second segment

            // We define the array containing the indexes of the vertexes
            belong_indexes[0] = belong_index_slave_start;
            belong_indexes[1] = belong_index_slave_end;
            belong_indexes[2] = belong_index_master_start + TNumNodes;
            belong_indexes[3] = belong_index_master_end   + TNumNodes;

            // We define the diffs between the extremes of segements
            noalias(diff1) = xs1 - xs2;
            noalias(diff2) = xe2 - xs2;
            noalias(diff3) = xe1 - xs1;

            // We compute the denominator and numerator of the clipping
            MathUtils<double>::CrossProduct(aux_num,   diff1, diff2);
            MathUtils<double>::CrossProduct(aux_denom, diff3, diff2);
            const double num   = inner_prod(aux_num,   rNormal);
            const double denom = inner_prod(aux_denom, rNormal);

            for (IndexType i_belong = 0; i_belong < 4; ++i_belong) {
                // The index of the node
                const IndexType belong_index = belong_indexes[i_belong];

                for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                    // We get the delta normal
                    if ((ConsiderNormalVariation == ELEMENTAL_DERIVATIVES || ConsiderNormalVariation == NODAL_ELEMENTAL_DERIVATIVES) && belong_index < TNumNodes) delta_normal = all_delta_normal[belong_index * TDim + i_dof] * (1.0/aux_nodes_coeff);
                    else delta_normal = zero_array;

                    auto& r_local_delta_vertex = rDerivativeData.DeltaCellVertex[belong_index * TDim + i_dof];

                    // Special cases (slave nodes)
                    if (i_belong == 0) { // First node of the slave
                        const double coeff = 1.0 + num/denom;
                        noalias(row(r_local_delta_vertex, i_triangle)) += LocalDeltaVertex( rNormal,  delta_normal, i_dof, belong_index, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry, coeff);
                    } else if (i_belong == 1) { // Second node of the slave
                        const double coeff = - num/denom;
                        noalias(row(r_local_delta_vertex, i_triangle)) += LocalDeltaVertex( rNormal,  delta_normal, i_dof, belong_index, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry, coeff);
                    }

                    // We define some auxiliar coefficients
                    const double coeff1 = - 1.0/denom;
                    const double coeff2 = num/std::pow(denom, 2);

                    // We add the part corresponding purely to delta normal
                    if (ConsiderNormalVariation == ELEMENTAL_DERIVATIVES || ConsiderNormalVariation == NODAL_ELEMENTAL_DERIVATIVES) {
                        noalias(row(r_local_delta_vertex, i_triangle)) += diff3 * coeff1 * inner_prod(aux_num,  delta_normal);
                        noalias(row(r_local_delta_vertex, i_triangle)) += diff3 * coeff2 * inner_prod(aux_denom, delta_normal);
                    }

                    // We compute the delta diffs
                    noalias(delta_diff1) = (i_belong == 0) ? LocalDeltaVertex(rNormal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry, 1.0) : (i_belong == 2) ? LocalDeltaVertex(rNormal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry, - 1.0) : zero_array;
                    noalias(delta_diff2) = (i_belong == 3) ? LocalDeltaVertex(rNormal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry, 1.0) : (i_belong == 2) ? LocalDeltaVertex(rNormal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry, - 1.0) : zero_array;
                    noalias(delta_diff3) = (i_belong == 1) ? LocalDeltaVertex(rNormal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry, 1.0) : (i_belong == 0) ? LocalDeltaVertex(rNormal, delta_normal, i_dof, belong_index, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry, - 1.0) : zero_array;

                    // We compute now the delta num and denom
                    MathUtils<double>::CrossProduct(aux_cross_product, delta_diff1, diff2);
                    double delta_num = inner_prod(aux_cross_product, rNormal);
                    MathUtils<double>::CrossProduct(aux_cross_product, diff1, delta_diff2);
                    delta_num += inner_prod(aux_cross_product, rNormal);

                    MathUtils<double>::CrossProduct(aux_cross_product, delta_diff3, diff2);
                    double delta_denom = inner_prod(aux_cross_product, rNormal);
                    MathUtils<double>::CrossProduct(aux_cross_product, diff3, delta_diff2);
                    delta_denom += inner_prod(aux_cross_product, rNormal);

                    // Finally we add the contributions of delta num and denom
                    noalias(row(r_local_delta_vertex, i_triangle)) += coeff1 * diff3 * delta_num;
                    noalias(row(r_local_delta_vertex, i_triangle)) += coeff2 * diff3 * delta_denom;
                }
            }
        } else { // It belongs to a master/slave node
            const IndexType belong_index = static_cast<IndexType>(rTheseBelongs[i_triangle]);

            for (unsigned i_dof = 0; i_dof < TDim; ++i_dof) {
                // We get the delta normal
                if ((ConsiderNormalVariation == ELEMENTAL_DERIVATIVES || ConsiderNormalVariation == NODAL_ELEMENTAL_DERIVATIVES) && belong_index < TNumNodes)
                    delta_normal = all_delta_normal[belong_index * TDim + i_dof] * (1.0/aux_nodes_coeff);
                else
                    delta_normal = zero_array;

                auto& r_local_delta_vertex = rDerivativeData.DeltaCellVertex[belong_index * TDim + i_dof];
                noalias(row(r_local_delta_vertex, i_triangle)) += LocalDeltaVertex( rNormal,  delta_normal, i_dof, belong_index, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaN1(
    const GeneralVariables& rVariables,
    DerivativeDataType& rDerivativeData,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    const array_1d<double, 3> rSlaveNormal,
    const DecompositionType& rDecompGeom,
    const PointType& rLocalPointDecomp,
    const PointType& rLocalPointParent,
    const NormalDerivativesComputation ConsiderNormalVariation
    )
{
    KRATOS_TRY

    // Auxiliar zero array
    const array_1d<double, 3> zero_array = ZeroVector(3);

    /* Shape functions */
    const Vector& r_N1 = rVariables.NSlave;

    /* Local gradients */
    const Matrix& r_DNDe1 = rVariables.DNDeSlave;

    // The Normal and delta Normal in the center of the element
    const array_1d<array_1d<double, 3>, TDim * TNumNodes> all_delta_normal = DeltaNormalCenter(rSlaveGeometry);

    /* Shape function decomposition */
    Vector N_decomp;
    rDecompGeom.ShapeFunctionsValues( N_decomp, rLocalPointDecomp.Coordinates() );

    if (TDim == 3) { // NOTE: This is not used in 2D
        array_1d<double, 3> delta_normal, aux_RHS1, aux_delta_node;
        array_1d<double, 2> aux_delta_coords1;
        for ( IndexType i_node = 0; i_node < (TNumNodesMaster + TNumNodes); ++i_node) {
            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                // We get the delta normal
                noalias(delta_normal) = ((ConsiderNormalVariation == ELEMENTAL_DERIVATIVES || ConsiderNormalVariation == NODAL_ELEMENTAL_DERIVATIVES) && i_node < TNumNodes) ? all_delta_normal[i_node * TDim + i_dof] : zero_array;

                // We compute the residuals
                noalias(aux_RHS1) = ZeroVector(3);

                // The vertex cell contribution
                const auto& r_local_delta_cell = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof];
                for(std::size_t i_belong = 0; i_belong < 3; ++i_belong) {
                    noalias(aux_RHS1) += N_decomp[i_belong] * row(r_local_delta_cell, i_belong);
                }

                // Local contribution
                noalias(aux_delta_node) = LocalDeltaVertex( rSlaveNormal, delta_normal, i_dof, i_node, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry );
                if (i_node < TNumNodes)
                    noalias(aux_RHS1) -= r_N1[i_node] * aux_delta_node;

                // We compute the delta coordinates
                DeltaPointLocalCoordinatesSlave(aux_delta_coords1, aux_RHS1, rVariables.DNDeSlave, rSlaveGeometry, rSlaveNormal);

                // Now we can compute the delta shape functions
                auto& r_delta_n1 = rDerivativeData.DeltaN1[i_node * TDim + i_dof];
                noalias(r_delta_n1) = (aux_delta_coords1[0] * column(r_DNDe1, 0) + aux_delta_coords1[1] * column(r_DNDe1, 1));
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaN(
    const GeneralVariables& rVariables,
    DerivativeDataType& rDerivativeData,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    const array_1d<double, 3>& rMasterNormal,
    const DecompositionType& rDecompGeom,
    const PointType& rLocalPointDecomp,
    const PointType& rLocalPointParent,
    const NormalDerivativesComputation ConsiderNormalVariation,
    const bool DualLM
    )
{
    KRATOS_TRY

    // Auxiliar zero array
    const array_1d<double, 3> zero_array = ZeroVector(3);

    /* Shape functions */
    const Vector& r_N1 = rVariables.NSlave;
    const Vector& r_N2 = rVariables.NMaster;

    /* Local gradients */
    const Matrix& r_DNDe1 = rVariables.DNDeSlave;
    const Matrix& r_DNDe2 = rVariables.DNDeMaster;

    // The Normal and delta Normal in the center of the element
    const array_1d<array_1d<double, 3>, TDim * TNumNodes> all_delta_normal = DeltaNormalCenter(rSlaveGeometry);

    /* Shape function decomposition */
    Vector N_decomp;
    rDecompGeom.ShapeFunctionsValues( N_decomp, rLocalPointDecomp.Coordinates() );

    if (TDim == 2) {
        array_1d<PointType, TNumNodes> projected_in_slave, projected_in_master;

        for (IndexType i_mortar_node = 0; i_mortar_node < TNumNodes; ++i_mortar_node) {
            // Projecting points in opposite geometry, defining mortar nodes
            GeometricalProjectionUtilities::FastProjectDirection( rSlaveGeometry,  rMasterGeometry[i_mortar_node], projected_in_slave[i_mortar_node],  rSlaveNormal, rMasterNormal );
            GeometricalProjectionUtilities::FastProjectDirection( rMasterGeometry, rSlaveGeometry[i_mortar_node],  projected_in_master[i_mortar_node], rMasterNormal, rSlaveNormal );
        }

        // Auxiliary point for local coordinates
        PointType aux_point_slave, aux_point_master;

        // Delta local coordinates
        array_1d<double, TNumNodes> DeltaXi_slave, DeltaXi_master;

        for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {

                double DeltaXi1 = 0.0, DeltaXi2 = 0.0;

                for (IndexType i_mortar_node = 0; i_mortar_node < TNumNodes; ++i_mortar_node) {
                    // Compute DeltaXi on the slave side if point is inside geometry (nodes of geometry excluded)
                    if (rSlaveGeometry.IsInside( projected_in_slave[i_mortar_node].Coordinates(), aux_point_slave.Coordinates() ) && (norm_2(projected_in_slave[i_mortar_node].Coordinates() - rSlaveGeometry[0].Coordinates()) > ZeroTolerance) && (norm_2(projected_in_slave[i_mortar_node].Coordinates() - rSlaveGeometry[1].Coordinates()) > ZeroTolerance))
                        DeltaXi_slave[i_mortar_node] = LocalDeltaSegmentN1( all_delta_normal, rSlaveNormal, rSlaveGeometry, rMasterGeometry, r_N1, r_DNDe1, i_mortar_node, i_node, i_dof, ConsiderNormalVariation );

                    // Compute DeltaXi on the master side if point is inside geometry (nodes of geometry excluded)
                    if (rMasterGeometry.IsInside( projected_in_master[i_mortar_node].Coordinates(), aux_point_master.Coordinates() ) && (norm_2(projected_in_master[i_mortar_node].Coordinates() - rMasterGeometry[0].Coordinates()) > ZeroTolerance) && (norm_2(projected_in_master[i_mortar_node].Coordinates() - rMasterGeometry[1].Coordinates()) > ZeroTolerance))
                        DeltaXi_master[i_mortar_node] = LocalDeltaSegmentN2( all_delta_normal, rSlaveNormal, rSlaveGeometry, rMasterGeometry, r_N2, r_DNDe2, i_mortar_node, i_node, i_dof, ConsiderNormalVariation );
                }

                // Evaluate DeltaXi1 expression
                DeltaXi1 = inner_prod(N_decomp, DeltaXi_slave);

                // Evaluate DeltaXi2 expression
                DeltaXi2 = inner_prod(N_decomp, DeltaXi_master);

                // Multiply for DNDe for obtaining DeltaN

                auto& r_delta_n1 = rDerivativeData.DeltaN1[i_node * TDim + i_dof];
                noalias(r_delta_n1) = DeltaXi1 * column(r_DNDe1, 0);

                auto& r_delta_n2 = rDerivativeData.DeltaN2[i_node * TDim + i_dof];
                noalias(r_delta_n2) = DeltaXi2 * column(r_DNDe2, 0);

                // The derivatives of the dual shape function
                auto& r_delta_phi = rDerivativeData.DeltaPhi[i_node * TDim + i_dof];
                if (DualLM) {
                    noalias(r_delta_phi) = prod(rDerivativeData.Ae, r_delta_n1);
                    if (i_node >= TNumNodes)
                        noalias(r_delta_phi) += prod(rDerivativeData.DeltaAe[i_node * TDim + i_dof], r_N1);
                } else {
                    noalias(r_delta_phi) = r_delta_n1;
                }
            }
        }
    } else {
        array_1d<double, 3> delta_normal, aux_delta_node, aux_RHS1, aux_RHS2;
        array_1d<double, 2> aux_delta_coords1, aux_delta_coords2;
        for ( IndexType i_node = 0; i_node < (TNumNodesMaster + TNumNodes); ++i_node) {
            for (IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
                // We get the delta normal
                noalias(delta_normal) = ((ConsiderNormalVariation == ELEMENTAL_DERIVATIVES || ConsiderNormalVariation == NODAL_ELEMENTAL_DERIVATIVES) && i_node < TNumNodes) ? all_delta_normal[i_node * TDim + i_dof] : zero_array;

                // We compute the residuals
                noalias(aux_RHS1) = ZeroVector(3);

                // The vertex cell contribution
                const auto& r_local_delta_cell = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof];
                for(std::size_t i_belong = 0; i_belong < 3; ++i_belong) {
                    noalias(aux_RHS1) += N_decomp[i_belong] * row(r_local_delta_cell, i_belong);
                }
                // Copy
                noalias(aux_RHS2) = aux_RHS1;

                // Local contribution
                noalias(aux_delta_node) = LocalDeltaVertex( rSlaveNormal, delta_normal, i_dof, i_node, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry );
                if (i_node < TNumNodes)
                    noalias(aux_RHS1) -= r_N1[i_node] * aux_delta_node;
                else
                    noalias(aux_RHS2) -= r_N2[i_node - TNumNodes] * aux_delta_node;

                // We compute the delta coordinates
                DeltaPointLocalCoordinatesSlave(aux_delta_coords1, aux_RHS1, rVariables.DNDeSlave, rSlaveGeometry, rSlaveNormal);
                DeltaPointLocalCoordinatesMaster(aux_delta_coords2, aux_RHS2, rVariables.DNDeMaster, rMasterGeometry, rSlaveNormal);

                // Now we can compute the delta shape functions
                auto& r_delta_n1 = rDerivativeData.DeltaN1[i_node * TDim + i_dof];
                noalias(r_delta_n1) = (aux_delta_coords1[0] * column(r_DNDe1, 0) + aux_delta_coords1[1] * column(r_DNDe1, 1));

                auto& r_delta_n2 = rDerivativeData.DeltaN2[i_node * TDim + i_dof];
                noalias(r_delta_n2) = (aux_delta_coords2[0] * column(r_DNDe2, 0) + aux_delta_coords2[1] * column(r_DNDe2, 1));

                // The derivatives of the dual shape function
                auto& r_delta_phi = rDerivativeData.DeltaPhi[i_node * TDim + i_dof];
                if (DualLM) {
                    noalias(r_delta_phi) = prod(rDerivativeData.DeltaAe[i_node * TDim + i_dof], r_N1);
                    noalias(r_delta_phi) += prod(rDerivativeData.Ae, r_delta_n1);
                } else {
                    noalias(r_delta_phi) = r_delta_n1;
                }
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
Matrix& DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaPosition(
    Matrix& rDeltaPosition,
    const GeometryType& rThisGeometry,
    const ConditionArrayType& rLocalCoordinates
    )
{
    KRATOS_TRY;

    // We check size
    if (rDeltaPosition.size1() != TDim || rDeltaPosition.size2() != TDim)
        rDeltaPosition.resize(TDim, TDim, false);

    // We initialize the values
    noalias(rDeltaPosition) = ZeroMatrix(TDim, TDim);

    Vector N;
    for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node ) {
        const array_1d<double, 3 > delta_displacement = rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - rThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( IndexType j_node = 0; j_node < TDim; ++j_node ) {
            rThisGeometry.ShapeFunctionsValues( N, rLocalCoordinates[j_node].Coordinates() );

            for ( IndexType j_dim = 0; j_dim < TDim; ++j_dim )
                rDeltaPosition(j_node, j_dim) += N[i_node] * delta_displacement[j_dim];
        }
    }

    return rDeltaPosition;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline Matrix& DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaPosition(
    Matrix& rDeltaPosition,
    const GeometryType& ThisGeometry
    )
{
    KRATOS_TRY;

    // We check size
    if (rDeltaPosition.size1() != TNumNodes || rDeltaPosition.size2() != TDim)
        rDeltaPosition.resize(TNumNodes, TDim, false);

    // We initialize the values
    noalias(rDeltaPosition) = ZeroMatrix(TNumNodes, TDim);

    for ( IndexType i_node = 0; i_node < TNumNodes; ++i_node ) {
        const array_1d<double, 3 > delta_displacement = ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( IndexType i_dim = 0; i_dim < TDim; ++i_dim )
            rDeltaPosition(i_node, i_dim) += delta_displacement[i_dim];
    }

    return rDeltaPosition;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaPosition(
    Vector& rDeltaPosition,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    const IndexType IndexNode
    )
{
    KRATOS_TRY;

    if (IndexNode < TNumNodes) {
        rDeltaPosition = rSlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - rSlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1);
    } else {
        const IndexType index_master = IndexNode - TNumNodes;
        rDeltaPosition = rMasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - rMasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaPosition(
    Vector& rDeltaPosition,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    const IndexType IndexNode,
    const IndexType iDoF
    )
{
    KRATOS_TRY;

    // We check size
    if (rDeltaPosition.size() != 3)
        rDeltaPosition.resize(3, false);

    // We initialize the values
    noalias(rDeltaPosition) = ZeroVector(3);

    if (IndexNode < TNumNodes) {
        rDeltaPosition[iDoF] = (rSlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - rSlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
    } else {
        const IndexType index_master = IndexNode - TNumNodes;
        rDeltaPosition[iDoF] = (rMasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - rMasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateDeltaPosition(
    double& rDeltaPosition,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    const IndexType IndexNode,
    const IndexType iDoF
    )
{
    KRATOS_TRY;

    if (IndexNode < TNumNodes) {
        rDeltaPosition = (rSlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - rSlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
    } else {
        const IndexType index_master = IndexNode - TNumNodes;
        rDeltaPosition = (rMasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - rMasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
bool DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateAeAndDeltaAe(
    const GeometryType& rSlaveGeometry,
    const array_1d<double, 3>& rSlaveNormal,
    const GeometryType& rMasterGeometry,
    DerivativeDataType& rDerivativeData,
    GeneralVariables& rVariables,
    const NormalDerivativesComputation ConsiderNormalVariation,
    ConditionArrayListType& rConditionsPointsSlave,
    GeometryData::IntegrationMethod ThisIntegrationMethod,
    const double AxiSymCoeff
    )
{
    KRATOS_TRY

    // We initilize the Ae components
    AeData ae_data;
    ae_data.Initialize();

    rDerivativeData.InitializeDeltaAeComponents();

    // Initialize general variables for the current master element
    rVariables.Initialize();

    for (IndexType i_geom = 0; i_geom < rConditionsPointsSlave.size(); ++i_geom) {
        PointerVector< PointType > points_array(TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        array_1d<BelongType, TDim> belong_array;
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            rSlaveGeometry.GlobalCoordinates(global_point, rConditionsPointsSlave[i_geom][i_node]);
            points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
            belong_array[i_node] = rConditionsPointsSlave[i_geom][i_node].GetBelong();
        }

        DecompositionType decomp_geom( points_array );

        const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, rSlaveGeometry.Length() * 1.0e-12) : MortarUtilities::HeronCheck(decomp_geom);

        if (!bad_shape) {
            const GeometryType::IntegrationPointsArrayType integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );

            // Integrating the mortar operators
            array_1d<double,3> gp_normal;
            PointType local_point_parent, gp_global, projected_gp_global;
            GeometryType::CoordinatesArrayType slave_gp_global, projected_gp_local;
            for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                // We reset the derivatives
                rDerivativeData.ResetDerivatives();

                // We compute the local coordinates
                const PointType local_point_decomp = PointType{integration_points_slave[point_number].Coordinates()};

                // Calculate the kinematic variables
                // We compute the current configuration
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                rSlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);

                // SLAVE KINEMATIC COMPUTATIONS
                rSlaveGeometry.ShapeFunctionsValues( rVariables.NSlave, local_point_parent.Coordinates() );
                rSlaveGeometry.ShapeFunctionsLocalGradients( rVariables.DNDeSlave, local_point_parent );
                rVariables.PhiLagrangeMultipliers = rVariables.NSlave;

                rVariables.jSlave = decomp_geom.Jacobian( rVariables.jSlave, local_point_decomp.Coordinates());
                rVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );

                // MASTER KINEMATIC COMPUTATIONS
                noalias(gp_normal) = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, rSlaveGeometry);

                rSlaveGeometry.GlobalCoordinates( slave_gp_global, local_point_parent );
                GeometricalProjectionUtilities::FastProjectDirection( rMasterGeometry, PointType{slave_gp_global}, projected_gp_global, rSlaveNormal, -gp_normal ); // The opposite direction

                rMasterGeometry.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                rMasterGeometry.ShapeFunctionsValues( rVariables.NMaster,    projected_gp_local );
                rMasterGeometry.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );
                rVariables.jMaster = rMasterGeometry.Jacobian( rVariables.jMaster, projected_gp_local);

                // Update the derivative of the integration vertex (just in 3D)
                if (TDim == 3) CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, ConsiderNormalVariation, rSlaveGeometry, rMasterGeometry, rSlaveNormal);

                // Update the derivative of DetJ
                CalculateDeltaDetjSlave(decomp_geom, rVariables, rDerivativeData);

                // Update the derivatives of the shape functions and the gap
                CalculateDeltaN1(rVariables, rDerivativeData, rSlaveGeometry, rMasterGeometry, rSlaveNormal, decomp_geom, local_point_decomp, local_point_parent, ConsiderNormalVariation);

                // Integrate
                const double integration_weight = AxiSymCoeff * integration_points_slave[point_number].Weight();

                ae_data.CalculateDeltaAeComponents(rVariables, rDerivativeData, integration_weight);
            }
        }
    }

    return ae_data.CalculateDeltaAe(rDerivativeData);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline array_1d<double, 3> DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::LocalDeltaVertex(
    const array_1d<double, 3>& rNormal,
    const array_1d<double, 3>& rDeltaNormal,
    const IndexType iDoF,
    const IndexType iBelong,
    const NormalDerivativesComputation ConsiderNormalVariation,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    double Coeff
    )
{
    KRATOS_TRY

    // We create the auxiliar array
    array_1d<double, 3> aux_delta_vertex = ZeroVector(3);

    // This is the coefficient of the center contribution
    const double auxiliar_coeff = 1.0/static_cast<double>(TNumNodes);

    //  We initialize some values
    const array_1d<double, 3> coords_center = rSlaveGeometry.Center().Coordinates();
    const array_1d<double, 3>& r_coords_node = (iBelong < TNumNodes) ? rSlaveGeometry[iBelong].Coordinates() : rMasterGeometry[iBelong - TNumNodes].Coordinates();

    // The corresponding part to the nodal coordinates
    array_1d<double, 3> aux_der = ZeroVector(3);
    aux_der[iDoF] = 1.0;
    aux_delta_vertex += aux_der;

    // The corresponding part to the normal
    const double coordsxdeltanormal = (ConsiderNormalVariation == ELEMENTAL_DERIVATIVES || ConsiderNormalVariation == NODAL_ELEMENTAL_DERIVATIVES) ? inner_prod(r_coords_node - coords_center, rDeltaNormal) : 0.0;

    const double factor_belong = (iBelong < TNumNodes) ? (1.0 - auxiliar_coeff) : 1.0;
    const double deltacoordsxnormal =  factor_belong * rNormal[iDoF];
    aux_delta_vertex += - rNormal * (deltacoordsxnormal + coordsxdeltanormal);

    // The corresponding part to delta normal
    const double coordsxnormal = - inner_prod(r_coords_node - coords_center, rNormal);
    if (ConsiderNormalVariation == ELEMENTAL_DERIVATIVES || ConsiderNormalVariation == NODAL_ELEMENTAL_DERIVATIVES)
        aux_delta_vertex += coordsxnormal * rDeltaNormal;

    return Coeff * aux_delta_vertex;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline BoundedMatrix<double, 3, 3> DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ComputeRenormalizerMatrix(
    const array_1d<double, 3>& rDiffVector,
    const array_1d<double, 3>& rDeltaNormal
    )
{
    KRATOS_TRY

    BoundedMatrix<double, 3, 3> aux_matrix;
    for (IndexType itry = 0; itry < 3; ++itry) {
        if (rDeltaNormal[itry] > ZeroTolerance) {

            const IndexType aux_index_1 = itry == 2 ? 0 : itry + 1;
            const IndexType aux_index_2 = itry == 2 ? 1 : (itry == 1 ? 0 : 2);

            const double diff = rDeltaNormal[aux_index_1] + rDeltaNormal[aux_index_2];
            const double coeff = rDeltaNormal[itry];

            aux_matrix(0, aux_index_1) = 1.0;
            aux_matrix(0, aux_index_2) = 1.0;
            aux_matrix(1, aux_index_1) = 1.0;
            aux_matrix(1, aux_index_2) = 1.0;
            aux_matrix(2, aux_index_1) = 1.0;
            aux_matrix(2, aux_index_2) = 1.0;

            aux_matrix(0, itry) = (rDiffVector[0] - diff)/coeff;
            aux_matrix(1, itry) = (rDiffVector[1] - diff)/coeff;
            aux_matrix(2, itry) = (rDiffVector[2] - diff)/coeff;

            return aux_matrix;
        }
    }

    return IdentityMatrix(3, 3);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline BoundedMatrix<double, 3, 3> DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ComputeRenormalizerMatrix(
    const BoundedMatrix<double, TNumNodes, TDim>& rDiffMatrix,
    const BoundedMatrix<double, TNumNodes, TDim>& rDeltaNormal,
    const IndexType iGeometry
    )
{
    KRATOS_TRY

    BoundedMatrix<double, 3, 3> aux_matrix;
    for (IndexType itry = 0; itry < 3; ++itry) {
        if (rDeltaNormal(iGeometry, itry) > ZeroTolerance) {

            const IndexType aux_index_1 = itry == 2 ? 0 : itry + 1;
            const IndexType aux_index_2 = itry == 2 ? 1 : (itry == 1 ? 0 : 2);

            const double diff = rDeltaNormal(iGeometry, aux_index_1) + rDeltaNormal(iGeometry, aux_index_2);
            const double coeff = rDeltaNormal(iGeometry, itry);

            aux_matrix(0, aux_index_1) = 1.0;
            aux_matrix(0, aux_index_2) = 1.0;
            aux_matrix(1, aux_index_1) = 1.0;
            aux_matrix(1, aux_index_2) = 1.0;
            aux_matrix(2, aux_index_1) = 1.0;
            aux_matrix(2, aux_index_2) = 1.0;

            aux_matrix(0, itry) = (rDiffMatrix(iGeometry, 0) - diff)/coeff;
            aux_matrix(1, itry) = (rDiffMatrix(iGeometry, 1) - diff)/coeff;
            aux_matrix(2, itry) = (rDiffMatrix(iGeometry, 2) - diff)/coeff;

            return aux_matrix;
        }
    }

    return IdentityMatrix(3, 3);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline array_1d<double, 3> DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::PreviousNormalGeometry(
    const GeometryType& rThisGeometry,
    const GeometryType::CoordinatesArrayType& rPointLocal
    )
{
    KRATOS_TRY

    // We compute the previous normal in the geometry
    Matrix previous_jacobian, delta_position;
    delta_position = CalculateDeltaPosition(delta_position, rThisGeometry);
    previous_jacobian = rThisGeometry.Jacobian( previous_jacobian, rPointLocal, delta_position);

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
    KRATOS_ERROR_IF(norm_normal < ZeroTolerance) << "ERROR: The normal norm is zero or almost zero. Norm. normal: " << norm_normal << std::endl;

    return previous_normal;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ConvertAuxHashIndex(
    const IndexType AuxIndex,
    IndexType& riBelongSlaveStart,
    IndexType& riBelongSlaveEnd,
    IndexType& riBelongMasterStart,
    IndexType& riBelongMasterEnd
    )
{
    KRATOS_TRY

    IndexType index_to_decompose = AuxIndex - (TNumNodesMaster + TNumNodes);

    riBelongMasterEnd = index_to_decompose/10000;
    index_to_decompose = std::fmod(index_to_decompose, 10000);
    riBelongMasterStart = index_to_decompose/1000;
    index_to_decompose = std::fmod(index_to_decompose, 1000);
    riBelongSlaveEnd = index_to_decompose/100;
    index_to_decompose = std::fmod(index_to_decompose, 100);
    riBelongSlaveStart = index_to_decompose/10;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::DeltaPointLocalCoordinatesSlave(
    array_1d<double, 2>& rResult,
    const array_1d<double, 3>& rDeltaPoint,
    const Matrix& rDNDe,
    const GeometryType& rThisGeometry,
    const array_1d<double, 3>& rThisNormal
    )
{
    KRATOS_TRY

    BoundedMatrix<double, 3, TNumNodes> X;
    for(IndexType i = 0; i < TNumNodes; ++i) {
        X(0, i) = rThisGeometry[i].X();
        X(1, i) = rThisGeometry[i].Y();
        X(2, i) = rThisGeometry[i].Z();
    }

    const BoundedMatrix<double, 3, 2> DN = prod(X, rDNDe);

    const BoundedMatrix<double, 2, 2> J = prod(trans(DN),DN);
    double det_j;
    BoundedMatrix<double, 2, 2> invJ;
    MathUtils<double>::InvertMatrix(J, invJ, det_j, -1.0);
    const bool good_condition_number = MathUtils<double>::CheckConditionNumber(J, invJ, std::numeric_limits<double>::epsilon(), false);
    if (!good_condition_number) // Reset in case of bad condition number
        noalias(invJ) = ZeroMatrix(2,2);

#ifdef KRATOS_DEBUG
    KRATOS_WARNING_IF("Jacobian invert", !good_condition_number) << "WARNING:: CANNOT INVERT JACOBIAN TO COMPUTE DELTA COORDINATES" << std::endl;
#endif

    const array_1d<double, 2> res = prod(trans(DN), rDeltaPoint);
    noalias(rResult) = prod(invJ, res);

//     BoundedMatrix<double, 3, 3> L;
//     for(IndexType i = 0; i < 3; ++i)
//     {
//         for(IndexType j = 0; j < 2; ++j) L(i, j) = DN(i, j);
//         L(i, 2) = ThisNormal[i];
//     }
//
//     double det_L;
//     BoundedMatrix<double, 3, 3> invL;
//     MathUtils<double>::InvertMatrix(L, invL, det_L, -1.0);
//     const bool good_condition_number = MathUtils<double>::CheckConditionNumber(L, invL, std::numeric_limits<double>::epsilon(), false);
//     if (!good_condition_number) // Reset in case of bad condition number
//         noalias(invL) = ZeroMatrix(3,3);
//     array_1d<double, 3> aux = prod(invL, DeltaPoint);
//     rResult[0] = aux[0];
//     rResult[1] = aux[1];
//     #ifdef KRATOS_DEBUG
//         KRATOS_WARNING_IF("Jacobian invert", !good_condition_number) << "WARNING: CANNOT INVERT JACOBIAN TO COMPUTE DELTA COORDINATES" << std::endl;
//     #endif

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline void DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::DeltaPointLocalCoordinatesMaster(
    array_1d<double, 2>& rResult,
    const array_1d<double, 3>& rDeltaPoint,
    const Matrix& rDNDe,
    const GeometryType& rThisGeometry,
    const array_1d<double, 3>& rThisNormal
    )
{
    KRATOS_TRY

    BoundedMatrix<double, 3, TNumNodesMaster> X;
    for(IndexType i = 0; i < TNumNodesMaster; ++i) {
        X(0, i) = rThisGeometry[i].X();
        X(1, i) = rThisGeometry[i].Y();
        X(2, i) = rThisGeometry[i].Z();
    }

    const BoundedMatrix<double, 3, 2> DN = prod(X, rDNDe);

    const BoundedMatrix<double, 2, 2> J = prod(trans(DN),DN);
    double det_j;
    BoundedMatrix<double, 2, 2> invJ;
    MathUtils<double>::InvertMatrix(J, invJ, det_j, -1.0);
    const bool good_condition_number = MathUtils<double>::CheckConditionNumber(J, invJ, std::numeric_limits<double>::epsilon(), false);
    if (!good_condition_number) // Reset in case of bad condition number
        noalias(invJ) = ZeroMatrix(2,2);

#ifdef KRATOS_DEBUG
    KRATOS_WARNING_IF("Jacobian invert", !good_condition_number) << "WARNING:: CANNOT INVERT JACOBIAN TO COMPUTE DELTA COORDINATES" << std::endl;
#endif

    const array_1d<double, 2> res = prod(trans(DN), rDeltaPoint);
    noalias(rResult) = prod(invJ, res);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline double DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::LocalDeltaSegmentN1(
    const array_1d<array_1d<double, 3>, TDim * TNumNodes>& rDeltaNormal,
    const array_1d<double, 3>& rSlaveNormal,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    const Vector& rN1,
    const Matrix& rDNDe1,
    const IndexType MortarNode,
    const IndexType iNode,
    const IndexType iDoF,
    const NormalDerivativesComputation ConsiderNormalVariation
    )
{
    KRATOS_TRY

    array_1d<double, TDim> Xa, DXa;
    BoundedMatrix<double, TDim, TNumNodes> X1, DX1;
    BoundedMatrix<double, 3, TNumNodes> n1, Dn1;

    // Projected node coordinates
    Xa[0] = rMasterGeometry[MortarNode].X();
    Xa[1] = rMasterGeometry[MortarNode].Y();

    // Projected node derivatives
    DXa = ZeroVector(TDim);
    DXa[iDoF] = 1.0;

    // Slave element nodes coordinates and normals with respective derivatives
    for(IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        X1(0, i_node) = rSlaveGeometry[i_node].X();
        X1(1, i_node) = rSlaveGeometry[i_node].Y();

        for(IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
            if(i_dof == iDoF)
                DX1(i_dof, i_node) = 1.0;
            else
                DX1(i_dof, i_node) = 0.0;
        }

        column(n1, i_node)  = rSlaveNormal;
        column(Dn1, i_node) = (ConsiderNormalVariation == ELEMENTAL_DERIVATIVES || ConsiderNormalVariation == NODAL_ELEMENTAL_DERIVATIVES) ? rDeltaNormal[i_node * TDim + iDoF] : ZeroVector(3);
    }

    // Computation of DeltaXi_a
    const double denom = -inner_prod(row(X1, 0), column(rDNDe1, 0))*(rN1[iNode]*n1(1,iNode)) + inner_prod(row(X1, 1), column(rDNDe1,0))*(rN1[iNode]*n1(0,iNode)) - (rN1[iNode]*X1(0,iNode) - Xa[0])*inner_prod(row(n1, 1), column(rDNDe1, 0)) + (rN1[iNode]*X1(1,iNode) - Xa[1])*inner_prod(row(n1, 0), column(rDNDe1, 0));
    const double num = (rN1[iNode]*DX1(0,iNode) - DXa[0])*(rN1[iNode]*n1(1,iNode)) - (rN1[iNode]*DX1(1,iNode) - DXa[1])*(rN1[iNode]*n1(0,iNode)) + (rN1[iNode]*X1(0,iNode) - Xa[0])*(rN1[iNode]*Dn1(1,iNode)) - (rN1[iNode]*X1(1,iNode) - Xa[1])*(rN1[iNode]*Dn1(0,iNode));

    return num/denom;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<const SizeType TDim, const SizeType TNumNodes, bool TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
inline double DerivativesUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::LocalDeltaSegmentN2(
    const array_1d<array_1d<double, 3>, TDim * TNumNodes>& rDeltaNormal,
    const array_1d<double, 3>& rSlaveNormal,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    const Vector& rN2,
    const Matrix& rDNDe2,
    const IndexType MortarNode,
    const IndexType iNode,
    const IndexType iDoF,
    const NormalDerivativesComputation ConsiderNormalVariation
    )
{
    KRATOS_TRY

    array_1d<double, TDim> Xa, DXa;
    array_1d<double, 3>    na, Dna;
    BoundedMatrix<double, TDim, TNumNodesMaster> X2, DX2;

    // Projected node coordinates
    Xa[0] = rSlaveGeometry[MortarNode].X();
    Xa[1] = rSlaveGeometry[MortarNode].Y();

    // Projected node derivatives
    DXa = ZeroVector(TDim);
    DXa[iDoF] = 1.0;

    // Projected normal and derivative
    na = rSlaveNormal;
    Dna = (ConsiderNormalVariation == ELEMENTAL_DERIVATIVES || ConsiderNormalVariation == NODAL_ELEMENTAL_DERIVATIVES) ? rDeltaNormal[MortarNode * TDim + iDoF]: ZeroVector(3);

    // Slave element nodes coordinates and derivatives
    for(IndexType i_node = 0; i_node < TNumNodesMaster; ++i_node) {
        X2(0, i_node) = rMasterGeometry[i_node].X();
        X2(1, i_node) = rMasterGeometry[i_node].Y();

        for(IndexType i_dof = 0; i_dof < TDim; ++i_dof) {
        if(i_dof == iDoF)
            DX2(i_dof, i_node) = 1.0;
        else
            DX2(i_dof, i_node) = 0.0;
        }
    }

    // Computation of DeltaXi_a
    const double lhs = -1.0/(inner_prod(row(X2, 0), column(rDNDe2,0))*na[1] - inner_prod(row(X2, 1), column(rDNDe2,0))*na[0]);
    const double rhs = (rN2[iNode]*DX2(0, iNode)-DXa[0])*na[1] - (rN2[iNode]*DX2(1, iNode)-DXa[1])*na[0] + (rN2[iNode]*X2(0, iNode)-Xa[0])*Dna[1] - (rN2[iNode]*X2(1, iNode)-Xa[1])*Dna[0];

    return lhs*rhs;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

// Frictionless cases
template class DerivativesUtilities<2, 2, false, false, 2>;
template class DerivativesUtilities<3, 3, false, false, 3>;
template class DerivativesUtilities<3, 4, false, false, 4>;
template class DerivativesUtilities<3, 3, false, false, 4>;
template class DerivativesUtilities<3, 4, false, false, 3>;
template class DerivativesUtilities<2, 2, false, true, 2>;
template class DerivativesUtilities<3, 3, false, true, 3>;
template class DerivativesUtilities<3, 4, false, true, 4>;
template class DerivativesUtilities<3, 3, false, true, 4>;
template class DerivativesUtilities<3, 4, false, true, 3>;

// Frictional cases
template class DerivativesUtilities<2, 2, true, false, 2>;
template class DerivativesUtilities<3, 3, true, false, 3>;
template class DerivativesUtilities<3, 4, true, false, 4>;
template class DerivativesUtilities<3, 3, true, false, 4>;
template class DerivativesUtilities<3, 4, true, false, 3>;
template class DerivativesUtilities<2, 2, true, true, 2>;
template class DerivativesUtilities<3, 3, true, true, 3>;
template class DerivativesUtilities<3, 4, true, true, 4>;
template class DerivativesUtilities<3, 3, true, true, 4>;
template class DerivativesUtilities<3, 4, true, true, 3>;

} // namespace Kratos
