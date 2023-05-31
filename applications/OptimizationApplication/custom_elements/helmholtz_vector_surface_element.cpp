//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/variables.h"

// Application incldues
#include "custom_utilities/entity_calculation_utils.h"
#include "optimization_application_variables.h"

// Include base h
#include "custom_elements/helmholtz_vector_surface_element.h"

namespace Kratos {

//************************************************************************************
//************************************************************************************

HelmholtzVectorSurfaceElement::HelmholtzVectorSurfaceElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
    mpSolidGeometry = EntityCalculationUtils::CreateSolidGeometry(this->GetGeometry());
}

//************************************************************************************
//************************************************************************************

HelmholtzVectorSurfaceElement::HelmholtzVectorSurfaceElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    // DO NOT ADD DOFS HERE!!!
    mpSolidGeometry = EntityCalculationUtils::CreateSolidGeometry(this->GetGeometry());
}

Element::Pointer HelmholtzVectorSurfaceElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzVectorSurfaceElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer HelmholtzVectorSurfaceElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzVectorSurfaceElement>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer HelmholtzVectorSurfaceElement::Clone(
    IndexType NewId,
    const NodesArrayType& rThisNodes) const
{
    KRATOS_TRY

    HelmholtzVectorSurfaceElement::Pointer p_new_cond =
        Kratos::make_intrusive<HelmholtzVectorSurfaceElement>(
            NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));

    return p_new_cond;

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
void HelmholtzVectorSurfaceElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); // resetting LHS

    // Resizing as needed the RHS
    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size, false);

    rRightHandSideVector = ZeroVector(mat_size); // resetting RHS

    MatrixType M;
    CalculateMassMatrix(M, rCurrentProcessInfo);

    MatrixType K;
    CalculateStiffnessMatrix(K, rCurrentProcessInfo);

    const bool is_inversed = rCurrentProcessInfo[COMPUTE_HELMHOLTZ_INVERSE];

    noalias(rLeftHandSideMatrix) += M;
    if (!is_inversed)
        noalias(rLeftHandSideMatrix) += K;

    const unsigned int number_of_points = r_geometry.size();
    Vector nodal_vals(number_of_points * 3);
    const bool is_integrated_field = rCurrentProcessInfo[HELMHOLTZ_INTEGRATED_FIELD];
    for (unsigned int node_element = 0; node_element < number_of_points; node_element++) {
        const auto& source = r_geometry[node_element].GetValue(HELMHOLTZ_VECTOR_SOURCE);
        auto node_weight = r_geometry[node_element].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        nodal_vals[3 * node_element + 0] = source[0];
        nodal_vals[3 * node_element + 1] = source[1];
        nodal_vals[3 * node_element + 2] = source[2];
        if (is_integrated_field) {
            nodal_vals[3 * node_element + 0] /= node_weight;
            nodal_vals[3 * node_element + 1] /= node_weight;
            nodal_vals[3 * node_element + 2] /= node_weight;
        }
    }

    if (is_integrated_field) {
        noalias(rRightHandSideVector) += nodal_vals;
    } else if (is_inversed) {
        noalias(rRightHandSideVector) += prod(K + M, nodal_vals);
    } else {
        noalias(rRightHandSideVector) += prod(M, nodal_vals);
    }

    // apply drichlet BC
    Vector temp;
    GetValuesVector(temp, 0);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp);

    KRATOS_CATCH("")
}

//******************************************************************************
//******************************************************************************
void HelmholtzVectorSurfaceElement::GetValuesVector(VectorType& rValues, int Step) const
{
    const GeometryType& rgeom = this->GetGeometry();
    const SizeType num_nodes = rgeom.PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int local_size = num_nodes * dimension;

    if (rValues.size() != local_size)
        rValues.resize(local_size, false);

    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        rValues[index++] =
            rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VECTOR_X, Step);
        rValues[index++] =
            rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VECTOR_Y, Step);
        rValues[index++] =
            rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VECTOR_Z, Step);
    }
}
//******************************************************************************
//******************************************************************************

void HelmholtzVectorSurfaceElement::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == ELEMENT_STRAIN_ENERGY) {
        MatrixType K;
        CalculateStiffnessMatrix(K, rCurrentProcessInfo);

        auto& r_geometry = this->GetGeometry();

        const unsigned int number_of_points = r_geometry.size();
        Vector nodal_vals(number_of_points * 3);
        for (unsigned int node_element = 0; node_element < number_of_points; node_element++) {
            nodal_vals[3 * node_element + 0] = r_geometry[node_element].X0();
            nodal_vals[3 * node_element + 1] = r_geometry[node_element].Y0();
            nodal_vals[3 * node_element + 2] = r_geometry[node_element].Z0();
        }

        rOutput = inner_prod(nodal_vals, prod(K, nodal_vals));
    }
    else {
        auto& parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
        parentElement[0].Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

//************************************************************************************
//************************************************************************************
void HelmholtzVectorSurfaceElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzVectorSurfaceElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0, 0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzVectorSurfaceElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes, false);

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(HELMHOLTZ_VECTOR_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const SizeType index = i * 3;
        rResult[index] = GetGeometry()[i].GetDof(HELMHOLTZ_VECTOR_X, pos).EquationId();
        rResult[index + 1] =
            GetGeometry()[i].GetDof(HELMHOLTZ_VECTOR_Y, pos + 1).EquationId();
        rResult[index + 2] =
            GetGeometry()[i].GetDof(HELMHOLTZ_VECTOR_Z, pos + 2).EquationId();
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void HelmholtzVectorSurfaceElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(HELMHOLTZ_VECTOR_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(HELMHOLTZ_VECTOR_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(HELMHOLTZ_VECTOR_Z));
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
int HelmholtzVectorSurfaceElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (IndexType i = 0; i < number_of_nodes; i++) {
        const NodeType& rnode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HELMHOLTZ_VECTOR, rnode)

        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VECTOR_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VECTOR_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VECTOR_Z, rnode)
    }

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(COMPUTE_HELMHOLTZ_INVERSE))
        << "COMPUTE_HELMHOLTZ_INVERSE not defined in the ProcessInfo!" << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_INTEGRATED_FIELD))
        << "HELMHOLTZ_INTEGRATED_FIELD not defined in the ProcessInfo!" << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_RADIUS))
        << "HELMHOLTZ_RADIUS not defined in the ProcessInfo!" << std::endl;

    return check;

    KRATOS_CATCH("");
}
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVectorSurfaceElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const auto& r_cond_geom = GetGeometry();
    SizeType dimension = r_cond_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_cond_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize(mat_size, mat_size, false);
    rMassMatrix = ZeroMatrix(mat_size, mat_size);

    IntegrationMethod integration_method = r_cond_geom.GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_cond_geom.IntegrationPoints(integration_method);

    VectorType Ws;
    MatrixType Ns;
    EntityCalculationUtils::CalculateSurfaceElementGaussPointData(Ws, Ns, r_cond_geom, integration_method);

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
        const double integration_weight = Ws[point_number];
        const Vector& rN = row(Ns, point_number);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index_i = i * dimension;

            for (IndexType j = 0; j < number_of_nodes; ++j) {
                const SizeType index_j = j * dimension;
                const double NiNj_weight = rN[i] * rN[j] * integration_weight;

                for (IndexType k = 0; k < dimension; ++k) {
                    rMassMatrix(index_i + k, index_j + k) += NiNj_weight;
                }
            }
        }
    }

    KRATOS_CATCH("");
}
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzVectorSurfaceElement::CalculateStiffnessMatrix(
    MatrixType& rStiffnessMatrix,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rStiffnessMatrix.size1() != mat_size || rStiffnessMatrix.size2() != mat_size)
        rStiffnessMatrix.resize(mat_size, mat_size, false);
    rStiffnessMatrix = ZeroMatrix(mat_size, mat_size);

    MatrixType A_dirc = ZeroMatrix(number_of_nodes, number_of_nodes);

    // reading integration points and local gradients
    const IntegrationMethod& integration_method = r_geom.GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geom.IntegrationPoints(integration_method);
    const unsigned int NumGauss = integration_points.size();
    Vector GaussPtsJDet = ZeroVector(NumGauss);
    r_geom.DeterminantOfJacobian(GaussPtsJDet, integration_method);

    VectorType n_surf;
    CalculateAvgSurfUnitNormal(n_surf);
    MatrixType id_matrix = IdentityMatrix(dimension, dimension);
    MatrixType tangent_projection_matrix = id_matrix - outer_prod(n_surf, n_surf);

    for (std::size_t i_point = 0; i_point < integration_points.size(); ++i_point) {

        Matrix DN_DX;
        EntityCalculationUtils::CalculateSurfaceElementShapeDerivatives(DN_DX, *mpSolidGeometry, r_geom, integration_method, i_point);

        const double IntToReferenceWeight =
            integration_points[i_point].Weight() * GaussPtsJDet[i_point];

        MatrixType DN_DX_t = prod(DN_DX, tangent_projection_matrix);

        const double r_helmholtz = rCurrentProcessInfo[HELMHOLTZ_RADIUS];
        noalias(A_dirc) += IntToReferenceWeight * r_helmholtz * r_helmholtz *
                           prod(DN_DX_t, trans(DN_DX_t));
    }

    // contruct the stifness matrix in all dims
    for (IndexType i = 0; i < number_of_nodes; i++) {
        for (IndexType j = 0; j < dimension; j++) {
            for (IndexType k = 0; k < number_of_nodes; k++) {
                rStiffnessMatrix(dimension * i + j, dimension * k + j) = A_dirc(i, k);
            }
        }
    }

    KRATOS_CATCH("");
}

void HelmholtzVectorSurfaceElement::CalculateAvgSurfUnitNormal(VectorType& rNormal) const
{
    const auto& r_geom = GetGeometry();
    const IntegrationMethod& integration_method = r_geom.GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points =
        r_geom.IntegrationPoints(integration_method);

    rNormal.resize(3);
    rNormal = ZeroVector(3);
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
        rNormal += r_geom.UnitNormal(point_number, integration_method);

    rNormal /= integration_points.size();
    rNormal /= MathUtils<double>::Norm3(rNormal);
}

} // Namespace Kratos
