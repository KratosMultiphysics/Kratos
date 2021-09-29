//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "wave_element.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/shallow_water_utilities.h"
#include "custom_friction_laws/friction_laws_factory.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TNumNodes>
int WaveElement<TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int err = Element::Check(rCurrentProcessInfo);
    if (err != 0) return err;

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry())
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TOPOGRAPHY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VERTICAL_VELOCITY, r_node)

        KRATOS_CHECK_DOF_IN_NODE(this->GetUnknownComponent(0), r_node)
        KRATOS_CHECK_DOF_IN_NODE(this->GetUnknownComponent(1), r_node)
        KRATOS_CHECK_DOF_IN_NODE(this->GetUnknownComponent(2), r_node)
    }

    return err;

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rResult.size() != mLocalSize)
        rResult.resize(mLocalSize, false); // False says not to preserve existing storage!!

    const GeometryType& r_geom = this->GetGeometry();
    int counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = r_geom[i].GetDof(this->GetUnknownComponent(0)).EquationId();
        rResult[counter++] = r_geom[i].GetDof(this->GetUnknownComponent(1)).EquationId();
        rResult[counter++] = r_geom[i].GetDof(this->GetUnknownComponent(2)).EquationId();
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rElementalDofList.size() != mLocalSize)
        rElementalDofList.resize(mLocalSize);

    const GeometryType& r_geom = this->GetGeometry();
    int counter=0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[counter++] = r_geom[i].pGetDof(this->GetUnknownComponent(0));
        rElementalDofList[counter++] = r_geom[i].pGetDof(this->GetUnknownComponent(1));
        rElementalDofList[counter++] = r_geom[i].pGetDof(this->GetUnknownComponent(2));
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(this->GetUnknownComponent(0), Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(this->GetUnknownComponent(1), Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(this->GetUnknownComponent(2), Step);
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(VERTICAL_VELOCITY, Step);
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_ERROR << "WaveElement::GetSecondDerivativesVector This method is not supported by the formulation" << std::endl;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    for (IndexType PointNumber = 0; PointNumber < 1; PointNumber++)
        rValues[PointNumber] = double(this->GetValue(rVariable));
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
}

template<std::size_t TNumNodes>
const Variable<double>& WaveElement<TNumNodes>::GetUnknownComponent(int Index) const
{
    switch (Index) {
        case 0: return VELOCITY_X;
        case 1: return VELOCITY_Y;
        case 2: return HEIGHT;
        default: KRATOS_ERROR << "WaveElement::GetUnknownComponent index out of bounds." << std::endl;
    }
}

template<std::size_t TNumNodes>
typename WaveElement<TNumNodes>::LocalVectorType WaveElement<TNumNodes>::GetUnknownVector(ElementData& rData)
{
    std::size_t index = 0;
    array_1d<double,mLocalSize> unknown;
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        unknown[index++] = rData.nodal_v[i][0];
        unknown[index++] = rData.nodal_v[i][1];
        unknown[index++] = rData.nodal_h[i];
    }
    return unknown;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::InitializeData(ElementData& rData, const ProcessInfo& rCurrentProcessInfo)
{
    rData.stab_factor = rCurrentProcessInfo[STABILIZATION_FACTOR];
    rData.shock_stab_factor = rCurrentProcessInfo[SHOCK_STABILIZATION_FACTOR];
    rData.relative_dry_height = rCurrentProcessInfo[RELATIVE_DRY_HEIGHT];
    rData.gravity = rCurrentProcessInfo[GRAVITY_Z];
    rData.p_bottom_friction = FrictionLawsFactory().CreateBottomFrictionLaw(
        this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetNodalData(ElementData& rData, const GeometryType& rGeometry, int Step)
{
    rData.length = rGeometry.Length();

    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rData.nodal_h[i] = rGeometry[i].FastGetSolutionStepValue(HEIGHT, Step);
        rData.nodal_z[i] = rGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY, Step);
        rData.nodal_v[i] = rGeometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        rData.nodal_q[i] = rGeometry[i].FastGetSolutionStepValue(MOMENTUM, Step);
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateGaussPointData(
    ElementData& rData,
    const array_1d<double,TNumNodes>& rN)
{
    const double h = inner_prod(rData.nodal_h, rN);
    const array_1d<double,3> v = VectorProduct(rData.nodal_v, rN);

    rData.height = h;
    rData.velocity = v;

    /**
     * A_1 = {{ 0   0   g },
     *        { 0   0   0 },
     *        { h   0   0 }}
     */
    rData.A1 = ZeroMatrix(3, 3);
    rData.A1(0,2) = rData.gravity;
    rData.A1(2,0) = h;

    /*
     * A_2 = {{ 0   0   0 },
     *        { 0   0   g },
     *        { 0   h   0 }}
     */
    rData.A2 = ZeroMatrix(3, 3);
    rData.A2(1,2) = rData.gravity;
    rData.A2(2,1) = h;

    /// b_1
    rData.b1 = ZeroVector(3);
    rData.b1[0] = rData.gravity;

    /// b_2
    rData.b2 = ZeroVector(3);
    rData.b2[1] = rData.gravity;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateArtificialViscosityData(
    ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX)
{
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateGeometryData(
    Vector &rGaussWeights,
    Matrix &rNContainer,
    ShapeFunctionsGradientsType &rDN_DXContainer) const
{
    Vector det_j_vector;
    const auto& r_geom = this->GetGeometry();
    const auto integration_method = r_geom.GetDefaultIntegrationMethod();
    r_geom.ShapeFunctionsIntegrationPointsGradients(rDN_DXContainer, det_j_vector, integration_method, rNContainer);

    const unsigned int number_of_gauss_points = r_geom.IntegrationPointsNumber(integration_method);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(integration_method);

    if (rGaussWeights.size() != number_of_gauss_points)
        rGaussWeights.resize(number_of_gauss_points, false);

    for (unsigned int g = 0; g < number_of_gauss_points; ++g)
        rGaussWeights[g] = det_j_vector[g] * integration_points[g].Weight();
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::AddWaveTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    const bool integrate_by_parts = false;
    const auto z = rData.nodal_z;
    const double l = StabilizationParameter(rData);

    const BoundedMatrix<double,3,3> AA11 = prod(rData.A1, rData.A1);
    const BoundedMatrix<double,3,3> AA22 = prod(rData.A2, rData.A2);
    const BoundedMatrix<double,3,3> AA12 = prod(rData.A1, rData.A2);
    const BoundedMatrix<double,3,3> AA21 = prod(rData.A2, rData.A1);
    const array_1d<double,3> Ab11 = prod(rData.A1, rData.b1);
    const array_1d<double,3> Ab22 = prod(rData.A2, rData.b2);
    const array_1d<double,3> Ab12 = prod(rData.A1, rData.b2);
    const array_1d<double,3> Ab21 = prod(rData.A2, rData.b1);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            double g1_ij;
            double g2_ij;
            double d_ij;
            if (integrate_by_parts) {
                g1_ij = -rDN_DX(i,0) * rN[j];
                g2_ij = -rDN_DX(i,1) * rN[j];
            } else {
                g1_ij = rN[i] * rDN_DX(j,0);
                g2_ij = rN[i] * rDN_DX(j,1);
            }

            /// First component
            MathUtils<double>::AddMatrix(rMatrix,  Weight*g1_ij*rData.A1, 3*i, 3*j);
            MathUtils<double>::AddVector(rVector, -Weight*g1_ij*rData.b1*z[j], 3*i);

            /// Second component
            MathUtils<double>::AddMatrix(rMatrix,  Weight*g2_ij*rData.A2, 3*i, 3*j);
            MathUtils<double>::AddVector(rVector, -Weight*g2_ij*rData.b2*z[j], 3*i);

            /// Stabilization x-x
            d_ij = rDN_DX(i,0) * rDN_DX(j,0);
            MathUtils<double>::AddMatrix(rMatrix,  Weight*l*d_ij*AA11, 3*i, 3*j);
            MathUtils<double>::AddVector(rVector, -Weight*l*d_ij*Ab11*z[j], 3*i);

            /// Stabilization y-y
            d_ij = rDN_DX(i,1) * rDN_DX(j,1);
            MathUtils<double>::AddMatrix(rMatrix,  Weight*l*d_ij*AA22, 3*i, 3*j);
            MathUtils<double>::AddVector(rVector, -Weight*l*d_ij*Ab22*z[j], 3*i);

            /// Stabilization x-y
            d_ij = rDN_DX(i,0) * rDN_DX(j,1);
            MathUtils<double>::AddMatrix(rMatrix,  Weight*l*d_ij*AA12, 3*i, 3*j);
            MathUtils<double>::AddVector(rVector, -Weight*l*d_ij*Ab12*z[j], 3*i);

            /// Stabilization y-x
            d_ij = rDN_DX(i,1) * rDN_DX(j,0);
            MathUtils<double>::AddMatrix(rMatrix,  Weight*l*d_ij*AA21, 3*i, 3*j);
            MathUtils<double>::AddVector(rVector, -Weight*l*d_ij*Ab21*z[j], 3*i);
        }
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::AddFrictionTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    const double s = rData.p_bottom_friction->CalculateLHS(rData.height, rData.velocity);
    const double lumping_factor = 1.0 / TNumNodes;
    const double l = StabilizationParameter(rData);
    BoundedMatrix<double,3,3> Sf = ZeroMatrix(3,3);
    Sf(0,0) = rData.gravity*s;
    Sf(1,1) = rData.gravity*s;

    BoundedMatrix<double,3,3> ASf1 = prod(rData.A1, Sf);
    BoundedMatrix<double,3,3> ASf2 = prod(rData.A2, Sf);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        MathUtils<double>::AddMatrix(rMatrix, Weight*lumping_factor*Sf, 3*i, 3*i);

        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            /// Stabilization x
            const double g1_ij = rDN_DX(i,0) * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*l*g1_ij*ASf1, 3*i, 3*j);

            /// Stabilization y
            const double g2_ij = rDN_DX(i,1) * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*l*g2_ij*ASf2, 3*i, 3*j);
        }
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::AddArtificialViscosityTerms(
    LocalMatrixType& rMatrix,
    const ElementData& rData,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{}

template<>
void WaveElement<3>::AddMassTerms(
    LocalMatrixType& rMatrix,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX,
    const double Weight)
{
    const double l = StabilizationParameter(rData);
    BoundedMatrix<double, 3, 3> M = IdentityMatrix(3);

    // Algebraic factor
    const double one_twelfth = 1.0 / 12.0;
    const double one_sixth = 1.0 / 6.0;

    for (IndexType i = 0; i < 3; ++i)
    {
        for (IndexType j = 0; j < 3; ++j)
        {
            // Consistent mass matrix
            const double n_ij = (i == j) ? one_sixth : one_twelfth;
            MathUtils<double>::AddMatrix(rMatrix, Weight*n_ij*M, 3*i, 3*j);

            /// Stabilization x
            const double g1_ij = rDN_DX(i,0) * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*l*g1_ij*rData.A1, 3*i, 3*j);

            /// Stabilization y
            const double g2_ij = rDN_DX(i,1) * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*l*g2_ij*rData.A2, 3*i, 3*j);
        }
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::AddMassTerms(
    LocalMatrixType& rMatrix,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    const double l = StabilizationParameter(rData);
    BoundedMatrix<double, 3, 3> M = IdentityMatrix(3);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            // Consistent mass matrix
            const double n_ij = rN[i] * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*n_ij*M, 3*i, 3*j);

            /// Stabilization x
            const double g1_ij = rDN_DX(i,0) * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*l*g1_ij*rData.A1, 3*i, 3*j);

            /// Stabilization y
            const double g2_ij = rDN_DX(i,1) * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*l*g2_ij*rData.A2, 3*i, 3*j);
        }
    }
}

template<std::size_t TNumNodes>
double WaveElement<TNumNodes>::StabilizationParameter(const ElementData& rData) const
{
    const double inv_c = std::sqrt(InverseHeight(rData) / rData.gravity);
    return rData.length * rData.stab_factor * inv_c;
}

template<std::size_t TNumNodes>
double WaveElement<TNumNodes>::InverseHeight(const ElementData& rData) const
{
    const double height = rData.height;
    const double epsilon = rData.relative_dry_height * rData.length;
    return ShallowWaterUtilities().InverseHeight(height, epsilon);
}

template<std::size_t TNumNodes>
const array_1d<double,3> WaveElement<TNumNodes>::VectorProduct(
    const array_1d<array_1d<double,3>,TNumNodes>& rV,
    const array_1d<double,TNumNodes>& rN) const
{
    array_1d<double,3> result = ZeroVector(3);
    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        result += rV[i] * rN[i];
    }
    return result;
}

template<>
void WaveElement<3>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != mLocalSize)
        rLeftHandSideMatrix.resize(mLocalSize, mLocalSize, false);

    if(rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);

    LocalMatrixType lhs = ZeroMatrix(mLocalSize, mLocalSize);
    LocalVectorType rhs = ZeroVector(mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, this->GetGeometry());

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, area);

    CalculateGaussPointData(data, N);

    AddWaveTerms(lhs, rhs, data, N, DN_DX);
    AddFrictionTerms(lhs, rhs, data, N, DN_DX);
    AddArtificialViscosityTerms(lhs, data, DN_DX);

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rhs) -= prod(lhs, this->GetUnknownVector(data));

    noalias(rLeftHandSideMatrix) = area * lhs;
    noalias(rRightHandSideVector) = area * rhs;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != mLocalSize)
        rLeftHandSideMatrix.resize(mLocalSize, mLocalSize, false);

    if(rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);

    LocalMatrixType lhs = ZeroMatrix(mLocalSize, mLocalSize);
    LocalVectorType rhs = ZeroVector(mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, this->GetGeometry());

    Vector weights;
    Matrix N_container;
    ShapeFunctionsGradientsType DN_DX_container;
    CalculateGeometryData(weights, N_container, DN_DX_container);
    const IndexType num_gauss_points = weights.size();

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);
        const BoundedMatrix<double,TNumNodes,2> DN_DX = DN_DX_container[g];

        CalculateGaussPointData(data, N);

        AddWaveTerms(lhs, rhs, data, N, DN_DX, weight);
        AddFrictionTerms(lhs, rhs, data, N, DN_DX, weight);
        AddArtificialViscosityTerms(lhs, data, DN_DX, weight);
    }

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rhs) -= prod(lhs, this->GetUnknownVector(data));

    noalias(rLeftHandSideMatrix) = lhs;
    noalias(rRightHandSideVector) = rhs;
}

template<>
void WaveElement<3>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if(rMassMatrix.size1() != mLocalSize)
        rMassMatrix.resize(mLocalSize, mLocalSize, false);

    LocalMatrixType mass_matrix = ZeroMatrix(mLocalSize, mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, this->GetGeometry());

    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, area);

    CalculateGaussPointData(data, N);

    AddMassTerms(mass_matrix, data, N, DN_DX);

    noalias(rMassMatrix) = area * mass_matrix;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if(rMassMatrix.size1() != mLocalSize)
        rMassMatrix.resize(mLocalSize, mLocalSize, false);

    LocalMatrixType mass_matrix = ZeroMatrix(mLocalSize, mLocalSize);

    // Struct to pass around the data
    ElementData data;
    InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, this->GetGeometry());

    Vector weights;
    Matrix N_container;
    ShapeFunctionsGradientsType DN_DX_container;
    CalculateGeometryData(weights, N_container, DN_DX_container);
    const IndexType num_gauss_points = weights.size();

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);
        const BoundedMatrix<double,TNumNodes,2> DN_DX = DN_DX_container[g];

        CalculateGaussPointData(data, N);

        AddMassTerms(mass_matrix, data, N, DN_DX, weight);
    }
    noalias(rMassMatrix) = mass_matrix;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

template class WaveElement<3>;
template class WaveElement<4>;
template class WaveElement<6>;
template class WaveElement<8>;
template class WaveElement<9>;

} // namespace Kratos
