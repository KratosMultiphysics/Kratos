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
#include "shallow_element.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/phase_function.h"
#include "shallow_water_application_variables.h"
#include "custom_friction_laws/friction_laws_factory.h"

namespace Kratos
{

template<class TElementData>
int ShallowElement<TElementData>::Check(const ProcessInfo& rCurrentProcessInfo) const
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

        KRATOS_CHECK_DOF_IN_NODE(TElementData::UnknownComponent(0), r_node)
        KRATOS_CHECK_DOF_IN_NODE(TElementData::UnknownComponent(1), r_node)
        KRATOS_CHECK_DOF_IN_NODE(TElementData::UnknownComponent(2), r_node)
    }

    return err;

    KRATOS_CATCH("")
}

template<class TElementData>
void ShallowElement<TElementData>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rResult.size() != TLocalSize)
        rResult.resize(TLocalSize, false); // False says not to preserve existing storage!!

    const GeometryType& r_geom = this->GetGeometry();
    int counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = r_geom[i].GetDof(TElementData::UnknownComponent(0)).EquationId();
        rResult[counter++] = r_geom[i].GetDof(TElementData::UnknownComponent(1)).EquationId();
        rResult[counter++] = r_geom[i].GetDof(TElementData::UnknownComponent(2)).EquationId();
    }

    KRATOS_CATCH("")
}

template<class TElementData>
void ShallowElement<TElementData>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rElementalDofList.size() != TLocalSize)
        rElementalDofList.resize(TLocalSize);

    const GeometryType& r_geom = this->GetGeometry();
    int counter=0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[counter++] = r_geom[i].pGetDof(TElementData::UnknownComponent(0));
        rElementalDofList[counter++] = r_geom[i].pGetDof(TElementData::UnknownComponent(1));
        rElementalDofList[counter++] = r_geom[i].pGetDof(TElementData::UnknownComponent(2));
    }

    KRATOS_CATCH("")
}

template<class TElementData>
void ShallowElement<TElementData>::GetValuesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != TLocalSize)
        rValues.resize(TLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(TElementData::UnknownComponent(0), Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(TElementData::UnknownComponent(1), Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(TElementData::UnknownComponent(2), Step);
    }
}

template<class TElementData>
void ShallowElement<TElementData>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != TLocalSize)
        rValues.resize(TLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(VERTICAL_VELOCITY, Step);
    }
}

template<class TElementData>
void ShallowElement<TElementData>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_ERROR << "ShallowElement::GetSecondDerivativesVector. This method is not supported by the formulation" << std::endl;
}

template<class TElementData>
void ShallowElement<TElementData>::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    for (IndexType PointNumber = 0; PointNumber < 1; PointNumber++)
        rValues[PointNumber] = double(this->GetValue(rVariable));
}

template<class TElementData>
void ShallowElement<TElementData>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
}

template<class TElementData>
void ShallowElement<TElementData>::CalculateArtificialViscosity(
    BoundedMatrix<double,3,3>& rViscosity,
    BoundedMatrix<double,2,2>& rDiffusion,
    const TElementData& rData,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX)
{
}

template<class TElementData>
void ShallowElement<TElementData>::CalculateArtificialDamping(
    BoundedMatrix<double,3,3>& rFriction,
    const TElementData& rData)
{
}

template<class TElementData>
void ShallowElement<TElementData>::AddWaveTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const TElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
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
            if (rData.integrate_by_parts) {
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

template<class TElementData>
void ShallowElement<TElementData>::AddFrictionTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const TElementData& rData,
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

    BoundedMatrix<double,3,3> art_s = ZeroMatrix(3,3);
    this->CalculateArtificialDamping(art_s, rData);
    Sf += art_s;

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

template<class TElementData>
void ShallowElement<TElementData>::AddDispersiveTerms(
    LocalVectorType& rVector,
    const TElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{}

template<class TElementData>
void ShallowElement<TElementData>::AddArtificialViscosityTerms(
    LocalMatrixType& rMatrix,
    const TElementData& rData,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    BoundedMatrix<double,3,3> D = ZeroMatrix(3,3);
    BoundedMatrix<double,2,2> C = ZeroMatrix(2,2);
    BoundedMatrix<double,2,3> biq = ZeroMatrix(2,3);
    BoundedMatrix<double,2,3> bjq = ZeroMatrix(2,3);
    BoundedMatrix<double,3,2> tmp = ZeroMatrix(3,2);
    array_1d<double,2> bih = ZeroVector(2);
    array_1d<double,2> bjh = ZeroVector(2);

    this->CalculateArtificialViscosity(D, C, rData, rDN_DX);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            biq(0,0) = rDN_DX(i,0);
            biq(0,2) = rDN_DX(i,1);
            biq(1,1) = rDN_DX(i,1);
            biq(1,2) = rDN_DX(i,0);

            bjq(0,0) = rDN_DX(j,0);
            bjq(0,2) = rDN_DX(j,1);
            bjq(1,1) = rDN_DX(j,1);
            bjq(1,2) = rDN_DX(j,0);

            bih[0] = rDN_DX(i,0);
            bih[1] = rDN_DX(i,1);

            bjh[0] = rDN_DX(j,0);
            bjh[1] = rDN_DX(j,1);

            tmp = prod(D, trans(bjq));
            MathUtils<double>::AddMatrix(rMatrix, Weight*prod(biq, tmp), 3*i, 3*j);
            rMatrix(3*i + 2, 3*j + 2) += inner_prod(bih, Weight*prod(C, bjh));
        }
    }
}

template<class TElementData>
void ShallowElement<TElementData>::AddMassTerms(
    LocalMatrixType& rMatrix,
    const TElementData& rData,
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
            const double n_ij = TElementData::ShapeFunctionProduct(rN, i, j);
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

template<class TElementData>
double ShallowElement<TElementData>::StabilizationParameter(const TElementData& rData) const
{
    const double inv_c = std::sqrt(InverseHeight(rData) / rData.gravity);
    return rData.length * rData.stab_factor * inv_c;
}


template<class TElementData>
double ShallowElement<TElementData>::InverseHeight(const TElementData& rData) const
{
    const double height = rData.height;
    const double epsilon = rData.relative_dry_height * rData.length;
    return PhaseFunction::InverseHeight(height, epsilon);
}


template<class TElementData>
void ShallowElement<TElementData>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != TLocalSize)
        rLeftHandSideMatrix.resize(TLocalSize, TLocalSize, false);

    if(rRightHandSideVector.size() != TLocalSize)
        rRightHandSideVector.resize(TLocalSize, false);

    LocalMatrixType lhs = ZeroMatrix(TLocalSize, TLocalSize);
    LocalVectorType rhs = ZeroVector(TLocalSize);

    // Struct to pass around the data
    TElementData data;
    data.InitializeData(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    data.CalculateGeometryData(this->GetGeometry());
    data.SetNodalData(this->GetGeometry());

    for (IndexType g = 0; g < data.NumberOfGaussPoints(); ++g)
    {
        const double weight = data.GetWeight(g);
        const array_1d<double,TNumNodes> N = data.GetShapeFunctions(g);
        const BoundedMatrix<double,TNumNodes,2> DN_DX = data.GetShapeFunctionDerivatives(g);

        data.UpdateGaussPointData(N);

        AddWaveTerms(lhs, rhs, data, N, DN_DX, weight);
        AddFrictionTerms(lhs, rhs, data, N, DN_DX, weight);
        AddDispersiveTerms(rhs, data, N, DN_DX, weight);
        AddArtificialViscosityTerms(lhs, data, DN_DX, weight);
    }

    // Substracting the Dirichlet term (since we use a residualbased approach)
    noalias(rhs) -= prod(lhs, data.GetUnknownVector());

    noalias(rLeftHandSideMatrix) = lhs;
    noalias(rRightHandSideVector) = rhs;
}


template<class TElementData>
void ShallowElement<TElementData>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if(rMassMatrix.size1() != TLocalSize)
        rMassMatrix.resize(TLocalSize, TLocalSize, false);

    LocalMatrixType mass_matrix = ZeroMatrix(TLocalSize, TLocalSize);

    // Struct to pass around the data
    TElementData data;
    data.InitializeData(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    data.CalculateGeometryData(this->GetGeometry());
    data.SetNodalData(this->GetGeometry());

    for (IndexType g = 0; g < data.NumberOfGaussPoints(); ++g)
    {
        const double weight = data.GetWeight(g);
        const array_1d<double,TNumNodes> N = data.GetShapeFunctions(g);
        const BoundedMatrix<double,TNumNodes,2> DN_DX = data.GetShapeFunctionDerivatives(g);

        data.UpdateGaussPointData(N);

        AddMassTerms(mass_matrix, data, N, DN_DX, weight);
    }
    noalias(rMassMatrix) = mass_matrix;
}


template<class TElementData>
void ShallowElement<TElementData>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

template class ShallowElement<ShallowElementData<3>>;
template class ShallowElement<ShallowElementData<4>>;
template class ShallowElement<ShallowElementData<6>>;
template class ShallowElement<ShallowElementData<8>>;
template class ShallowElement<ShallowElementData<9>>;

} // namespace Kratos
