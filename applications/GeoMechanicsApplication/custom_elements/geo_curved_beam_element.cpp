// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// Application includes
#include "custom_elements/geo_curved_beam_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

#include <cmath>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoCurvedBeamElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                               NodesArrayType const& ThisNodes,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new GeoCurvedBeamElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoCurvedBeamElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                               GeometryType::Pointer pGeom,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new GeoCurvedBeamElement(NewId, pGeom, pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
int GeoCurvedBeamElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive area and Id > 0
    int ierr = GeoStructuralBaseElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& rProp = this->GetProperties();

    if (I33.Key() == 0 || rProp.Has(I33) == false || rProp[I33] < 0.0)
        KRATOS_ERROR << "I33 has Key zero, is not defined or has an invalid "
                        "value at element: "
                     << this->Id() << std::endl;

    if (CROSS_AREA.Key() == 0 || rProp.Has(CROSS_AREA) == false || rProp[CROSS_AREA] < 0.0)
        KRATOS_ERROR << "CROSS_AREA has Key zero, is not defined or has an "
                        "invalid value at element: "
                     << this->Id() << std::endl;

    if constexpr (TDim > 2) {
        if (TORSIONAL_INERTIA.Key() == 0 || rProp.Has(TORSIONAL_INERTIA) == false || rProp[TORSIONAL_INERTIA] < 0.0)
            KRATOS_ERROR << "TORSIONAL_INERTIA has Key zero, is not defined or "
                            "has an invalid value at element: "
                         << this->Id() << std::endl;

        if (I22.Key() == 0 || rProp.Has(I22) == false || rProp[I22] < 0.0)
            KRATOS_ERROR << "I22 has Key zero, is not defined or has an "
                            "invalid value at element: "
                         << this->Id() << std::endl;
    }

    return 0;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <>
void GeoCurvedBeamElement<3, 3>::SetRotationalInertiaVector(const PropertiesType& rProp, Vector& rRotationalInertia) const
{
    KRATOS_TRY

    if (rRotationalInertia.size() != N_DOF_NODE_ROT)
        rRotationalInertia.resize(N_DOF_NODE_ROT, false);

    unsigned int index          = 0;
    rRotationalInertia[index++] = rProp[TORSIONAL_INERTIA];
    rRotationalInertia[index++] = rProp[I22];
    rRotationalInertia[index++] = rProp[I33];

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <>
void GeoCurvedBeamElement<2, 3>::SetRotationalInertiaVector(const PropertiesType& rProp, Vector& rRotationalInertia) const
{
    KRATOS_TRY

    if (rRotationalInertia.size() != N_DOF_NODE_ROT)
        rRotationalInertia.resize(N_DOF_NODE_ROT, false);

    rRotationalInertia[0] = rProp[I33];

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resizing mass matrix
    if (rMassMatrix.size1() != N_DOF_ELEMENT)
        rMassMatrix.resize(N_DOF_ELEMENT, N_DOF_ELEMENT, false);
    noalias(rMassMatrix) = ZeroMatrix(N_DOF_ELEMENT, N_DOF_ELEMENT);

    const PropertiesType&                           rProp = this->GetProperties();
    const GeometryType&                             rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(mThisIntegrationMethod);

    // Defining shape functions and the determinant of the jacobian at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);
    Vector        detJContainer(IntegrationPoints.size());
    rGeom.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);

    // Defining necessary variables
    const double Density = rProp[DENSITY];

    Vector RotationalInertia;
    this->SetRotationalInertiaVector(rProp, RotationalInertia);

    unsigned int index = 0;
    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); GPoint++) {
        // calculating weighting coefficient for integration
        double IntegrationCoefficient = IntegrationPoints[GPoint].Weight() * detJContainer[GPoint];

        // Adding contribution to Mass matrix
        //  loop over nodes
        for (unsigned int node = 0; node < TNumNodes; ++node) {
            // displacement degrees of freedom
            for (unsigned int dof = 0; dof < N_DOF_NODE_DISP; ++dof) {
                const unsigned int i = index++;
                rMassMatrix(i, i) += Density * NContainer(GPoint, node) * IntegrationCoefficient;
            }

            // rotational degrees of freedom
            for (unsigned int dof = 0; dof < N_DOF_NODE_ROT; ++dof) {
                const unsigned int i = index++;
                rMassMatrix(i, i) += RotationalInertia[dof] * NContainer(GPoint, node) * IntegrationCoefficient;
            }
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                         VectorType&        rRightHandSideVector,
                                                         const ProcessInfo& rCurrentProcessInfo,
                                                         const bool CalculateStiffnessMatrixFlag,
                                                         const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           rProp = this->GetProperties();
    const GeometryType&                             rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPointsAlong =
        rGeom.IntegrationPoints(mThisIntegrationMethod);

    // Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);

    // calculating the local gradients
    const ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
    if (CalculateStiffnessMatrixFlag)
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, ConstitutiveParameters, rGeom, rProp, rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPointsAlong.size(); ++GPointAlong) {
        // Compute Nu, GradNe, B and StrainVector
        noalias(Variables.GradNe) = DN_De[GPointAlong];
        noalias(Variables.Nu)     = row(NContainer, GPointAlong);

        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.NuTot, NContainer, GPointAlong);

        this->CalculateTransformationMatrix(Variables.TransformationMatrix, Variables.GradNe);

        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.GaussVolumeAcceleration, NContainer, Variables.NodalVolumeAcceleration, GPointAlong);

        for (unsigned int GPointCross = 0; GPointCross < GetCrossNumberIntegrationPoints(); ++GPointCross) {
            int GPoint = GPointAlong * GetCrossNumberIntegrationPoints() + GPointCross;

            BoundedMatrix<double, TDim, TDim> JacobianMatrix;
            this->CalculateJacobianMatrix(GPointCross, Variables, JacobianMatrix);

            double                            detJacobian;
            BoundedMatrix<double, TDim, TDim> InvertJacobianMatrix;
            MathUtils<double>::InvertMatrix(JacobianMatrix, InvertJacobianMatrix, detJacobian);

            this->CalculateBMatrix(Variables.B, GPointCross, InvertJacobianMatrix, Variables);

            this->CalculateStrainVector(Variables);

            // Compute constitutive tensor
            ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            // Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(
                GPointCross, detJacobian, IntegrationPointsAlong[GPointAlong].Weight());

            // Contributions to the left hand side
            if (CalculateStiffnessMatrixFlag)
                this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

            // Contributions to the right hand side
            if (CalculateResidualVectorFlag)
                this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::InitializeElementVariables(ElementVariables& rVariables,
                                                                       ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                                                       const GeometryType&   rGeom,
                                                                       const PropertiesType& rProp,
                                                                       const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Properties variables
    rVariables.HalfThickness = 0.5 * std::sqrt(12.0 * (rProp[I33] / rProp[CROSS_AREA]));

    // ProcessInfo variables

    // Nodal Variables
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.DisplacementVector, rGeom, DISPLACEMENT);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VelocityVector, rGeom, VELOCITY);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.NodalVolumeAcceleration,
                                                                 rGeom, VOLUME_ACCELERATION);

    rVariables.DofValuesVector.resize(N_DOF_ELEMENT);
    GetNodalDofValuesVector(rVariables.DofValuesVector, rGeom);

    rVariables.NodalCrossDirection.resize(TNumNodes, TDim);
    CalculateNodalCrossDirection(rVariables.NodalCrossDirection);

    // Variables computed at each GP
    rVariables.B.resize(VoigtSize, N_DOF_ELEMENT, false);

    noalias(rVariables.NuTot) = ZeroMatrix(TDim, TNumNodes * TDim);

    rVariables.TransformationMatrix.resize(VoigtSize, VoigtSize, false);
    rVariables.UVoigtMatrix.resize(N_DOF_ELEMENT, VoigtSize, false);

    // Constitutive Law parameters
    rVariables.StrainVector.resize(VoigtSize, false);
    rVariables.StressVector.resize(VoigtSize, false);
    rVariables.ConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);
    rVariables.Nu.resize(TNumNodes, false);
    rVariables.GradNe.resize(TNumNodes, 1, false);

    rVariables.F.resize(TDim, TDim, false);
    rVariables.detF = 1.0;
    rConstitutiveParameters.SetStrainVector(rVariables.StrainVector);
    rConstitutiveParameters.SetStressVector(rVariables.StressVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rConstitutiveParameters.SetShapeFunctionsValues(rVariables.Nu);
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rVariables.GradNe);
    rConstitutiveParameters.SetDeformationGradientF(rVariables.F);
    rConstitutiveParameters.SetDeterminantF(rVariables.detF);
    // Auxiliary variables

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                               ElementVariables& rVariables) const
{
    KRATOS_TRY

    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B), rVariables.ConstitutiveMatrix);
    rLeftHandSideMatrix += prod(rVariables.UVoigtMatrix, rVariables.B) * rVariables.IntegrationCoefficient;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                               ElementVariables& rVariables,
                                                               unsigned int      GPoint) const
{
    KRATOS_TRY

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables, GPoint);
    this->CalculateAndAddBodyForce(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateAndAddBodyForce(VectorType& rRightHandSideVector,
                                                                     ElementVariables& rVariables) const
{
    KRATOS_TRY

    const PropertiesType& rProp   = this->GetProperties();
    const double&         density = rProp[DENSITY];

    // Distribute body force block vector into elemental vector
    noalias(rVariables.UVector) =
        density * prod(trans(rVariables.NuTot), rVariables.GaussVolumeAcceleration) * rVariables.IntegrationCoefficient;

    // Distribute body force block vector into elemental vector
    GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, rVariables.UVector);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector,
                                                                          ElementVariables& rVariables,
                                                                          unsigned int GPoint) const
{
    KRATOS_TRY

    // Distribute stiffness block vector into elemental vector
    rRightHandSideVector -= prod(trans(rVariables.B), mStressVector[GPoint]) * rVariables.IntegrationCoefficient;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
double GeoCurvedBeamElement<TDim, TNumNodes>::CalculateAngleAtNode(
    unsigned int GPoint, const BoundedMatrix<double, TNumNodes, TNumNodes>& DN_DeContainer) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    double dx = 0;
    double dy = 0;

    // loop over nodes
    for (unsigned int node = 0; node < TNumNodes; ++node) {
        dx += DN_DeContainer(GPoint, node) * rGeom[node].X0();
        dy += DN_DeContainer(GPoint, node) * rGeom[node].Y0();
    }

    return atan2(dx, -dy);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
double GeoCurvedBeamElement<TDim, TNumNodes>::CalculateAngleAtGaussPoint(const Matrix& GradNe) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();

    double dx = 0;
    double dy = 0;
    // loop over nodes
    for (unsigned int node = 0; node < TNumNodes; ++node) {
        dx += GradNe(node, 0) * rGeom[node].X0();
        dy += GradNe(node, 0) * rGeom[node].Y0();
    }

    return atan2(dy, dx);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <>
void GeoCurvedBeamElement<2, 3>::CalculateTransformationMatrix(Matrix&       TransformationMatrix,
                                                               const Matrix& GradNe) const
{
    KRATOS_TRY

    // details can be found in
    // "1. Geometrically non-linear formulation for the axisymmetric shell
    // elements" Eq.(8)

    const double phi = CalculateAngleAtGaussPoint(GradNe);

    const double cosPhi = cos(phi);
    const double sinPhi = sin(phi);

    TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_XX) = cosPhi * cosPhi;
    TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_YY) = sinPhi * sinPhi;
    TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_XY) = -2.0 * sinPhi * cosPhi;

    TransformationMatrix(INDEX_2D_BEAM_YY, INDEX_2D_BEAM_XX) =
        TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_YY);
    TransformationMatrix(INDEX_2D_BEAM_YY, INDEX_2D_BEAM_YY) =
        TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_XX);
    TransformationMatrix(INDEX_2D_BEAM_YY, INDEX_2D_BEAM_XY) =
        -TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_XY);

    TransformationMatrix(INDEX_2D_BEAM_XY, INDEX_2D_BEAM_XX) = sinPhi * cosPhi;
    TransformationMatrix(INDEX_2D_BEAM_XY, INDEX_2D_BEAM_YY) =
        -TransformationMatrix(INDEX_2D_BEAM_XY, INDEX_2D_BEAM_XX);
    TransformationMatrix(INDEX_2D_BEAM_XY, INDEX_2D_BEAM_XY) =
        TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_XX) -
        TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_YY);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <>
void GeoCurvedBeamElement<3, 3>::CalculateTransformationMatrix(Matrix&       TransformationMatrix,
                                                               const Matrix& GradNe) const
{
    KRATOS_TRY

    KRATOS_ERROR << "Undefined dimension in CalculateTransformationMatrix" << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateNodalCrossDirection(Matrix& NodalCrossDirection) const
{
    KRATOS_TRY
    // details can be found in
    // "1. Geometrically non-linear formulation for the axisymmetric shell
    // elements" Eq.(1)

    BoundedMatrix<double, TNumNodes, TNumNodes> DN_DeContainer;
    GeoElementUtilities::CalculateNewtonCotesLocalShapeFunctionsGradients(DN_DeContainer);

    for (unsigned int node = 0; node < TNumNodes; ++node) {
        double phi = CalculateAngleAtNode(node, DN_DeContainer);

        NodalCrossDirection(node, INDEX_X) = cos(phi);
        NodalCrossDirection(node, INDEX_Y) = sin(phi);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateJacobianMatrix(unsigned int GPointCross,
                                                                    const ElementVariables& rVariables,
                                                                    BoundedMatrix<double, TDim, TDim>& JacobianMatrix) const
{
    KRATOS_TRY
    // details can be found in
    // "1. Geometrically non-linear formulation for the axisymmetric shell
    // elements" Eqs.(1), (20) - (22)

    const GeometryType&       rGeom = this->GetGeometry();
    const double&             t     = rVariables.HalfThickness;
    const std::vector<double> CrossEta{-1.0 / sqrt(3), 1.0 / sqrt(3)};

    noalias(JacobianMatrix) = ZeroMatrix(TDim, TDim);

    for (unsigned int node = 0; node < TNumNodes; ++node) {
        const double& Vx = rVariables.NodalCrossDirection(node, INDEX_X);
        const double& Vy = rVariables.NodalCrossDirection(node, INDEX_Y);

        // dx/d_xi
        JacobianMatrix(0, 0) +=
            rVariables.GradNe(node, 0) * (rGeom[node].X0() + t * CrossEta[GPointCross] * Vx);

        // dy/d_xi
        JacobianMatrix(0, 1) +=
            rVariables.GradNe(node, 0) * (rGeom[node].Y0() + t * CrossEta[GPointCross] * Vy);

        // dx/d_eta
        JacobianMatrix(1, 0) += rVariables.Nu(node) * t * Vx;

        // dy/d_eta
        JacobianMatrix(1, 1) += rVariables.Nu(node) * t * Vy;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateBMatrix(Matrix&      BTransformed,
                                                             unsigned int GPointCross,
                                                             const BoundedMatrix<double, TDim, TDim>& InvJ,
                                                             ElementVariables& rVariables) const
{
    KRATOS_TRY

    Matrix B;
    this->CalculateLocalBMatrix(B, GPointCross, InvJ, rVariables);

    if (BTransformed.size1() != rVariables.TransformationMatrix.size1() || BTransformed.size2() != B.size2())
        BTransformed.resize(rVariables.TransformationMatrix.size1(), B.size2(), false);

    noalias(BTransformed) = prod(trans(rVariables.TransformationMatrix), B);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateLocalBMatrix(Matrix&      B,
                                                                  unsigned int GPointCross,
                                                                  const BoundedMatrix<double, TDim, TDim>& InvJ,
                                                                  ElementVariables& rVariables) const
{
    KRATOS_TRY
    // Details of derivation of linear part of B-Matrix can be found in:
    // "1. Geometrically non-linear formulation for the axisymmetric shell elements"
    // Note: In order to find B-Matrix, substitute Eqs.(23) & (26) into Eq.(11) to obtain Eq.(10)

    if (B.size1() != VoigtSize || B.size2() != N_DOF_ELEMENT)
        B.resize(VoigtSize, N_DOF_ELEMENT, false);

    const double&             t = rVariables.HalfThickness;
    const std::vector<double> CrossEta{-1.0 / sqrt(3), 1.0 / sqrt(3)};

    noalias(B)        = ZeroMatrix(VoigtSize, N_DOF_ELEMENT);
    const double& A11 = InvJ(0, 0);
    const double& A21 = InvJ(1, 0);
    const double& A12 = InvJ(0, 1);
    const double& A22 = InvJ(1, 1);
    const double& eta = CrossEta[GPointCross];

    for (unsigned int node = 0; node < TNumNodes; ++node) {
        const unsigned int index   = N_DOF_NODE * node;
        const unsigned int index_x = index + INDEX_2D_BEAM_X;
        const unsigned int index_y = index + INDEX_2D_BEAM_Y;
        const double&      Fjx     = -rVariables.NodalCrossDirection(node, INDEX_Y);
        const double&      Fjy     = rVariables.NodalCrossDirection(node, INDEX_X);

        // Parts of B-matrix due to u_x and u_y (u, v)
        B(INDEX_2D_BEAM_XX, index_x) = rVariables.GradNe(node, 0) * A11;
        B(INDEX_2D_BEAM_YY, index_y) = rVariables.GradNe(node, 0) * A21;
        B(INDEX_2D_BEAM_XY, index_x) = B(INDEX_2D_BEAM_YY, index_y);
        B(INDEX_2D_BEAM_XY, index_y) = B(INDEX_2D_BEAM_XX, index_x);

        // Parts of B-matrix due to rotation (theta)
        const unsigned int index_t = index + INDEX_2D_BEAM_T;

        // eps_xx due to theta
        const double term_xx = A12 * rVariables.Nu(node) + eta * rVariables.GradNe(node, 0) * A11;
        B(INDEX_2D_BEAM_XX, index_t) = t * Fjx * term_xx;

        // eps_yy due to theta
        const double term_yy = A22 * rVariables.Nu(node) + eta * rVariables.GradNe(node, 0) * A21;
        B(INDEX_2D_BEAM_YY, index_t) = t * Fjy * term_yy;

        // eps_xy due to theta
        B(INDEX_2D_BEAM_XY, index_t) = t * Fjx * term_yy + t * Fjy * term_xx;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateStrainVector(ElementVariables& rVariables) const
{
    KRATOS_TRY

    noalias(rVariables.StrainVector) = prod(rVariables.B, rVariables.DofValuesVector);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
double GeoCurvedBeamElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(unsigned int GPointCross,
                                                                              double detJ,
                                                                              double weight) const
{
    const std::vector<double> CrossWeight{1.0, 1.0};

    return weight * CrossWeight[GPointCross] * detJ;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
SizeType GeoCurvedBeamElement<TDim, TNumNodes>::GetAlongNumberIntegrationPoints() const
{
    return this->GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
SizeType GeoCurvedBeamElement<TDim, TNumNodes>::GetCrossNumberIntegrationPoints() const
{
    return N_POINT_CROSS;
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateLocalInternalForce(VectorType& rInternalForceVector,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           rProp = this->GetProperties();
    const GeometryType&                             rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPointsAlong =
        rGeom.IntegrationPoints(mThisIntegrationMethod);

    // Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);

    // calculating the local gradients
    const ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    // Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, ConstitutiveParameters, rGeom, rProp, rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPointsAlong.size(); ++GPointAlong) {
        // Compute Nu, GradNe, B and StrainVector
        noalias(Variables.GradNe) = DN_De[GPointAlong];
        noalias(Variables.Nu)     = row(NContainer, GPointAlong);

        this->CalculateTransformationMatrix(Variables.TransformationMatrix, Variables.GradNe);

        for (unsigned int GPointCross = 0; GPointCross < GetCrossNumberIntegrationPoints(); ++GPointCross) {
            int GPoint = GPointAlong * GetCrossNumberIntegrationPoints() + GPointCross;

            BoundedMatrix<double, TDim, TDim> JacobianMatrix;
            this->CalculateJacobianMatrix(GPointCross, Variables, JacobianMatrix);

            double                            detJacobian;
            BoundedMatrix<double, TDim, TDim> InvertJacobianMatrix;
            MathUtils<double>::InvertMatrix(JacobianMatrix, InvertJacobianMatrix, detJacobian);

            this->CalculateLocalBMatrix(Variables.B, GPointCross, InvertJacobianMatrix, Variables);

            this->CalculateStrainVector(Variables);

            // Compute constitutive tensor
            ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            // Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(
                GPointCross, detJacobian, IntegrationPointsAlong[GPointAlong].Weight());

            // Contributions to the right hand side
            this->CalculateAndAddStiffnessForce(rInternalForceVector, Variables, GPoint);
        }
    }

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        for (unsigned int j = 0; j < N_DOF_NODE; ++j) {
            const int index = i * TNumNodes + j;
            if (i == 0) rInternalForceVector[index] = -rInternalForceVector[index];
            if (i == (TNumNodes - 1)) {
                rInternalForceVector[index] =
                    0.5 * (rInternalForceVector[j] + rInternalForceVector[TNumNodes + j]);
            }
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                         std::vector<Matrix>& rOutput,
                                                                         const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Defining necessary variables
    const GeometryType& rGeom      = this->GetGeometry();
    const unsigned int  NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);

    if (rOutput.size() != NumGPoints) rOutput.resize(NumGPoints);

    if (rVariable == CAUCHY_STRESS_TENSOR) {
        // Previous definitions
        const PropertiesType&                           rProp = this->GetProperties();
        const GeometryType::IntegrationPointsArrayType& IntegrationPointsAlong =
            rGeom.IntegrationPoints(mThisIntegrationMethod);

        // Containers of variables at all integration points
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);

        // calculating the local gradients
        const ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);

        // Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        // Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, ConstitutiveParameters, rGeom, rProp, rCurrentProcessInfo);

        Matrix AverageStresses = ZeroMatrix(VoigtSize, NumGPoints);

        // Loop over integration points
        for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPointsAlong.size(); ++GPointAlong) {
            // Compute Nu, GradNe, B and StrainVector
            noalias(Variables.GradNe) = DN_De[GPointAlong];
            noalias(Variables.Nu)     = row(NContainer, GPointAlong);

            this->CalculateTransformationMatrix(Variables.TransformationMatrix, Variables.GradNe);

            Vector AverageStressVector = ZeroVector(VoigtSize);

            for (unsigned int GPointCross = 0; GPointCross < GetCrossNumberIntegrationPoints(); ++GPointCross) {
                int GPoint = GPointAlong * GetCrossNumberIntegrationPoints() + GPointCross;

                BoundedMatrix<double, TDim, TDim> JacobianMatrix;
                this->CalculateJacobianMatrix(GPointCross, Variables, JacobianMatrix);

                double                            detJacobian;
                BoundedMatrix<double, TDim, TDim> InvertJacobianMatrix;
                MathUtils<double>::InvertMatrix(JacobianMatrix, InvertJacobianMatrix, detJacobian);

                this->CalculateBMatrix(Variables.B, GPointCross, InvertJacobianMatrix, Variables);

                this->CalculateStrainVector(Variables);
                // Compute constitutive tensor
                ConstitutiveParameters.SetStressVector(Variables.StressVector);
                mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
                AverageStressVector += Variables.StressVector;
            }
            AverageStressVector /= GetCrossNumberIntegrationPoints();
            for (unsigned i = 0; i < AverageStresses.size1(); ++i) {
                AverageStresses(i, GPointAlong) = AverageStressVector[i];
            }
        }

        this->InterpolateOnOutputPoints(AverageStresses);

        // Loop over integration points
        for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPointsAlong.size(); ++GPointAlong) {
            if (rOutput[GPointAlong].size2() != TDim)
                rOutput[GPointAlong].resize(TDim, TDim, false);

            Vector AverageStressVector = ZeroVector(VoigtSize);
            AverageStressVector        = column(AverageStresses, GPointAlong);

            rOutput[GPointAlong] = MathUtils<double>::StrainVectorToTensor(AverageStressVector);
        }
    } else if (rVariable == ENGINEERING_STRAIN_TENSOR) {
        // Previous definitions
        const PropertiesType&                           rProp = this->GetProperties();
        const GeometryType::IntegrationPointsArrayType& IntegrationPointsAlong =
            rGeom.IntegrationPoints(mThisIntegrationMethod);

        // Containers of variables at all integration points
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(mThisIntegrationMethod);

        // calculating the local gradients
        const ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);

        // Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        // Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, ConstitutiveParameters, rGeom, rProp, rCurrentProcessInfo);

        Matrix AverageStrains = ZeroMatrix(VoigtSize, NumGPoints);

        // Loop over integration points
        for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPointsAlong.size(); ++GPointAlong) {
            // Compute Nu, GradNe, B and StrainVector
            noalias(Variables.GradNe) = DN_De[GPointAlong];
            noalias(Variables.Nu)     = row(NContainer, GPointAlong);

            this->CalculateTransformationMatrix(Variables.TransformationMatrix, Variables.GradNe);

            Vector AverageStrainVector = ZeroVector(VoigtSize);

            for (unsigned int GPointCross = 0; GPointCross < GetCrossNumberIntegrationPoints(); ++GPointCross) {
                BoundedMatrix<double, TDim, TDim> JacobianMatrix;
                this->CalculateJacobianMatrix(GPointCross, Variables, JacobianMatrix);

                double                            detJacobian;
                BoundedMatrix<double, TDim, TDim> InvertJacobianMatrix;
                MathUtils<double>::InvertMatrix(JacobianMatrix, InvertJacobianMatrix, detJacobian);

                this->CalculateBMatrix(Variables.B, GPointCross, InvertJacobianMatrix, Variables);

                this->CalculateStrainVector(Variables);
                AverageStrainVector += Variables.StrainVector;
            }
            AverageStrainVector /= GetCrossNumberIntegrationPoints();
            for (unsigned i = 0; i < AverageStrains.size1(); ++i) {
                AverageStrains(i, GPointAlong) = AverageStrainVector[i];
            }
        }

        this->InterpolateOnOutputPoints(AverageStrains);

        // Loop over integration points
        for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPointsAlong.size(); ++GPointAlong) {
            if (rOutput[GPointAlong].size2() != TDim)
                rOutput[GPointAlong].resize(TDim, TDim, false);

            Vector AverageStrainVector = ZeroVector(VoigtSize);
            AverageStrainVector        = column(AverageStrains, GPointAlong);

            rOutput[GPointAlong] = MathUtils<double>::StrainVectorToTensor(AverageStrainVector);
        }
    } else {
        if (rOutput.size() != mConstitutiveLawVector.size())
            rOutput.resize(mConstitutiveLawVector.size());

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            rOutput[i].resize(TDim, TDim, false);
            noalias(rOutput[i]) = ZeroMatrix(TDim, TDim);
            rOutput[i]          = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::InterpolateOnOutputPoints(Matrix& ValuesMatrix) const
{
    KRATOS_TRY

    Vector ValuesVector = ZeroVector(ValuesMatrix.size2());

    for (unsigned int i = 0; i < ValuesMatrix.size1(); ++i) {
        noalias(ValuesVector) = row(ValuesMatrix, i);
        this->InterpolateOnOutputPoints(ValuesVector);
        for (unsigned j = 0; j < ValuesMatrix.size2(); ++j) {
            ValuesMatrix(i, j) = ValuesVector[j];
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::InterpolateOnOutputPoints(Vector& Values) const
{
    KRATOS_TRY
    // solve: y = m * x + b

    const std::vector<double> XiOutput{-1.0 / 3.0, 1.0 / 3.0};
    const std::vector<double> Xi{-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};

    KRATOS_ERROR_IF_NOT(Values.size() == Xi.size()) << "Wrong size for interpolation" << std::endl;

    const double m = (Values[1] - Values[0]) / (Xi[1] - Xi[0]);
    const double b = Values[0] - m * Xi[0];

    for (unsigned int i = 0; i < Values.size(); ++i)
        Values[i] = m * XiOutput[i] + b;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCurvedBeamElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                                         std::vector<array_1d<double, 3>>& rOutput,
                                                                         const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom      = this->GetGeometry();
    const unsigned int  NumGPoints = rGeom.IntegrationPointsNumber(mThisIntegrationMethod);
    if (rOutput.size() != NumGPoints) rOutput.resize(NumGPoints);

    if (rVariable == FORCE || rVariable == MOMENT) {
        // Containers of variables at all integration points
        //  GiD accepts equally distributed points
        Matrix NContainer;
        GeoElementUtilities::CalculateEquallyDistributedPointsLineShapeFunctions3N(NContainer);

        Vector InternalForce = ZeroVector(N_DOF_ELEMENT);
        this->CalculateLocalInternalForce(InternalForce, rCurrentProcessInfo);

        // Loop over integration points
        for (unsigned int GPointAlong = 0; GPointAlong < NumGPoints; ++GPointAlong) {
            Vector GaussInternalForce = ZeroVector(N_DOF_NODE);

            GeoElementUtilities::InterpolateVariableWithComponents<N_DOF_NODE, TNumNodes>(
                GaussInternalForce, NContainer, InternalForce, GPointAlong);

            Vector Force = ZeroVector(3);
            if (rVariable == FORCE) {
                Force[INDEX_2D_BEAM_X] = GaussInternalForce[INDEX_2D_BEAM_X];
                Force[INDEX_2D_BEAM_Y] = GaussInternalForce[INDEX_2D_BEAM_Y];
            } else if (rVariable == MOMENT) {
                Force[INDEX_2D_BEAM_T] = GaussInternalForce[INDEX_2D_BEAM_T];
            }

            if (rOutput[GPointAlong].size() != 3) rOutput[GPointAlong].resize(3, false);

            noalias(rOutput[GPointAlong]) = Force;
        }
    } else {
        if (rOutput.size() != mConstitutiveLawVector.size())
            rOutput.resize(mConstitutiveLawVector.size());

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            noalias(rOutput[i]) = ZeroVector(3);
            rOutput[i]          = mConstitutiveLawVector[i]->GetValue(rVariable, rOutput[i]);
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class GeoCurvedBeamElement<2, 3>;
template class GeoCurvedBeamElement<3, 3>;

} // Namespace Kratos
