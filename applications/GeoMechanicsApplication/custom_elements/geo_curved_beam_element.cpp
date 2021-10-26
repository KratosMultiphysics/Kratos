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
#include "geo_mechanics_application_variables.h"
#include "custom_elements/geo_curved_beam_element.hpp"
#include "custom_utilities/element_utilities.hpp"

#include <math.h>

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer GeoCurvedBeamElement<TDim,TNumNodes>::
    Create( IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties ) const
{
    // KRATOS_INFO("0-GeoCurvedBeamElement::Create") << this->Id() << std::endl;

    return Element::Pointer( new GeoCurvedBeamElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer GeoCurvedBeamElement<TDim,TNumNodes>::
    Create( IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties ) const
{
    // KRATOS_INFO("0-GeoCurvedBeamElement::Create") << this->Id() << std::endl;

    return Element::Pointer( new GeoCurvedBeamElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
                           const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_ERROR << "GeoCurvedBeamElement::CalculateLeftHandSide not implemented" << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
int GeoCurvedBeamElement<TDim,TNumNodes>::
    Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    // KRATOS_INFO("0-GeoCurvedBeamElement::Check") << this->Id() << std::endl;

    // Base class checks for positive area and Id > 0
    int ierr = GeoStructuralBaseElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& rProp = this->GetProperties();

    if ( I33.Key() == 0 ||
         rProp.Has( I33 ) == false ||
         rProp[I33] < 0.0 )
        KRATOS_ERROR << "I33 has Key zero, is not defined or has an invalid value at element: "
                     << this->Id()
                     << std::endl;

    if ( CROSS_AREA.Key() == 0 ||
         rProp.Has( CROSS_AREA ) == false ||
         rProp[CROSS_AREA] < 0.0 )
        KRATOS_ERROR << "CROSS_AREA has Key zero, is not defined or has an invalid value at element: "
                     << this->Id()
                     << std::endl;

    if (TDim > 2) {
        if ( TORSIONAL_INERTIA.Key() == 0 ||
             rProp.Has( TORSIONAL_INERTIA ) == false ||
             rProp[TORSIONAL_INERTIA] < 0.0 )
             KRATOS_ERROR << "TORSIONAL_INERTIA has Key zero, is not defined or has an invalid value at element: "
                          << this->Id()
                          << std::endl;

        if ( I22.Key() == 0 ||
             rProp.Has( I22 ) == false ||
             rProp[I22] < 0.0 )
             KRATOS_ERROR << "I22 has Key zero, is not defined or has an invalid value at element: "
                          << this->Id()
                          << std::endl;
    }

    // KRATOS_INFO("1-GeoCurvedBeamElement::Check") << std::endl;

    return 0;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<3,3>::
    SetRotationalInertiaVector(const PropertiesType& rProp, Vector& rRotationalInertia) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::SetRotationalInertiaVector") << std::endl;

    if ( rRotationalInertia.size() != N_DOF_NODE_ROT )
        rRotationalInertia.resize( N_DOF_NODE_ROT, false );

    unsigned int index = 0;
    rRotationalInertia[index++] = rProp[TORSIONAL_INERTIA];
    rRotationalInertia[index++] = rProp[I22];
    rRotationalInertia[index++] = rProp[I33];

    // KRATOS_INFO("1-GeoCurvedBeamElement::SetRotationalInertiaVector") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<2,3>::
    SetRotationalInertiaVector(const PropertiesType& rProp, Vector& rRotationalInertia) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::SetRotationalInertiaVector") << std::endl;

    if ( rRotationalInertia.size() != N_DOF_NODE_ROT )
        rRotationalInertia.resize( N_DOF_NODE_ROT, false );

    rRotationalInertia[0] = rProp[I33];

    // KRATOS_INFO("1-GeoCurvedBeamElement::SetRotationalInertiaVector") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateMassMatrix( MatrixType& rMassMatrix,
                         const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateMassMatrix") << std::endl;

    //Resizing mass matrix
    if ( rMassMatrix.size1() != N_DOF_ELEMENT )
        rMassMatrix.resize( N_DOF_ELEMENT, N_DOF_ELEMENT, false );
    noalias( rMassMatrix ) = ZeroMatrix( N_DOF_ELEMENT, N_DOF_ELEMENT );

    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( mThisIntegrationMethod );

    //Defining shape functions and the determinant of the jacobian at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );
    Vector detJContainer(IntegrationPoints.size());
    rGeom.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);

    //Defining necessary variables
    const double Density = rProp[DENSITY];

    Vector RotationalInertia;
    this->SetRotationalInertiaVector(rProp, RotationalInertia);

     unsigned int index = 0;
    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); GPoint++ ) {
        //calculating weighting coefficient for integration
        double IntegrationCoefficient = IntegrationPoints[GPoint].Weight() * detJContainer[GPoint];

        //Adding contribution to Mass matrix
        // loop over nodes
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

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateMassMatrix") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix,
                              const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateStiffnessMatrix") << std::endl;

    //Resizing mass matrix
    if ( rStiffnessMatrix.size1() != N_DOF_ELEMENT )
        rStiffnessMatrix.resize( N_DOF_ELEMENT, N_DOF_ELEMENT, false );
    noalias( rStiffnessMatrix ) = ZeroMatrix( N_DOF_ELEMENT, N_DOF_ELEMENT );

    //Previous definitions
    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( mThisIntegrationMethod );

    //Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );

    //calculating the local gradients
    const ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,
                                     ConstitutiveParameters,
                                     rGeom,
                                     rProp,
                                     CurrentProcessInfo);

    //Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPoints.size(); ++GPointAlong) {
        //Compute Nu, GradNe, B and StrainVector
        noalias(Variables.Nu)     = row(NContainer, GPointAlong);
        noalias(Variables.GradNe) = DN_De[GPointAlong];

        this->CalculateTransformationMatrix( Variables.TransformationMatrix,
                                             Variables.GradNe );

        for (unsigned int GPointCross = 0; GPointCross < N_POINT_CROSS; ++GPointCross) {
            int GPoint= GPointAlong * N_POINT_CROSS + GPointCross;

            BoundedMatrix<double,TDim, TDim> JacobianMatrix;
            this->CalculateJacobianMatrix( GPointCross,
                                           Variables,
                                           JacobianMatrix );

            double detJacobian;
            BoundedMatrix<double,TDim, TDim> InvertJacobianMatrix;
            MathUtils<double>::InvertMatrix(JacobianMatrix,
                                            InvertJacobianMatrix,
                                            detJacobian);

            this->CalculateBMatrix(Variables.B, GPointCross, InvertJacobianMatrix, Variables);

            this->CalculateStrainVector(Variables);

            //Compute constitutive tensor
            ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            //Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(GPointCross,
                                                                                     detJacobian,
                                                                                     IntegrationPoints[GPointAlong].Weight());

            //Compute stiffness matrix
            this->CalculateAndAddLHS(rStiffnessMatrix, Variables);
        }
    }

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateStiffnessMatrix") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAll") << this->Id() << std::endl;

    //Previous definitions
    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType&
        IntegrationPointsAlong = rGeom.IntegrationPoints( mThisIntegrationMethod );

    //Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );

    //calculating the local gradients
    const ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,rProp,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables( Variables,
                                      ConstitutiveParameters,
                                      rGeom,
                                      rProp,
                                      CurrentProcessInfo );

    //Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPointsAlong.size(); ++GPointAlong) {

        //Compute Nu, GradNe, B and StrainVector
        noalias(Variables.GradNe) = DN_De[GPointAlong];
        noalias(Variables.Nu)     = row(NContainer, GPointAlong);

        GeoElementUtilities::
            CalculateNuMatrix<TDim, TNumNodes>(Variables.NuTot,
                                               NContainer,
                                               GPointAlong);

        this->CalculateTransformationMatrix( Variables.TransformationMatrix,
                                             Variables.GradNe );

        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim,TNumNodes>( Variables.BodyAcceleration,
                                                               NContainer,
                                                               Variables.VolumeAcceleration,
                                                               GPointAlong );

        for (unsigned int GPointCross = 0; GPointCross < GetCrossNumberIntegrationPoints(); ++GPointCross) {
            int GPoint = GPointAlong * GetCrossNumberIntegrationPoints() + GPointCross;

            BoundedMatrix<double,TDim, TDim> JacobianMatrix;
            this->CalculateJacobianMatrix( GPointCross,
                                           Variables,
                                           JacobianMatrix );

            double detJacobian;
            BoundedMatrix<double,TDim, TDim> InvertJacobianMatrix;
            MathUtils<double>::InvertMatrix(JacobianMatrix,
                                            InvertJacobianMatrix,
                                            detJacobian);

            this->CalculateBMatrix(Variables.B, GPointCross, InvertJacobianMatrix, Variables);

            this->CalculateStrainVector(Variables);

            //Compute constitutive tensor
            ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            //Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(GPointCross,
                                                                                     detJacobian,
                                                                                     IntegrationPointsAlong[GPointAlong].Weight());

            //Contributions to the left hand side
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

            //Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
        }
    }

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAll") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    InitializeElementVariables( ElementVariables& rVariables,
                                ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                const GeometryType& rGeom,
                                const PropertiesType& rProp,
                                const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY

    // KRATOS_INFO("0-GeoCurvedBeamElement::InitializeElementVariables") << std::endl;

    //Properties variables
    rVariables.HalfThickness = 0.5 * std::sqrt(12.0 * (rProp[I33] / rProp[CROSS_AREA]));

    //ProcessInfo variables

    //Nodal Variables
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.DisplacementVector, rGeom, DISPLACEMENT);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VelocityVector,     rGeom, VELOCITY);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VolumeAcceleration, rGeom, VOLUME_ACCELERATION);

    rVariables.DofValuesVector.resize(N_DOF_ELEMENT);
    GetNodalDofValuesVector(rVariables.DofValuesVector, rGeom);

    rVariables.NodalCrossDirection.resize(TNumNodes, TDim);
    CalculateNodalCrossDirection(rVariables.NodalCrossDirection);

    //Variables computed at each GP
    rVariables.B.resize(VoigtSize, N_DOF_ELEMENT, false);

    noalias(rVariables.NuTot) = ZeroMatrix(TDim, TNumNodes*TDim);

    rVariables.TransformationMatrix.resize(VoigtSize, VoigtSize, false);
    rVariables.UVoigtMatrix.resize(N_DOF_ELEMENT, VoigtSize, false);

    //Constitutive Law parameters
    rVariables.StrainVector.resize(VoigtSize, false);
    rVariables.StressVector.resize(VoigtSize, false);
    rVariables.ConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);
    rVariables.Nu.resize(TNumNodes, false);
    rVariables.GradNe.resize(TNumNodes, 1, false);

    rVariables.F.resize(TDim,TDim,false);
    rVariables.detF = 1.0;
    rConstitutiveParameters.SetStrainVector(rVariables.StrainVector);
    rConstitutiveParameters.SetStressVector(rVariables.StressVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rConstitutiveParameters.SetShapeFunctionsValues(rVariables.Nu);
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rVariables.GradNe);
    rConstitutiveParameters.SetDeformationGradientF(rVariables.F);
    rConstitutiveParameters.SetDeterminantF(rVariables.detF);
    //Auxiliary variables

    // KRATOS_INFO("1-GeoCurvedBeamElement::InitializeElementVariables") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateRHS( VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateRHS") << std::endl;

    //Previous definitions
    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( mThisIntegrationMethod );

    //Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );

    //calculating the local gradients
    const ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,rProp,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables( Variables,
                                      ConstitutiveParameters,
                                      rGeom,
                                      rProp,
                                      CurrentProcessInfo );

    //Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPoints.size(); ++GPointAlong) {
        //Compute Nu, GradNe, B and StrainVector
        noalias(Variables.GradNe) = DN_De[GPointAlong];
        noalias(Variables.Nu)     = row(NContainer, GPointAlong);
        GeoElementUtilities::
            CalculateNuMatrix<TDim, TNumNodes>(Variables.NuTot,
                                               NContainer,
                                               GPointAlong);

        this->CalculateTransformationMatrix( Variables.TransformationMatrix,
                                             Variables.GradNe );

        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim,TNumNodes>( Variables.BodyAcceleration,
                                                               NContainer,
                                                               Variables.VolumeAcceleration,
                                                               GPointAlong );

        for (unsigned int GPointCross = 0; GPointCross < N_POINT_CROSS; ++GPointCross) {
            int GPoint= GPointAlong * N_POINT_CROSS + GPointCross;

            BoundedMatrix<double,TDim, TDim> JacobianMatrix;
            this->CalculateJacobianMatrix(GPointCross,
                                          Variables,
                                          JacobianMatrix);

            double detJacobian;
            BoundedMatrix<double,TDim, TDim> InvertJacobianMatrix;
            MathUtils<double>::InvertMatrix(JacobianMatrix,
                                            InvertJacobianMatrix,
                                            detJacobian);

            this->CalculateBMatrix(Variables.B, GPointCross, InvertJacobianMatrix, Variables);

            this->CalculateStrainVector(Variables);

            //Compute constitutive tensor
            ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            //Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(GPointCross,
                                                                                     detJacobian,
                                                                                     IntegrationPoints[GPointAlong].Weight());

            //Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
        }
    }

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateRHS") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, 
                       ElementVariables& rVariables ) const
{
    KRATOS_TRY

    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAndAddLHS") << std::endl;

    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B), rVariables.ConstitutiveMatrix);
    rLeftHandSideMatrix +=  prod(rVariables.UVoigtMatrix, rVariables.B)
                          * rVariables.IntegrationCoefficient;

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAndAddLHS") << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddRHS(VectorType& rRightHandSideVector,
                       ElementVariables& rVariables,
                       unsigned int GPoint ) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAndAddRHS") << std::endl;

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables, GPoint);

    this->CalculateAndAddBodyForce(rRightHandSideVector, rVariables);

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAndAddRHS") << std::endl;
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddBodyForce(VectorType& rRightHandSideVector,
                             ElementVariables& rVariables) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAndAddBodyForce") << std::endl;

    const PropertiesType& rProp = this->GetProperties();
    const double &density = rProp[DENSITY];

    //Distribute body force block vector into elemental vector
    noalias(rVariables.UVector) = density
                                 * prod(trans(rVariables.NuTot), rVariables.BodyAcceleration)
                                 * rVariables.IntegrationCoefficient;

    //Distribute body force block vector into elemental vector
    GeoElementUtilities::AssembleUBlockVector<TDim, TNumNodes>(rRightHandSideVector,rVariables.UVector);

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAndAddBodyForce") << std::endl;
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector,
                                  ElementVariables& rVariables,
                                  unsigned int GPoint ) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAndAddStiffnessForce") << std::endl;

    //Distribute stiffness block vector into elemental vector
    rRightHandSideVector -=  prod(trans(rVariables.B),mStressVector[GPoint])
                           * rVariables.IntegrationCoefficient;

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAndAddStiffnessForce") << force << std::endl;
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAngleAtNode(unsigned int GPoint,
                         const BoundedMatrix<double,TNumNodes, TNumNodes> & DN_DeContainer) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAngleAtNode") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();

    double dx = 0;
    double dy = 0;

    // loop over nodes
    for (unsigned int node = 0; node < TNumNodes; ++node) {
        dx += DN_DeContainer(GPoint, node) * rGeom[node].X0();
        dy += DN_DeContainer(GPoint, node) * rGeom[node].Y0();
    }

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAngleAtNode, dx") << dx << ", dy" << dy << std::endl;

    return atan2(dx, -dy);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAngleAtGaussPoint(const Matrix &GradNe) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAngleAtGaussPoint") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();

    double dx = 0;
    double dy = 0;
    // loop over nodes
    for (unsigned int node = 0; node < TNumNodes; ++node) {
        dx += GradNe(node, 0) * rGeom[node].X0();
        dy += GradNe(node, 0) * rGeom[node].Y0();
    }

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAngleAtGaussPoint, dy, dx") << dy << ", " << dx << std::endl;

    return atan2(dy, dx);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<2,3>::
    CalculateTransformationMatrix( Matrix &TransformationMatrix,
                                   const Matrix &GradNe ) const
{
    KRATOS_TRY

    // details can be found in
    // "1. Geometrically non-linear formulation for the axisymmetric shell elements"
    // Eq.(8)

    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateTransformationMatrix") << std::endl;

    const double phi = CalculateAngleAtGaussPoint(GradNe);

    const double cosPhi = cos(phi);
    const double sinPhi = sin(phi);

    TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_XX) = cosPhi * cosPhi;
    TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_YY) = sinPhi * sinPhi;
    TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_XY) = -2.0 * sinPhi * cosPhi;

    TransformationMatrix(INDEX_2D_BEAM_YY, INDEX_2D_BEAM_XX) = TransformationMatrix(INDEX_2D_BEAM_XX,INDEX_2D_BEAM_YY);
    TransformationMatrix(INDEX_2D_BEAM_YY, INDEX_2D_BEAM_YY) = TransformationMatrix(INDEX_2D_BEAM_XX,INDEX_2D_BEAM_XX);
    TransformationMatrix(INDEX_2D_BEAM_YY, INDEX_2D_BEAM_XY) = -TransformationMatrix(INDEX_2D_BEAM_XX,INDEX_2D_BEAM_XY);

    TransformationMatrix(INDEX_2D_BEAM_XY, INDEX_2D_BEAM_XX) = sinPhi * cosPhi;
    TransformationMatrix(INDEX_2D_BEAM_XY, INDEX_2D_BEAM_YY) = -TransformationMatrix(INDEX_2D_BEAM_XY, INDEX_2D_BEAM_XX);
    TransformationMatrix(INDEX_2D_BEAM_XY, INDEX_2D_BEAM_XY) =  TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_XX)
                                                              - TransformationMatrix(INDEX_2D_BEAM_XX, INDEX_2D_BEAM_YY);

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateTransformationMatrix") << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<3,3>::
    CalculateTransformationMatrix( Matrix &TransformationMatrix,
                                   const Matrix &GradNe ) const
{
    KRATOS_TRY

    KRATOS_ERROR << "Undefined dimension in CalculateTransformationMatrix" << std::endl;;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateNodalCrossDirection(Matrix& NodalCrossDirection) const
{
    KRATOS_TRY
    // details can be found in
    // "1. Geometrically non-linear formulation for the axisymmetric shell elements"
    // Eq.(1)

    //if (this->Id()==9) KRATOS_INFO("0-GeoCurvedBeamElement::CalculateNodalCrossDirection") << std::endl;

    BoundedMatrix<double, TNumNodes, TNumNodes> DN_DeContainer;
    GeoElementUtilities::CalculateNewtonCotesLocalShapeFunctionsGradients(DN_DeContainer);
    // KRATOS_INFO("NewtonCotes_DN_DeContainer") << DN_DeContainer << std::endl;

    for (unsigned int node = 0; node < TNumNodes; ++node) {
        double phi = CalculateAngleAtNode(node, DN_DeContainer);
        // KRATOS_INFO("phi") << phi << std::endl;

        NodalCrossDirection(node, INDEX_X) = cos(phi);
        NodalCrossDirection(node, INDEX_Y) = sin(phi);
    }

    //if (this->Id()==9) KRATOS_INFO("1-GeoCurvedBeamElement::CalculateNodalCrossDirection") << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
     CalculateJacobianMatrix(unsigned int GPointCross,
                             const ElementVariables &rVariables,
                             BoundedMatrix<double,TDim, TDim> &JacobianMatrix) const
{
    KRATOS_TRY
    // details can be found in
    // "1. Geometrically non-linear formulation for the axisymmetric shell elements"
    // Eqs.(1), (20) - (22)
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateJacobianMatrix") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();
    const double &t = rVariables.HalfThickness;
    const std::vector<double> CrossEta{-1.0/sqrt(3), 1.0/sqrt(3)};

    noalias(JacobianMatrix) = ZeroMatrix(TDim, TDim);

    // KRATOS_INFO("rVariables.GradNe")         << rVariables.GradNe << std::endl;
    // KRATOS_INFO("rVariables.NodalCrossDirection") << rVariables.NodalCrossDirection << std::endl;
    // KRATOS_INFO("CrossEta")         << CrossEta << std::endl;
    // KRATOS_INFO("t")                << t << std::endl;

    for (unsigned int node=0; node < TNumNodes; ++node) {
        const double &Vx = rVariables.NodalCrossDirection(node, INDEX_X);
        const double &Vy = rVariables.NodalCrossDirection(node, INDEX_Y);

        // dx/d_xi
        JacobianMatrix(0, 0) += rVariables.GradNe(node, 0) * (rGeom[node].X0() + t*CrossEta[GPointCross]*Vx);

        // dy/d_xi
        JacobianMatrix(0, 1) += rVariables.GradNe(node, 0) * (rGeom[node].Y0() + t*CrossEta[GPointCross]*Vy);

        // dx/d_eta
        JacobianMatrix(1, 0) += rVariables.Nu(node) * t*Vx;

        // dy/d_eta
        JacobianMatrix(1, 1) += rVariables.Nu(node) * t*Vy;
    }

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateJacobianMatrix") << JacobianMatrix << std::endl;
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateBMatrix( Matrix &BTransformed,
                      unsigned int GPointCross,
                      const BoundedMatrix<double,TDim, TDim> &InvJ,
                      ElementVariables &rVariables ) const
{
    KRATOS_TRY
    // Details of derivation of linear part of B-Matrix can be found in:
    // "1. Geometrically non-linear formulation for the axisymmetric shell elements"
    // Note: In order to find B-Matrix, substitute Eqs.(23) & (26) into Eq.(11) to obtain Eq.(10)

    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateBMatrix") << std::endl;

    const double &t = rVariables.HalfThickness;
    const std::vector<double> CrossEta{-1.0/sqrt(3), 1.0/sqrt(3)};

    Matrix B;
    B.resize(VoigtSize, N_DOF_ELEMENT, false);
    noalias(B) = ZeroMatrix(VoigtSize, N_DOF_ELEMENT);
    const double &A11 = InvJ(0, 0);
    const double &A21 = InvJ(1, 0);
    const double &A12 = InvJ(0, 1);
    const double &A22 = InvJ(1, 1);
    const double &eta = CrossEta[GPointCross];

    for (unsigned int node=0; node < TNumNodes; ++node) {
        const unsigned int index = N_DOF_NODE * node;
        const unsigned int index_x = index + INDEX_2D_BEAM_X;
        const unsigned int index_y = index + INDEX_2D_BEAM_Y;
        const double& Fjx = - rVariables.NodalCrossDirection(node, INDEX_Y);
        const double& Fjy =   rVariables.NodalCrossDirection(node, INDEX_X);

        // Parts of B-matrix due to u_x and u_y (u, v)
        B(INDEX_2D_BEAM_XX, index_x) = rVariables.GradNe(node, 0) * A11;
        B(INDEX_2D_BEAM_YY, index_y) = rVariables.GradNe(node, 0) * A21;
        B(INDEX_2D_BEAM_XY, index_x) = B(INDEX_2D_BEAM_YY, index_y);
        B(INDEX_2D_BEAM_XY, index_y) = B(INDEX_2D_BEAM_XX, index_x);

        // Parts of B-matrix due to rotation (theta)
        const unsigned int index_t = index + INDEX_2D_BEAM_T;

        // eps_xx due to theta
        const double term_xx = A12*rVariables.Nu(node) + eta*rVariables.GradNe(node, 0)*A11;
        B(INDEX_2D_BEAM_XX, index_t) = t * Fjx * term_xx;

        // eps_yy due to theta
        const double term_yy = A22*rVariables.Nu(node) + eta*rVariables.GradNe(node, 0)*A21;
        B(INDEX_2D_BEAM_YY, index_t) =  t * Fjy * term_yy;

        // eps_xy due to theta
        B(INDEX_2D_BEAM_XY, index_t) =   t * Fjx * term_yy
                                       + t * Fjy * term_xx;
    }

    // KRATOS_INFO("B") << B << std::endl;
    // KRATOS_INFO("T") << rVariables.TransformationMatrix << std::endl;

    if (BTransformed.size1() != rVariables.TransformationMatrix.size1() ||
        BTransformed.size2() != B.size2()     )
        BTransformed.resize( rVariables.TransformationMatrix.size1(), B.size2(), false );

    // KRATOS_INFO("0-BTransformed") << BTransformed << std::endl;

    noalias(BTransformed) = prod(trans(rVariables.TransformationMatrix), B);

    // KRATOS_INFO("1-BTransformed") << BTransformed << std::endl;

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateBMatrix") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
     CalculateStrainVector(ElementVariables &rVariables) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateStrainVector") << std::endl;
    
    // KRATOS_INFO("rVariables.DofValuesVector") << rVariables.DofValuesVector << std::endl;

    noalias(rVariables.StrainVector) = prod(rVariables.B, rVariables.DofValuesVector);

    // KRATOS_INFO("1-GeoCurvedBeamElement::CalculateStrainVector") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateIntegrationCoefficient(unsigned int GPointCross,
                                    double detJ,
                                    double weight) const
{
    // KRATOS_INFO("0-GeoCurvedBeamElement::CalculateIntegrationCoefficient") << std::endl;

    const std::vector<double> CrossWeight{1.0, 1.0};

    // KRATOS_INFO("weight") << weight << std::endl;
    // KRATOS_INFO("CrossWeight[GPointCross]") << CrossWeight[GPointCross] << std::endl;
    // KRATOS_INFO("detJ") << detJ << std::endl;
    // KRATOS_INFO("IntegrationCoefficient") << weight * CrossWeight[GPointCross] * detJ << std::endl;

    return weight * CrossWeight[GPointCross] * detJ;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
SizeType GeoCurvedBeamElement<TDim,TNumNodes>::
    GetAlongNumberIntegrationPoints() const
{
    return this->GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
SizeType GeoCurvedBeamElement<TDim,TNumNodes>::
    GetCrossNumberIntegrationPoints() const
{
    return N_POINT_CROSS;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class GeoCurvedBeamElement<2,3>;
template class GeoCurvedBeamElement<3,3>;

} // Namespace Kratos
