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
    KRATOS_ERROR << "calling the default Create method for a particular element ... illegal operation!!" << std::endl;

    return Element::Pointer( new GeoCurvedBeamElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer GeoCurvedBeamElement<TDim,TNumNodes>::
    Create( IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties ) const
{
    KRATOS_ERROR << "calling the default Create method for a particular element ... illegal operation!!" << std::endl;

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

    // Base class checks for positive area and Id > 0
    int ierr = GeoStructuralBaseElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& Prop = this->GetProperties();

    if ( I33.Key() == 0 ||
         Prop.Has( I33 ) == false ||
         Prop[I33] < 0.0 )
        KRATOS_ERROR << "I33 has Key zero, is not defined or has an invalid value at element: "
                     << this->Id()
                     << std::endl;

    if ( CROSS_AREA.Key() == 0 ||
         Prop.Has( CROSS_AREA ) == false ||
         Prop[CROSS_AREA] < 0.0 )
        KRATOS_ERROR << "CROSS_AREA has Key zero, is not defined or has an invalid value at element: "
                     << this->Id()
                     << std::endl;

    if (TDim > 2) {
        if ( TORSIONAL_INERTIA.Key() == 0 ||
             Prop.Has( TORSIONAL_INERTIA ) == false ||
             Prop[TORSIONAL_INERTIA] < 0.0 )
             KRATOS_ERROR << "TORSIONAL_INERTIA has Key zero, is not defined or has an invalid value at element: "
                          << this->Id()
                          << std::endl;

        if ( I22.Key() == 0 ||
             Prop.Has( I22 ) == false ||
             Prop[I22] < 0.0 )
             KRATOS_ERROR << "I22 has Key zero, is not defined or has an invalid value at element: "
                          << this->Id()
                          << std::endl;
    }

    return 0;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<3,3>::
    SetRotationalInertiaVector(const PropertiesType& Prop, Vector& rRotationalInertia) const
{
    KRATOS_TRY

    if ( rRotationalInertia.size() != N_DOF_NODE_ROT )
        rRotationalInertia.resize( N_DOF_NODE_ROT, false );

    unsigned int index = 0;
    rRotationalInertia[index++] = Prop[TORSIONAL_INERTIA];
    rRotationalInertia[index++] = Prop[I22];
    rRotationalInertia[index++] = Prop[I33];

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<2,3>::
    SetRotationalInertiaVector(const PropertiesType& Prop, Vector& rRotationalInertia) const
{
    KRATOS_TRY

    if ( rRotationalInertia.size() != N_DOF_NODE_ROT )
        rRotationalInertia.resize( N_DOF_NODE_ROT, false );

    rRotationalInertia[0] = Prop[I33];

    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateMassMatrix( MatrixType& rMassMatrix,
                         const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //Resizing mass matrix
    if ( rMassMatrix.size1() != N_DOF_ELEMENT )
        rMassMatrix.resize( N_DOF_ELEMENT, N_DOF_ELEMENT, false );
    noalias( rMassMatrix ) = ZeroMatrix( N_DOF_ELEMENT, N_DOF_ELEMENT );

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( mThisIntegrationMethod );

    //Defining shape functions and the determinant of the jacobian at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    Vector detJContainer(IntegrationPoints.size());
    Geom.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);

    //Defining necessary variables
    const double Density = Prop[DENSITY];

    Vector RotationalInertia;
    this->SetRotationalInertiaVector(Prop, RotationalInertia);

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

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix,
                              const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    //Resizing mass matrix
    if ( rStiffnessMatrix.size1() != N_DOF_ELEMENT )
        rStiffnessMatrix.resize( N_DOF_ELEMENT, N_DOF_ELEMENT, false );
    noalias( rStiffnessMatrix ) = ZeroMatrix( N_DOF_ELEMENT, N_DOF_ELEMENT );

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( mThisIntegrationMethod );

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(IntegrationPoints.size());
    Vector detJContainer(IntegrationPoints.size());
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, detJContainer, mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom, Prop, CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,
                                     ConstitutiveParameters,
                                     Geom,
                                     Prop,
                                     CurrentProcessInfo);

    //Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPoints.size(); ++GPointAlong) {
        //Compute Np, GradNpT, B and StrainVector
        noalias(Variables.Np)      = row(NContainer, GPointAlong);
        noalias(Variables.GradNpT) = DN_DXContainer[GPointAlong];

        this->CalculateTransformationMatrix( Variables.TransformationMatrix,
                                             Variables.GradNpT );

        for (unsigned int GPointCross = 0; GPointCross < N_POINT_CROSS; ++GPointCross) {
            int GPoint= GPointAlong * N_POINT_CROSS + GPointCross;

            BoundedMatrix<double,TDim, TDim> DetJacobianMatrix;

            this->CalculateDeterminantJacobian( GPointCross,
                                                Variables,
                                                DetJacobianMatrix );

            double detJacobian;
            BoundedMatrix<double,TDim, TDim> InvertDetJacobianMatrix;
            MathUtils<double>::InvertMatrix(DetJacobianMatrix,
                                            InvertDetJacobianMatrix,
                                            detJacobian);

            this->CalculateBMatrix(Variables.B, GPointCross, InvertDetJacobianMatrix, Variables);

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

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType&
        IntegrationPoints = Geom.IntegrationPoints( mThisIntegrationMethod );

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(IntegrationPoints.size());
    Vector detJContainer(IntegrationPoints.size());
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,
                                                  detJContainer,
                                                  mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables( Variables,
                                      ConstitutiveParameters,
                                      Geom,
                                      Prop,
                                      CurrentProcessInfo );

    //Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPoints.size(); ++GPointAlong) {
        //Compute Np, GradNpT, B and StrainVector
        noalias(Variables.GradNpT) = DN_DXContainer[GPointAlong];
        noalias(Variables.Np)      = row(NContainer, GPointAlong);
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, NContainer, GPointAlong);

        this->CalculateTransformationMatrix( Variables.TransformationMatrix,
                                             Variables.GradNpT );

        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim,TNumNodes>( Variables.BodyAcceleration,
                                                               NContainer,
                                                               Variables.VolumeAcceleration,
                                                               GPointAlong );

        for (unsigned int GPointCross = 0; GPointCross < N_POINT_CROSS; ++GPointCross) {
            int GPoint= GPointAlong * N_POINT_CROSS + GPointCross;

            BoundedMatrix<double,TDim, TDim> DetJacobianMatrix;

            this->CalculateDeterminantJacobian( GPointCross,
                                                Variables,
                                                DetJacobianMatrix );

            double detJacobian;
            BoundedMatrix<double,TDim, TDim> InvertDetJacobianMatrix;
            MathUtils<double>::InvertMatrix(DetJacobianMatrix,
                                            InvertDetJacobianMatrix,
                                            detJacobian);

            this->CalculateBMatrix(Variables.B, GPointCross, InvertDetJacobianMatrix, Variables);

            this->CalculateStrainVector(Variables);

            //Compute constitutive tensor
            ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            //Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(GPointCross,
                                                                                     detJacobian,
                                                                                     IntegrationPoints[GPointAlong].Weight());
            //Contributions to the left hand side
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

            //Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    InitializeElementVariables( ElementVariables& rVariables,
                                ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                const GeometryType& Geom,
                                const PropertiesType& Prop,
                                const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY
    GeoStructuralBaseElement<TDim,TNumNodes>::
        InitializeElementVariables( rVariables,
                                    rConstitutiveParameters,
                                    Geom,
                                    Prop,
                                    CurrentProcessInfo );

    rVariables.thickness = std::sqrt(12.0 * (Prop[I33] / Prop[CROSS_AREA]));

    KRATOS_CATCH( "" )
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateRHS( VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( mThisIntegrationMethod );

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(IntegrationPoints.size());
    Vector detJContainer(IntegrationPoints.size());
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, detJContainer, mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables( Variables,
                                      ConstitutiveParameters,
                                      Geom,
                                      Prop,
                                      CurrentProcessInfo );

    //Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPoints.size(); ++GPointAlong) {
        //Compute Np, GradNpT, B and StrainVector
        noalias(Variables.GradNpT) = DN_DXContainer[GPointAlong];
        noalias(Variables.Np)      = row(NContainer, GPointAlong);
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, NContainer, GPointAlong);

        this->CalculateTransformationMatrix( Variables.TransformationMatrix,
                                             Variables.GradNpT );

        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim,TNumNodes>( Variables.BodyAcceleration,
                                                               NContainer,
                                                               Variables.VolumeAcceleration,
                                                               GPointAlong );

        for (unsigned int GPointCross = 0; GPointCross < N_POINT_CROSS; ++GPointCross) {
            int GPoint= GPointAlong * N_POINT_CROSS + GPointCross;

            BoundedMatrix<double,TDim, TDim> DetJacobianMatrix;

            this->CalculateDeterminantJacobian( GPointCross,
                                                Variables,
                                                DetJacobianMatrix );

            double detJacobian;
            BoundedMatrix<double,TDim, TDim> InvertDetJacobianMatrix;
            MathUtils<double>::InvertMatrix(DetJacobianMatrix,
                                            InvertDetJacobianMatrix,
                                            detJacobian);

            this->CalculateBMatrix(Variables.B, GPointCross, InvertDetJacobianMatrix, Variables);

            this->CalculateStrainVector(Variables);

            //Compute constitutive tensor
            ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            //Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(GPointCross,
                                                                                     detJacobian,
                                                                                     IntegrationPoints[GPointAlong].Weight());

            //Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, 
                       ElementVariables& rVariables ) const
{
    KRATOS_TRY;

    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B),rVariables.ConstitutiveMatrix);
    rLeftHandSideMatrix =  prod(rVariables.UVoigtMatrix,rVariables.B)
                         * rVariables.IntegrationCoefficient;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddRHS(VectorType& rRightHandSideVector,
                       ElementVariables& rVariables) const
{
    KRATOS_TRY;

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddBodyForce(rRightHandSideVector, rVariables);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddBodyForce(VectorType& rRightHandSideVector,
                             ElementVariables& rVariables) const
{
    KRATOS_TRY;

    const PropertiesType& Prop = this->GetProperties();
    const double density = Prop[DENSITY];

    //Distribute body force block vector into elemental vector
    rRightHandSideVector +=  density
                           * prod(trans(rVariables.Nu), rVariables.BodyAcceleration)
                           * rVariables.IntegrationCoefficient;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddStiffnessForce( VectorType& rRightHandSideVector,
                                   ElementVariables& rVariables ) const
{
    KRATOS_TRY;

    //Distribute stiffness block vector into elemental vector
    rRightHandSideVector -=  prod(trans(rVariables.B),rVariables.StressVector)
                           * rVariables.IntegrationCoefficient;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateElementCrossAngle(unsigned int GPoint,
                               const BoundedMatrix<double,TNumNodes, TNumNodes> & DN_DXContainer) const
{
    KRATOS_TRY;

    const GeometryType& Geom = this->GetGeometry();

    double dx = 0;
    double dy = 0;

    // loop over nodes
    for (unsigned int node = 0; node < TNumNodes; ++node) {
        dx += DN_DXContainer(GPoint, node) * Geom[node].X0();
        dy += DN_DXContainer(GPoint, node) * Geom[node].Y0();
    }

    return atan2(dx, -dy);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
double GeoCurvedBeamElement<2,3>::
    CalculateElementAngle(unsigned int GPoint,
                          const BoundedMatrix<double,TNumNodes, TNumNodes> &DN_DXContainer) const
{
    KRATOS_TRY;

    const GeometryType& Geom = this->GetGeometry();

    double dx = 0;
    double dy = 0;

    // loop over nodes
    for (unsigned int node = 0; node < TNumNodes; ++node) {
        dx += DN_DXContainer(GPoint, node) * Geom[node].X0();
        dy += DN_DXContainer(GPoint, node) * Geom[node].Y0();
    }

    return atan2(dy, dx);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
double GeoCurvedBeamElement<2, 3>::
    CalculateElementAngle(const Matrix &GradNpT) const
{
    KRATOS_TRY;

    const GeometryType& Geom = this->GetGeometry();

    double dx = 0;
    double dy = 0;

    // loop over nodes
    for (unsigned int node = 0; node < TNumNodes; ++node) {
        dx += GradNpT(node, 0) * Geom[node].X0();
        dy += GradNpT(node, 0) * Geom[node].Y0();
    }

    return atan2(dy, dx);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<2,3>::
    CalculateTransformationMatrix( Matrix &TransformationMatrix,
                                   const Matrix &GradNpT ) const
{
    KRATOS_TRY;

    const double phi = CalculateElementAngle(GradNpT);
    const double cosPhi = cos(phi);
    const double sinPhi = sin(phi);

    TransformationMatrix(0,0) = cosPhi * cosPhi;
    TransformationMatrix(1,1) = TransformationMatrix(0,0)

    TransformationMatrix(0,1) = sinPhi * sinPhi;
    TransformationMatrix(1,0) = TransformationMatrix(0,1)

    TransformationMatrix(2,2) =   TransformationMatrix(0,0) - TransformationMatrix(0,1);
    TransformationMatrix(0,2) =  -2.0 * sinPhi * cosPhi;
    TransformationMatrix(1,2) =  -TransformationMatrix(0,2);

    TransformationMatrix(0,2) = sinPhi * cosPhi;
    TransformationMatrix(1,2) = - TransformationMatrix(0,2);


    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<3,3>::
    CalculateTransformationMatrix( Matrix &TransformationMatrix,
                                   const Matrix &GradNpT ) const
{
    KRATOS_TRY;

    KRATOS_ERROR << "Undefined dimension in CalculateTransformationMatrix" << std::endl;;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateCrossDirection( Matrix &CrossDirection ) const
{
    KRATOS_TRY;

    BoundedMatrix<double, TNumNodes, TNumNodes> DN_DXContainer;
    GeoElementUtilities::CalculateShapeFunctionsNodesGradients(DN_DXContainer);

    for (unsigned int IntegrationNode = 0; IntegrationNode < TNumNodes; ++IntegrationNode) {
        double phi = CalculateElementCrossAngle(IntegrationNode, DN_DXContainer);
        CrossDirection(0, IntegrationNode) = cos(phi);
        CrossDirection(1, IntegrationNode) = sin(phi);
    }

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
     CalculateDeterminantJacobian(unsigned int GPointCross,
                                  const ElementVariables &rVariables,
                                  BoundedMatrix<double,TDim, TDim> &DeterminantJacobian) const
{
    KRATOS_TRY;
    const GeometryType& Geom = this->GetGeometry();

    const double thick = 0.5 * rVariables.thickness;
    const std::vector<double> CrossXi{-1.0/sqrt(3), 1.0/sqrt(3)};

    noalias(DeterminantJacobian) = ZeroMatrix(TDim, TDim);

    for (unsigned int node=0; node < TNumNodes; ++node) {
        DeterminantJacobian(0, 0) += rVariables.GradNpT(node, 0) * Geom[node].X0() + thick * CrossXi[GPointCross] * rVariables.CrossDirection(0, node);
        DeterminantJacobian(0, 1) += rVariables.GradNpT(node, 0) * Geom[node].Y0() + thick * CrossXi[GPointCross] * rVariables.CrossDirection(1, node);

        DeterminantJacobian(1, 0) += rVariables.Np(node) * thick * rVariables.CrossDirection(0, node);
        DeterminantJacobian(1, 1) += rVariables.Np(node) * thick * rVariables.CrossDirection(1, node);
    }

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateBMatrix( Matrix &BTransformed,
                      unsigned int GPointCross,
                      const BoundedMatrix<double,TDim, TDim> &InvertDetJacobian,
                      ElementVariables &rVariables ) const
{
    KRATOS_TRY

    const double thick = 0.5 * rVariables.thickness;
    const std::vector<double> CrossXi{-1.0/sqrt(3), 1.0/sqrt(3)};

    Matrix B;
    B.resize(VoigtSize, N_DOF_ELEMENT, false);
    noalias(B) = ZeroMatrix(VoigtSize, N_DOF_ELEMENT);

    for (unsigned int node=0; node < TNumNodes; ++node)
    {
        const unsigned int index_0 = N_DOF_NODE * node;
        const unsigned int index_1 = index_0 + 1;
        const unsigned int index_2 = index_1 + 1;

        B(0, index_0) = rVariables.GradNpT(node, 0) * InvertDetJacobian(0, 0);
        B(1, index_1) = rVariables.GradNpT(node, 0) * InvertDetJacobian(1, 0);
        B(2, index_0) = B(1, index_1);
        B(2, index_1) = B(0, index_0);


        double term_1 = InvertDetJacobian(0,1)*rVariables.Np(node) + CrossXi[GPointCross] * B(0, index_0);
        B(0, index_2) = -thick * rVariables.CrossDirection(1, node) * term_1;

        double term_2 = InvertDetJacobian(1,1)*rVariables.Np(node) + CrossXi[GPointCross] * B(1,index_1);
        B(1, index_2) =  thick * rVariables.CrossDirection(0, node) * term_2;

        term_1 = InvertDetJacobian(1,1)*rVariables.Np(node) + CrossXi[GPointCross] * B(1, index_1);
        term_2 = InvertDetJacobian(0,1)*rVariables.Np(node) + CrossXi[GPointCross] * B(0, index_0);

        B(2, index_2) = - thick * rVariables.CrossDirection(1, node) * term_1
                        + thick * rVariables.CrossDirection(0, node) * term_2;
    }

    //noalias(BTransformed) = prod(rVariables.TransformationMatrix, B);
    noalias(BTransformed) = ZeroMatrix(BTransformed.size1(), BTransformed.size2());
    for (unsigned int idim=0; idim < TDim; ++idim) {
        for (unsigned int idof=0; idof < N_DOF_ELEMENT; ++idof) {
            for (unsigned int jvoigt=0; jvoigt < VoigtSize; ++jvoigt) {
                BTransformed(idim, idof) += rVariables.TransformationMatrix(idim, jvoigt) * B(jvoigt,idof);
            }
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
     CalculateStrainVector(ElementVariables &rVariables) const
{
    KRATOS_TRY

    noalias(rVariables.StrainVector) = prod(rVariables.B, rVariables.DofValuesVector);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateIntegrationCoefficient(unsigned int GPointCross,
                                    double detJ,
                                    double weight) const
{
    return weight * CrossWeight[GPointCross] * detJ;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class GeoCurvedBeamElement<2,3>;
template class GeoCurvedBeamElement<3,3>;

} // Namespace Kratos
