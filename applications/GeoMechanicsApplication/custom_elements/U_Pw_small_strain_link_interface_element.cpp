// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// Application includes
#include "custom_elements/U_Pw_small_strain_link_interface_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwSmallStrainLinkInterfaceElement<TDim,TNumNodes>::
    Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UPwSmallStrainLinkInterfaceElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainLinkInterfaceElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<array_1d<double,3>>& rVariable,
                                  std::vector<array_1d<double,3>>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if (rVariable == FLUID_FLUX_VECTOR)
    {
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );

        //Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( this->GetIntegrationMethod() );
        GeometryType::JacobiansType JContainer(NumGPoints);
        Geom.Jacobian( JContainer, this->GetIntegrationMethod() );

        //Defining necessary variables
        array_1d<double,TNumNodes> PressureVector;
        for (unsigned int i=0; i<TNumNodes; i++)
            PressureVector[i] = Geom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);
        array_1d<double,TNumNodes*TDim> VolumeAcceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(VolumeAcceleration,Geom,VOLUME_ACCELERATION);
        array_1d<double,TDim> BodyAcceleration;
        BoundedMatrix<double,TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix,Geom);
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
        array_1d<double,TDim> LocalRelDispVector;
        array_1d<double,TDim> RelDispVector;
        const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
        BoundedMatrix<double,TNumNodes, TDim> GradNpT;
        double JointWidth;
        BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim,TDim);
        const double& DynamicViscosityInverse = 1.0/Prop[DYNAMIC_VISCOSITY];
        const double& FluidDensity = Prop[DENSITY_WATER];
        array_1d<double,TDim> LocalFluidFlux;
        array_1d<double,TDim> GradPressureTerm;
        array_1d<double,TDim> FluidFlux;
        SFGradAuxVariables SFGradAuxVars;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

            noalias(RelDispVector) = prod(Nu,DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);

            this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], MinimumJointWidth,GPoint);

            this->template CalculateShapeFunctionsGradients< BoundedMatrix<double,TNumNodes,TDim> >(GradNpT,SFGradAuxVars,JContainer[GPoint],RotationMatrix,
                                                                                                                DN_DeContainer[GPoint],NContainer,JointWidth,GPoint);

            GeoElementUtilities::
                InterpolateVariableWithComponents<TDim, TNumNodes>( BodyAcceleration,
                                                                    NContainer,
                                                                    VolumeAcceleration,
                                                                    GPoint );

            InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(LocalPermeabilityMatrix,JointWidth);

            noalias(GradPressureTerm) = prod(trans(GradNpT),PressureVector);
            noalias(GradPressureTerm) += PORE_PRESSURE_SIGN_FACTOR * FluidDensity*BodyAcceleration;

            noalias(LocalFluidFlux) = -DynamicViscosityInverse*prod(LocalPermeabilityMatrix,GradPressureTerm);

            noalias(FluidFlux) = prod(trans(RotationMatrix),LocalFluidFlux);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint],FluidFlux);
        }
    }
    else if (rVariable == LOCAL_STRESS_VECTOR)
    {
        //Defining necessary variables
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);
        BoundedMatrix<double,TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix,Geom);
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
        array_1d<double,TDim> RelDispVector;
        const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
        double JointWidth;
        array_1d<double,TDim> LocalStressVector;

        //Create constitutive law parameters:
        Vector StrainVector(TDim);
        Vector StressVectorDynamic(TDim);
        Matrix ConstitutiveMatrix(TDim,TDim);
        Vector Np(TNumNodes);
        Matrix GradNpT(TNumNodes,TDim);
        Matrix F = identity_matrix<double>(TDim);
        double detF = 1.0;
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveParameters.SetConstitutiveMatrix(ConstitutiveMatrix);
        ConstitutiveParameters.SetStressVector(StressVectorDynamic);
        ConstitutiveParameters.SetStrainVector(StrainVector);
        ConstitutiveParameters.SetShapeFunctionsValues(Np);
        ConstitutiveParameters.SetShapeFunctionsDerivatives(GradNpT);
        ConstitutiveParameters.SetDeterminantF(detF);
        ConstitutiveParameters.SetDeformationGradientF(F);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); GPoint++ )
        {
            InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

            noalias(RelDispVector) = prod(Nu,DisplacementVector);

            noalias(StrainVector) = prod(RotationMatrix,RelDispVector);

            this->CheckAndCalculateJointWidth(JointWidth, ConstitutiveParameters, StrainVector[TDim-1], MinimumJointWidth, GPoint);

            noalias(Np) = row(NContainer,GPoint);

            //compute constitutive tensor and/or stresses
            UPwSmallStrainInterfaceElement<TDim, TNumNodes>::UpdateElementalVariableStressVector(StressVectorDynamic, GPoint);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            noalias(LocalStressVector) = StressVectorDynamic;

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint],LocalStressVector);
        }
    }
    else if (rVariable == LOCAL_RELATIVE_DISPLACEMENT_VECTOR)
    {
        //Defining necessary variables
        const GeometryType& Geom = this->GetGeometry();
        const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);
        BoundedMatrix<double,TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix,Geom);
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
        array_1d<double,TDim> LocalRelDispVector;
        array_1d<double,TDim> RelDispVector;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); GPoint++ )
        {
            InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

            noalias(RelDispVector) = prod(Nu,DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint],LocalRelDispVector);
        }
    }
    else if (rVariable == LOCAL_FLUID_FLUX_VECTOR)
    {
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );

        //Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( this->GetIntegrationMethod() );
        GeometryType::JacobiansType JContainer(NumGPoints);
        Geom.Jacobian( JContainer, this->GetIntegrationMethod() );

        //Defining necessary variables
        array_1d<double,TNumNodes> PressureVector;
        for (unsigned int i=0; i<TNumNodes; i++)
            PressureVector[i] = Geom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);
        array_1d<double,TNumNodes*TDim> VolumeAcceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(VolumeAcceleration,Geom,VOLUME_ACCELERATION);
        array_1d<double,TDim> BodyAcceleration;
        BoundedMatrix<double,TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix,Geom);
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
        array_1d<double,TDim> LocalRelDispVector;
        array_1d<double,TDim> RelDispVector;
        const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
        BoundedMatrix<double,TNumNodes, TDim> GradNpT;
        double JointWidth;
        BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim,TDim);
        const double& DynamicViscosityInverse = 1.0/Prop[DYNAMIC_VISCOSITY];
        const double& FluidDensity = Prop[DENSITY_WATER];
        array_1d<double,TDim> LocalFluidFlux;
        array_1d<double,TDim> GradPressureTerm;
        SFGradAuxVariables SFGradAuxVars;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

            noalias(RelDispVector) = prod(Nu,DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);

            this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], MinimumJointWidth,GPoint);

            this->template CalculateShapeFunctionsGradients< BoundedMatrix<double,TNumNodes,TDim> >(GradNpT,SFGradAuxVars,JContainer[GPoint],RotationMatrix,
                                                                                                                DN_DeContainer[GPoint],NContainer,JointWidth,GPoint);

            GeoElementUtilities::
                InterpolateVariableWithComponents<TDim, TNumNodes>( BodyAcceleration,
                                                                    NContainer,
                                                                    VolumeAcceleration,
                                                                    GPoint );

            InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(LocalPermeabilityMatrix,JointWidth);

            noalias(GradPressureTerm) = prod(trans(GradNpT),PressureVector);
            noalias(GradPressureTerm) += PORE_PRESSURE_SIGN_FACTOR * FluidDensity*BodyAcceleration;

            noalias(LocalFluidFlux) = -DynamicViscosityInverse*prod(LocalPermeabilityMatrix,GradPressureTerm);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint],LocalFluidFlux);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainLinkInterfaceElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable,
                                  std::vector<Matrix>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if (rVariable == PERMEABILITY_MATRIX)
    {
        const GeometryType& Geom = this->GetGeometry();
        const PropertiesType& Prop = this->GetProperties();

        //Defining the shape functions container
        const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );

        //Defining necessary variables
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);
        BoundedMatrix<double,TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix,Geom);
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
        array_1d<double,TDim> LocalRelDispVector;
        array_1d<double,TDim> RelDispVector;
        const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
        double JointWidth;
        BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim,TDim);
        BoundedMatrix<double,TDim, TDim> PermeabilityMatrix;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); GPoint++ )
        {
            InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

            noalias(RelDispVector) = prod(Nu,DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);

            this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], MinimumJointWidth,GPoint);

            InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(LocalPermeabilityMatrix,JointWidth);

            noalias(PermeabilityMatrix) = prod(trans(RotationMatrix),BoundedMatrix<double,TDim, TDim>(prod(LocalPermeabilityMatrix,RotationMatrix)));

            rOutput[GPoint].resize(TDim,TDim,false);
            noalias(rOutput[GPoint]) = PermeabilityMatrix;
        }
    }
    else if (rVariable == LOCAL_PERMEABILITY_MATRIX)
    {
        const GeometryType& Geom = this->GetGeometry();
        const PropertiesType& Prop = this->GetProperties();

        //Defining the shape functions container
        const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );

        //Defining necessary variables
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);
        BoundedMatrix<double,TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix,Geom);
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
        array_1d<double,TDim> LocalRelDispVector;
        array_1d<double,TDim> RelDispVector;
        const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
        double JointWidth;
        BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim,TDim);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); GPoint++ )
        {
            InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

            noalias(RelDispVector) = prod(Nu,DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);

            this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], MinimumJointWidth,GPoint);

            InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(LocalPermeabilityMatrix,JointWidth);

            rOutput[GPoint].resize(TDim,TDim,false);
            noalias(rOutput[GPoint]) = LocalPermeabilityMatrix;
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainLinkInterfaceElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo,
                  const bool CalculateStiffnessMatrixFlag,
                  const bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( this->GetIntegrationMethod() );
    const unsigned int NumGPoints = IntegrationPoints.size();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( this->GetIntegrationMethod() );
    GeometryType::JacobiansType JContainer(NumGPoints);
    Geom.Jacobian( JContainer, this->GetIntegrationMethod() );
    Vector detJContainer(NumGPoints);
    Geom.DeterminantOfJacobian(detJContainer,this->GetIntegrationMethod());

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,CurrentProcessInfo);
    // if (CalculateStiffnessMatrixFlag) ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // Stiffness matrix is always needed to calculate Biot coefficient
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag)  ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);


    //Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables, ConstitutiveParameters, Geom, Prop, CurrentProcessInfo);

    //Auxiliary variables
    const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
    array_1d<double,TDim> RelDispVector;
    SFGradAuxVariables SFGradAuxVars;

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(Geom, this->GetProperties(), CurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute Np, StrainVector, JointWidth, GradNpT
        noalias(Variables.Np) = row(NContainer,GPoint);
        InterfaceElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        noalias(RelDispVector) = prod(Variables.Nu,Variables.DisplacementVector);
        noalias(Variables.StrainVector) = prod(Variables.RotationMatrix,RelDispVector);
        this->CheckAndCalculateJointWidth(Variables.JointWidth,
                                          ConstitutiveParameters,
                                          Variables.StrainVector[TDim-1],
                                          MinimumJointWidth,
                                          GPoint);

        this->template CalculateShapeFunctionsGradients< Matrix >(Variables.GradNpT,
                                                                  SFGradAuxVars,
                                                                  JContainer[GPoint],
                                                                  Variables.RotationMatrix,
                                                                  DN_DeContainer[GPoint],
                                                                  NContainer,
                                                                  Variables.JointWidth,
                                                                  GPoint);

        //Compute BodyAcceleration and Permeability Matrix
        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim, TNumNodes>( Variables.BodyAcceleration,
                                                                NContainer,
                                                                Variables.VolumeAcceleration,
                                                                GPoint );

        InterfaceElementUtilities::CalculateLinkPermeabilityMatrix( Variables.LocalPermeabilityMatrix,
                                                                    Variables.JointWidth );

        //Compute constitutive tensor and stresses
        UPwSmallStrainInterfaceElement<TDim, TNumNodes>::UpdateElementalVariableStressVector(Variables, GPoint);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
        UPwSmallStrainInterfaceElement<TDim, TNumNodes>::UpdateStressVector(Variables, GPoint);

        this->CalculateRetentionResponse( Variables,
                                          RetentionParameters,
                                          GPoint );

        // calculate Bulk modulus from stiffness matrix
        const double BulkModulus = this->CalculateBulkModulus(Variables.ConstitutiveMatrix);
        this->InitializeBiotCoefficients(Variables, BulkModulus);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient,
                                              detJContainer[GPoint],
                                              IntegrationPoints[GPoint].Weight() );

        //Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if (CalculateResidualVectorFlag)  this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------

template class UPwSmallStrainLinkInterfaceElement<2,4>;
template class UPwSmallStrainLinkInterfaceElement<3,6>;
template class UPwSmallStrainLinkInterfaceElement<3,8>;

} // Namespace Kratos
