//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


// Application includes
#include "custom_elements/one-phase_flow/U_Pl_small_strain_link_interface_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPlSmallStrainLinkInterfaceElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UPlSmallStrainLinkInterfaceElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainLinkInterfaceElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<array_1d<double,3>>& rVariable, 
                                                                                std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if(rVariable == LIQUID_FLUX_VECTOR || rVariable == CONTACT_STRESS_VECTOR || rVariable == LOCAL_STRESS_VECTOR || rVariable == LOCAL_RELATIVE_DISPLACEMENT_VECTOR || rVariable == LOCAL_LIQUID_FLUX_VECTOR)
    {
        //Variables computed on Lobatto points
        const GeometryType& Geom = this->GetGeometry();
        std::vector<array_1d<double,3>> GPValues(Geom.IntegrationPointsNumber( mThisIntegrationMethod ));

        if(rVariable == LIQUID_FLUX_VECTOR)
        {
            const PropertiesType& Prop = this->GetProperties();
            const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

            //Defining the shape functions, the jacobian and the shape functions local gradients Containers
            const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
            const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );
            GeometryType::JacobiansType JContainer(NumGPoints);
            Geom.Jacobian( JContainer, mThisIntegrationMethod );
            
            //Defining necessary variables
            array_1d<double,TNumNodes> PressureVector;
            for(unsigned int i=0; i<TNumNodes; i++)
                PressureVector[i] = Geom[i].FastGetSolutionStepValue(LIQUID_PRESSURE);
            array_1d<double,TNumNodes*TDim> DisplacementVector;
            PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
            array_1d<double,TNumNodes*TDim> VolumeAcceleration;
            PoroElementUtilities::GetNodalVariableVector(VolumeAcceleration,Geom,VOLUME_ACCELERATION);
            array_1d<double,TDim> BodyAcceleration;
            BoundedMatrix<double,TDim, TDim> RotationMatrix;
            this->CalculateRotationMatrix(RotationMatrix,Geom);
            BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
            array_1d<double,TDim> LocalRelDispVector;
            array_1d<double,TDim> RelDispVector;
            const double& InitialJointWidth = Prop[INITIAL_JOINT_WIDTH];
            BoundedMatrix<double,TNumNodes, TDim> GradNpT; 
            double JointWidth;
            BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim,TDim);
            const double& DynamicViscosityInverse = 1.0/Prop[DYNAMIC_VISCOSITY_LIQUID];
            const double& LiquidDensity = Prop[DENSITY_LIQUID];
            array_1d<double,TDim> LocalLiquidFlux;
            array_1d<double,TDim> GradPressureTerm;
            array_1d<double,TDim> LiquidFlux;
            SFGradAuxVariables SFGradAuxVars;

            //Loop over integration points
            for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
            {
                InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

                noalias(RelDispVector) = prod(Nu,DisplacementVector);
                
                noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);
                
                this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], InitialJointWidth,GPoint);

                this->template CalculateShapeFunctionsGradients< BoundedMatrix<double,TNumNodes,TDim> >(GradNpT,SFGradAuxVars,JContainer[GPoint],RotationMatrix,
                                                                                                                    DN_DeContainer[GPoint],NContainer,JointWidth,GPoint);
                
                PoroElementUtilities::InterpolateVariableWithComponents(BodyAcceleration,NContainer,VolumeAcceleration,GPoint);

                InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(LocalPermeabilityMatrix,JointWidth);
                
                noalias(GradPressureTerm) = prod(trans(GradNpT),PressureVector);
                noalias(GradPressureTerm) += -LiquidDensity*BodyAcceleration;
                
                noalias(LocalLiquidFlux) = -DynamicViscosityInverse*prod(LocalPermeabilityMatrix,GradPressureTerm);
                
                noalias(LiquidFlux) = prod(trans(RotationMatrix),LocalLiquidFlux);
                
                PoroElementUtilities::FillArray1dOutput(GPValues[GPoint],LiquidFlux);
            }
        }
        else if(rVariable == CONTACT_STRESS_VECTOR)
        {
            //Defining necessary variables
            const PropertiesType& Prop = this->GetProperties();
            const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
            array_1d<double,TNumNodes*TDim> DisplacementVector;
            PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
            BoundedMatrix<double,TDim, TDim> RotationMatrix;
            this->CalculateRotationMatrix(RotationMatrix,Geom);
            BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
            array_1d<double,TDim> RelDispVector;
            const double& InitialJointWidth = Prop[INITIAL_JOINT_WIDTH];
            double JointWidth;
            array_1d<double,TDim> LocalStressVector;
            array_1d<double,TDim> ContactStressVector;
            
            //Create constitutive law parameters:
            Vector StrainVector(TDim);
            Vector StressVectorDynamic(TDim);
            Matrix ConstitutiveMatrix(TDim,TDim);
            Vector Np(TNumNodes);
            Matrix GradNpT(TNumNodes,TDim);
            Matrix F = identity_matrix<double>(TDim);
            double detF = 1.0;
            ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
            ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
            ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
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
                
                this->CheckAndCalculateJointWidth(JointWidth, ConstitutiveParameters, StrainVector[TDim-1], InitialJointWidth, GPoint);
                
                noalias(Np) = row(NContainer,GPoint);
                
                //compute constitutive tensor and/or stresses
                mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
                
                noalias(LocalStressVector) = StressVectorDynamic;

                noalias(ContactStressVector) = prod(trans(RotationMatrix),LocalStressVector);
                
                PoroElementUtilities::FillArray1dOutput(GPValues[GPoint],ContactStressVector);
            }
        }
        else if(rVariable == LOCAL_STRESS_VECTOR)
        {
            //Defining necessary variables
            const PropertiesType& Prop = this->GetProperties();
            const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
            array_1d<double,TNumNodes*TDim> DisplacementVector;
            PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
            BoundedMatrix<double,TDim, TDim> RotationMatrix;
            this->CalculateRotationMatrix(RotationMatrix,Geom);
            BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
            array_1d<double,TDim> RelDispVector;
            const double& InitialJointWidth = Prop[INITIAL_JOINT_WIDTH];
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
            ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
            ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
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
                
                this->CheckAndCalculateJointWidth(JointWidth, ConstitutiveParameters, StrainVector[TDim-1], InitialJointWidth, GPoint);
                
                noalias(Np) = row(NContainer,GPoint);
                
                //compute constitutive tensor and/or stresses
                mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
                
                noalias(LocalStressVector) = StressVectorDynamic;
                
                PoroElementUtilities::FillArray1dOutput(GPValues[GPoint],LocalStressVector);
            }
        }
        else if(rVariable == LOCAL_RELATIVE_DISPLACEMENT_VECTOR)
        {
            //Defining necessary variables
            const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
            array_1d<double,TNumNodes*TDim> DisplacementVector;
            PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
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
                            
                PoroElementUtilities::FillArray1dOutput(GPValues[GPoint],LocalRelDispVector);
            }
        }
        else if(rVariable == LOCAL_LIQUID_FLUX_VECTOR)
        {
            const PropertiesType& Prop = this->GetProperties();
            const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

            //Defining the shape functions, the jacobian and the shape functions local gradients Containers
            const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
            const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );
            GeometryType::JacobiansType JContainer(NumGPoints);
            Geom.Jacobian( JContainer, mThisIntegrationMethod );
            
            //Defining necessary variables
            array_1d<double,TNumNodes> PressureVector;
            for(unsigned int i=0; i<TNumNodes; i++)
                PressureVector[i] = Geom[i].FastGetSolutionStepValue(LIQUID_PRESSURE);
            array_1d<double,TNumNodes*TDim> DisplacementVector;
            PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
            array_1d<double,TNumNodes*TDim> VolumeAcceleration;
            PoroElementUtilities::GetNodalVariableVector(VolumeAcceleration,Geom,VOLUME_ACCELERATION);
            array_1d<double,TDim> BodyAcceleration;
            BoundedMatrix<double,TDim, TDim> RotationMatrix;
            this->CalculateRotationMatrix(RotationMatrix,Geom);
            BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
            array_1d<double,TDim> LocalRelDispVector;
            array_1d<double,TDim> RelDispVector;
            const double& InitialJointWidth = Prop[INITIAL_JOINT_WIDTH];
            BoundedMatrix<double,TNumNodes, TDim> GradNpT; 
            double JointWidth;
            BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim,TDim);
            const double& DynamicViscosityInverse = 1.0/Prop[DYNAMIC_VISCOSITY_LIQUID];
            const double& LiquidDensity = Prop[DENSITY_LIQUID];
            array_1d<double,TDim> LocalLiquidFlux;
            array_1d<double,TDim> GradPressureTerm;
            SFGradAuxVariables SFGradAuxVars;
            
            //Loop over integration points
            for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
            {
                InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

                noalias(RelDispVector) = prod(Nu,DisplacementVector);
                
                noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);
                
                this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], InitialJointWidth,GPoint);

                this->template CalculateShapeFunctionsGradients< BoundedMatrix<double,TNumNodes,TDim> >(GradNpT,SFGradAuxVars,JContainer[GPoint],RotationMatrix,
                                                                                                                    DN_DeContainer[GPoint],NContainer,JointWidth,GPoint);
                
                PoroElementUtilities::InterpolateVariableWithComponents(BodyAcceleration,NContainer,VolumeAcceleration,GPoint);

                InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(LocalPermeabilityMatrix,JointWidth);
                
                noalias(GradPressureTerm) = prod(trans(GradNpT),PressureVector);
                noalias(GradPressureTerm) += -LiquidDensity*BodyAcceleration;
                
                noalias(LocalLiquidFlux) = -DynamicViscosityInverse*prod(LocalPermeabilityMatrix,GradPressureTerm);
                
                PoroElementUtilities::FillArray1dOutput(GPValues[GPoint],LocalLiquidFlux);
            }
        }

        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );
        if ( rOutput.size() != OutputGPoints )
            rOutput.resize( OutputGPoints );

        this->template InterpolateOutputValues< array_1d<double,3> >(rOutput,GPValues);
    }
    else
    {
        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = this->GetGeometry().IntegrationPointsNumber( this->GetIntegrationMethod() );
        if ( rOutput.size() != OutputGPoints )
            rOutput.resize( OutputGPoints );

        for(unsigned int i=0; i < OutputGPoints; i++)
        {
            noalias(rOutput[i]) = ZeroVector(3);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainLinkInterfaceElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, 
                                                                                    const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if(rVariable == PERMEABILITY_MATRIX || rVariable == LOCAL_PERMEABILITY_MATRIX)
    {
        //Variables computed on Lobatto points
        const GeometryType& Geom = this->GetGeometry();
        std::vector<Matrix> GPValues(Geom.IntegrationPointsNumber( mThisIntegrationMethod ));

        if(rVariable == PERMEABILITY_MATRIX)
        {
            const PropertiesType& Prop = this->GetProperties();

            //Defining the shape functions container
            const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
            
            //Defining necessary variables
            array_1d<double,TNumNodes*TDim> DisplacementVector;
            PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
            BoundedMatrix<double,TDim, TDim> RotationMatrix;
            this->CalculateRotationMatrix(RotationMatrix,Geom);
            BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
            array_1d<double,TDim> LocalRelDispVector;
            array_1d<double,TDim> RelDispVector;
            const double& InitialJointWidth = Prop[INITIAL_JOINT_WIDTH];
            double JointWidth;
            BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim,TDim);
            BoundedMatrix<double,TDim, TDim> PermeabilityMatrix;
        
            //Loop over integration points
            for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); GPoint++ )
            {
                InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

                noalias(RelDispVector) = prod(Nu,DisplacementVector);
                
                noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);
                
                this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], InitialJointWidth,GPoint);

                InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(LocalPermeabilityMatrix,JointWidth);
                
                noalias(PermeabilityMatrix) = prod(trans(RotationMatrix),BoundedMatrix<double,TDim, TDim>(prod(LocalPermeabilityMatrix,RotationMatrix)));
                
                GPValues[GPoint].resize(TDim,TDim,false);
                noalias(GPValues[GPoint]) = PermeabilityMatrix;
            }
        }
        else if(rVariable == LOCAL_PERMEABILITY_MATRIX)
        {
            const PropertiesType& Prop = this->GetProperties();

            //Defining the shape functions container
            const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
            
            //Defining necessary variables
            array_1d<double,TNumNodes*TDim> DisplacementVector;
            PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
            BoundedMatrix<double,TDim, TDim> RotationMatrix;
            this->CalculateRotationMatrix(RotationMatrix,Geom);
            BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
            array_1d<double,TDim> LocalRelDispVector;
            array_1d<double,TDim> RelDispVector;
            const double& InitialJointWidth = Prop[INITIAL_JOINT_WIDTH];
            double JointWidth;
            BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim,TDim);
        
            //Loop over integration points
            for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); GPoint++ )
            {
                InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

                noalias(RelDispVector) = prod(Nu,DisplacementVector);
                
                noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);
                
                this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], InitialJointWidth,GPoint);

                InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(LocalPermeabilityMatrix,JointWidth);

                GPValues[GPoint].resize(TDim,TDim,false);
                noalias(GPValues[GPoint]) = LocalPermeabilityMatrix;
            }
        }

        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );
        if ( rOutput.size() != OutputGPoints )
            rOutput.resize( OutputGPoints );

        for(unsigned int GPoint=0; GPoint<OutputGPoints; GPoint++)
            rOutput[GPoint].resize(TDim,TDim,false);

        this->template InterpolateOutputValues< Matrix >(rOutput,GPValues);
    }
    else
    {
        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = this->GetGeometry().IntegrationPointsNumber( this->GetIntegrationMethod() );
        if ( rOutput.size() != OutputGPoints )
            rOutput.resize( OutputGPoints );

        for(unsigned int i=0; i < OutputGPoints; i++)
        {
            rOutput[i].resize(TDim,TDim,false);
            noalias(rOutput[i]) = ZeroMatrix(TDim,TDim);
        }
    }
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainLinkInterfaceElement<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{    
    KRATOS_TRY
    
    //Previous definitions 
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer(NumGPoints);
    Geom.Jacobian( JContainer, mThisIntegrationMethod );
    Vector detJContainer(NumGPoints);
    Geom.DeterminantOfJacobian(detJContainer,mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    
    //Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,rCurrentProcessInfo);
    
    //Auxiliary variables
    const double& InitialJointWidth = Prop[INITIAL_JOINT_WIDTH];
    array_1d<double,TDim> RelDispVector;
    SFGradAuxVariables SFGradAuxVars;
    
    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute Np, StrainVector, JointWidth, GradNpT
        noalias(Variables.Np) = row(NContainer,GPoint);
        InterfaceElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        noalias(RelDispVector) = prod(Variables.Nu,Variables.DisplacementVector);
        noalias(Variables.StrainVector) = prod(Variables.RotationMatrix,RelDispVector);        
        this->CheckAndCalculateJointWidth(Variables.JointWidth,ConstitutiveParameters,Variables.StrainVector[TDim-1], InitialJointWidth, GPoint);
        this->template CalculateShapeFunctionsGradients< Matrix >(Variables.GradNpT,SFGradAuxVars,JContainer[GPoint],Variables.RotationMatrix,
                                                        DN_DeContainer[GPoint],NContainer,Variables.JointWidth,GPoint);
        
        //Compute BodyAcceleration and Permeability Matrix
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);
        InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(Variables.LocalPermeabilityMatrix,Variables.JointWidth);
        
        //Compute constitutive tensor and stresses
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );
        
        //Contributions to the left hand side
        this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
        
        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainLinkInterfaceElement<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{     
    KRATOS_TRY
       
    //Previous definitions 
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer(NumGPoints);
    Geom.Jacobian( JContainer, mThisIntegrationMethod );
    Vector detJContainer(NumGPoints);
    Geom.DeterminantOfJacobian(detJContainer,mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    
    //Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,rCurrentProcessInfo);
    
    //Auxiliary variables
    const double& InitialJointWidth = Prop[INITIAL_JOINT_WIDTH];
    array_1d<double,TDim> RelDispVector;
    SFGradAuxVariables SFGradAuxVars;
    
    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute Np, StrainVector, JointWidth, GradNpT
        noalias(Variables.Np) = row(NContainer,GPoint);        
        InterfaceElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        noalias(RelDispVector) = prod(Variables.Nu,Variables.DisplacementVector);
        noalias(Variables.StrainVector) = prod(Variables.RotationMatrix,RelDispVector);
        this->CheckAndCalculateJointWidth(Variables.JointWidth,ConstitutiveParameters,Variables.StrainVector[TDim-1], InitialJointWidth, GPoint);
        this->template CalculateShapeFunctionsGradients< Matrix >(Variables.GradNpT,SFGradAuxVars,JContainer[GPoint],Variables.RotationMatrix,
                                                        DN_DeContainer[GPoint],NContainer,Variables.JointWidth,GPoint);

        //Compute BodyAcceleration and Permeability Matrix
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);
        InterfaceElementUtilities::CalculateLinkPermeabilityMatrix(Variables.LocalPermeabilityMatrix,Variables.JointWidth);
        
        //Compute constitutive tensor and stresses
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );
        
        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlSmallStrainLinkInterfaceElement<2,4>;
template class UPlSmallStrainLinkInterfaceElement<3,6>;
template class UPlSmallStrainLinkInterfaceElement<3,8>;

} // Namespace Kratos
