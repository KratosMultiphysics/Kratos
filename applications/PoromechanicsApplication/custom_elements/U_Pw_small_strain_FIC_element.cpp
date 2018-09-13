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
#include "custom_elements/U_Pw_small_strain_FIC_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwSmallStrainFICElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UPwSmallStrainFICElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwSmallStrainFICElement<TDim,TNumNodes>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwSmallStrainFICElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::Initialize()
{
    KRATOS_TRY
    
    UPwElement<TDim,TNumNodes>::Initialize();
    
    unsigned int VoigtSize = 6;
    if(TDim == 2) VoigtSize = 3;
        
    for(unsigned int i = 0; i < TDim; i++)
    {
        mNodalConstitutiveTensor[i].resize(VoigtSize);
        
        for(unsigned int j = 0; j < VoigtSize; j++)
        {
            for(unsigned int k = 0; k < TNumNodes; k++)
                mNodalConstitutiveTensor[i][j][k] = 0.0;
        }
    }
    
    for(unsigned int i = 0; i < TDim; i++)
    {
        for(unsigned int j = 0; j < TNumNodes; j++)
            mNodalDtStress[i][j] = 0.0;
    }
        
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);
    
    unsigned int VoigtSize = 6;
    if(TDim == 2) VoigtSize = 3;
    Matrix B(VoigtSize,TNumNodes*TDim);
    noalias(B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
    array_1d<double,TNumNodes*TDim> VelocityVector;
    PoroElementUtilities::GetNodalVariableVector(VelocityVector,Geom,VELOCITY);

    //Create constitutive law parameters:
    Vector StrainVector(VoigtSize);
    Vector StressVector(VoigtSize);
    Matrix ConstitutiveMatrix(VoigtSize,VoigtSize);
    Vector Np(TNumNodes);
    Matrix GradNpT(TNumNodes,TDim);
    Matrix F = identity_matrix<double>(TDim);
    double detF = 1.0;
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,this->GetProperties(),rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); //Note: this is for nonlocal damage
    ConstitutiveParameters.SetConstitutiveMatrix(ConstitutiveMatrix);
    ConstitutiveParameters.SetStressVector(StressVector);
    ConstitutiveParameters.SetStrainVector(StrainVector);
    ConstitutiveParameters.SetShapeFunctionsValues(Np);
    ConstitutiveParameters.SetShapeFunctionsDerivatives(GradNpT);
    ConstitutiveParameters.SetDeformationGradientF(F);
    ConstitutiveParameters.SetDeterminantF(detF);
    
    //Extrapolation variables
    array_1d<Matrix,TDim> ConstitutiveTensorContainer;
    for(unsigned int i = 0; i < TDim; i++)
    {
        ConstitutiveTensorContainer[i].resize(NumGPoints,VoigtSize,false);
    }
    Matrix DtStressContainer(NumGPoints,TDim);
    
    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        noalias(Np) = row(NContainer,GPoint);
        noalias(GradNpT) = DN_DXContainer[GPoint];
        this->CalculateBMatrix(B, GradNpT);
        
        // Compute ConstitutiveTensor
        noalias(StrainVector) = prod(B,DisplacementVector);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
        this->SaveGPConstitutiveTensor(ConstitutiveTensorContainer,ConstitutiveMatrix,GPoint);
        
        // Compute DtStress
        noalias(StrainVector) = prod(B,VelocityVector);
        noalias(StressVector) = prod(ConstitutiveMatrix,StrainVector);
        this->SaveGPDtStress(DtStressContainer,StressVector,GPoint);
    }
    
    this->ExtrapolateGPConstitutiveTensor(ConstitutiveTensorContainer);
    this->ExtrapolateGPDtStress(DtStressContainer);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);
    
    unsigned int VoigtSize = 6;
    if(TDim == 2) VoigtSize = 3;
    Matrix B(VoigtSize,TNumNodes*TDim);
    noalias(B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
    array_1d<double,TNumNodes*TDim> VelocityVector;
    PoroElementUtilities::GetNodalVariableVector(VelocityVector,Geom,VELOCITY);

    //Create constitutive law parameters:
    Vector StrainVector(VoigtSize);
    Vector StressVector(VoigtSize);
    Matrix ConstitutiveMatrix(VoigtSize,VoigtSize);
    Vector Np(TNumNodes);
    Matrix GradNpT(TNumNodes,TDim);
    Matrix F = identity_matrix<double>(TDim);
    double detF = 1.0;
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,this->GetProperties(),rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); //Note: this is for nonlocal damage
    ConstitutiveParameters.SetConstitutiveMatrix(ConstitutiveMatrix);
    ConstitutiveParameters.SetStressVector(StressVector);
    ConstitutiveParameters.SetStrainVector(StrainVector);
    ConstitutiveParameters.SetShapeFunctionsValues(Np);
    ConstitutiveParameters.SetShapeFunctionsDerivatives(GradNpT);
    ConstitutiveParameters.SetDeformationGradientF(F);
    ConstitutiveParameters.SetDeterminantF(detF);
    
    //Containers for extrapolation variables
    Matrix DtStressContainer(NumGPoints,TDim);
        
    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        noalias(Np) = row(NContainer,GPoint);
        noalias(GradNpT) = DN_DXContainer[GPoint];
        this->CalculateBMatrix(B, GradNpT);
        
        // Compute ConstitutiveTensor
        noalias(StrainVector) = prod(B,DisplacementVector);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
        
        // Compute DtStress
        noalias(StrainVector) = prod(B,VelocityVector);
        noalias(StressVector) = prod(ConstitutiveMatrix,StrainVector);
        this->SaveGPDtStress(DtStressContainer,StressVector,GPoint);
    }
    
    this->ExtrapolateGPDtStress(DtStressContainer);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::SaveGPConstitutiveTensor(array_1d<Matrix,TDim>& rConstitutiveTensorContainer, const Matrix& ConstitutiveMatrix, const unsigned int& GPoint)
{
    for(unsigned int i = 0; i < TDim; i++)
    {
        for(unsigned int j = 0; j < ConstitutiveMatrix.size1(); j++) //VoigtSize
        {
            rConstitutiveTensorContainer[i](GPoint,j) = ConstitutiveMatrix(i,j);
        }
    }

    /* INFO: (Quadrilateral_2D_4 with GI_GAUSS_2)
     * 
     *                                ( |D00-0 D01-0 D02-0|   |D10-0 D11-0 D12-0| )
     * rConstitutiveTensorContainer = ( |D00-1 D01-1 D02-1|   |D10-1 D11-1 D12-1| )
     *                                ( |D00-2 D01-2 D02-2|   |D10-2 D11-2 D12-2| )
     *                                ( |D00-3 D01-3 D02-3| , |D10-3 D11-3 D12-3| )
     * 
     * D00-0 = D(0,0) at GP 0
    */
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::SaveGPDtStress(Matrix& rDtStressContainer, const Vector& StressVector, const unsigned int& GPoint)
{
    for(unsigned int i = 0; i < TDim; i++)
    {
        rDtStressContainer(GPoint,i) = StressVector[i];
    }

    /* INFO: (Quadrilateral_2D_4 with GI_GAUSS_2)
     * 
     *                      |S0-0 S1-0|
     * rDtStressContainer = |S0-1 S1-1|
     *                      |S0-2 S1-2|
     *                      |S0-3 S1-3|
     * 
     * S0-0 = S[0] at GP 0
    */
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,3>::ExtrapolateGPConstitutiveTensor(const array_1d<Matrix,2>& ConstitutiveTensorContainer)
{
    // Triangle_2d_3 with GI_GAUSS_1
    
    for(unsigned int i = 0; i < 2; i++) //TDim
    {
        for(unsigned int j = 0; j < 3; j++) // VoigtSize
        {
            for(unsigned int k = 0; k < 3; k++) // TNumNodes
                mNodalConstitutiveTensor[i][j][k] = ConstitutiveTensorContainer[i](0,j);
        }
    }
    
    /* INFO:
     * 
     *                            [ ( |D00-0|   |D01-0|   |D02-0| )   ( |D10-0|   |D11-0|   |D12-0| ) ]
     * mNodalConstitutiveTensor = [ ( |D00-1|   |D01-1|   |D02-1| )   ( |D10-1|   |D11-1|   |D12-1| ) ]
     *                            [ ( |D00-2| , |D01-2| , |D02-2| ) , ( |D10-2| , |D11-2| , |D12-2| ) ]
     * 
     * D00-0 = D(0,0) at node 0
    */
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,4>::ExtrapolateGPConstitutiveTensor(const array_1d<Matrix,2>& ConstitutiveTensorContainer)
{
    // Quadrilateral_2d_4 with GI_GAUSS_2
    
    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    PoroElementUtilities::CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,4,3> AuxNodalConstitutiveTensor;
    
    for(unsigned int i = 0; i < 2; i++) //TDim
    {
        noalias(AuxNodalConstitutiveTensor) = prod(ExtrapolationMatrix,ConstitutiveTensorContainer[i]);
        
        for(unsigned int j = 0; j < 3; j++) // VoigtSize
            noalias(mNodalConstitutiveTensor[i][j]) = column(AuxNodalConstitutiveTensor,j);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,4>::ExtrapolateGPConstitutiveTensor(const array_1d<Matrix,3>& ConstitutiveTensorContainer)
{
    // Tetrahedra_3d_4 with GI_GAUSS_1
    
    for(unsigned int i = 0; i < 3; i++) //TDim
    {
        for(unsigned int j = 0; j < 6; j++) // VoigtSize
        {
            for(unsigned int k = 0; k < 4; k++) // TNumNodes
                mNodalConstitutiveTensor[i][j][k] = ConstitutiveTensorContainer[i](0,j);
        }
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,8>::ExtrapolateGPConstitutiveTensor(const array_1d<Matrix,3>& ConstitutiveTensorContainer)
{
    // Hexahedra_3d_8 with GI_GAUSS_2
    
    BoundedMatrix<double,8,8> ExtrapolationMatrix;
    PoroElementUtilities::CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,8,6> AuxNodalConstitutiveTensor;
    
    for(unsigned int i = 0; i < 3; i++) //TDim
    {
        noalias(AuxNodalConstitutiveTensor) = prod(ExtrapolationMatrix,ConstitutiveTensorContainer[i]);
        
        for(unsigned int j = 0; j < 6; j++) // VoigtSize
            noalias(mNodalConstitutiveTensor[i][j]) = column(AuxNodalConstitutiveTensor,j);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,3>::ExtrapolateGPDtStress(const Matrix& DtStressContainer)
{
    // Triangle_2d_3 with GI_GAUSS_1
        
    for(unsigned int i = 0; i < 2; i++) //TDim
    {
        for(unsigned int j = 0; j < 3; j++) // TNumNodes
            mNodalDtStress[i][j] = DtStressContainer(0,i);
    }
    
    /* INFO:
     * 
     *                  ( |S0-0|   |S1-0| )
     * mNodalDtStress = ( |S0-1|   |S1-1| )
     *                  ( |S0-2| , |S1-2| )
     * 
     * S0-0 = S[0] at node 0
    */
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,4>::ExtrapolateGPDtStress(const Matrix& DtStressContainer)
{
    // Quadrilateral_2d_4 with GI_GAUSS_2

    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    PoroElementUtilities::CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,4,2> AuxNodalDtStress;
    noalias(AuxNodalDtStress) = prod(ExtrapolationMatrix,DtStressContainer);
    
    for(unsigned int i = 0; i < 2; i++) // TDim
        noalias(mNodalDtStress[i]) = column(AuxNodalDtStress,i);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,4>::ExtrapolateGPDtStress(const Matrix& DtStressContainer)
{
    // Tetrahedra_3d_4 with GI_GAUSS_1
    
    for(unsigned int i = 0; i < 3; i++) //TDim
    {
        for(unsigned int j = 0; j < 4; j++) // TNumNodes
            mNodalDtStress[i][j] = DtStressContainer(0,i);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,8>::ExtrapolateGPDtStress(const Matrix& DtStressContainer)
{
    // Hexahedra_3d_8 with GI_GAUSS_2
    
    BoundedMatrix<double,8,8> ExtrapolationMatrix;
    PoroElementUtilities::CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,8,3> AuxNodalDtStress;
    noalias(AuxNodalDtStress) = prod(ExtrapolationMatrix,DtStressContainer);
    
    for(unsigned int i = 0; i < 3; i++) // TDim
        noalias(mNodalDtStress[i]) = column(AuxNodalDtStress,i);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    
    //Previous definitions 
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Vector detJContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    
    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,CurrentProcessInfo);
    
    FICElementVariables FICVariables;
    this->InitializeFICElementVariables(FICVariables,DN_DXContainer,Geom,Prop,CurrentProcessInfo);
    
    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        noalias(Variables.GradNpT) = DN_DXContainer[GPoint];
        this->CalculateBMatrix(Variables.B, Variables.GradNpT);
        noalias(Variables.StrainVector) = prod(Variables.B,Variables.DisplacementVector);
        
        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

        //Compute ShapeFunctionsSecondOrderGradients
        this->CalculateShapeFunctionsSecondOrderGradients(FICVariables,Variables);

        //Compute constitutive tensor and stresses
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
        
        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient( Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );
        
        //Contributions to the left hand side
        this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
        
        this->CalculateAndAddLHSStabilization(rLeftHandSideMatrix, Variables, FICVariables);
    
        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    
        this->CalculateAndAddRHSStabilization(rRightHandSideVector, Variables, FICVariables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{     
    KRATOS_TRY
       
    //Previous definitions 
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Vector detJContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    
    //Element variables
    ElementVariables Variables; 
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,CurrentProcessInfo);

    FICElementVariables FICVariables;
    this->InitializeFICElementVariables(FICVariables,DN_DXContainer,Geom,Prop,CurrentProcessInfo);
    
    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        noalias(Variables.GradNpT) = DN_DXContainer[GPoint];
        this->CalculateBMatrix(Variables.B, Variables.GradNpT);
        noalias(Variables.StrainVector) = prod(Variables.B,Variables.DisplacementVector);
        
        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

        //Compute ShapeFunctionsSecondOrderGradients
        this->CalculateShapeFunctionsSecondOrderGradients(FICVariables,Variables);
        
        //Compute stresses
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );
                
        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
        
        this->CalculateAndAddRHSStabilization(rRightHandSideVector, Variables, FICVariables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::InitializeFICElementVariables(FICElementVariables& rFICVariables,const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer,
                                                                                const GeometryType& Geom,const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo)
{   
    KRATOS_TRY
    
    //Properties variables    
    rFICVariables.ShearModulus = Prop[YOUNG_MODULUS]/(2.0*(1.0+Prop[POISSON_RATIO]));   
    
    //Nodal Variables
    this->ExtrapolateShapeFunctionsGradients(rFICVariables.NodalShapeFunctionsGradients,DN_DXContainer);
    
    //General Variables
    this->CalculateElementLength(rFICVariables.ElementLength,Geom);
    
    //Variables computed at each GP
    this->InitializeSecondOrderTerms(rFICVariables);
        
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,3>::ExtrapolateShapeFunctionsGradients(array_1d< array_1d<double,6> , 3 >& rNodalShapeFunctionsGradients,  
                                                                        const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer)
{
    // Triangle_2d_3 with GI_GAUSS_1
    // No necessary
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,4>::ExtrapolateShapeFunctionsGradients(array_1d< array_1d<double,8> , 4 >& rNodalShapeFunctionsGradients, 
                                                                        const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer)
{
    // Quadrilateral_2d_4 with GI_GAUSS_2

    BoundedMatrix<double,4,8> ShapeFunctionsGradientsContainer; //NumGPoints X TDim*TNumNodes
    unsigned int index;
    
    for(unsigned int i = 0; i < 4; i++) //NumGPoints
    {
        for(unsigned int j = 0; j < 4; j++) //TNumNodes
        {
            index = j*2;
            
            ShapeFunctionsGradientsContainer(i,index) = DN_DXContainer[i](j,0);
            ShapeFunctionsGradientsContainer(i,index+1) = DN_DXContainer[i](j,1);
        }
    }
    
    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    PoroElementUtilities::CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,4,8> AuxNodalShapeFunctionsGradients;
    noalias(AuxNodalShapeFunctionsGradients) = prod(ExtrapolationMatrix,ShapeFunctionsGradientsContainer);
    
    for(unsigned int i = 0; i < 4; i++) //TNumNodes
    {
        index = i*2;
        
        rNodalShapeFunctionsGradients[i][0] = AuxNodalShapeFunctionsGradients(0,index);
        rNodalShapeFunctionsGradients[i][1] = AuxNodalShapeFunctionsGradients(0,index+1);
        rNodalShapeFunctionsGradients[i][2] = AuxNodalShapeFunctionsGradients(1,index);
        rNodalShapeFunctionsGradients[i][3] = AuxNodalShapeFunctionsGradients(1,index+1);
        rNodalShapeFunctionsGradients[i][4] = AuxNodalShapeFunctionsGradients(2,index);
        rNodalShapeFunctionsGradients[i][5] = AuxNodalShapeFunctionsGradients(2,index+1);
        rNodalShapeFunctionsGradients[i][6] = AuxNodalShapeFunctionsGradients(3,index);
        rNodalShapeFunctionsGradients[i][7] = AuxNodalShapeFunctionsGradients(3,index+1);
    }

    /* INFO:
     * 
     *                                 ( |N0x-0|   |N1x-0|   |N2x-0|   |N3x-0| )
     *                                 ( |N0y-0|   |N1y-0|   |N2y-0|   |N3y-0| )
     *                                 ( |N0x-1|   |N1x-1|   |N2x-1|   |N3x-1| )
     * rNodalShapeFunctionsGradients = ( |N0y-1|   |N1y-1|   |N2y-1|   |N3y-1| )
     *                                 ( |N0x-2|   |N1x-2|   |N2x-2|   |N3x-2| )
     *                                 ( |N0y-2|   |N1y-2|   |N2y-2|   |N3y-2| )
     *                                 ( |N0x-3|   |N1x-3|   |N2x-3|   |N3x-3| )
     *                                 ( |N0y-3| , |N1y-3| , |N2y-3| , |N3y-3| )
     * 
     * N0x-0 = aN0/ax at node 0
    */
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,4>::ExtrapolateShapeFunctionsGradients(array_1d< array_1d<double,12> , 4 >& rNodalShapeFunctionsGradients, 
                                                                        const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer)
{
    // Tetrahedra_3d_4 with GI_GAUSS_1
    // No necessary
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,8>::ExtrapolateShapeFunctionsGradients(array_1d< array_1d<double,24> , 8 >& rNodalShapeFunctionsGradients, 
                                                                        const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer)
{
    // Hexahedra_3d_8 with GI_GAUSS_2

    BoundedMatrix<double,8,24> ShapeFunctionsGradientsContainer; //NumGPoints X TDim*TNumNodes
    unsigned int index;
    
    for(unsigned int i = 0; i < 8; i++) //NumGPoints
    {
        for(unsigned int j = 0; j < 8; j++) //TNumNodes
        {
            index = j*3;
            
            ShapeFunctionsGradientsContainer(i,index) = DN_DXContainer[i](j,0);
            ShapeFunctionsGradientsContainer(i,index+1) = DN_DXContainer[i](j,1);
            ShapeFunctionsGradientsContainer(i,index+2) = DN_DXContainer[i](j,2);
        }
    }
    
    BoundedMatrix<double,8,8> ExtrapolationMatrix;
    PoroElementUtilities::CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,8,24> AuxNodalShapeFunctionsGradients;
    noalias(AuxNodalShapeFunctionsGradients) = prod(ExtrapolationMatrix,ShapeFunctionsGradientsContainer);
    
    for(unsigned int i = 0; i < 8; i++) //TNumNodes
    {
        index = i*3;
        
        rNodalShapeFunctionsGradients[i][0] = AuxNodalShapeFunctionsGradients(0,index);
        rNodalShapeFunctionsGradients[i][1] = AuxNodalShapeFunctionsGradients(0,index+1);
        rNodalShapeFunctionsGradients[i][2] = AuxNodalShapeFunctionsGradients(0,index+2);
        rNodalShapeFunctionsGradients[i][3] = AuxNodalShapeFunctionsGradients(1,index);
        rNodalShapeFunctionsGradients[i][4] = AuxNodalShapeFunctionsGradients(1,index+1);
        rNodalShapeFunctionsGradients[i][5] = AuxNodalShapeFunctionsGradients(1,index+2);
        rNodalShapeFunctionsGradients[i][6] = AuxNodalShapeFunctionsGradients(2,index);
        rNodalShapeFunctionsGradients[i][7] = AuxNodalShapeFunctionsGradients(2,index+1);
        rNodalShapeFunctionsGradients[i][8] = AuxNodalShapeFunctionsGradients(2,index+2);
        rNodalShapeFunctionsGradients[i][9] = AuxNodalShapeFunctionsGradients(3,index);
        rNodalShapeFunctionsGradients[i][10] = AuxNodalShapeFunctionsGradients(3,index+1);
        rNodalShapeFunctionsGradients[i][11] = AuxNodalShapeFunctionsGradients(3,index+2);
        rNodalShapeFunctionsGradients[i][12] = AuxNodalShapeFunctionsGradients(4,index);
        rNodalShapeFunctionsGradients[i][13] = AuxNodalShapeFunctionsGradients(4,index+1);
        rNodalShapeFunctionsGradients[i][14] = AuxNodalShapeFunctionsGradients(4,index+2);
        rNodalShapeFunctionsGradients[i][15] = AuxNodalShapeFunctionsGradients(5,index);
        rNodalShapeFunctionsGradients[i][16] = AuxNodalShapeFunctionsGradients(5,index+1);
        rNodalShapeFunctionsGradients[i][17] = AuxNodalShapeFunctionsGradients(5,index+2);
        rNodalShapeFunctionsGradients[i][18] = AuxNodalShapeFunctionsGradients(6,index);
        rNodalShapeFunctionsGradients[i][19] = AuxNodalShapeFunctionsGradients(6,index+1);
        rNodalShapeFunctionsGradients[i][20] = AuxNodalShapeFunctionsGradients(6,index+2);
        rNodalShapeFunctionsGradients[i][21] = AuxNodalShapeFunctionsGradients(7,index);
        rNodalShapeFunctionsGradients[i][22] = AuxNodalShapeFunctionsGradients(7,index+1);
        rNodalShapeFunctionsGradients[i][23] = AuxNodalShapeFunctionsGradients(7,index+2);
    }
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,3>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = sqrt(4.0*Geom.Area()/Globals::Pi);
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,4>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = sqrt(4.0*Geom.Area()/Globals::Pi);
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,4>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = pow( (6.0*Geom.Volume()/Globals::Pi) , (1.0/3.0) );
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,8>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = pow( (6.0*Geom.Volume()/Globals::Pi) , (1.0/3.0) );
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,3>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    for(unsigned int i = 0; i < 2; i++) //TDim
        rFICVariables.ConstitutiveTensorGradients[i].resize(3); //VoigtSize
    
    rFICVariables.DimVoigtMatrix.resize(2,3,false); //TDim X VoigtSize
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,4>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    //Voigt identity matrix
    rFICVariables.VoigtMatrix.resize(3,3,false); //VoigtSize X VoigtSize
    noalias(rFICVariables.VoigtMatrix) = ZeroMatrix(3,3);
    rFICVariables.VoigtMatrix(0,0) = 1.0;
    rFICVariables.VoigtMatrix(1,1) = 1.0;
    rFICVariables.VoigtMatrix(2,2) = 0.5;

    for(unsigned int i = 0; i < 4; i++) //TNumNodes
        rFICVariables.ShapeFunctionsSecondOrderGradients[i].resize(3,false); //VoigtSize

    for(unsigned int i = 0; i < 2; i++) //TDim
        rFICVariables.ConstitutiveTensorGradients[i].resize(3); //VoigtSize
    
    rFICVariables.DimVoigtMatrix.resize(2,3,false); //TDim X VoigtSize
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,4>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    for(unsigned int i = 0; i < 3; i++) //TDim
        rFICVariables.ConstitutiveTensorGradients[i].resize(6); //VoigtSize
    
    rFICVariables.DimVoigtMatrix.resize(3,6,false); //TDim X VoigtSize
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,8>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    //Voigt identity matrix
    rFICVariables.VoigtMatrix.resize(6,6,false); //VoigtSize X VoigtSize
    noalias(rFICVariables.VoigtMatrix) = ZeroMatrix(6,6);
    rFICVariables.VoigtMatrix(0,0) = 1.0;
    rFICVariables.VoigtMatrix(1,1) = 1.0;
    rFICVariables.VoigtMatrix(2,2) = 1.0;
    rFICVariables.VoigtMatrix(3,3) = 0.5;
    rFICVariables.VoigtMatrix(4,4) = 0.5;
    rFICVariables.VoigtMatrix(5,5) = 0.5;

    for(unsigned int i = 0; i < 8; i++) //TNumNodes
        rFICVariables.ShapeFunctionsSecondOrderGradients[i].resize(6,false); //VoigtSize

    for(unsigned int i = 0; i < 3; i++) //TDim
        rFICVariables.ConstitutiveTensorGradients[i].resize(6); //VoigtSize
    
    rFICVariables.DimVoigtMatrix.resize(3,6,false); //TDim X VoigtSize
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,3>::CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables, ElementVariables& rVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,4>::CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables, ElementVariables& rVariables)
{
    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B),rFICVariables.VoigtMatrix);
    unsigned int index;
    for(unsigned int i = 0; i < 4; i++) //TNumNodes
    {
        index = 2*i;
        
        noalias(rFICVariables.ShapeFunctionsSecondOrderGradients[i]) = prod(trans(rVariables.UVoigtMatrix),rFICVariables.NodalShapeFunctionsGradients[i]);
        
        rFICVariables.StrainGradients(0, index)   = rFICVariables.ShapeFunctionsSecondOrderGradients[i][0]+0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][1];
        rFICVariables.StrainGradients(1, index+1) = 0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][0]+rFICVariables.ShapeFunctionsSecondOrderGradients[i][1];
        rFICVariables.StrainGradients(0, index+1) = 0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][2];
        rFICVariables.StrainGradients(1, index)   = rFICVariables.StrainGradients(0, index+1);
    }

    /* INFO:
     * 
     *                                                    ( |N0xx|   |N1xx|   |N2xx|   |N3xx| ) 
     * rFICVariables.ShapeFunctionsSecondOrderGradients = ( |N0yy|   |N1yy|   |N2yy|   |N3yy| )
     *                                                    ( |N0xy| , |N1xy| , |N2xy| , |N3xy| )
     * 
     * N0xx = a2N0/ax2 at current GP
    */
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,4>::CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables, ElementVariables& rVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,8>::CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables, ElementVariables& rVariables)
{
    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B),rFICVariables.VoigtMatrix);
    unsigned int index;
    for(unsigned int i = 0; i < 8; i++) //TNumNodes
    {
        index = 3*i;
        
        noalias(rFICVariables.ShapeFunctionsSecondOrderGradients[i]) = prod(trans(rVariables.UVoigtMatrix),rFICVariables.NodalShapeFunctionsGradients[i]);
        
        rFICVariables.StrainGradients(0, index)   = rFICVariables.ShapeFunctionsSecondOrderGradients[i][0]+0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][1]+0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][2];
        rFICVariables.StrainGradients(1, index+1) = 0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][0]+rFICVariables.ShapeFunctionsSecondOrderGradients[i][1]+0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][2];
        rFICVariables.StrainGradients(2, index+2) = 0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][0]+0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][1]+rFICVariables.ShapeFunctionsSecondOrderGradients[i][2];
        rFICVariables.StrainGradients(0, index+1) = 0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][3];
        rFICVariables.StrainGradients(1, index)   = rFICVariables.StrainGradients(0, index+1);
        rFICVariables.StrainGradients(1, index+2) = 0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][4];
        rFICVariables.StrainGradients(2, index+1) = rFICVariables.StrainGradients(1, index+2);
        rFICVariables.StrainGradients(0, index+2) = 0.5*rFICVariables.ShapeFunctionsSecondOrderGradients[i][5];
        rFICVariables.StrainGradients(2, index)   = rFICVariables.StrainGradients(0, index+2);
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::CalculateAndAddLHSStabilization(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    this->CalculateAndAddStrainGradientMatrix(rLeftHandSideMatrix,rVariables,rFICVariables);

    this->CalculateAndAddDtStressGradientMatrix(rLeftHandSideMatrix,rVariables,rFICVariables);

    this->CalculateAndAddPressureGradientMatrix(rLeftHandSideMatrix,rVariables,rFICVariables);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,3>::CalculateAndAddStrainGradientMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,4>::CalculateAndAddStrainGradientMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    noalias(rVariables.PUMatrix) = -rVariables.VelocityCoefficient*0.25*rFICVariables.ElementLength*rFICVariables.ElementLength*rVariables.BiotCoefficient*
                                        prod(rVariables.GradNpT,rFICVariables.StrainGradients)*rVariables.IntegrationCoefficient;

    //Distribute strain gradient matrix into the elemental matrix
    PoroElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix,rVariables.PUMatrix);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,4>::CalculateAndAddStrainGradientMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,8>::CalculateAndAddStrainGradientMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    noalias(rVariables.PUMatrix) = -rVariables.VelocityCoefficient*0.25*rFICVariables.ElementLength*rFICVariables.ElementLength*rVariables.BiotCoefficient*
                                        prod(rVariables.GradNpT,rFICVariables.StrainGradients)*rVariables.IntegrationCoefficient;

    //Distribute strain gradient matrix into the elemental matrix
    PoroElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix,rVariables.PUMatrix);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::CalculateAndAddDtStressGradientMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    this->CalculateConstitutiveTensorGradients(rFICVariables,rVariables);
    
    double StabilizationParameter = rFICVariables.ElementLength*rFICVariables.ElementLength*rVariables.BiotCoefficient/(8.0*rFICVariables.ShearModulus);
    
    noalias(rVariables.PUMatrix) = -rVariables.VelocityCoefficient*StabilizationParameter/3.0*prod(rVariables.GradNpT,rFICVariables.DimUMatrix)*rVariables.IntegrationCoefficient;

    //Distribute DtStressGradient Matrix into the elemental matrix
    PoroElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix,rVariables.PUMatrix);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,3>::CalculateConstitutiveTensorGradients(FICElementVariables& rFICVariables, const ElementVariables& Variables)
{    
    for(unsigned int i = 0; i < 2; i++) //TDim
    {
        for(unsigned int j = 0; j < 3; j++) //VoigtSize
        {            
            for(unsigned int k = 0; k < 2; k++) //TDim
            {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;
                
                for(unsigned int l = 0; l < 3; l++) //TNumNodes
                    (rFICVariables.ConstitutiveTensorGradients[i][j][k]) += Variables.GradNpT(l,k)*(mNodalConstitutiveTensor[i][j][l]);
            }
        }
    }

    for(unsigned int i = 0; i < 2; i++) //TDim
    {
        for(unsigned int j = 0; j < 3; j++) //VoigtSize
        {
            rFICVariables.DimVoigtMatrix(i,j) = 0.0;
            
            for(unsigned int k = 0; k < 2; k++) //TDim
                rFICVariables.DimVoigtMatrix(i,j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }

    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix,Variables.B);

    /* INFO:
     * 
     * rFICVariables.ConstitutiveTensorGradients = [ ( |D00x|   |D01x|   |D02x| )   ( |D10x|   |D11x|   |D12x| ) ]
     *                                             [ ( |D00y| , |D01y| , |D02y| ) , ( |D10y| , |D11y| , |D12y| ) ]
     * 
     * rFICVariables.DimVoigtMatrix = |D00x+D10x D01x+D11x D02x+D12x|
     *                                |D00y+D10y D01y+D11y D02y+D12y|
     * 
     * D00x = aD(0,0)/ax at current GP
    */
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,4>::CalculateConstitutiveTensorGradients(FICElementVariables& rFICVariables, const ElementVariables& Variables)
{    
    for(unsigned int i = 0; i < 2; i++) //TDim
    {
        for(unsigned int j = 0; j < 3; j++) //VoigtSize
        {            
            for(unsigned int k = 0; k < 2; k++) //TDim
            {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;
                
                for(unsigned int l = 0; l < 4; l++) //TNumNodes
                    rFICVariables.ConstitutiveTensorGradients[i][j][k] += Variables.GradNpT(l,k)*mNodalConstitutiveTensor[i][j][l];
            }
        }
    }

    for(unsigned int i = 0; i < 2; i++) //TDim
    {
        for(unsigned int j = 0; j < 3; j++) //VoigtSize
        {
            rFICVariables.DimVoigtMatrix(i,j) = 0.0;
            
            for(unsigned int k = 0; k < 2; k++) //TDim
                rFICVariables.DimVoigtMatrix(i,j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }

    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix,Variables.B);
    
    // Adding ShapeFunctionsSecondOrderGradients terms
    unsigned int index;
    
    for(unsigned int i = 0; i < 4; i++) //TNumNodes
    {
        index = 2*i;
        
        rFICVariables.DimUMatrix(0,index)   += rFICVariables.ShapeFunctionsSecondOrderGradients[i][0]*( Variables.ConstitutiveMatrix(0,0) + Variables.ConstitutiveMatrix(1,0) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][2]*( Variables.ConstitutiveMatrix(0,2) + Variables.ConstitutiveMatrix(1,2) );
        rFICVariables.DimUMatrix(0,index+1) += rFICVariables.ShapeFunctionsSecondOrderGradients[i][0]*( Variables.ConstitutiveMatrix(0,2) + Variables.ConstitutiveMatrix(1,2) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][2]*( Variables.ConstitutiveMatrix(0,1) + Variables.ConstitutiveMatrix(1,1) );
        rFICVariables.DimUMatrix(1,index)   += rFICVariables.ShapeFunctionsSecondOrderGradients[i][1]*( Variables.ConstitutiveMatrix(0,2) + Variables.ConstitutiveMatrix(1,2) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][2]*( Variables.ConstitutiveMatrix(0,0) + Variables.ConstitutiveMatrix(1,0) );
        rFICVariables.DimUMatrix(1,index+1) += rFICVariables.ShapeFunctionsSecondOrderGradients[i][1]*( Variables.ConstitutiveMatrix(0,1) + Variables.ConstitutiveMatrix(1,1) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][2]*( Variables.ConstitutiveMatrix(0,2) + Variables.ConstitutiveMatrix(1,2) );
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,4>::CalculateConstitutiveTensorGradients(FICElementVariables& rFICVariables, const ElementVariables& Variables)
{    
    for(unsigned int i = 0; i < 3; i++) //TDim
    {
        for(unsigned int j = 0; j < 6; j++) //VoigtSize
        {            
            for(unsigned int k = 0; k < 3; k++) //TDim
            {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;
                
                for(unsigned int l = 0; l < 4; l++) //TNumNodes
                    rFICVariables.ConstitutiveTensorGradients[i][j][k] += Variables.GradNpT(l,k)*mNodalConstitutiveTensor[i][j][l];
            }
        }
    }

    for(unsigned int i = 0; i < 3; i++) //TDim
    {
        for(unsigned int j = 0; j < 6; j++) //VoigtSize
        {
            rFICVariables.DimVoigtMatrix(i,j) = 0.0;
            
            for(unsigned int k = 0; k < 3; k++) //TDim
                rFICVariables.DimVoigtMatrix(i,j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }
    
    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix,Variables.B);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,8>::CalculateConstitutiveTensorGradients(FICElementVariables& rFICVariables, const ElementVariables& Variables)
{    
    for(unsigned int i = 0; i < 3; i++) //TDim
    {
        for(unsigned int j = 0; j < 6; j++) //VoigtSize
        {            
            for(unsigned int k = 0; k < 3; k++) //TDim
            {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;
                
                for(unsigned int l = 0; l < 8; l++) //TNumNodes
                    rFICVariables.ConstitutiveTensorGradients[i][j][k] += Variables.GradNpT(l,k)*mNodalConstitutiveTensor[i][j][l];
            }
        }
    }

    for(unsigned int i = 0; i < 3; i++) //TDim
    {
        for(unsigned int j = 0; j < 6; j++) //VoigtSize
        {
            rFICVariables.DimVoigtMatrix(i,j) = 0.0;
            
            for(unsigned int k = 0; k < 3; k++) //TDim
                rFICVariables.DimVoigtMatrix(i,j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }

    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix,Variables.B);
    
    // Adding ShapeFunctionsSecondOrderGradients terms
    unsigned int index;
    
    for(unsigned int i = 0; i < 8; i++) //TNumNodes
    {
        index = 3*i;
        
        rFICVariables.DimUMatrix(0,index)   += rFICVariables.ShapeFunctionsSecondOrderGradients[i][0]*( Variables.ConstitutiveMatrix(0,0) + Variables.ConstitutiveMatrix(1,0) + Variables.ConstitutiveMatrix(2,0) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][3]*( Variables.ConstitutiveMatrix(0,3) + Variables.ConstitutiveMatrix(1,3) + Variables.ConstitutiveMatrix(2,3) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][5]*( Variables.ConstitutiveMatrix(0,5) + Variables.ConstitutiveMatrix(1,5) + Variables.ConstitutiveMatrix(2,5) );
        rFICVariables.DimUMatrix(0,index+1) += rFICVariables.ShapeFunctionsSecondOrderGradients[i][0]*( Variables.ConstitutiveMatrix(0,3) + Variables.ConstitutiveMatrix(1,3) + Variables.ConstitutiveMatrix(2,3) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][3]*( Variables.ConstitutiveMatrix(0,1) + Variables.ConstitutiveMatrix(1,1) + Variables.ConstitutiveMatrix(2,1) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][5]*( Variables.ConstitutiveMatrix(0,4) + Variables.ConstitutiveMatrix(1,4) + Variables.ConstitutiveMatrix(2,4) );
        rFICVariables.DimUMatrix(0,index+2) += rFICVariables.ShapeFunctionsSecondOrderGradients[i][0]*( Variables.ConstitutiveMatrix(0,5) + Variables.ConstitutiveMatrix(1,5) + Variables.ConstitutiveMatrix(2,5) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][3]*( Variables.ConstitutiveMatrix(0,4) + Variables.ConstitutiveMatrix(1,4) + Variables.ConstitutiveMatrix(2,4) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][5]*( Variables.ConstitutiveMatrix(0,2) + Variables.ConstitutiveMatrix(1,2) + Variables.ConstitutiveMatrix(2,2) );
        rFICVariables.DimUMatrix(1,index)   += rFICVariables.ShapeFunctionsSecondOrderGradients[i][1]*( Variables.ConstitutiveMatrix(0,3) + Variables.ConstitutiveMatrix(1,3) + Variables.ConstitutiveMatrix(2,3) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][3]*( Variables.ConstitutiveMatrix(0,0) + Variables.ConstitutiveMatrix(1,0) + Variables.ConstitutiveMatrix(2,0) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][4]*( Variables.ConstitutiveMatrix(0,5) + Variables.ConstitutiveMatrix(1,5) + Variables.ConstitutiveMatrix(2,5) );
        rFICVariables.DimUMatrix(1,index+1) += rFICVariables.ShapeFunctionsSecondOrderGradients[i][1]*( Variables.ConstitutiveMatrix(0,1) + Variables.ConstitutiveMatrix(1,1) + Variables.ConstitutiveMatrix(2,1) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][3]*( Variables.ConstitutiveMatrix(0,3) + Variables.ConstitutiveMatrix(1,3) + Variables.ConstitutiveMatrix(2,3) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][4]*( Variables.ConstitutiveMatrix(0,4) + Variables.ConstitutiveMatrix(1,4) + Variables.ConstitutiveMatrix(2,4) );
        rFICVariables.DimUMatrix(1,index+2) += rFICVariables.ShapeFunctionsSecondOrderGradients[i][1]*( Variables.ConstitutiveMatrix(0,4) + Variables.ConstitutiveMatrix(1,4) + Variables.ConstitutiveMatrix(2,4) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][3]*( Variables.ConstitutiveMatrix(0,5) + Variables.ConstitutiveMatrix(1,5) + Variables.ConstitutiveMatrix(2,5) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][4]*( Variables.ConstitutiveMatrix(0,2) + Variables.ConstitutiveMatrix(1,2) + Variables.ConstitutiveMatrix(2,2) );
        rFICVariables.DimUMatrix(2,index)   += rFICVariables.ShapeFunctionsSecondOrderGradients[i][2]*( Variables.ConstitutiveMatrix(0,5) + Variables.ConstitutiveMatrix(1,5) + Variables.ConstitutiveMatrix(2,5) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][4]*( Variables.ConstitutiveMatrix(0,3) + Variables.ConstitutiveMatrix(1,3) + Variables.ConstitutiveMatrix(2,3) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][5]*( Variables.ConstitutiveMatrix(0,0) + Variables.ConstitutiveMatrix(1,0) + Variables.ConstitutiveMatrix(2,0) );
        rFICVariables.DimUMatrix(2,index+1) += rFICVariables.ShapeFunctionsSecondOrderGradients[i][2]*( Variables.ConstitutiveMatrix(0,4) + Variables.ConstitutiveMatrix(1,4) + Variables.ConstitutiveMatrix(2,4) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][4]*( Variables.ConstitutiveMatrix(0,1) + Variables.ConstitutiveMatrix(1,1) + Variables.ConstitutiveMatrix(2,1) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][5]*( Variables.ConstitutiveMatrix(0,3) + Variables.ConstitutiveMatrix(1,3) + Variables.ConstitutiveMatrix(2,3) );
        rFICVariables.DimUMatrix(2,index+2) += rFICVariables.ShapeFunctionsSecondOrderGradients[i][2]*( Variables.ConstitutiveMatrix(0,2) + Variables.ConstitutiveMatrix(1,2) + Variables.ConstitutiveMatrix(2,2) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][4]*( Variables.ConstitutiveMatrix(0,4) + Variables.ConstitutiveMatrix(1,4) + Variables.ConstitutiveMatrix(2,4) ) +
                                               rFICVariables.ShapeFunctionsSecondOrderGradients[i][5]*( Variables.ConstitutiveMatrix(0,5) + Variables.ConstitutiveMatrix(1,5) + Variables.ConstitutiveMatrix(2,5) );
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::CalculateAndAddPressureGradientMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    double StabilizationParameter = rFICVariables.ElementLength*rFICVariables.ElementLength*rVariables.BiotCoefficient/(8.0*rFICVariables.ShearModulus);
    
    noalias(rVariables.PMatrix) = rVariables.DtPressureCoefficient*StabilizationParameter*(rVariables.BiotCoefficient-2.0*rFICVariables.ShearModulus*rVariables.BiotModulusInverse/(3.0*rVariables.BiotCoefficient))*
                                      prod(rVariables.GradNpT,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;
    
    //Distribute pressure gradient block matrix into the elemental matrix
    PoroElementUtilities::AssemblePBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::CalculateAndAddRHSStabilization(VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    this->CalculateAndAddStrainGradientFlow(rRightHandSideVector, rVariables, rFICVariables);

    this->CalculateAndAddDtStressGradientFlow(rRightHandSideVector, rVariables, rFICVariables);

    this->CalculateAndAddPressureGradientFlow(rRightHandSideVector, rVariables, rFICVariables);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,3>::CalculateAndAddStrainGradientFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,4>::CalculateAndAddStrainGradientFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    noalias(rVariables.PUMatrix) = 0.25*rFICVariables.ElementLength*rFICVariables.ElementLength*rVariables.BiotCoefficient*
                                    prod(rVariables.GradNpT,rFICVariables.StrainGradients)*rVariables.IntegrationCoefficient;
    
    noalias(rVariables.PVector) = prod(rVariables.PUMatrix,rVariables.VelocityVector);
    
    //Distribute Strain Gradient vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,4> >(rRightHandSideVector,rVariables.PVector,2,4);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,4>::CalculateAndAddStrainGradientFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,8>::CalculateAndAddStrainGradientFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    noalias(rVariables.PUMatrix) = 0.25*rFICVariables.ElementLength*rFICVariables.ElementLength*rVariables.BiotCoefficient*
                                    prod(rVariables.GradNpT,rFICVariables.StrainGradients)*rVariables.IntegrationCoefficient;
    
    noalias(rVariables.PVector) = prod(rVariables.PUMatrix,rVariables.VelocityVector);
    
    //Distribute Strain Gradient vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,8> >(rRightHandSideVector,rVariables.PVector,3,8);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::CalculateAndAddDtStressGradientFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{
    this->CalculateDtStressGradients(rFICVariables,rVariables);
    
    double StabilizationParameter = rFICVariables.ElementLength*rFICVariables.ElementLength*rVariables.BiotCoefficient/(8.0*rFICVariables.ShearModulus);
    
    noalias(rVariables.PVector) = StabilizationParameter/3.0*prod(rVariables.GradNpT,rFICVariables.DimVector)*rVariables.IntegrationCoefficient;
    
    //Distribute DtStressGradient block vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::CalculateDtStressGradients(FICElementVariables& rFICVariables, const ElementVariables& Variables)
{
    for(unsigned int i = 0; i < TDim; i++)
    {
        for(unsigned int j = 0; j < TDim; j++)
        {
            rFICVariables.DtStressGradients[i][j] = 0.0;
            
            for(unsigned int k = 0; k < TNumNodes; k++)
                rFICVariables.DtStressGradients[i][j] += Variables.GradNpT(k,j)*mNodalDtStress[i][k];
        }
    }
    
    for(unsigned int i = 0; i < TDim; i++)
    {
        rFICVariables.DimVector[i] = 0.0;
        
        for(unsigned int j = 0; j < TDim; j++)
            rFICVariables.DimVector[i] += rFICVariables.DtStressGradients[j][i];
    }
    
    /* INFO: (2D)
     * 
     * rFICVariables.DtStressGradients = ( |S0x|   |S1x| )
     *                                   ( |S0y| , |S1y| )
     * 
     * rFICVariables.DimVector = |S0x+S1x|
     *                           |S0y+S1y|
     * 
     * S0x = aS[0]/ax at current GP
    */
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::CalculateAndAddPressureGradientFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables, FICElementVariables& rFICVariables)
{    
    double StabilizationParameter = rFICVariables.ElementLength*rFICVariables.ElementLength*rVariables.BiotCoefficient/(8.0*rFICVariables.ShearModulus);
    
    noalias(rVariables.PMatrix) = StabilizationParameter*(rVariables.BiotCoefficient-2.0*rFICVariables.ShearModulus*rVariables.BiotModulusInverse/(3.0*rVariables.BiotCoefficient))*
                                      prod(rVariables.GradNpT,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;
    
    noalias(rVariables.PVector) = -1.0*prod(rVariables.PMatrix,rVariables.DtPressureVector);

    //Distribute PressureGradient block vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwSmallStrainFICElement<2,3>;
template class UPwSmallStrainFICElement<2,4>;
template class UPwSmallStrainFICElement<3,4>;
template class UPwSmallStrainFICElement<3,8>;

} // Namespace Kratos
