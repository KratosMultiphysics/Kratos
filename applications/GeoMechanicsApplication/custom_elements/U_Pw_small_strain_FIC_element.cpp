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
void UPwSmallStrainFICElement<TDim,TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    UPwBaseElement<TDim,TNumNodes>::Initialize(rCurrentProcessInfo);

    unsigned int VoigtSize = VOIGT_SIZE_3D;
    if (TDim == 2) VoigtSize = VOIGT_SIZE_2D_PLANE_STRESS;

    for (unsigned int i = 0; i < TDim; i++)
    {
        mNodalConstitutiveTensor[i].resize(VoigtSize);

        for (unsigned int j = 0; j < VoigtSize; j++)
        {
            for (unsigned int k = 0; k < TNumNodes; k++)
                mNodalConstitutiveTensor[i][j][k] = 0.0;
        }
    }

    for (unsigned int i = 0; i < TDim; i++)
    {
        for (unsigned int j = 0; j < TNumNodes; j++)
            mNodalDtStress[i][j] = 0.0;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    UPwSmallStrainElement<TDim,TNumNodes>::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    UPwSmallStrainElement<TDim,TNumNodes>::FinalizeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,this->GetProperties(),rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); //Note: this is for nonlocal damage

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,
                                     rCurrentProcessInfo);

    //Extrapolation variables
    array_1d<Matrix,TDim> ConstitutiveTensorContainer;
    for (unsigned int i = 0; i < TDim; i++)
    {
        ConstitutiveTensorContainer[i].resize(NumGPoints, Variables.ConstitutiveMatrix.size1(), false);
    }
    Matrix DtStressContainer(NumGPoints,TDim);

    Vector StressVector;
    if (TDim > 2)
        StressVector.resize(VOIGT_SIZE_3D);
    else
        StressVector.resize(VOIGT_SIZE_2D_PLANE_STRESS);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        this->CalculateKinematics(Variables, GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        UPwSmallStrainElement<TDim,TNumNodes>::UpdateElementalVariableStressVector(Variables, GPoint);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
        this->SaveGPConstitutiveTensor( ConstitutiveTensorContainer,
                                        Variables.ConstitutiveMatrix,
                                        GPoint );

        // Compute DtStress
        noalias(Variables.StrainVector) = prod(Variables.B, Variables.VelocityVector);

        noalias(StressVector) = prod(Variables.ConstitutiveMatrix, Variables.StrainVector);
        this->SaveGPDtStress(DtStressContainer, StressVector, GPoint);
    }

    this->ExtrapolateGPConstitutiveTensor(ConstitutiveTensorContainer);
    this->ExtrapolateGPDtStress(DtStressContainer);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );

    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,this->GetProperties(),rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); //Note: this is for nonlocal damage

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,
                                     rCurrentProcessInfo);

    //Containers for extrapolation variables
    Matrix DtStressContainer(NumGPoints,TDim);

    Vector StressVector;
    if (TDim > 2)
        StressVector.resize(VOIGT_SIZE_3D);
    else
        StressVector.resize(VOIGT_SIZE_2D_PLANE_STRESS);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        this->CalculateKinematics(Variables, GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Compute ConstitutiveTensor
        UPwSmallStrainElement<TDim,TNumNodes>::UpdateElementalVariableStressVector(Variables, GPoint);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        // Compute DtStress
        noalias(Variables.StrainVector) = prod(Variables.B, Variables.VelocityVector);
        noalias(StressVector) = prod(Variables.ConstitutiveMatrix, Variables.StrainVector);
        this->SaveGPDtStress(DtStressContainer, StressVector, GPoint);
    }

    this->ExtrapolateGPDtStress(DtStressContainer);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    SaveGPConstitutiveTensor( array_1d<Matrix,TDim>& rConstitutiveTensorContainer,
                              const Matrix& ConstitutiveMatrix,
                              const unsigned int& GPoint )
{
    for (unsigned int i = 0; i < TDim; i++)
    {
        for (unsigned int j = 0; j < ConstitutiveMatrix.size1(); j++) //VoigtSize
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
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    SaveGPDtStress( Matrix& rDtStressContainer,
                    const Vector& StressVector,
                    const unsigned int& GPoint )
{
    for (unsigned int i = 0; i < TDim; i++)
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
void UPwSmallStrainFICElement<2,3>::
    ExtrapolateGPConstitutiveTensor(const array_1d<Matrix,2>& ConstitutiveTensorContainer)
{
    // Triangle_2d_3 with GI_GAUSS_2
    
    BoundedMatrix<double,3,3> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,3,3> AuxNodalConstitutiveTensor;
    
    for (unsigned int i = 0; i < 2; i++) //TDim
    {
        noalias(AuxNodalConstitutiveTensor) = prod(ExtrapolationMatrix,ConstitutiveTensorContainer[i]);
        
        for (unsigned int j = 0; j < VOIGT_SIZE_2D_PLANE_STRESS; j++) // VoigtSize
            noalias(mNodalConstitutiveTensor[i][j]) = column(AuxNodalConstitutiveTensor,j);
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
void UPwSmallStrainFICElement<2,4>::
    ExtrapolateGPConstitutiveTensor(const array_1d<Matrix,2>& ConstitutiveTensorContainer)
{
    // Quadrilateral_2d_4 with GI_GAUSS_2

    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,4,3> AuxNodalConstitutiveTensor;

    for (unsigned int i = 0; i < 2; i++) //TDim
    {
        noalias(AuxNodalConstitutiveTensor) = prod(ExtrapolationMatrix,ConstitutiveTensorContainer[i]);

        for (unsigned int j = 0; j < VOIGT_SIZE_2D_PLANE_STRESS; j++) // VoigtSize
            noalias(mNodalConstitutiveTensor[i][j]) = column(AuxNodalConstitutiveTensor,j);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,4>::
    ExtrapolateGPConstitutiveTensor(const array_1d<Matrix,3>& ConstitutiveTensorContainer)
{
    // Tetrahedra_3d_4 with GI_GAUSS_2
    
    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,4,6> AuxNodalConstitutiveTensor;
    
    for (unsigned int i = 0; i < 3; i++) //TDim
    {
        noalias(AuxNodalConstitutiveTensor) = prod(ExtrapolationMatrix,ConstitutiveTensorContainer[i]);
        
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; j++) // VoigtSize
            noalias(mNodalConstitutiveTensor[i][j]) = column(AuxNodalConstitutiveTensor,j);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,8>::
    ExtrapolateGPConstitutiveTensor(const array_1d<Matrix,3>& ConstitutiveTensorContainer)
{
    // Hexahedra_3d_8 with GI_GAUSS_2

    BoundedMatrix<double,8,8> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double,8,6> AuxNodalConstitutiveTensor;

    for (unsigned int i = 0; i < 3; i++) //TDim
    {
        noalias(AuxNodalConstitutiveTensor) = prod(ExtrapolationMatrix,ConstitutiveTensorContainer[i]);

        for (unsigned int j = 0; j < VOIGT_SIZE_3D; j++) // VoigtSize
            noalias(mNodalConstitutiveTensor[i][j]) = column(AuxNodalConstitutiveTensor,j);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,3>::ExtrapolateGPDtStress(const Matrix& DtStressContainer)
{
    // Triangle_2d_3 with GI_GAUSS_2

    BoundedMatrix<double,3,3> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,3,2> AuxNodalDtStress;
    noalias(AuxNodalDtStress) = prod(ExtrapolationMatrix,DtStressContainer);
    
    for (unsigned int i = 0; i < 2; i++) // TDim
        noalias(mNodalDtStress[i]) = column(AuxNodalDtStress,i);

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
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,4,2> AuxNodalDtStress;
    noalias(AuxNodalDtStress) = prod(ExtrapolationMatrix,DtStressContainer);

    for (unsigned int i = 0; i < 2; i++) // TDim
        noalias(mNodalDtStress[i]) = column(AuxNodalDtStress,i);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,4>::ExtrapolateGPDtStress(const Matrix& DtStressContainer)
{
    // Tetrahedra_3d_4 with GI_GAUSS_2
    
    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);
    
    BoundedMatrix<double,4,3> AuxNodalDtStress;
    noalias(AuxNodalDtStress) = prod(ExtrapolationMatrix,DtStressContainer);
    
    for (unsigned int i = 0; i < 3; i++) // TDim
        noalias(mNodalDtStress[i]) = column(AuxNodalDtStress,i);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,8>::ExtrapolateGPDtStress(const Matrix& DtStressContainer)
{
    // Hexahedra_3d_8 with GI_GAUSS_2

    BoundedMatrix<double,8,8> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double,8,3> AuxNodalDtStress;
    noalias(AuxNodalDtStress) = prod(ExtrapolationMatrix,DtStressContainer);

    for (unsigned int i = 0; i < 3; i++) // TDim
        noalias(mNodalDtStress[i]) = column(AuxNodalDtStress,i);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo,
                  const bool CalculateStiffnessMatrixFlag,
                  const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( this->GetIntegrationMethod() );
    const unsigned int NumGPoints = IntegrationPoints.size();

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom, Prop, CurrentProcessInfo);
    if (CalculateStiffnessMatrixFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag)  ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, CurrentProcessInfo);

    FICElementVariables FICVariables;
    this->InitializeFICElementVariables(FICVariables, Variables.DN_DXContainer, Geom, Prop, CurrentProcessInfo);

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(Geom, this->GetProperties(), CurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute Np, GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        GeoElementUtilities::
            CalculateNuMatrix<TDim, TNumNodes>( Variables.Nu,
                                                Variables.NContainer,
                                                GPoint );
        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim, TNumNodes>( Variables.BodyAcceleration,
                                                                Variables.NContainer,
                                                                Variables.VolumeAcceleration,
                                                                GPoint );

        //Compute ShapeFunctionsSecondOrderGradients
        this->CalculateShapeFunctionsSecondOrderGradients(FICVariables,Variables);

        //Compute constitutive tensor and stresses
        UPwSmallStrainElement<TDim,TNumNodes>::UpdateElementalVariableStressVector(Variables, GPoint);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
        UPwSmallStrainElement<TDim,TNumNodes>::UpdateStressVector(Variables, GPoint);

        this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        // set shear modulus from stiffness matrix
        FICVariables.ShearModulus = CalculateShearModulus(Variables.ConstitutiveMatrix);

        // calculate Bulk modulus from stiffness matrix
        const double BulkModulus = CalculateBulkModulus(Variables.ConstitutiveMatrix);
        this->InitializeBiotCoefficients(Variables, BulkModulus);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient,
                                              Variables.detJ0,
                                              IntegrationPoints[GPoint].Weight());

        if (CalculateStiffnessMatrixFlag)
        {
            //Contributions to the left hand side
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
            this->CalculateAndAddLHSStabilization(rLeftHandSideMatrix, Variables, FICVariables);
        }

        if (CalculateResidualVectorFlag)
        {
            //Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
            this->CalculateAndAddRHSStabilization(rRightHandSideVector, Variables, FICVariables);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwSmallStrainFICElement<TDim,TNumNodes>::
    CalculateShearModulus(const Matrix &ConstitutiveMatrix) const
{
    const int IndexG = ConstitutiveMatrix.size1() - 1;
    return ConstitutiveMatrix(IndexG, IndexG);
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    InitializeFICElementVariables( FICElementVariables& rFICVariables,
                                   const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer,
                                   const GeometryType& Geom,
                                   const PropertiesType& Prop,
                                   const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

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
void UPwSmallStrainFICElement<2,3>::
    ExtrapolateShapeFunctionsGradients( array_1d< array_1d<double,6> , 3 >& rNodalShapeFunctionsGradients,
                                        const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer )
{
    // Triangle_2d_3 with GI_GAUSS_2
    // No necessary
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,4>::
    ExtrapolateShapeFunctionsGradients( array_1d< array_1d<double,8> , 4 >& rNodalShapeFunctionsGradients,
                                        const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer )
{
    // Quadrilateral_2d_4 with GI_GAUSS_2

    BoundedMatrix<double,4,8> ShapeFunctionsGradientsContainer; //NumGPoints X TDim*TNumNodes
    unsigned int index;

    for (unsigned int i = 0; i < 4; i++) //NumGPoints
    {
        for (unsigned int j = 0; j < 4; j++) //TNumNodes
        {
            index = j*2;

            ShapeFunctionsGradientsContainer(i,index) = DN_DXContainer[i](j,0);
            ShapeFunctionsGradientsContainer(i,index+1) = DN_DXContainer[i](j,1);
        }
    }

    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double,4,8> AuxNodalShapeFunctionsGradients;
    noalias(AuxNodalShapeFunctionsGradients) = prod(ExtrapolationMatrix,ShapeFunctionsGradientsContainer);

    for (unsigned int i = 0; i < 4; i++) //TNumNodes
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
void UPwSmallStrainFICElement<3,4>::
    ExtrapolateShapeFunctionsGradients( array_1d< array_1d<double,12> , 4 >& rNodalShapeFunctionsGradients,
                                        const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer )
{
    // Tetrahedra_3d_4 with GI_GAUSS_2
    // No necessary
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,8>::
    ExtrapolateShapeFunctionsGradients( array_1d< array_1d<double,24> , 8 >& rNodalShapeFunctionsGradients,
                                        const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer )
{
    // Hexahedra_3d_8 with GI_GAUSS_2

    BoundedMatrix<double,8,24> ShapeFunctionsGradientsContainer; //NumGPoints X TDim*TNumNodes
    unsigned int index;

    for (unsigned int i = 0; i < 8; i++) //NumGPoints
    {
        for (unsigned int j = 0; j < 8; j++) //TNumNodes
        {
            index = j*3;

            ShapeFunctionsGradientsContainer(i,index) = DN_DXContainer[i](j,0);
            ShapeFunctionsGradientsContainer(i,index+1) = DN_DXContainer[i](j,1);
            ShapeFunctionsGradientsContainer(i,index+2) = DN_DXContainer[i](j,2);
        }
    }

    BoundedMatrix<double,8,8> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double,8,24> AuxNodalShapeFunctionsGradients;
    noalias(AuxNodalShapeFunctionsGradients) = prod(ExtrapolationMatrix,ShapeFunctionsGradientsContainer);

    for (unsigned int i = 0; i < 8; i++) //TNumNodes
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
    for (unsigned int i = 0; i < 2; i++) //TDim
        rFICVariables.ConstitutiveTensorGradients[i].resize(VOIGT_SIZE_2D_PLANE_STRESS); //VoigtSize

    rFICVariables.DimVoigtMatrix.resize(2,3,false); //TDim X VoigtSize
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,4>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    //Voigt identity matrix
    rFICVariables.VoigtMatrix.resize(VOIGT_SIZE_2D_PLANE_STRESS, VOIGT_SIZE_2D_PLANE_STRESS,false); //VoigtSize X VoigtSize
    noalias(rFICVariables.VoigtMatrix) = ZeroMatrix(VOIGT_SIZE_2D_PLANE_STRESS, VOIGT_SIZE_2D_PLANE_STRESS);
    rFICVariables.VoigtMatrix(0,0) = 1.0;
    rFICVariables.VoigtMatrix(1,1) = 1.0;
    rFICVariables.VoigtMatrix(2,2) = 0.5;

    for (unsigned int i = 0; i < 4; i++) //TNumNodes
        rFICVariables.ShapeFunctionsSecondOrderGradients[i].resize(VOIGT_SIZE_2D_PLANE_STRESS,false); //VoigtSize

    for (unsigned int i = 0; i < 2; i++) //TDim
        rFICVariables.ConstitutiveTensorGradients[i].resize(VOIGT_SIZE_2D_PLANE_STRESS); //VoigtSize

    rFICVariables.DimVoigtMatrix.resize(2,VOIGT_SIZE_2D_PLANE_STRESS,false); //TDim X VoigtSize
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,4>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    for (unsigned int i = 0; i < 3; i++) //TDim
        rFICVariables.ConstitutiveTensorGradients[i].resize(VOIGT_SIZE_3D); //VoigtSize

    rFICVariables.DimVoigtMatrix.resize(3,VOIGT_SIZE_3D,false); //TDim X VoigtSize
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,8>::InitializeSecondOrderTerms(FICElementVariables& rFICVariables)
{
    //Voigt identity matrix
    rFICVariables.VoigtMatrix.resize(VOIGT_SIZE_3D, VOIGT_SIZE_3D,false); //VoigtSize X VoigtSize
    noalias(rFICVariables.VoigtMatrix) = ZeroMatrix(VOIGT_SIZE_3D,VOIGT_SIZE_3D);
    rFICVariables.VoigtMatrix(0,0) = 1.0;
    rFICVariables.VoigtMatrix(1,1) = 1.0;
    rFICVariables.VoigtMatrix(2,2) = 1.0;
    rFICVariables.VoigtMatrix(3,3) = 0.5;
    rFICVariables.VoigtMatrix(4,4) = 0.5;
    rFICVariables.VoigtMatrix(5,5) = 0.5;

    for (unsigned int i = 0; i < 8; i++) //TNumNodes
        rFICVariables.ShapeFunctionsSecondOrderGradients[i].resize(VOIGT_SIZE_3D,false); //VoigtSize

    for (unsigned int i = 0; i < 3; i++) //TDim
        rFICVariables.ConstitutiveTensorGradients[i].resize(VOIGT_SIZE_3D); //VoigtSize

    rFICVariables.DimVoigtMatrix.resize(3,VOIGT_SIZE_3D,false); //TDim X VoigtSize
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,3>::
    CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables, ElementVariables& rVariables)
{
    // Not necessary
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<2,4>::
    CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables, ElementVariables& rVariables)
{
    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B),rFICVariables.VoigtMatrix);
    unsigned int index;
    for (unsigned int i = 0; i < 4; i++) //TNumNodes
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
void UPwSmallStrainFICElement<3,4>::
    CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables, ElementVariables& rVariables)
{
    // No necessary
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainFICElement<3,8>::
    CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables, ElementVariables& rVariables)
{
    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B),rFICVariables.VoigtMatrix);
    unsigned int index;
    for (unsigned int i = 0; i < 8; i++) //TNumNodes
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
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    CalculateAndAddLHSStabilization( MatrixType& rLeftHandSideMatrix,
                                     ElementVariables& rVariables,
                                     FICElementVariables& rFICVariables )
{
    KRATOS_TRY;

    this->CalculateAndAddStrainGradientMatrix(rLeftHandSideMatrix,rVariables,rFICVariables);

    this->CalculateAndAddDtStressGradientMatrix(rLeftHandSideMatrix,rVariables,rFICVariables);

    this->CalculateAndAddPressureGradientMatrix(rLeftHandSideMatrix,rVariables,rFICVariables);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,3>::
    CalculateAndAddStrainGradientMatrix( MatrixType& rLeftHandSideMatrix,
                                         ElementVariables& rVariables,
                                         FICElementVariables& rFICVariables )
{
    // No necessary
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,4>::
    CalculateAndAddStrainGradientMatrix( MatrixType& rLeftHandSideMatrix,
                                         ElementVariables& rVariables,
                                         FICElementVariables& rFICVariables )
{
    noalias(rVariables.PUMatrix) = PORE_PRESSURE_SIGN_FACTOR * rVariables.VelocityCoefficient
                                  * 0.25      * rFICVariables.ElementLength 
                                              * rFICVariables.ElementLength
                                              * rVariables.BiotCoefficient
                                              * prod(rVariables.GradNpT, rFICVariables.StrainGradients)
                                              * rVariables.IntegrationCoefficient;

    //Distribute strain gradient matrix into the elemental matrix
    GeoElementUtilities::AssemblePUBlockMatrix<2, 4>(rLeftHandSideMatrix,rVariables.PUMatrix);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,4>::
    CalculateAndAddStrainGradientMatrix( MatrixType& rLeftHandSideMatrix,
                                         ElementVariables& rVariables,
                                         FICElementVariables& rFICVariables )
{
    // No necessary
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,8>::
    CalculateAndAddStrainGradientMatrix( MatrixType& rLeftHandSideMatrix,
                                         ElementVariables& rVariables,
                                         FICElementVariables& rFICVariables )
{
    noalias(rVariables.PUMatrix) =  PORE_PRESSURE_SIGN_FACTOR * rVariables.VelocityCoefficient
                                  * 0.25       * rFICVariables.ElementLength
                                               * rFICVariables.ElementLength
                                               * rVariables.BiotCoefficient
                                               * prod(rVariables.GradNpT, rFICVariables.StrainGradients)
                                               * rVariables.IntegrationCoefficient;

    //Distribute strain gradient matrix into the elemental matrix
    GeoElementUtilities::AssemblePUBlockMatrix<3, 8>(rLeftHandSideMatrix,rVariables.PUMatrix);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    CalculateAndAddDtStressGradientMatrix( MatrixType& rLeftHandSideMatrix,
                                           ElementVariables& rVariables,
                                           FICElementVariables& rFICVariables )
{
    this->CalculateConstitutiveTensorGradients(rFICVariables,rVariables);

    double StabilizationParameter = - PORE_PRESSURE_SIGN_FACTOR * rFICVariables.ElementLength
                                                 * rFICVariables.ElementLength
                                                 * rVariables.BiotCoefficient/(8.0*rFICVariables.ShearModulus);

    noalias(rVariables.PUMatrix) = -  rVariables.VelocityCoefficient
                                    * StabilizationParameter / 3.0 * prod(rVariables.GradNpT, rFICVariables.DimUMatrix)
                                                                   * rVariables.IntegrationCoefficient;

    //Distribute DtStressGradient Matrix into the elemental matrix
    GeoElementUtilities::AssemblePUBlockMatrix<TDim, TNumNodes>(rLeftHandSideMatrix,rVariables.PUMatrix);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,3>::
    CalculateConstitutiveTensorGradients( FICElementVariables& rFICVariables,
                                          const ElementVariables& Variables )
{
    for (unsigned int i = 0; i < 2; i++) //TDim
    {
        for (unsigned int j = 0; j < VOIGT_SIZE_2D_PLANE_STRESS; j++) //VoigtSize
        {
            for (unsigned int k = 0; k < 2; k++) //TDim
            {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;

                for (unsigned int l = 0; l < 3; l++) //TNumNodes
                    (rFICVariables.ConstitutiveTensorGradients[i][j][k]) += Variables.GradNpT(l,k)*(mNodalConstitutiveTensor[i][j][l]);
            }
        }
    }

    for (unsigned int i = 0; i < 2; i++) //TDim
    {
        for (unsigned int j = 0; j < VOIGT_SIZE_2D_PLANE_STRESS; j++) //VoigtSize
        {
            rFICVariables.DimVoigtMatrix(i,j) = 0.0;

            for (unsigned int k = 0; k < 2; k++) //TDim
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
void UPwSmallStrainFICElement<2,4>::
    CalculateConstitutiveTensorGradients( FICElementVariables& rFICVariables,
                                          const ElementVariables& Variables )
{
    for (unsigned int i = 0; i < 2; i++) //TDim
    {
        for (unsigned int j = 0; j < VOIGT_SIZE_2D_PLANE_STRESS; j++) //VoigtSize
        {
            for (unsigned int k = 0; k < 2; k++) //TDim
            {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;

                for (unsigned int l = 0; l < 4; l++) //TNumNodes
                    rFICVariables.ConstitutiveTensorGradients[i][j][k] += Variables.GradNpT(l,k)*mNodalConstitutiveTensor[i][j][l];
            }
        }
    }

    for (unsigned int i = 0; i < 2; i++) //TDim
    {
        for (unsigned int j = 0; j < VOIGT_SIZE_2D_PLANE_STRESS; j++) //VoigtSize
        {
            rFICVariables.DimVoigtMatrix(i,j) = 0.0;

            for (unsigned int k = 0; k < 2; k++) //TDim
                rFICVariables.DimVoigtMatrix(i,j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }

    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix,Variables.B);

    // Adding ShapeFunctionsSecondOrderGradients terms
    unsigned int index;

    for (unsigned int i = 0; i < 4; i++) //TNumNodes
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
void UPwSmallStrainFICElement<3,4>::
    CalculateConstitutiveTensorGradients( FICElementVariables& rFICVariables,
                                          const ElementVariables& Variables )
{
    for (unsigned int i = 0; i < 3; i++) //TDim
    {
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; j++) //VoigtSize
        {
            for (unsigned int k = 0; k < 3; k++) //TDim
            {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;

                for (unsigned int l = 0; l < 4; l++) //TNumNodes
                    rFICVariables.ConstitutiveTensorGradients[i][j][k] += Variables.GradNpT(l,k)*mNodalConstitutiveTensor[i][j][l];
            }
        }
    }

    for (unsigned int i = 0; i < 3; i++) //TDim
    {
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; j++) //VoigtSize
        {
            rFICVariables.DimVoigtMatrix(i,j) = 0.0;

            for (unsigned int k = 0; k < 3; k++) //TDim
                rFICVariables.DimVoigtMatrix(i,j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }

    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix,Variables.B);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,8>::
    CalculateConstitutiveTensorGradients( FICElementVariables& rFICVariables,
                                          const ElementVariables& Variables )
{
    for (unsigned int i = 0; i < 3; i++) //TDim
    {
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; j++) //VoigtSize
        {
            for (unsigned int k = 0; k < 3; k++) //TDim
            {
                rFICVariables.ConstitutiveTensorGradients[i][j][k] = 0.0;

                for (unsigned int l = 0; l < 8; l++) //TNumNodes
                    rFICVariables.ConstitutiveTensorGradients[i][j][k] += Variables.GradNpT(l,k)*mNodalConstitutiveTensor[i][j][l];
            }
        }
    }

    for (unsigned int i = 0; i < 3; i++) //TDim
    {
        for (unsigned int j = 0; j < VOIGT_SIZE_3D; j++) //VoigtSize
        {
            rFICVariables.DimVoigtMatrix(i,j) = 0.0;

            for (unsigned int k = 0; k < 3; k++) //TDim
                rFICVariables.DimVoigtMatrix(i,j) += rFICVariables.ConstitutiveTensorGradients[k][j][i];
        }
    }

    noalias(rFICVariables.DimUMatrix) = prod(rFICVariables.DimVoigtMatrix,Variables.B);

    // Adding ShapeFunctionsSecondOrderGradients terms
    unsigned int index;

    for (unsigned int i = 0; i < 8; i++) //TNumNodes
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
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    CalculateAndAddPressureGradientMatrix( MatrixType& rLeftHandSideMatrix,
                                           ElementVariables& rVariables,
                                           FICElementVariables& rFICVariables)
{
    const double SignBiotCoefficient = - PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotCoefficient;

    const double StabilizationParameter =  rFICVariables.ElementLength
                                         * rFICVariables.ElementLength
                                         * SignBiotCoefficient/(8.0*rFICVariables.ShearModulus);

    noalias(rVariables.PMatrix) =   rVariables.DtPressureCoefficient * StabilizationParameter
                                  * ( SignBiotCoefficient - 2.0* rFICVariables.ShearModulus
                                                               * rVariables.BiotModulusInverse/(3.0*SignBiotCoefficient) )
                                  * prod(rVariables.GradNpT,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;

    //Distribute pressure gradient block matrix into the elemental matrix
    GeoElementUtilities::AssemblePBlockMatrix< TDim, TNumNodes >(rLeftHandSideMatrix,rVariables.PMatrix);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    CalculateAndAddRHSStabilization( VectorType& rRightHandSideVector,
                                     ElementVariables& rVariables,
                                     FICElementVariables& rFICVariables )
{
    KRATOS_TRY;

    this->CalculateAndAddStrainGradientFlow(rRightHandSideVector, rVariables, rFICVariables);

    this->CalculateAndAddDtStressGradientFlow(rRightHandSideVector, rVariables, rFICVariables);

    this->CalculateAndAddPressureGradientFlow(rRightHandSideVector, rVariables, rFICVariables);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,3>::
    CalculateAndAddStrainGradientFlow( VectorType& rRightHandSideVector,
                                       ElementVariables& rVariables,
                                       FICElementVariables& rFICVariables )
{
    // No necessary
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<2,4>::
    CalculateAndAddStrainGradientFlow( VectorType& rRightHandSideVector,
                                       ElementVariables& rVariables,
                                       FICElementVariables& rFICVariables )
{
    noalias(rVariables.PUMatrix) = 0.25 * rFICVariables.ElementLength * rFICVariables.ElementLength
                                        * rVariables.BiotCoefficient  * (-PORE_PRESSURE_SIGN_FACTOR)
                                        * prod(rVariables.GradNpT, rFICVariables.StrainGradients)
                                        * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = prod(rVariables.PUMatrix,rVariables.VelocityVector);

    //Distribute Strain Gradient vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector< 2, 4 >(rRightHandSideVector,rVariables.PVector);
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,4>::
    CalculateAndAddStrainGradientFlow( VectorType& rRightHandSideVector,
                                       ElementVariables& rVariables,
                                       FICElementVariables& rFICVariables)
{
    // Not necessary
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainFICElement<3,8>::
    CalculateAndAddStrainGradientFlow( VectorType& rRightHandSideVector,
                                       ElementVariables& rVariables,
                                       FICElementVariables& rFICVariables )
{
    noalias(rVariables.PUMatrix) = 0.25 * rFICVariables.ElementLength * rFICVariables.ElementLength
                                        * rVariables.BiotCoefficient  * (-PORE_PRESSURE_SIGN_FACTOR)
                                        * prod(rVariables.GradNpT, rFICVariables.StrainGradients)
                                        * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = prod(rVariables.PUMatrix,rVariables.VelocityVector);

    //Distribute Strain Gradient vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector< 3, 8 >(rRightHandSideVector,rVariables.PVector);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    CalculateAndAddDtStressGradientFlow( VectorType& rRightHandSideVector,
                                         ElementVariables& rVariables,
                                         FICElementVariables& rFICVariables )
{
    this->CalculateDtStressGradients(rFICVariables,rVariables);

    double StabilizationParameter =   rFICVariables.ElementLength * rFICVariables.ElementLength
                                    * (rVariables.BiotCoefficient  * (-PORE_PRESSURE_SIGN_FACTOR)) /(8.0*rFICVariables.ShearModulus);

    noalias(rVariables.PVector) = StabilizationParameter/3.0 * prod(rVariables.GradNpT,rFICVariables.DimVector)
                                                             * rVariables.IntegrationCoefficient;

    //Distribute DtStressGradient block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector< TDim, TNumNodes >(rRightHandSideVector,rVariables.PVector);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    CalculateDtStressGradients( FICElementVariables& rFICVariables,
                                const ElementVariables& Variables )
{
    for (unsigned int i = 0; i < TDim; i++)
    {
        for (unsigned int j = 0; j < TDim; j++)
        {
            rFICVariables.DtStressGradients[i][j] = 0.0;

            for (unsigned int k = 0; k < TNumNodes; k++)
                rFICVariables.DtStressGradients[i][j] += Variables.GradNpT(k,j)*mNodalDtStress[i][k];
        }
    }

    for (unsigned int i = 0; i < TDim; i++)
    {
        rFICVariables.DimVector[i] = 0.0;

        for (unsigned int j = 0; j < TDim; j++)
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
void UPwSmallStrainFICElement<TDim,TNumNodes>::
    CalculateAndAddPressureGradientFlow( VectorType& rRightHandSideVector,
                                         ElementVariables& rVariables,
                                         FICElementVariables& rFICVariables )
{
    double SignBiotCoefficient = -PORE_PRESSURE_SIGN_FACTOR * rVariables.BiotCoefficient;
    double StabilizationParameter =  rFICVariables.ElementLength * rFICVariables.ElementLength
                                    *SignBiotCoefficient / (8.0*rFICVariables.ShearModulus);

    noalias(rVariables.PMatrix) =  StabilizationParameter
                                 * (SignBiotCoefficient -2.0* rFICVariables.ShearModulus
                                                            * rVariables.BiotModulusInverse/(3.0*SignBiotCoefficient))
                                 * prod(rVariables.GradNpT, trans(rVariables.GradNpT))
                                 * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = -1.0*prod(rVariables.PMatrix,rVariables.DtPressureVector);

    //Distribute PressureGradient block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector< TDim, TNumNodes >(rRightHandSideVector,rVariables.PVector);
}

//----------------------------------------------------------------------------------------------------

template class UPwSmallStrainFICElement<2,3>;
template class UPwSmallStrainFICElement<2,4>;
template class UPwSmallStrainFICElement<3,4>;
template class UPwSmallStrainFICElement<3,8>;

} // Namespace Kratos
