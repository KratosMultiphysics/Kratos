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
#include "custom_elements/U_Pw_small_strain_element.hpp"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwSmallStrainElement<TDim,TNumNodes>::Create(IndexType NewId,
                                                               NodesArrayType const& ThisNodes,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwSmallStrainElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwSmallStrainElement<TDim,TNumNodes>::Create(IndexType NewId,
                                                               GeometryType::Pointer pGeom,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwSmallStrainElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
int UPwSmallStrainElement<TDim,TNumNodes>::
    Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::Check()") << std::endl;

    // Base class checks for positive area and Id > 0
    // Verify generic variables
    int ierr = UPwBaseElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();

    if (Geom.DomainSize() < 1.0e-15)
        KRATOS_THROW_ERROR( std::logic_error, "DomainSize < 1.0e-15 for the element ", this->Id() )

    // Verify specific properties
    bool IgnoreUndrained = false;
    if (Prop.Has(IGNORE_UNDRAINED))
        IgnoreUndrained = Prop[IGNORE_UNDRAINED];

    if (!IgnoreUndrained)
    {
        if ( BULK_MODULUS_FLUID.Key() == 0 || Prop.Has( BULK_MODULUS_FLUID ) == false || Prop[BULK_MODULUS_FLUID] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "BULK_MODULUS_FLUID has Key zero, is not defined or has an invalid value at element",
                                this->Id() )

        if ( DYNAMIC_VISCOSITY.Key() == 0 || Prop.Has( DYNAMIC_VISCOSITY ) == false || Prop[DYNAMIC_VISCOSITY] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "DYNAMIC_VISCOSITY has Key zero, is not defined or has an invalid value at element",
                                this->Id() )

        if ( PERMEABILITY_XX.Key() == 0 || Prop.Has( PERMEABILITY_XX ) == false || Prop[PERMEABILITY_XX] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "PERMEABILITY_XX has Key zero, is not defined or has an invalid value at element",
                                this->Id() )

        if ( PERMEABILITY_YY.Key() == 0 || Prop.Has( PERMEABILITY_YY ) == false || Prop[PERMEABILITY_YY] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "PERMEABILITY_YY has Key zero, is not defined or has an invalid value at element",
                                this->Id() )

        if ( PERMEABILITY_XY.Key() == 0 || Prop.Has( PERMEABILITY_XY ) == false || Prop[PERMEABILITY_XY] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "PERMEABILITY_XY has Key zero, is not defined or has an invalid value at element",
                                this->Id() )
        if (TDim > 2)
        {
            if ( PERMEABILITY_ZZ.Key() == 0 || Prop.Has( PERMEABILITY_ZZ ) == false || Prop[PERMEABILITY_ZZ] < 0.0 )
                KRATOS_THROW_ERROR( std::invalid_argument,
                                    "PERMEABILITY_ZZ has Key zero, is not defined or has an invalid value at element", this->Id() )

            if ( PERMEABILITY_YZ.Key() == 0 || Prop.Has( PERMEABILITY_YZ ) == false || Prop[PERMEABILITY_YZ] < 0.0 )
                KRATOS_THROW_ERROR( std::invalid_argument,
                                    "PERMEABILITY_YZ has Key zero, is not defined or has an invalid value at element", this->Id() )

            if ( PERMEABILITY_ZX.Key() == 0 || Prop.Has( PERMEABILITY_ZX ) == false || Prop[PERMEABILITY_ZX] < 0.0 )
                KRATOS_THROW_ERROR( std::invalid_argument,
                                    "PERMEABILITY_ZX has Key zero, is not defined or has an invalid value at element", this->Id() )
        }
    }


    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW )) 
        << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strainSize = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    if ( TDim == 2 ) {
        KRATOS_ERROR_IF( strainSize < 3 || strainSize > 4) 
        << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 (el id = ) "
        << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strainSize == 6)
        << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) "
        <<  this->Id() << std::endl;
    }

    // Check constitutive law
    if ( mConstitutiveLawVector.size() > 0 ) {
        return mConstitutiveLawVector[0]->Check( Prop, Geom, rCurrentProcessInfo );
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::Check()") << std::endl;

    return ierr;

    KRATOS_CATCH( "" );
}

//-------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    UpdateElementalVariableStressVector(ElementVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::UpdateElementalVariableStressVector()") << std::endl;

    for (unsigned int i=0; i < rVariables.StressVector.size(); ++i)
    {
        rVariables.StressVector(i) = mStressVector[PointNumber][i];
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::UpdateElementalVariableStressVector()") << std::endl;

    KRATOS_CATCH("");
}

//-------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    UpdateElementalVariableStressVector(Vector &StressVector, unsigned int PointNumber)
{
    KRATOS_TRY;
    // KRATOS_INFO("01-UPwSmallStrainElement::UpdateElementalVariableStressVector()") << std::endl;

    for (unsigned int i=0; i < StressVector.size(); ++i)
    {
        StressVector(i) = mStressVector[PointNumber][i];
    }
    // KRATOS_INFO("11-UPwSmallStrainElement::UpdateElementalVariableStressVector()") << std::endl;

    KRATOS_CATCH("");
}

//-------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    UpdateStressVector(const ElementVariables &rVariables, unsigned int PointNumber)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::UpdateStressVector()") << std::endl;

    for (unsigned int i=0; i < mStressVector[PointNumber].size(); ++i)
    {
        mStressVector[PointNumber][i] = rVariables.StressVector(i);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::UpdateStressVector()") << std::endl;

    KRATOS_CATCH("");
}

//-------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    UpdateStressVector(const Vector &StressVector, unsigned int PointNumber)
{
    KRATOS_TRY;
    // KRATOS_INFO("01-UPwSmallStrainElement::UpdateStressVector()") << std::endl;

    for (unsigned int i=0; i < mStressVector[PointNumber].size(); ++i)
    {
        mStressVector[PointNumber][i] = StressVector(i);
    }

    // KRATOS_INFO("11-UPwSmallStrainElement::UpdateStressVector()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeSolutionStep()") << std::endl;

    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );
    const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,this->GetIntegrationMethod());

    unsigned int VoigtSize = VOIGT_SIZE_3D;
    if (TDim == 2) VoigtSize = VOIGT_SIZE_2D_PLANE_STRESS;
    Matrix B(VoigtSize,TNumNodes*TDim);
    noalias(B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);

    //Create constitutive law parameters:
    Vector StrainVector(VoigtSize);
    Matrix ConstitutiveMatrix(VoigtSize,VoigtSize);

    unsigned int VoigtSizeStress = VOIGT_SIZE_3D;
    if (TDim == 2) VoigtSizeStress = VOIGT_SIZE_2D_PLANE_STRAIN;
    Vector StressVector(VoigtSizeStress);

    Vector Np(TNumNodes);
    Matrix GradNpT(TNumNodes,TDim);
    Matrix F = identity_matrix<double>(TDim);
    double detF = 1.0;
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,this->GetProperties(),rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); //Note: this is for nonlocal damage
    ConstitutiveParameters.SetConstitutiveMatrix(ConstitutiveMatrix);
    ConstitutiveParameters.SetStressVector(StressVector);
    ConstitutiveParameters.SetStrainVector(StrainVector);
    ConstitutiveParameters.SetShapeFunctionsValues(Np);
    ConstitutiveParameters.SetShapeFunctionsDerivatives(GradNpT);
    ConstitutiveParameters.SetDeformationGradientF(F);
    ConstitutiveParameters.SetDeterminantF(detF);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        noalias(Np) = row(NContainer,GPoint);
        noalias(GradNpT) = DN_DXContainer[GPoint];
        this->CalculateBMatrix(B, GradNpT);

        // Compute Stress
        noalias(StrainVector) = prod(B, DisplacementVector);
        UpdateElementalVariableStressVector(StressVector, GPoint);
        mConstitutiveLawVector[GPoint]->InitializeMaterialResponseCauchy(ConstitutiveParameters);
    }
    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeSolutionStep()") << std::endl;
    KRATOS_CATCH("");
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeNonLinearIteration()") << std::endl;

    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );
    const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,this->GetIntegrationMethod());

    unsigned int VoigtSize = VOIGT_SIZE_3D;
    if (TDim == 2) VoigtSize = VOIGT_SIZE_2D_PLANE_STRESS;
    Matrix B(VoigtSize,TNumNodes*TDim);
    noalias(B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector, Geom, DISPLACEMENT);

    //Create constitutive law parameters:
    Vector StrainVector(VoigtSize);
    Matrix ConstitutiveMatrix(VoigtSize,VoigtSize);

    unsigned int VoigtSizeStress = VOIGT_SIZE_3D;
    if (TDim == 2) VoigtSizeStress = VOIGT_SIZE_2D_PLANE_STRAIN;
    Vector StressVector(VoigtSizeStress);
    Vector Np(TNumNodes);
    Matrix GradNpT(TNumNodes,TDim);
    Matrix F = identity_matrix<double>(TDim);
    double detF = 1.0;
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,this->GetProperties(),rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); //Note: this is for nonlocal damage
    ConstitutiveParameters.SetConstitutiveMatrix(ConstitutiveMatrix);
    ConstitutiveParameters.SetStressVector(StressVector);
    ConstitutiveParameters.SetStrainVector(StrainVector);
    ConstitutiveParameters.SetShapeFunctionsValues(Np);
    ConstitutiveParameters.SetShapeFunctionsDerivatives(GradNpT);
    ConstitutiveParameters.SetDeformationGradientF(F);
    ConstitutiveParameters.SetDeterminantF(detF);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        noalias(Np) = row(NContainer,GPoint);
        noalias(GradNpT) = DN_DXContainer[GPoint];
        this->CalculateBMatrix(B, GradNpT);

        // Compute Stress
        noalias(StrainVector) = prod(B, DisplacementVector);
        UpdateElementalVariableStressVector(StressVector, GPoint);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
        UpdateStressVector(StressVector, GPoint);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeNonLinearIteration()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::FinalizeNonLinearIteration()") << std::endl;

    this->InitializeNonLinearIteration(rCurrentProcessInfo);

    // KRATOS_INFO("1-UPwSmallStrainElement::FinalizeNonLinearIteration()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::FinalizeSolutionStep()") << std::endl;

    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );
    const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,this->GetIntegrationMethod());

    unsigned int VoigtSize = VOIGT_SIZE_3D;
    if (TDim == 2) VoigtSize = VOIGT_SIZE_2D_PLANE_STRESS;
    Matrix B(VoigtSize,TNumNodes*TDim);
    noalias(B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector, Geom, DISPLACEMENT);

    //Create constitutive law parameters:
    Vector StrainVector(VoigtSize);
    Matrix ConstitutiveMatrix(VoigtSize,VoigtSize);

    unsigned int VoigtSizeStress = VOIGT_SIZE_3D;
    if (TDim == 2) VoigtSizeStress = VOIGT_SIZE_2D_PLANE_STRESS;
    Vector StressVector(VoigtSizeStress);

    Vector Np(TNumNodes);
    Matrix GradNpT(TNumNodes,TDim);
    Matrix F = identity_matrix<double>(TDim);
    double detF = 1.0;
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,this->GetProperties(),rCurrentProcessInfo);
    ConstitutiveParameters.SetConstitutiveMatrix(ConstitutiveMatrix);
    ConstitutiveParameters.SetStressVector(StressVector);
    ConstitutiveParameters.SetStrainVector(StrainVector);
    ConstitutiveParameters.SetShapeFunctionsValues(Np);
    ConstitutiveParameters.SetShapeFunctionsDerivatives(GradNpT);
    ConstitutiveParameters.SetDeterminantF(detF);
    ConstitutiveParameters.SetDeformationGradientF(F);
    //ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);

    if (rCurrentProcessInfo[NODAL_SMOOTHING] == true)
    {
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        Matrix StressContainer(NumGPoints, StressVector.size());

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            this->CalculateBMatrix(B, GradNpT);

            noalias(StrainVector) = prod(B, DisplacementVector);

            noalias(Np) = row(NContainer,GPoint);

            //compute constitutive tensor and/or stresses
            UpdateElementalVariableStressVector(StressVector, GPoint);
            mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
            mStateVariablesFinalized[GPoint] = 
                mConstitutiveLawVector[GPoint]->GetValue( STATE_VARIABLES,
                                                          mStateVariablesFinalized[GPoint] );

            this->SaveGPStress(StressContainer,StressVector,StressVector.size(),GPoint);
        }
        this->ExtrapolateGPValues(StressContainer,StressVector.size());
    }
    else
    {
        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            this->CalculateBMatrix(B, GradNpT);

            noalias(StrainVector) = prod(B, DisplacementVector);

            noalias(Np) = row(NContainer,GPoint);

            //compute constitutive tensor and/or stresses
            UpdateElementalVariableStressVector(StressVector, GPoint);
            mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
            mStateVariablesFinalized[GPoint] = 
                mConstitutiveLawVector[GPoint]->GetValue( STATE_VARIABLES,
                                                          mStateVariablesFinalized[GPoint] );
        }
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::FinalizeSolutionStep()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::SaveGPStress(Matrix& rStressContainer,
                                                         const Vector& StressVector,
                                                         const unsigned int& VoigtSize,
                                                         const unsigned int& GPoint)
{

    KRATOS_TRY

    // KRATOS_INFO("0-UPwSmallStrainElement::SaveGPStress()") << std::endl;

    for (unsigned int i = 0; i < VoigtSize; i++)
    {
        rStressContainer(GPoint,i) = StressVector[i];
    }

    /* INFO: (Quadrilateral_2D_4 with GI_GAUSS_2)
     *
     *                      |S0-0 S1-0 S2-0|
     * rStressContainer =   |S0-1 S1-1 S2-1|
     *                      |S0-2 S1-2 S2-2|
     *                      |S0-3 S1-3 S2-3|
     *
     * S1-0 = S[1] at GP 0
    */
    // KRATOS_INFO("1-UPwSmallStrainElement::SaveGPStress()") << std::endl;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim, TNumNodes>::
    ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "ExtrapolateGPValues is not implemented for this element", "" )

    KRATOS_CATCH( "" )
}

/*
//----------------------------------------------------------------------------------------
template< >
void UPwSmallStrainElement<2,8>::
    ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "ExtrapolateGPValues is not implemented for 2D8N elemnents", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< >
void UPwSmallStrainElement<2,9>::
    ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "ExtrapolateGPValues is not implemented for 2D9N elemnents", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< >
void UPwSmallStrainElement<3,10>::
    ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "ExtrapolateGPValues is not implemented for 3D10N elemnents", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< >
void UPwSmallStrainElement<3,20>::
    ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "ExtrapolateGPValues is not implemented for 3D20N elemnents", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< >
void UPwSmallStrainElement<3,27>::
    ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "ExtrapolateGPValues is not implemented for 3D27N elemnents", "" )

    KRATOS_CATCH( "" )
}
*/

//----------------------------------------------------------------------------------------
template< >
void UPwSmallStrainElement<2,3>::
    ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement<2,3>::ExtrapolateGPValues()") << std::endl;

    // Triangle_2d_3 with GI_GAUSS_2

    array_1d<double,3> DamageContainer; // 3 GPoints

    for ( unsigned int i = 0;  i < 3; i++ ) //NumGPoints
    {
        DamageContainer[i] = 0.0;
        DamageContainer[i] = mConstitutiveLawVector[i]->GetValue( DAMAGE_VARIABLE, DamageContainer[i] );
    }

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area();
    array_1d<Vector,3> NodalStressVector; //List with stresses at each node
    array_1d<Matrix,3> NodalStressTensor;

    for (unsigned int Node = 0; Node < 3; Node++)
    {
        NodalStressVector[Node].resize(VoigtSize);
        NodalStressTensor[Node].resize(2,2);
    }

    BoundedMatrix<double,3,3> ExtrapolationMatrix;
    GeoElementUtilities::Calculate2DExtrapolationMatrix(ExtrapolationMatrix);


    BoundedMatrix<double, 3, 3> AuxNodalStress;
    noalias(AuxNodalStress) = prod(ExtrapolationMatrix, StressContainer);

    /* INFO:
        *
        *                  |S0-0 S1-0 S2-0|
        * AuxNodalStress = |S0-1 S1-1 S2-1|
        *                  |S0-2 S1-2 S2-2|
        *
        * S1-0 = S[1] at node 0
    */

    array_1d<double,3> NodalDamage;
    noalias(NodalDamage) = prod(ExtrapolationMatrix,DamageContainer);

    for (unsigned int i = 0; i < 3; i++) //TNumNodes
    {
        noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) += NodalDamage[i]*Area;
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    // KRATOS_INFO("1-UPwSmallStrainElement<2,3>::ExtrapolateGPValues()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainElement<2,4>::
    ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement<2,4>::ExtrapolateGPValues()") << std::endl;

    // Quadrilateral_2d_4 with GI_GAUSS_2

    array_1d<double,4> DamageContainer; // 4 GPoints

    for ( unsigned int i = 0;  i < 4; i++ ) //NumGPoints
    {
        DamageContainer[i] = 0.0;
        DamageContainer[i] = mConstitutiveLawVector[i]->GetValue( DAMAGE_VARIABLE, DamageContainer[i] );
    }

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area();
    array_1d<Vector,4> NodalStressVector; //List with stresses at each node
    array_1d<Matrix,4> NodalStressTensor;

    for (unsigned int Node = 0; Node < 4; Node ++)
    {
        NodalStressVector[Node].resize(VoigtSize);
        NodalStressTensor[Node].resize(2,2);
    }

    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    GeoElementUtilities::Calculate2DExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double,4,3> AuxNodalStress;
    noalias(AuxNodalStress) = prod(ExtrapolationMatrix,StressContainer);

    /* INFO:
        *
        *                  |S0-0 S1-0 S2-0|
        * AuxNodalStress = |S0-1 S1-1 S2-1|
        *                  |S0-2 S1-2 S2-2|
        *                  |S0-3 S1-3 S2-3|
        *
        * S1-0 = S[1] at node 0
    */

    array_1d<double,4> NodalDamage; //List with stresses at each node
    noalias(NodalDamage) = prod(ExtrapolationMatrix,DamageContainer);

    for (unsigned int i = 0; i < 4; i++) //TNumNodes
    {
        noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) += NodalDamage[i]*Area;
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    // KRATOS_INFO("1-UPwSmallStrainElement<2,4>::ExtrapolateGPValues()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainElement<3,4>::
    ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement<3,4>::ExtrapolateGPValues()") << std::endl;

    // Tetrahedra_3d_4 with GI_GAUSS_2

    array_1d<double,4> DamageContainer; // 4 GPoints

    for ( unsigned int i = 0;  i < 4; i++ ) //NumGPoints
    {
        DamageContainer[i] = 0.0;
        DamageContainer[i] = mConstitutiveLawVector[i]->GetValue( DAMAGE_VARIABLE, DamageContainer[i] );
    }

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area(); // In 3D this is Volume
    array_1d<Vector,4> NodalStressVector; //List with stresses at each node
    array_1d<Matrix,4> NodalStressTensor;

    for (unsigned int Node = 0; Node < 4; Node ++)
    {
        NodalStressVector[Node].resize(VoigtSize);
        NodalStressTensor[Node].resize(3,3);
    }

    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    GeoElementUtilities::Calculate3DExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double,4,6> AuxNodalStress;
    noalias(AuxNodalStress) = prod(ExtrapolationMatrix,StressContainer);

    array_1d<double,4> NodalDamage;
    noalias(NodalDamage) = prod(ExtrapolationMatrix,DamageContainer);

    for (unsigned int i = 0; i < 4; i++) //TNumNodes
    {
        noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) += NodalDamage[i]*Area;
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    // KRATOS_INFO("1-UPwSmallStrainElement<3,4>::ExtrapolateGPValues()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainElement<3,8>::
    ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement<3,8>::ExtrapolateGPValues()") << std::endl;

    // Hexahedra_3d_8 with GI_GAUSS_2

    array_1d<double,8> DamageContainer; // 8 GPoints

    for ( unsigned int i = 0;  i < 8; i++ ) //NumGPoints
    {
        DamageContainer[i] = 0.0;
        DamageContainer[i] = mConstitutiveLawVector[i]->GetValue( DAMAGE_VARIABLE, DamageContainer[i] );
    }

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area(); // In 3D this is Volume
    array_1d<Vector,8> NodalStressVector; //List with stresses at each node
    array_1d<Matrix,8> NodalStressTensor;

    for (unsigned int Node = 0; Node < 8; Node ++)
    {
        NodalStressVector[Node].resize(VoigtSize);
        NodalStressTensor[Node].resize(3,3);
    }

    BoundedMatrix<double,8,8> ExtrapolationMatrix;
    GeoElementUtilities::Calculate3DExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double,8,6> AuxNodalStress;
    noalias(AuxNodalStress) = prod(ExtrapolationMatrix,StressContainer);

    array_1d<double,8> NodalDamage; //List with stresses at each node
    noalias(NodalDamage) = prod(ExtrapolationMatrix,DamageContainer);

    for (unsigned int i = 0; i < 8; i++) //TNumNodes
    {
        noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) += NodalDamage[i]*Area;
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    // KRATOS_INFO("1-UPwSmallStrainElement<3,8>::ExtrapolateGPValues()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<double>& rVariable,
                                  std::vector<double>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateOnIntegrationPoints()") << std::endl;

    if (rVariable == VON_MISES_STRESS)
    {
        //Defining necessary variables
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );
        const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,this->GetIntegrationMethod());

        unsigned int VoigtSize = VOIGT_SIZE_3D;
        if (TDim == 2) VoigtSize = VOIGT_SIZE_2D_PLANE_STRESS;
        Matrix B(VoigtSize,TNumNodes*TDim);
        noalias(B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);

        //Create constitutive law parameters:
        Vector StrainVector(VoigtSize);
        Matrix ConstitutiveMatrix(VoigtSize,VoigtSize);

        unsigned int VoigtSizeStress = VOIGT_SIZE_3D;
        if (TDim == 2) VoigtSizeStress = VOIGT_SIZE_2D_PLANE_STRAIN;
        Vector StressVector(VoigtSizeStress);

        Vector Np(TNumNodes);
        Matrix GradNpT(TNumNodes,TDim);
        Matrix F = identity_matrix<double>(TDim);
        double detF = 1.0;
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,this->GetProperties(),rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveParameters.SetConstitutiveMatrix(ConstitutiveMatrix);
        ConstitutiveParameters.SetStressVector(StressVector);
        ConstitutiveParameters.SetStrainVector(StrainVector);
        ConstitutiveParameters.SetShapeFunctionsValues(Np);
        ConstitutiveParameters.SetShapeFunctionsDerivatives(GradNpT);
        ConstitutiveParameters.SetDeterminantF(detF);
        ConstitutiveParameters.SetDeformationGradientF(F);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            this->CalculateBMatrix(B, GradNpT);

            noalias(StrainVector) = prod(B,DisplacementVector);

            noalias(Np) = row(NContainer,GPoint);

            //compute constitutive tensor and/or stresses
            UpdateElementalVariableStressVector(StressVector, GPoint);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
            UpdateStressVector(StressVector, GPoint);

            ComparisonUtilities EquivalentStress;
            rOutput[GPoint] = EquivalentStress.CalculateVonMises(StressVector);
        }
    }
    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<array_1d<double,3>>& rVariable,
                                  std::vector<array_1d<double,3>>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("01-UPwSmallStrainElement::CalculateOnIntegrationPoints()") << std::endl;

    if (rVariable == FLUID_FLUX_VECTOR)
    {
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );

        //Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,this->GetIntegrationMethod());

        //Defining necessary variables
        array_1d<double,TNumNodes> PressureVector;
        for (unsigned int i=0; i<TNumNodes; i++)
            PressureVector[i] = Geom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        array_1d<double,TNumNodes*TDim> VolumeAcceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(VolumeAcceleration,Geom,VOLUME_ACCELERATION);
        array_1d<double,TDim> BodyAcceleration;
        BoundedMatrix<double,TNumNodes, TDim> GradNpT;
        BoundedMatrix<double,TDim, TDim> PermeabilityMatrix;
        GeoElementUtilities::CalculatePermeabilityMatrix(PermeabilityMatrix,Prop);
        const double& DynamicViscosityInverse = 1.0/Prop[DYNAMIC_VISCOSITY];
        const double& FluidDensity = Prop[DENSITY_WATER];
        array_1d<double,TDim> GradPressureTerm;
        array_1d<double,TDim> FluidFlux;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            GeoElementUtilities::
                InterpolateVariableWithComponents<TDim, TNumNodes>( BodyAcceleration,
                                                                    NContainer,
                                                                    VolumeAcceleration,
                                                                    GPoint );

            noalias(GradPressureTerm) = prod(trans(GradNpT),PressureVector);
            noalias(GradPressureTerm) += -FluidDensity*BodyAcceleration;

            noalias(FluidFlux) = -DynamicViscosityInverse*prod(PermeabilityMatrix,GradPressureTerm);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint],FluidFlux);
        }
    }
    // KRATOS_INFO("11-UPwSmallStrainElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                 std::vector<Matrix>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateOnIntegrationPoints()") << std::endl;

    if (rVariable == CAUCHY_STRESS_TENSOR)
    {
        //Defining necessary variables
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );
        const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,this->GetIntegrationMethod());

        unsigned int VoigtSize = VOIGT_SIZE_3D;
        if (TDim == 2) VoigtSize = VOIGT_SIZE_2D_PLANE_STRESS;
        Matrix B(VoigtSize,TNumNodes*TDim);
        noalias(B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);

        //Create constitutive law parameters:
        Vector StrainVector(VoigtSize);
        Matrix ConstitutiveMatrix(VoigtSize,VoigtSize);

        unsigned int VoigtSizeStress = VOIGT_SIZE_3D;
        if (TDim == 2) VoigtSizeStress = VOIGT_SIZE_2D_PLANE_STRAIN;
        Vector StressVector(VoigtSizeStress);

        Vector Np(TNumNodes);
        Matrix GradNpT(TNumNodes,TDim);
        Matrix F = identity_matrix<double>(TDim);
        double detF = 1.0;
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,this->GetProperties(),rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveParameters.SetConstitutiveMatrix(ConstitutiveMatrix);
        ConstitutiveParameters.SetStressVector(StressVector);
        ConstitutiveParameters.SetStrainVector(StrainVector);
        ConstitutiveParameters.SetShapeFunctionsValues(Np);
        ConstitutiveParameters.SetShapeFunctionsDerivatives(GradNpT);
        ConstitutiveParameters.SetDeterminantF(detF);
        ConstitutiveParameters.SetDeformationGradientF(F);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            this->CalculateBMatrix(B, GradNpT);

            noalias(StrainVector) = prod(B,DisplacementVector);

            noalias(Np) = row(NContainer,GPoint);

            //compute constitutive tensor and/or stresses
            UpdateElementalVariableStressVector(StressVector, GPoint);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
            UpdateStressVector(StressVector, GPoint);

            rOutput[GPoint].resize(TDim,TDim,false );
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector);
        }
    }
    else if (rVariable == TOTAL_STRESS_TENSOR)
    {
        //Defining necessary variables
        const GeometryType& Geom = this->GetGeometry();
        const PropertiesType& Prop = this->GetProperties();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );
        const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,this->GetIntegrationMethod());

        unsigned int VoigtSizeStress = VOIGT_SIZE_3D;
        if (TDim == 2) VoigtSizeStress = VOIGT_SIZE_2D_PLANE_STRAIN;

        unsigned int VoigtSize = VOIGT_SIZE_3D;
        if (TDim == 2) VoigtSize = VOIGT_SIZE_2D_PLANE_STRESS;

        Vector VoigtVector = ZeroVector(VoigtSizeStress);

        if (VoigtVector.size() == VOIGT_SIZE_3D)
        {
            VoigtVector[INDEX_3D_XX] = 1.0;
            VoigtVector[INDEX_3D_YY] = 1.0;
            VoigtVector[INDEX_3D_ZZ] = 1.0;
        }
        else if (VoigtVector.size() == VOIGT_SIZE_2D_PLANE_STRAIN)
        {
            VoigtVector[INDEX_2D_PLANE_STRAIN_XX] = 1.0;
            VoigtVector[INDEX_2D_PLANE_STRAIN_YY] = 1.0;
            VoigtVector[INDEX_2D_PLANE_STRAIN_ZZ] = 1.0;
        }
        else
        {
            VoigtVector[INDEX_2D_PLANE_STRESS_XX] = 1.0;
            VoigtVector[INDEX_2D_PLANE_STRESS_YY] = 1.0;
        }


        Matrix B(VoigtSize,TNumNodes*TDim);
        noalias(B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);

        array_1d<double,TNumNodes> PressureVector;
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            PressureVector[i] = Geom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        }

        //Create constitutive law parameters:
        Vector StrainVector(VoigtSize);
        Vector StressVector(VoigtSizeStress);
        Matrix ConstitutiveMatrix(VoigtSize,VoigtSize);
        Vector Np(TNumNodes);
        Matrix GradNpT(TNumNodes,TDim);
        Matrix F = identity_matrix<double>(TDim);
        double detF = 1.0;
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveParameters.SetConstitutiveMatrix(ConstitutiveMatrix);
        ConstitutiveParameters.SetStressVector(StressVector);
        ConstitutiveParameters.SetStrainVector(StrainVector);
        ConstitutiveParameters.SetShapeFunctionsValues(Np);
        ConstitutiveParameters.SetShapeFunctionsDerivatives(GradNpT);
        ConstitutiveParameters.SetDeterminantF(detF);
        ConstitutiveParameters.SetDeformationGradientF(F);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            this->CalculateBMatrix(B, GradNpT);

            noalias(StrainVector) = prod(B,DisplacementVector);

            noalias(Np) = row(NContainer,GPoint);

            //compute constitutive tensor and/or stresses
            UpdateElementalVariableStressVector(StressVector, GPoint);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
            UpdateStressVector(StressVector, GPoint);

            const double BulkModulus = CalculateBulkModulus(ConstitutiveMatrix);
            const double BiotCoefficient = 1.0 - BulkModulus / Prop[BULK_MODULUS_SOLID];

            double Pressure = 0.0;
            for (unsigned int node = 0; node < TNumNodes; ++node)
            {
                Pressure += NContainer(GPoint, node)*PressureVector[node];
            }

            noalias(StressVector) += signFactor * BiotCoefficient*Pressure*VoigtVector;

            rOutput[GPoint].resize(TDim,TDim,false );
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector);
        }
    }
    else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
    {
        //Defining necessary variables
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( this->GetIntegrationMethod() );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,this->GetIntegrationMethod());

        unsigned int VoigtSize = VOIGT_SIZE_3D;
        if (TDim == 2) VoigtSize = VOIGT_SIZE_2D_PLANE_STRESS;
        Matrix B(VoigtSize,TNumNodes*TDim);
        noalias(B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(DisplacementVector,Geom,DISPLACEMENT);
        Vector StrainVector(VoigtSize);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            this->CalculateBMatrix(B, DN_DXContainer[GPoint]);

            noalias(StrainVector) = prod(B,DisplacementVector);

            rOutput[GPoint].resize(TDim,TDim,false );
            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector);
        }
    }
    else if (rVariable == PERMEABILITY_MATRIX)
    {
        const unsigned int NumGPoints = this->GetGeometry().IntegrationPointsNumber( this->GetIntegrationMethod() );

        //If the permeability of the element is a given property
        BoundedMatrix<double,TDim,TDim> PermeabilityMatrix;
        GeoElementUtilities::CalculatePermeabilityMatrix(PermeabilityMatrix,this->GetProperties());

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            rOutput[GPoint].resize(TDim,TDim,false);
            noalias(rOutput[GPoint]) = PermeabilityMatrix;
        }
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateMaterialStiffnessMatrix( MatrixType& rStiffnessMatrix,
                                      const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateMaterialStiffnessMatrix()") << std::endl;

    const unsigned int N_DOF = TNumNodes * (TDim + 1);

    //Resizing mass matrix
    if ( rStiffnessMatrix.size1() != N_DOF )
        rStiffnessMatrix.resize( N_DOF, N_DOF, false );
    noalias( rStiffnessMatrix ) = ZeroMatrix( N_DOF, N_DOF );

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( this->GetIntegrationMethod() );
    const unsigned int NumGPoints = IntegrationPoints.size();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Vector detJContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, detJContainer, this->GetIntegrationMethod());

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,CurrentProcessInfo);
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
    for (unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute Np, GradNpT, B and StrainVector
        noalias(Variables.Np) = row(NContainer,GPoint);
        noalias(Variables.GradNpT) = DN_DXContainer[GPoint];
        this->CalculateBMatrix(Variables.B, Variables.GradNpT);
        noalias(Variables.StrainVector) = prod(Variables.B, Variables.DisplacementVector);

        //Compute constitutive tensor
        UpdateElementalVariableStressVector(Variables, GPoint);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        // calculate Bulk modulus from stiffness matrix
        const double BulkModulus = CalculateBulkModulus(Variables.ConstitutiveMatrix);
        this->InitializeBiotCoefficients(Variables, Prop, BulkModulus);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient,
                                              detJContainer[GPoint],
                                              IntegrationPoints[GPoint].Weight());

        //Compute stiffness matrix
        this->CalculateAndAddStiffnessMatrix(rStiffnessMatrix, Variables);
    }
    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateMaterialStiffnessMatrix()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo,
                  const bool CalculateStiffnessMatrixFlag,
                  const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAll()") << std::endl;

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( this->GetIntegrationMethod() );
    const unsigned int NumGPoints = IntegrationPoints.size();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( this->GetIntegrationMethod() );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Vector detJContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients( DN_DXContainer,
                                                   detJContainer,
                                                   this->GetIntegrationMethod() );

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom, Prop, CurrentProcessInfo);
    if (CalculateStiffnessMatrixFlag) ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag)  ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables( Variables,
                                      ConstitutiveParameters,
                                      Geom,
                                      Prop,
                                      CurrentProcessInfo );

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        noalias(Variables.GradNpT) = DN_DXContainer[GPoint];

        this->CalculateBMatrix(Variables.B, Variables.GradNpT);

        noalias(Variables.StrainVector) = prod(Variables.B, Variables.DisplacementVector);

        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, NContainer, GPoint);
        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim, TNumNodes>( Variables.BodyAcceleration,
                                                                NContainer,
                                                                Variables.VolumeAcceleration,
                                                                GPoint );


        //Compute constitutive tensor and stresses
        UpdateElementalVariableStressVector(Variables, GPoint);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
        UpdateStressVector(Variables, GPoint);

        // calculate Bulk modulus from stiffness matrix
        const double BulkModulus = CalculateBulkModulus(Variables.ConstitutiveMatrix);
        this->InitializeBiotCoefficients(Variables, Prop, BulkModulus);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient,
                                              detJContainer[GPoint],
                                              IntegrationPoints[GPoint].Weight());

        //Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if (CalculateResidualVectorFlag)  this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAll()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    InitializeBiotCoefficients( ElementVariables& rVariables,
                               const PropertiesType& Prop,
                               const double &BulkModulus )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeBiotCoefficients()") << std::endl;

    //Properties variables
    rVariables.BiotCoefficient = 1.0 - BulkModulus / Prop[BULK_MODULUS_SOLID];
    rVariables.BiotModulusInverse =  (rVariables.BiotCoefficient - Prop[POROSITY])/Prop[BULK_MODULUS_SOLID]
                                   + Prop[POROSITY]/Prop[BULK_MODULUS_FLUID];

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeBiotCoefficients()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateBulkModulus(const Matrix &ConstitutiveMatrix)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateBulkModulus()") << std::endl;

    const int IndexG = ConstitutiveMatrix.size1() - 1;

    const double M = ConstitutiveMatrix(0, 0);
    const double G = ConstitutiveMatrix(IndexG, IndexG);

    const double BulkModulus = M - (4.0/3.0)*G;

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateBulkModulus()") << std::endl;

    return BulkModulus;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    InitializeElementVariables( ElementVariables& rVariables,
                                ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                const GeometryType& Geom,
                                const PropertiesType& Prop,
                                const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeElementVariables()") << std::endl;

    //Properties variables
    rVariables.IgnoreUndrained = false;
    if (Prop.Has(IGNORE_UNDRAINED))
        rVariables.IgnoreUndrained = Prop[IGNORE_UNDRAINED];

    rVariables.DynamicViscosityInverse = 1.0/Prop[DYNAMIC_VISCOSITY];
    rVariables.FluidDensity = Prop[DENSITY_WATER];
    rVariables.Density = Prop[POROSITY]*rVariables.FluidDensity + (1.0-Prop[POROSITY])*Prop[DENSITY_SOLID];
    GeoElementUtilities::CalculatePermeabilityMatrix(rVariables.PermeabilityMatrix, Prop);

    //ProcessInfo variables
    rVariables.VelocityCoefficient = CurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtPressureCoefficient = CurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    //Nodal Variables
    for (unsigned int i=0; i<TNumNodes; i++)
    {
        rVariables.PressureVector[i] = Geom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.DtPressureVector[i] = Geom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.DisplacementVector, Geom, DISPLACEMENT);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VelocityVector,     Geom, VELOCITY);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VolumeAcceleration, Geom, VOLUME_ACCELERATION);

    //General Variables
    unsigned int VoigtSizeStress;
    unsigned int VoigtSize;
    if (TDim > 2)
    {
        VoigtSize = VOIGT_SIZE_3D;
        VoigtSizeStress = VOIGT_SIZE_3D;
        rVariables.VoigtVector.resize(VoigtSize);
        rVariables.VoigtVector[INDEX_3D_XX] = 1.0;
        rVariables.VoigtVector[INDEX_3D_YY] = 1.0;
        rVariables.VoigtVector[INDEX_3D_ZZ] = 1.0;
        rVariables.VoigtVector[INDEX_3D_XY] = 0.0;
        rVariables.VoigtVector[INDEX_3D_YZ] = 0.0;
        rVariables.VoigtVector[INDEX_3D_XZ] = 0.0;
    }
    else
    {
        VoigtSize       = VOIGT_SIZE_2D_PLANE_STRESS;
        VoigtSizeStress = VOIGT_SIZE_2D_PLANE_STRAIN;
        rVariables.VoigtVector.resize(VoigtSize);
        rVariables.VoigtVector[INDEX_2D_PLANE_STRESS_XX] = 1.0;
        rVariables.VoigtVector[INDEX_2D_PLANE_STRESS_YY] = 1.0;
        rVariables.VoigtVector[INDEX_2D_PLANE_STRESS_XY] = 0.0;
    }

    //Variables computed at each GP
    rVariables.B.resize(VoigtSize,TNumNodes*TDim,false);
    noalias(rVariables.B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);
    noalias(rVariables.Nu) = ZeroMatrix(TDim, TNumNodes*TDim);
    //Constitutive Law parameters
    rVariables.StrainVector.resize(VoigtSize,false);
    rVariables.StressVector.resize(VoigtSizeStress,false);
    rVariables.ConstitutiveMatrix.resize(VoigtSize,VoigtSize,false);
    rVariables.Np.resize(TNumNodes,false);
    rVariables.GradNpT.resize(TNumNodes,TDim,false);
    rVariables.F.resize(TDim,TDim,false);
    rVariables.detF = 1.0;
    rConstitutiveParameters.SetStrainVector(rVariables.StrainVector);
    rConstitutiveParameters.SetStressVector(rVariables.StressVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rConstitutiveParameters.SetShapeFunctionsValues(rVariables.Np);
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rVariables.GradNpT);
    rConstitutiveParameters.SetDeformationGradientF(rVariables.F);
    rConstitutiveParameters.SetDeterminantF(rVariables.detF);
    //Auxiliary variables
    rVariables.UVoigtMatrix.resize(TNumNodes*TDim,VoigtSize,false);

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeElementVariables()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateBMatrix(Matrix& rB, const Matrix& GradNpT)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateBMatrix()") << std::endl;

    unsigned int index;

    if (TDim > 2)
    {
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            index = TDim * i;

            rB( 0, index + 0 ) = GradNpT( i, 0 );
            rB( 1, index + 1 ) = GradNpT( i, 1 );
            rB( 2, index + 2 ) = GradNpT( i, 2 );
            rB( 3, index + 0 ) = GradNpT( i, 1 );
            rB( 3, index + 1 ) = GradNpT( i, 0 );
            rB( 4, index + 1 ) = GradNpT( i, 2 );
            rB( 4, index + 2 ) = GradNpT( i, 1 );
            rB( 5, index + 0 ) = GradNpT( i, 2 );
            rB( 5, index + 2 ) = GradNpT( i, 0 );
        }
    }
    else
    {
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            index = TDim * i;

            rB( 0, index + 0 ) = GradNpT( i, 0 );
            rB( 1, index + 1 ) = GradNpT( i, 1 );
            rB( 2, index + 0 ) = GradNpT( i, 1 );
            rB( 2, index + 1 ) = GradNpT( i, 0 );
        }
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateBMatrix()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddLHS()") << std::endl;

    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix,rVariables);

    if (!rVariables.IgnoreUndrained)
    {
        this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix,rVariables);

        this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

        this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddLHS()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddStiffnessMatrix()") << std::endl;

    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B), rVariables.ConstitutiveMatrix);
    noalias(rVariables.UMatrix) = prod(rVariables.UVoigtMatrix, rVariables.B)*rVariables.IntegrationCoefficient;

    //Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUBlockMatrix<TDim, TNumNodes>(rLeftHandSideMatrix, rVariables.UMatrix);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddStiffnessMatrix()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddCouplingMatrix()") << std::endl;

    noalias(rVariables.UVector) = prod(trans(rVariables.B),rVariables.VoigtVector);

    noalias(rVariables.UPMatrix) = signFactor * rVariables.BiotCoefficient
                                              * outer_prod(rVariables.UVector,rVariables.Np)
                                              * rVariables.IntegrationCoefficient;

    //Distribute coupling block matrix into the elemental matrix
    GeoElementUtilities::AssembleUPBlockMatrix<TDim, TNumNodes>(rLeftHandSideMatrix,rVariables.UPMatrix);

    noalias(rVariables.PUMatrix) = signFactor * rVariables.VelocityCoefficient*trans(rVariables.UPMatrix);

    //Distribute transposed coupling block matrix into the elemental matrix
    GeoElementUtilities::AssemblePUBlockMatrix<TDim, TNumNodes>(rLeftHandSideMatrix,rVariables.PUMatrix);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddCouplingMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddCompressibilityMatrix()") << std::endl;

    noalias(rVariables.PMatrix) = - signFactor * rVariables.DtPressureCoefficient
                                               * rVariables.BiotModulusInverse
                                               * outer_prod(rVariables.Np, rVariables.Np)
                                               * rVariables.IntegrationCoefficient;

    //Distribute compressibility block matrix into the elemental matrix
    GeoElementUtilities::
        AssemblePBlockMatrix< TDim, TNumNodes >(rLeftHandSideMatrix,
                                                rVariables.PMatrix);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddCompressibilityMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddPermeabilityMatrix()") << std::endl;

    noalias(rVariables.PDimMatrix) = - signFactor * prod(rVariables.GradNpT,rVariables.PermeabilityMatrix);

    noalias(rVariables.PMatrix) =  rVariables.DynamicViscosityInverse
                                 * prod(rVariables.PDimMatrix,trans(rVariables.GradNpT))
                                 * rVariables.IntegrationCoefficient;

    //Distribute permeability block matrix into the elemental matrix
    GeoElementUtilities::AssemblePBlockMatrix< TDim, TNumNodes >(rLeftHandSideMatrix, rVariables.PMatrix);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddPermeabilityMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddRHS()") << std::endl;

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    if (!rVariables.IgnoreUndrained)
    {
        this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

        this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

        this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddRHS()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddStiffnessForce( VectorType& rRightHandSideVector,
                                   ElementVariables& rVariables )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddStiffnessForce()") << std::endl;

    if ( TDim > 2 )
    {
        noalias(rVariables.UVector) = -1.0 * prod(trans(rVariables.B),rVariables.StressVector)
                                           * rVariables.IntegrationCoefficient;
    }
    else
    {
        unsigned int VoigtSize = VOIGT_SIZE_2D_PLANE_STRESS;
        Vector StressVector;
        StressVector.resize(VoigtSize);

        StressVector[INDEX_2D_PLANE_STRESS_XX] = rVariables.StressVector(INDEX_2D_PLANE_STRAIN_XX);
        StressVector[INDEX_2D_PLANE_STRESS_YY] = rVariables.StressVector(INDEX_2D_PLANE_STRAIN_YY);
        StressVector[INDEX_2D_PLANE_STRESS_XY] = rVariables.StressVector(INDEX_2D_PLANE_STRAIN_XY);

        noalias(rVariables.UVector) = -1.0 * prod(trans(rVariables.B), StressVector)
                                           * rVariables.IntegrationCoefficient;
    }

    //Distribute stiffness block vector into elemental vector
    GeoElementUtilities::AssembleUBlockVector<TDim, TNumNodes>(rRightHandSideVector,rVariables.UVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddStiffnessForce()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddMixBodyForce()") << std::endl;

    noalias(rVariables.UVector) = rVariables.Density * prod(trans(rVariables.Nu),rVariables.BodyAcceleration)
                                                     * rVariables.IntegrationCoefficient;

    //Distribute body force block vector into elemental vector
    GeoElementUtilities::AssembleUBlockVector<TDim, TNumNodes>(rRightHandSideVector,rVariables.UVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddMixBodyForce()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddCouplingTerms()") << std::endl;

    noalias(rVariables.UVector) = prod(trans(rVariables.B),rVariables.VoigtVector);

    noalias(rVariables.UPMatrix) = - signFactor * rVariables.BiotCoefficient
                                                * outer_prod(rVariables.UVector,rVariables.Np)
                                                * rVariables.IntegrationCoefficient;

    noalias(rVariables.UVector) = prod(rVariables.UPMatrix,rVariables.PressureVector);

    //Distribute coupling block vector 1 into elemental vector
    GeoElementUtilities::AssembleUBlockVector<TDim, TNumNodes>(rRightHandSideVector, rVariables.UVector);

    if (!rVariables.IgnoreUndrained)
    {
        noalias(rVariables.PVector) = signFactor * prod(trans(rVariables.UPMatrix),
                                                      rVariables.VelocityVector);

        //Distribute coupling block vector 2 into elemental vector
        GeoElementUtilities::AssemblePBlockVector< TDim, TNumNodes >(rRightHandSideVector, rVariables.PVector);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddCouplingTerms()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddCompressibilityFlow( VectorType& rRightHandSideVector,
                                        ElementVariables& rVariables )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddCompressibilityFlow()") << std::endl;

    noalias(rVariables.PMatrix) = - signFactor * rVariables.BiotModulusInverse
                                               * outer_prod(rVariables.Np,rVariables.Np)
                                               * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = -1.0*prod(rVariables.PMatrix,rVariables.DtPressureVector);

    //Distribute compressibility block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<TDim, TNumNodes>(rRightHandSideVector,rVariables.PVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddCompressibilityFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddPermeabilityFlow( VectorType& rRightHandSideVector,
                                     ElementVariables& rVariables )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddPermeabilityFlow()") << std::endl;

    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT,rVariables.PermeabilityMatrix);

    noalias(rVariables.PMatrix) = - signFactor * rVariables.DynamicViscosityInverse
                                               * prod(rVariables.PDimMatrix,trans(rVariables.GradNpT))
                                               * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = -1.0*prod(rVariables.PMatrix,rVariables.PressureVector);

    //Distribute permeability block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<TDim, TNumNodes>(rRightHandSideVector,rVariables.PVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddPermeabilityFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                 ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddFluidBodyFlow()") << std::endl;

    noalias(rVariables.PDimMatrix) =  prod(rVariables.GradNpT,rVariables.PermeabilityMatrix)
                                    * rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) =  rVariables.DynamicViscosityInverse
                                 * rVariables.FluidDensity
                                 * prod(rVariables.PDimMatrix,rVariables.BodyAcceleration);

    //Distribute fluid body flow block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<TDim, TNumNodes>(rRightHandSideVector,rVariables.PVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddFluidBodyFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwSmallStrainElement<2,3>;
template class UPwSmallStrainElement<2,4>;
template class UPwSmallStrainElement<3,4>;
template class UPwSmallStrainElement<3,8>;

template class UPwSmallStrainElement<2,6>;
template class UPwSmallStrainElement<2,8>;
template class UPwSmallStrainElement<2,9>;
template class UPwSmallStrainElement<3,10>;
template class UPwSmallStrainElement<3,20>;
template class UPwSmallStrainElement<3,27>;

} // Namespace Kratos
