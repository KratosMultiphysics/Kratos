//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Danilo Cavalcanti
//                   Lorena Casallas
//                   Ignasi de Pouplana
//


// Application includes
#include "custom_elements/two-phase_flow/U_Pl_Pg_small_strain_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPlPgSmallStrainElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UPlPgSmallStrainElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPlPgSmallStrainElement<TDim,TNumNodes>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPlPgSmallStrainElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
int UPlPgSmallStrainElement<TDim,TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();

    // Base class checks for positive area and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if(Geom.DomainSize() < 1.0e-15)
        KRATOS_THROW_ERROR( std::logic_error, "DomainSize < 1.0e-15 for the element ", this->Id() )

    // Verify generic variables
    ierr = UPlPgElement<TDim,TNumNodes>::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Verify specific properties
    if ( PERMEABILITY_XX.Key() == 0 || Prop.Has( PERMEABILITY_XX ) == false || Prop[PERMEABILITY_XX] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_XX has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( PERMEABILITY_YY.Key() == 0 || Prop.Has( PERMEABILITY_YY ) == false || Prop[PERMEABILITY_YY] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_YY has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( PERMEABILITY_XY.Key() == 0 || Prop.Has( PERMEABILITY_XY ) == false || Prop[PERMEABILITY_XY] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_XY has Key zero, is not defined or has an invalid value at element", this->Id() )
    if constexpr (TDim > 2)
    {
        if ( PERMEABILITY_ZZ.Key() == 0 || Prop.Has( PERMEABILITY_ZZ ) == false || Prop[PERMEABILITY_ZZ] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_ZZ has Key zero, is not defined or has an invalid value at element", this->Id() )
        if ( PERMEABILITY_YZ.Key() == 0 || Prop.Has( PERMEABILITY_YZ ) == false || Prop[PERMEABILITY_YZ] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_YZ has Key zero, is not defined or has an invalid value at element", this->Id() )
        if ( PERMEABILITY_ZX.Key() == 0 || Prop.Has( PERMEABILITY_ZX ) == false || Prop[PERMEABILITY_ZX] < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_ZX has Key zero, is not defined or has an invalid value at element", this->Id() )
    }
    if ( BIOT_COEFFICIENT.Key() == 0 || Prop.Has( BIOT_COEFFICIENT ) == false || Prop[BIOT_COEFFICIENT] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"BIOT_COEFFICIENT has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( GAS_HENRY_SOLUBILITY_COEFF.Key() == 0 || Prop.Has( GAS_HENRY_SOLUBILITY_COEFF ) == false || Prop[GAS_HENRY_SOLUBILITY_COEFF] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"GAS_HENRY_SOLUBILITY_COEFF has Key zero, is not defined or has an invalid value at element", this->Id() )

    // Verify the constitutive law
    if ( CONSTITUTIVE_LAW.Key() == 0 || Prop.Has( CONSTITUTIVE_LAW ) == false )
        KRATOS_THROW_ERROR( std::invalid_argument, "CONSTITUTIVE_LAW has Key zero or is not defined at element ", this->Id() )
    if ( Prop[CONSTITUTIVE_LAW] != NULL )
    {
        // Verify compatibility of the element with the constitutive law
        ConstitutiveLaw::Features LawFeatures;
        Prop[CONSTITUTIVE_LAW]->GetLawFeatures(LawFeatures);
        bool correct_strain_measure = false;
        for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
        {
            if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal)
                correct_strain_measure = true;
        }
        if( correct_strain_measure == false )
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the element type", " StrainMeasure_Infinitesimal " );

        // Check constitutive law
        ierr = Prop[CONSTITUTIVE_LAW]->Check( Prop, Geom, rCurrentProcessInfo );
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "A constitutive law needs to be specified for the element ", this->Id() )

    // Verify the saturation law
    if ( SATURATION_LAW_NAME.Key() == 0 || Prop.Has( SATURATION_LAW_NAME ) == false ) {
        KRATOS_THROW_ERROR( std::invalid_argument, "SATURATION_LAW_NAME has Key zero or is not defined at element ", this->Id() )
    } else {
        // Check saturation law
        ierr = mSaturationLawVector[0]->Check( Prop, Geom, rCurrentProcessInfo );
    }

    return ierr;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

    const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    Matrix B(strain_size,TNumNodes*TDim);
    noalias(B) = ZeroMatrix(strain_size,TNumNodes*TDim);
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);

    //Create constitutive law parameters:
    Vector StrainVector(strain_size);
    Vector StressVector(strain_size);
    Matrix ConstitutiveMatrix(strain_size,strain_size);
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
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        noalias(Np) = row(NContainer,GPoint);

        this->CalculateKinematics(GradNpT,B,StrainVector,DN_DXContainer,DisplacementVector,GPoint);

        // Compute Stress
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    this->InitializeNonLinearIteration(rCurrentProcessInfo);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //Defining necessary variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

    const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    Matrix B(strain_size,TNumNodes*TDim);
    noalias(B) = ZeroMatrix(strain_size,TNumNodes*TDim);
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);

    //Create constitutive law parameters:
    Vector StrainVector(strain_size);
    Vector StressVector(strain_size);
    Matrix ConstitutiveMatrix(strain_size,strain_size);
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
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    if(rCurrentProcessInfo[NODAL_SMOOTHING] == true)
    {
        Matrix StressContainer(NumGPoints,strain_size);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            this->CalculateKinematics(GradNpT,B,StrainVector,DN_DXContainer,DisplacementVector,GPoint);

            noalias(Np) = row(NContainer,GPoint);

            //compute constitutive tensor and/or stresses
            mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
            this->SaveGPStress(StressContainer,StressVector,strain_size,GPoint);
        }
        this->ExtrapolateGPValues(StressContainer,strain_size);
    }
    else
    {
        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            this->CalculateKinematics(GradNpT,B,StrainVector,DN_DXContainer,DisplacementVector,GPoint);

            noalias(Np) = row(NContainer,GPoint);

            //compute constitutive tensor and/or stresses
            mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rOutput,
                                                                            const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != NumGPoints )
        rOutput.resize( NumGPoints, false );

    if (rVariable == VON_MISES_STRESS) {
        //Defining necessary variables
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);
        const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

        Matrix B(strain_size,TNumNodes*TDim);
        noalias(B) = ZeroMatrix(strain_size,TNumNodes*TDim);
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);

        //Create constitutive law parameters:
        Vector StrainVector(strain_size);
        Vector StressVector(strain_size);
        Matrix ConstitutiveMatrix(strain_size,strain_size);
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
            this->CalculateKinematics(GradNpT,B,StrainVector,DN_DXContainer,DisplacementVector,GPoint);

            noalias(Np) = row(NContainer,GPoint);

            //compute constitutive tensor and/or stresses
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            rOutput[GPoint] = ElementUtilities::CalculateVonMises(StressVector);
        }

    } else if (rVariable == LIQUID_SATURATION_DEGREE) {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);

        //Saturation Law paremeters
        SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            //Compute GradNpT (passed to the SaturationLaw)
            this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

            //Compute Np
            noalias(Variables.Np) = row(NContainer,GPoint);

            //Compute Saturation Law variables: Sl, dSldPc
            mSaturationLawVector[GPoint]->CalculateSaturation(SaturationParameters);

            rOutput[GPoint] = Variables.Sl;
        }

    } else if (rVariable == GAS_SATURATION_DEGREE) {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);

        //Saturation Law paremeters
        SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            //Compute GradNpT (passed to the SaturationLaw)
            this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

            //Compute Np
            noalias(Variables.Np) = row(NContainer,GPoint);

            //Compute Saturation Law variables: Sl, dSldPc
            mSaturationLawVector[GPoint]->CalculateSaturation(SaturationParameters);

            rOutput[GPoint] = 1.0-Variables.Sl;
        }
    
    } else if (rVariable == LIQUID_RELATIVE_PERMEABILITY) {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);

        //Saturation Law paremeters
        SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            //Compute GradNpT (passed to the SaturationLaw)
            this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

            //Compute Np
            noalias(Variables.Np) = row(NContainer,GPoint);

            //Compute Saturation Law variables: Sl, dSldPc, krl and krg
            mSaturationLawVector[GPoint]->CalculateMaterialResponse(SaturationParameters);

            rOutput[GPoint] = Variables.krl;
        }

    } else if (rVariable == GAS_RELATIVE_PERMEABILITY) {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);

        //Saturation Law paremeters
        SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            //Compute GradNpT (passed to the SaturationLaw)
            this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

            //Compute Np
            noalias(Variables.Np) = row(NContainer,GPoint);

            //Compute Saturation Law variables: Sl, dSldPc, krl and krg
            mSaturationLawVector[GPoint]->CalculateMaterialResponse(SaturationParameters);

            rOutput[GPoint] = Variables.krg;
        }

    } else {
        UPlPgElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(rVariable,rOutput,rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<array_1d<double,3>>& rVariable,std::vector<array_1d<double,3>>& rOutput,
                                                                            const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != NumGPoints )
        rOutput.resize( NumGPoints );

    if (rVariable == LIQUID_FLUX_VECTOR) {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);

        //Saturation Law paremeters
        SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

        // Necessary variables
        array_1d<double,TDim> GradLiquidPressureTerm;
        array_1d<double,TDim> LiquidFlux;

        //Loop over integration points
        for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
        {
            //Compute GradNpT, B and StrainVector
            this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

            //Compute Np, Nu and BodyAcceleration
            noalias(Variables.Np) = row(NContainer,GPoint);
            PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

            //Compute Saturation Law variables: Sl, dSldPc, krl and krg
            mSaturationLawVector[GPoint]->CalculateMaterialResponse(SaturationParameters);

            noalias(GradLiquidPressureTerm) = prod(trans(Variables.GradNpT),Variables.LiquidPressureVector);
            noalias(GradLiquidPressureTerm) += -Variables.LiquidDensity*Variables.BodyAcceleration;

            noalias(LiquidFlux) = -Variables.krl/(Variables.Porosity*Variables.Sl*Variables.LiquidDynamicViscosity)*
                                    prod(mIntrinsicPermeability,GradLiquidPressureTerm);

            PoroElementUtilities::FillArray1dOutput(rOutput[GPoint],LiquidFlux);
        }
        
    } else if (rVariable == GAS_FLUX_VECTOR) {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);

        //Saturation Law paremeters
        SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

        // Necessary variables
        array_1d<double,TDim> GradGasPressureTerm;
        array_1d<double,TDim> GasFlux;

        //Loop over integration points
        for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
        {
            //Compute GradNpT, B and StrainVector
            this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

            //Compute Np, Nu and BodyAcceleration
            noalias(Variables.Np) = row(NContainer,GPoint);
            PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

            //Compute Saturation Law variables: Sl, dSldPc, krl and krg
            mSaturationLawVector[GPoint]->CalculateMaterialResponse(SaturationParameters);

            // Update the density based on the new saturation degree
            this->UpdateDensity(Variables);

            noalias(GradGasPressureTerm) = prod(trans(Variables.GradNpT),Variables.GasPressureVector);
            noalias(GradGasPressureTerm) += -Variables.GasDensity*Variables.BodyAcceleration;

            noalias(GasFlux) = -Variables.krg/(Variables.Porosity*(1.0-Variables.Sl)*Variables.GasDynamicViscosity)*
                                    prod(mIntrinsicPermeability,GradGasPressureTerm);

            PoroElementUtilities::FillArray1dOutput(rOutput[GPoint],GasFlux);
        }

    } else if (rVariable == LIQUID_PRESSURE_GRADIENT) {

        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Defining necessary variables
        array_1d<double,TNumNodes> LiquidPressureVector;
        for(unsigned int i=0; i<TNumNodes; i++)
            LiquidPressureVector[i] = Geom[i].FastGetSolutionStepValue(LIQUID_PRESSURE);
        BoundedMatrix<double,TNumNodes, TDim> GradNpT;
        array_1d<double,TDim> GradLiquidPressure;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            noalias(GradLiquidPressure) = prod(trans(GradNpT),LiquidPressureVector);

            PoroElementUtilities::FillArray1dOutput(rOutput[GPoint],GradLiquidPressure);
        }

    } else if (rVariable == GAS_PRESSURE_GRADIENT) {
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Defining necessary variables
        array_1d<double,TNumNodes> GasPressureVector;
        for(unsigned int i=0; i<TNumNodes; i++)
            GasPressureVector[i] = Geom[i].FastGetSolutionStepValue(GAS_PRESSURE);
        BoundedMatrix<double,TNumNodes, TDim> GradNpT;
        array_1d<double,TDim> GradGasPressure;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            noalias(GradGasPressure) = prod(trans(GradNpT),GasPressureVector);

            PoroElementUtilities::FillArray1dOutput(rOutput[GPoint],GradGasPressure);
        }

    } else {
        UPlPgElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(rVariable,rOutput,rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput,
                                                                                const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != NumGPoints )
        rOutput.resize( NumGPoints );

    if(rVariable == EFFECTIVE_STRESS_TENSOR)
    {
        //Defining necessary variables
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);
        const unsigned int cl_dimension = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();
        const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

        Matrix B(strain_size,TNumNodes*TDim);
        noalias(B) = ZeroMatrix(strain_size,TNumNodes*TDim);
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);

        //Create constitutive law parameters:
        Vector StrainVector(strain_size);
        Vector StressVector(strain_size);
        Matrix ConstitutiveMatrix(strain_size,strain_size);
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
            this->CalculateKinematics(GradNpT,B,StrainVector,DN_DXContainer,DisplacementVector,GPoint);

            noalias(Np) = row(NContainer,GPoint);

            //compute constitutive tensor and/or stresses
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            rOutput[GPoint].resize(cl_dimension,cl_dimension,false );
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector);
        }
    }
    else if(rVariable == TOTAL_STRESS_TENSOR)
    {
        //Defining necessary variables
        const PropertiesType& Prop = this->GetProperties();
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);
        const unsigned int cl_dimension = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();

        //Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        //Saturation Law paremeters
        SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

        // Total pressure (of the solid phase)
        double TotalPressure = 0.0;
        double LiquidPressure;
        double GasPressure;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            //Compute GradNpT, B and StrainVector
            this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

            //Compute Np, Nu and BodyAcceleration
            noalias(Variables.Np) = row(NContainer,GPoint);
            PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
            PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

            //Compute Saturation Law variables: Sl, dSldPc, krl and krg
            mSaturationLawVector[GPoint]->CalculateMaterialResponse(SaturationParameters);

            // Update the density based on the new saturation degree
            this->UpdateDensity(Variables);

            //Compute constitutive tensor and stresses
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            LiquidPressure = 0.0;
            GasPressure = 0.0;
            for(unsigned int i = 0; i < TNumNodes; i++)
            {
                LiquidPressure += NContainer(GPoint,i)*Variables.LiquidPressureVector[i];
                GasPressure += NContainer(GPoint,i)*Variables.GasPressureVector[i];
            }
            TotalPressure = Variables.Sl*LiquidPressure + (1.0-Variables.Sl)*GasPressure;

            noalias(Variables.StressVector) += -Variables.BiotCoefficient*TotalPressure*Variables.VoigtVector;

            rOutput[GPoint].resize(cl_dimension,cl_dimension,false );
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(Variables.StressVector);
        }
    }
    else if(rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
    {
        //Defining necessary variables
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);
        const unsigned int cl_dimension = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();
        const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

        Matrix B(strain_size,TNumNodes*TDim);
        noalias(B) = ZeroMatrix(strain_size,TNumNodes*TDim);
        Matrix GradNpT(TNumNodes,TDim);
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
        Vector StrainVector(strain_size);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            this->CalculateKinematics(GradNpT,B,StrainVector,DN_DXContainer,DisplacementVector,GPoint);

            rOutput[GPoint].resize(cl_dimension,cl_dimension,false );
            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector);
        }
    }
    else if(rVariable == PERMEABILITY_MATRIX)
    {
        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            rOutput[GPoint].resize(TDim,TDim,false);
            noalias(rOutput[GPoint]) = mIntrinsicPermeability;
        }
    } else if(rVariable == LIQUID_PERMEABILITY_MATRIX)
    {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);

        //Saturation Law paremeters
        SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            //Compute GradNpT (passed to the SaturationLaw)
            this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

            //Compute Np
            noalias(Variables.Np) = row(NContainer,GPoint);

            //Compute Saturation Law variables: Sl, dSldPc, krl and krg
            mSaturationLawVector[GPoint]->CalculateMaterialResponse(SaturationParameters);

            rOutput[GPoint].resize(TDim,TDim,false);
            noalias(rOutput[GPoint]) = Variables.krl*mIntrinsicPermeability;
        }
    } else if(rVariable == GAS_PERMEABILITY_MATRIX)
    {
        //Previous definitions
        const PropertiesType& Prop = this->GetProperties();

        //Containers of variables at all integration points
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Constitutive Law parameters
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);

        //Saturation Law paremeters
        SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            //Compute GradNpT (passed to the SaturationLaw)
            this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

            //Compute Np
            noalias(Variables.Np) = row(NContainer,GPoint);

            //Compute Saturation Law variables: Sl, dSldPc, krl and krg
            mSaturationLawVector[GPoint]->CalculateMaterialResponse(SaturationParameters);

            rOutput[GPoint].resize(TDim,TDim,false);
            noalias(rOutput[GPoint]) = Variables.krg*mIntrinsicPermeability;
        }
    } else {
        UPlPgElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(rVariable,rOutput,rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 2);

    //Resizing mass matrix
    if ( rStiffnessMatrix.size1() != element_size )
        rStiffnessMatrix.resize( element_size, element_size, false );
    noalias( rStiffnessMatrix ) = ZeroMatrix( element_size, element_size );

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
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Saturation Law paremeters
    SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute Np, GradNpT, B and StrainVector
        noalias(Variables.Np) = row(NContainer,GPoint);
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute constitutive tensor
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Compute stiffness matrix
        this->CalculateAndAddStiffnessMatrix(rStiffnessMatrix, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 2);

    //Resizing mass matrix
    if ( rMassMatrix.size1() != element_size )
        rMassMatrix.resize( element_size, element_size, false );
    noalias( rMassMatrix ) = ZeroMatrix( element_size, element_size );

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
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);

    //Saturation Law paremeters
    SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
    {
        //Compute GradNpT (passed to the SaturationLaw)
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute Np and Nu
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);

        //Compute Saturation Law variables: Sl, dSldPc
        mSaturationLawVector[GPoint]->CalculateSaturation(SaturationParameters);

        // Update the density based on the new saturation degree
        this->UpdateDensity(Variables);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Adding contribution to Mass matrix
        this->CalculateAndAddMassMatrix(rMassMatrix, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateLumpedMassMatrix( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 2);

    //Resizing mass matrix
    if ( rLeftHandSideMatrix.size1() != element_size )
        rLeftHandSideMatrix.resize( element_size, element_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( element_size, element_size );

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();

    // Initialize elemenatal total mass
    double total_mass = 0.0;

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Vector detJContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);

    //Saturation Law paremeters
    SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
    {
        //Compute GradNpT (passed to the SaturationLaw)
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute Np
        noalias(Variables.Np) = row(NContainer,GPoint);

        //Compute Saturation Law variables: Sl, dSldPc
        mSaturationLawVector[GPoint]->CalculateSaturation(SaturationParameters);

        // Update the density based on the new saturation degree
        this->UpdateDensity(Variables);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        total_mass += Variables.Density * Variables.IntegrationCoefficient;
    }

    // LUMPED MASS MATRIX
    Vector lumping_factors;
    lumping_factors = Geom.LumpingFactors( lumping_factors );

    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        const double temp = lumping_factors[i] * total_mass;
        for ( IndexType j = 0; j < TDim; ++j ) {
            IndexType index = i * (TDim + 1) + j;
            rLeftHandSideMatrix(index,index) = temp;
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddMassMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    noalias(rVariables.UMatrix) = rVariables.Density*prod(trans(rVariables.Nu),rVariables.Nu)*rVariables.IntegrationCoefficient;

    //Distribute stiffness block matrix into the elemental matrix
    PoroElementUtilities::AssembleUBlockTwoPhaseFlowMatrix(rLeftHandSideMatrix,rVariables.UMatrix);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Saturation Law paremeters
    SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

        //Compute Saturation Law variables: Sl, dSldPc, krl and krg
        mSaturationLawVector[GPoint]->CalculateMaterialResponse(SaturationParameters);

        // Update the density based on the new saturation degree
        this->UpdateDensity(Variables);

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
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Saturation Law paremeters
    SaturationLaw::Parameters SaturationParameters(Geom,Prop,rCurrentProcessInfo);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,SaturationParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

        //Compute Saturation Law variables: Sl, dSldPc, krl and krg
        mSaturationLawVector[GPoint]->CalculateMaterialResponse(SaturationParameters);
        
        // Update the density based on the new saturation degree
        this->UpdateDensity(Variables);

        //Compute stresses
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddHydromechanicalCouplingMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B),rVariables.ConstitutiveMatrix);
    noalias(rVariables.UMatrix) = prod(rVariables.UVoigtMatrix,rVariables.B)*rVariables.IntegrationCoefficient;

    //Distribute stiffness block matrix into the elemental matrix
    PoroElementUtilities::AssembleUBlockTwoPhaseFlowMatrix(rLeftHandSideMatrix,rVariables.UMatrix);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddHydromechanicalCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    // Compute intermediate auxiliary vector and matrix
    noalias(rVariables.UVector) = prod(trans(rVariables.B),rVariables.VoigtVector);
    noalias(rVariables.UPMatrixAux) = -outer_prod(rVariables.UVector,rVariables.Np)*rVariables.IntegrationCoefficient;

    // Get hydromechanical coupling coefficients
    double Clu, Cgu;
    this->GetHydromechanicalCouplingCoefficients(Clu, Cgu, rVariables);

    // Compute the liquid element coupling matrices
    noalias(rVariables.UPMatrix) = Clu*rVariables.UPMatrixAux;
    PoroElementUtilities::AssembleUPlBlockMatrix(rLeftHandSideMatrix,rVariables.UPMatrix);
    noalias(rVariables.PUMatrix) = -rVariables.VelocityCoefficient*trans(rVariables.UPMatrix);
    PoroElementUtilities::AssemblePlUBlockMatrix(rLeftHandSideMatrix,rVariables.PUMatrix);
    
    // Compute the gas element coupling matrices
    noalias(rVariables.UPMatrix) = Cgu*rVariables.UPMatrixAux;
    PoroElementUtilities::AssembleUPgBlockMatrix(rLeftHandSideMatrix,rVariables.UPMatrix);
    noalias(rVariables.PUMatrix) = -rVariables.VelocityCoefficient*trans(rVariables.UPMatrix);
    PoroElementUtilities::AssemblePgUBlockMatrix(rLeftHandSideMatrix,rVariables.PUMatrix);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    // Compute auxiliary matrix: Np^T * Np * w * detJ
    noalias(rVariables.PMatrixAux) = outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;

    // Get compressibility coefficients
    double Cll, Clg, Cgl, Cgg;
    this->GetCompressibilityCoefficients(Cll, Clg, Cgl, Cgg, rVariables);

    // Compute the element compressibility sub-matrices
    noalias(rVariables.PMatrix) = rVariables.DtLiquidPressureCoefficient * Cll * rVariables.PMatrixAux;
    PoroElementUtilities::AssemblePlPlBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);
    noalias(rVariables.PMatrix) = rVariables.DtGasPressureCoefficient * Clg * rVariables.PMatrixAux;
    PoroElementUtilities::AssemblePlPgBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);
    noalias(rVariables.PMatrix) = rVariables.DtLiquidPressureCoefficient * Cgl * rVariables.PMatrixAux;
    PoroElementUtilities::AssemblePgPlBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);
    noalias(rVariables.PMatrix) = rVariables.DtGasPressureCoefficient * Cgg * rVariables.PMatrixAux;
    PoroElementUtilities::AssemblePgPgBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{

    // Compute the auxiliary matrix of the product between: gradNp^T * K * gradNp
    // where gradNp is the gradient of the shape function vector associated with the pressure fields,
    // and K is the permeability matrix.
    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT,mIntrinsicPermeability);
    noalias(rVariables.PMatrixAux) = prod(rVariables.PDimMatrix,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;

    // Compute the fluid-flow sub-matrices
    // Liquid (Hl)
    noalias(rVariables.PMatrix) = rVariables.krl/rVariables.LiquidDynamicViscosity * rVariables.PMatrixAux;
    PoroElementUtilities::AssemblePlPlBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);
    // Gas (Hg)
    noalias(rVariables.PMatrix) = rVariables.krg/rVariables.GasDynamicViscosity * rVariables.PMatrixAux;
    PoroElementUtilities::AssemblePgPgBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);

    //TODO(U_Pl_pg)
    // Add the permeability matrices associated with the diffusion of the gas into the liquid phase (Hgl and one term of Hg)
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddHydromechanicalCouplingTerms(rRightHandSideVector, rVariables);

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
}

//----------------------------------------------------------------------------------------

// Integral of (B^T * sigma)
template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rVariables.UVector) = -1.0*prod(trans(rVariables.B),rVariables.StressVector)*rVariables.IntegrationCoefficient;

    //Distribute stiffness block vector into elemental vector
    PoroElementUtilities::AssembleUBlockTwoPhaseFlowVector(rRightHandSideVector,rVariables.UVector);
}

//----------------------------------------------------------------------------------------

// Integral of (N^T * rho * b)
template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rVariables.UVector) = rVariables.Density*prod(trans(rVariables.Nu),rVariables.BodyAcceleration)*rVariables.IntegrationCoefficient;

    //Distribute body force block vector into elemental vector
    PoroElementUtilities::AssembleUBlockTwoPhaseFlowVector(rRightHandSideVector,rVariables.UVector);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddHydromechanicalCouplingTerms(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{

    // Compute intermediate auxiliary vector and matrix
    noalias(rVariables.UVector) = prod(trans(rVariables.B),rVariables.VoigtVector);
    noalias(rVariables.UPMatrixAux) = outer_prod(rVariables.UVector,rVariables.Np)*rVariables.IntegrationCoefficient;

    // Get compressibility coefficients
    double Clu, Cgu;
    this->GetHydromechanicalCouplingCoefficients(Clu, Cgu, rVariables);

    // Compute the element coupling liquid hydromechanical matrices
    noalias(rVariables.UPMatrix) = Clu*rVariables.UPMatrixAux;
    //Add the contribution of the coupling between the displacement and the liquid pressure in the residual of the momentum equation: Qw * pl
    noalias(rVariables.UVector) = prod(rVariables.UPMatrix,rVariables.LiquidPressureVector);
    PoroElementUtilities::AssembleUBlockTwoPhaseFlowVector(rRightHandSideVector,rVariables.UVector);
    //Add the contribution of the coupling between the displacement and the liquid pressure in the residual of the mass balance of liquid equation: -Qw^T * v
    noalias(rVariables.PVector) = -1.0*prod(trans(rVariables.UPMatrix),rVariables.VelocityVector);
    PoroElementUtilities::AssemblePlBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);

    // Compute the element coupling gas hydromechanical matrices
    noalias(rVariables.UPMatrix) = Cgu*rVariables.UPMatrixAux;
    //Add the contribution of the coupling between the displacement and the gas pressure in the residual of the momentum equation: Qg * pg
    noalias(rVariables.UVector) = prod(rVariables.UPMatrix,rVariables.GasPressureVector);
    PoroElementUtilities::AssembleUBlockTwoPhaseFlowVector(rRightHandSideVector,rVariables.UVector);
    //Add the contribution of the coupling between the displacement and the gas pressure in the residual of the mass balance of gas equation: -Qg^T * v
    noalias(rVariables.PVector) = -1.0*prod(trans(rVariables.UPMatrix),rVariables.VelocityVector);
    PoroElementUtilities::AssemblePgBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    // Compute auxiliary matrix: Np^T * Np * w * detJ
    noalias(rVariables.PMatrixAux) = outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;

    // Get compressibility coefficients
    double Cll, Clg, Cgl, Cgg;
    this->GetCompressibilityCoefficients(Cll, Clg, Cgl, Cgg, rVariables);

    //Add the contribution of the compressibility terms in the residual of the mass balance of liquid equation: -Sll * dpldt - Slg * dpgdt
    noalias(rVariables.PVector) = -Cll*prod(rVariables.PMatrixAux,rVariables.DtLiquidPressureVector) - Clg*prod(rVariables.PMatrixAux,rVariables.DtGasPressureVector);
    PoroElementUtilities::AssemblePlBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);

    //Add the contribution of the compressibility terms in the residual of the mass balance of gas equation: -Sgl * dpldt - Sgg * dpgdt
    noalias(rVariables.PVector) = -Cgl*prod(rVariables.PMatrixAux,rVariables.DtLiquidPressureVector) - Cgg*prod(rVariables.PMatrixAux,rVariables.DtGasPressureVector);
    PoroElementUtilities::AssemblePgBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    // Compute the auxiliary matrix of the product between: gradNp^T * K * gradNp
    // where gradNp is the gradient of the shape function vector associated with the pressure fields,
    // and K is the permeability matrix.
    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT,mIntrinsicPermeability);
    noalias(rVariables.PMatrixAux) = prod(rVariables.PDimMatrix,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;

    noalias(rVariables.PMatrix) = rVariables.krl/rVariables.LiquidDynamicViscosity * rVariables.PMatrixAux;
    //Add the contribution of the permeability term in the residual of the mass balance of liquid equation: -Hll * pl
    noalias(rVariables.PVector) = -1.0*prod(rVariables.PMatrix,rVariables.LiquidPressureVector);
    PoroElementUtilities::AssemblePlBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);

    noalias(rVariables.PMatrix) = rVariables.krg/rVariables.GasDynamicViscosity * rVariables.PMatrixAux;
    //Add the contribution of the permeability term in the residual of the mass balance of liquid equation: -Hgg * pg
    noalias(rVariables.PVector) = -1.0*prod(rVariables.PMatrix,rVariables.GasPressureVector);
    PoroElementUtilities::AssemblePgBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);

    //TODO(U_Pl_pg)
    // Add the permeability flow associated with the diffusion of the gas into the liquid phase (Hgl*Pl and one term of Hg*Pg)
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    // Compute the auxiliary matrix of the product between: gradNp^T * K
    // where gradNp is the gradient of the shape function vector associated with the pressure fields,
    // and K is the permeability matrix.
    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT,mIntrinsicPermeability);
    noalias(rVariables.PVectorAux) = prod(rVariables.PDimMatrix,rVariables.BodyAcceleration)*rVariables.IntegrationCoefficient;

    //NOTE. Maybe it would be cleaner if these two terms were calculated in different methods

    noalias(rVariables.PVector) = (rVariables.krl * rVariables.LiquidDensity / rVariables.LiquidDynamicViscosity) * rVariables.PVectorAux;
    //Add the contribution of the fluid body flow block vector into the residual of the mass balance equation of liquid
    PoroElementUtilities::AssemblePlBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);

    noalias(rVariables.PVector) = (rVariables.krg * rVariables.GasDensity   / rVariables.GasDynamicViscosity  ) * rVariables.PVectorAux;
    //Add the contribution of the fluid body flow block vector into the residual of the mass balance equation of gas
    PoroElementUtilities::AssemblePgBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);

    //TODO(U_Pl_pg)
    // Add the term associated with the diffusion of the gas into the liquid phase
    // We are also missing the conditions associated with the diffusion of gas into liquid
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::SaveGPStress(Matrix& rStressContainer, const Vector& StressVector, const unsigned int& VoigtSize, const unsigned int& GPoint)
{
    for(unsigned int i = 0; i < VoigtSize; i++)
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
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgSmallStrainElement<2,3>::ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    // Triangle_2d_3 with GI_GAUSS_2

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area();
    array_1d<array_1d<double,2>,3> NodalGradPressure; //List with 2D GradP at each node
    array_1d<Vector,3> NodalStressVector; //List with stresses at each node
    array_1d<Matrix,3> NodalStressTensor;
    const unsigned int cl_dimension = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();

    for(unsigned int iNode = 0; iNode < 3; iNode++)
    {
        NodalStressVector[iNode].resize(VoigtSize);
        NodalStressTensor[iNode].resize(cl_dimension,cl_dimension);
    }

    BoundedMatrix<double,3,3> ExtrapolationMatrix;
    PoroElementUtilities::Calculate2DExtrapolationMatrix(ExtrapolationMatrix);

    Matrix AuxNodalStress(3,VoigtSize);
    noalias(AuxNodalStress) = prod(ExtrapolationMatrix,StressContainer);

    /* INFO:
        *
        *                  |S0-0 S1-0 S2-0|
        * AuxNodalStress = |S0-1 S1-1 S2-1|
        *                  |S0-2 S1-2 S2-2|
        *
        * S1-0 = S[1] at node 0
    */

    for(unsigned int i = 0; i < 3; i++) //TNumNodes
    {
        noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        Matrix& rNodalStress = rGeom[i].FastGetSolutionStepValue(NODAL_EFFECTIVE_STRESS_TENSOR);
        for(unsigned int j = 0; j < cl_dimension; j++){
            for(unsigned int k = 0; k < cl_dimension; k++){
                rNodalStress(j,k) += NodalStressTensor[i](j,k);
            }
        }
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgSmallStrainElement<2,4>::ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    // Quadrilateral_2d_4 with GI_GAUSS_2

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area();
    array_1d<Vector,4> NodalStressVector; //List with stresses at each node
    array_1d<Matrix,4> NodalStressTensor;
    const unsigned int cl_dimension = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();

    for(unsigned int iNode = 0; iNode < 4; iNode ++)
    {
        NodalStressVector[iNode].resize(VoigtSize);
        NodalStressTensor[iNode].resize(cl_dimension,cl_dimension);
    }

    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    PoroElementUtilities::Calculate2DExtrapolationMatrix(ExtrapolationMatrix);

    Matrix AuxNodalStress(4,VoigtSize);
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

    for(unsigned int i = 0; i < 4; i++) //TNumNodes
    {
        noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        Matrix& rNodalStress = rGeom[i].FastGetSolutionStepValue(NODAL_EFFECTIVE_STRESS_TENSOR);
        for(unsigned int j = 0; j < cl_dimension; j++){
            for(unsigned int k = 0; k < cl_dimension; k++){
                rNodalStress(j,k) += NodalStressTensor[i](j,k);
            }
        }
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgSmallStrainElement<3,4>::ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    // Tetrahedra_3d_4 with GI_GAUSS_2

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area(); // In 3D this is Volume
    array_1d<Vector,4> NodalStressVector; //List with stresses at each node
    array_1d<Matrix,4> NodalStressTensor;

    for(unsigned int iNode = 0; iNode < 4; iNode ++)
    {
        NodalStressVector[iNode].resize(VoigtSize);
        NodalStressTensor[iNode].resize(3,3);
    }

    BoundedMatrix<double,4,4> ExtrapolationMatrix;
    PoroElementUtilities::Calculate3DExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double,4,6> AuxNodalStress;
    noalias(AuxNodalStress) = prod(ExtrapolationMatrix,StressContainer);

    for(unsigned int i = 0; i < 4; i++) //TNumNodes
    {
        noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        noalias(rGeom[i].FastGetSolutionStepValue(NODAL_EFFECTIVE_STRESS_TENSOR)) += NodalStressTensor[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgSmallStrainElement<3,8>::ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    // Hexahedra_3d_8 with GI_GAUSS_2

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area(); // In 3D this is Volume
    array_1d<Vector,8> NodalStressVector; //List with stresses at each node
    array_1d<Matrix,8> NodalStressTensor;

    for(unsigned int iNode = 0; iNode < 8; iNode ++)
    {
        NodalStressVector[iNode].resize(VoigtSize);
        NodalStressTensor[iNode].resize(3,3);
    }

    BoundedMatrix<double,8,8> ExtrapolationMatrix;
    PoroElementUtilities::Calculate3DExtrapolationMatrix(ExtrapolationMatrix);

    BoundedMatrix<double,8,6> AuxNodalStress;
    noalias(AuxNodalStress) = prod(ExtrapolationMatrix,StressContainer);

    for(unsigned int i = 0; i < 8; i++) //TNumNodes
    {
        noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        noalias(rGeom[i].FastGetSolutionStepValue(NODAL_EFFECTIVE_STRESS_TENSOR)) += NodalStressTensor[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::InitializeElementVariables(ElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                                                            SaturationLaw::Parameters& rSaturationParameters,
                                                                            const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //Properties variables
    rVariables.Porosity                   = Prop[POROSITY];
    rVariables.LiquidDynamicViscosity     = Prop[DYNAMIC_VISCOSITY_LIQUID];
    rVariables.GasDynamicViscosity        = Prop[DYNAMIC_VISCOSITY_GAS];
    rVariables.SolidDensity               = Prop[DENSITY_SOLID];
    rVariables.LiquidDensity              = Prop[DENSITY_LIQUID];
    rVariables.GasDensity                 = Prop[DENSITY_GAS];
    rVariables.GasHenrySolubilityCoeff    = Prop[GAS_HENRY_SOLUBILITY_COEFF];
    rVariables.GasMolarWeight             = Prop[GAS_MOLAR_WEIGHT];
    rVariables.Temperature                = Prop[TEMPERATURE];
    rVariables.TortuosityCoeff            = Prop[TORTUOSITY_COEFF];
    rVariables.GasDiffusionCoeff          = Prop[GAS_DIFFUSION_COEFF];
    rVariables.LiquidDensityRefPressure   = Prop[LIQUID_DENSITY_REF_PRESSURE];
    rVariables.BiotCoefficient            = Prop[BIOT_COEFFICIENT];
    rVariables.SolidCompressibilityCoeff  = (rVariables.BiotCoefficient-rVariables.Porosity)/Prop[BULK_MODULUS_SOLID];
    rVariables.LiquidCompressibilityCoeff = rVariables.Porosity/Prop[BULK_MODULUS_LIQUID];
    rVariables.GasCompressibilityCoeff    = rVariables.Porosity/Prop[BULK_MODULUS_GAS];
    rVariables.BiotModulusInverse         = rVariables.SolidCompressibilityCoeff + rVariables.LiquidCompressibilityCoeff;

    //ProcessInfo variables
    rVariables.VelocityCoefficient         = rCurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtLiquidPressureCoefficient = rCurrentProcessInfo[DT_LIQUID_PRESSURE_COEFFICIENT];
    rVariables.DtGasPressureCoefficient    = rCurrentProcessInfo[DT_GAS_PRESSURE_COEFFICIENT];

    //Nodal Variables
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        rVariables.LiquidPressureVector[i]   = Geom[i].FastGetSolutionStepValue(LIQUID_PRESSURE);
        rVariables.GasPressureVector[i]      = Geom[i].FastGetSolutionStepValue(GAS_PRESSURE);
        rVariables.DtLiquidPressureVector[i] = Geom[i].FastGetSolutionStepValue(DT_LIQUID_PRESSURE);
        rVariables.DtGasPressureVector[i]    = Geom[i].FastGetSolutionStepValue(DT_GAS_PRESSURE);
    }
    PoroElementUtilities::GetNodalVariableVector(rVariables.DisplacementVector,Geom,DISPLACEMENT);
    PoroElementUtilities::GetNodalVariableVector(rVariables.VelocityVector,Geom,VELOCITY);
    PoroElementUtilities::GetNodalVariableVector(rVariables.VolumeAcceleration,Geom,VOLUME_ACCELERATION);

    //General Variables
    const unsigned int strain_size = Prop.GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    const unsigned int cl_dimension = Prop.GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();

    rVariables.VoigtVector.resize(strain_size);
    noalias(rVariables.VoigtVector) = ZeroVector(strain_size);
    if(cl_dimension == 3) {
        rVariables.VoigtVector[0] = 1.0;
        rVariables.VoigtVector[1] = 1.0;
        rVariables.VoigtVector[2] = 1.0;
    } else {
        rVariables.VoigtVector[0] = 1.0;
        rVariables.VoigtVector[1] = 1.0;
    }

    //Variables computed at each GP
    rVariables.B.resize(strain_size,TNumNodes*TDim,false);
    noalias(rVariables.B) = ZeroMatrix(strain_size,TNumNodes*TDim);
    noalias(rVariables.Nu) = ZeroMatrix(TDim, TNumNodes*TDim);
    //Constitutive Law parameters
    rVariables.StrainVector.resize(strain_size,false);
    rVariables.StressVector.resize(strain_size,false);
    rVariables.ConstitutiveMatrix.resize(strain_size,strain_size,false);
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
    //Saturation Law parameters
    rSaturationParameters.SetSl(rVariables.Sl);
    rSaturationParameters.SetdSldPc(rVariables.dSldPc);
    rSaturationParameters.Setkrl(rVariables.krl);
    rSaturationParameters.Setkrg(rVariables.krg);
    rSaturationParameters.SetShapeFunctionsValues(rVariables.Np);
    rSaturationParameters.SetShapeFunctionsDerivatives(rVariables.GradNpT);
    //Auxiliary variables
    rVariables.UVoigtMatrix.resize(TNumNodes*TDim,strain_size,false);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::CalculateKinematics(Matrix& rGradNpT,
                                        Matrix& rB,
                                        Vector& rStrainVector,
                                        const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer,
                                        const array_1d<double,TNumNodes*TDim>& DisplacementVector,
                                        const unsigned int& GPoint)
{
    KRATOS_TRY

    noalias(rGradNpT) = DN_DXContainer[GPoint];
    this->CalculateBMatrix(rB, rGradNpT);
    noalias(rStrainVector) = prod(rB,DisplacementVector);

    // 2.5D element (2D Geometry with 3D ConstitutiveLaw)
    const unsigned int cl_dimension = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->WorkingSpaceDimension();
    if (cl_dimension > TDim) {

        // StrainVector must have the shape of a 3D element
        rStrainVector[3] = rStrainVector[2];
        rStrainVector[2] = mImposedZStrainVector[GPoint];

        unsigned int index;
        // B matrix must have the shape of a 3D element
        for ( unsigned int i = 0; i < TNumNodes; i++ )
        {
            index = 2 * i;

            rB( 3, index + 0 ) = rB( 2, index + 0 );
            rB( 3, index + 1 ) = rB( 2, index + 1 );
            rB( 2, index + 0 ) = 0.0;
            rB( 2, index + 1 ) = 0.0;
        }
    }
    
    KRATOS_CATCH( "" )

}
//----------------------------------------------------------------------------------------

template< >
void UPlPgSmallStrainElement<2,3>::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT)
{
    KRATOS_TRY

    unsigned int index;

    for ( unsigned int i = 0; i < 3; i++ )
    {
        index = 2 * i;

        rB( 0, index + 0 ) = GradNpT( i, 0 );
        rB( 1, index + 1 ) = GradNpT( i, 1 );
        rB( 2, index + 0 ) = GradNpT( i, 1 );
        rB( 2, index + 1 ) = GradNpT( i, 0 );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgSmallStrainElement<2,4>::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT)
{
    KRATOS_TRY

    unsigned int index;

    for ( unsigned int i = 0; i < 4; i++ )
    {
        index = 2 * i;

        rB( 0, index + 0 ) = GradNpT( i, 0 );
        rB( 1, index + 1 ) = GradNpT( i, 1 );
        rB( 2, index + 0 ) = GradNpT( i, 1 );
        rB( 2, index + 1 ) = GradNpT( i, 0 );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgSmallStrainElement<3,4>::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT)
{
    KRATOS_TRY

    unsigned int index;

    for ( unsigned int i = 0; i < 4; i++ )
    {
        index = 3 * i;

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

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgSmallStrainElement<3,8>::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT)
{
    KRATOS_TRY

    unsigned int index;

    for ( unsigned int i = 0; i < 8; i++ )
    {
        index = 3 * i;

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

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::GetHydromechanicalCouplingCoefficients(double& Clu, double& Cgu, ElementVariables& rVariables)
{
    Clu = rVariables.BiotCoefficient*rVariables.Sl;
    Cgu = rVariables.BiotCoefficient*(1.0 - rVariables.Sl) + 
            rVariables.GasHenrySolubilityCoeff * rVariables.Sl * (2.0 * rVariables.Porosity - rVariables.BiotCoefficient);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::GetCompressibilityCoefficients(double& Cll, double& Clg, double& Cgl, double& Cgg, const ElementVariables& Variables)
{

    // Get variables
    double Sl     = Variables.Sl;
    double dSldPc = Variables.dSldPc;
    double pc     = inner_prod(Variables.Np,(Variables.GasPressureVector - Variables.LiquidPressureVector)); // Capillary pore-pressure: pc = pg - pl
    double Ms     = Variables.SolidCompressibilityCoeff;
    double Ml     = Variables.LiquidCompressibilityCoeff;
    double Mg     = Variables.GasCompressibilityCoeff;
    double Phi    = Variables.Porosity;
    double Sg     = 1.0 - Sl;

    // Compute the compressibility coefficients
    Cll = Ms * Sl * (Sl + pc * dSldPc) - dSldPc * Phi + Sl * Ml;
    Clg = Ms * Sl * (Sg - pc * dSldPc) + dSldPc * Phi;
    Cgl = Ms * (Sg-Variables.GasHenrySolubilityCoeff*Sl) * (Sl + pc * dSldPc) + (1.0-Variables.GasHenrySolubilityCoeff) * dSldPc * Phi;
    Cgg = Ms * (Sg-Variables.GasHenrySolubilityCoeff*Sl) * (Sg - pc * dSldPc) + Mg * (Sg + Sl * Variables.GasHenrySolubilityCoeff)
            - (1.0-Variables.GasHenrySolubilityCoeff) * dSldPc * Phi;
    // TODO. Check whether Cgg should be like this
    // Cgg = Ms * (Sg-Variables.GasHenrySolubilityCoeff*Sl) * (Sg - pc * dSldPc) + Mg * (Sg + Sl * Variables.GasHenrySolubilityCoeff)
    //         - (1.0+Variables.GasHenrySolubilityCoeff) * dSldPc * Phi;
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgSmallStrainElement<TDim,TNumNodes>::UpdateDensity(ElementVariables& rVariables)
{
    //TODO. Ignasi
    // Gas density defined from the ideal gas law (Liakopoulos test in OGS)
    const double ideal_gas_constant = 8.31446261815324;
    const double pg = inner_prod(rVariables.Np,rVariables.GasPressureVector);
    rVariables.GasDensity = pg*rVariables.GasMolarWeight/(ideal_gas_constant*rVariables.Temperature);

    // Update the density based on the new saturation degree
    rVariables.Density = rVariables.Porosity * (rVariables.Sl * rVariables.LiquidDensity + (1.0 - rVariables.Sl) * rVariables.GasDensity) + 
                            (1.0 - rVariables.Porosity) * rVariables.SolidDensity;

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlPgSmallStrainElement<2,3>;
template class UPlPgSmallStrainElement<2,4>;
template class UPlPgSmallStrainElement<3,4>;
template class UPlPgSmallStrainElement<3,8>;

} // Namespace Kratos
