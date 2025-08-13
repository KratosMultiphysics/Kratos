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
#include "custom_elements/one-phase_flow/U_Pl_small_strain_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPlSmallStrainElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UPlSmallStrainElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPlSmallStrainElement<TDim,TNumNodes>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPlSmallStrainElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
int UPlSmallStrainElement<TDim,TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo ) const
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
    ierr = UPlElement<TDim,TNumNodes>::Check(rCurrentProcessInfo);
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

    return ierr;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
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
void UPlSmallStrainElement<TDim,TNumNodes>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    this->InitializeNonLinearIteration(rCurrentProcessInfo);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
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
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rOutput,
                                                                            const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != NumGPoints )
        rOutput.resize( NumGPoints, false );

    if(rVariable == VON_MISES_STRESS)
    {
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
    } else {
        UPlElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(rVariable,rOutput,rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<array_1d<double,3>>& rVariable,std::vector<array_1d<double,3>>& rOutput,
                                                                            const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != NumGPoints )
        rOutput.resize( NumGPoints );

    if(rVariable == LIQUID_FLUX_VECTOR) {
        const PropertiesType& Prop = this->GetProperties();

        //Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Defining necessary variables
        array_1d<double,TNumNodes> PressureVector;
        for(unsigned int i=0; i<TNumNodes; i++)
            PressureVector[i] = Geom[i].FastGetSolutionStepValue(LIQUID_PRESSURE);
        array_1d<double,TNumNodes*TDim> VolumeAcceleration;
        PoroElementUtilities::GetNodalVariableVector(VolumeAcceleration,Geom,VOLUME_ACCELERATION);
        array_1d<double,TDim> BodyAcceleration;
        BoundedMatrix<double,TNumNodes, TDim> GradNpT;
        const double& DynamicViscosityInverse = 1.0/Prop[DYNAMIC_VISCOSITY_LIQUID];
        const double& LiquidDensity = Prop[DENSITY_LIQUID];
        array_1d<double,TDim> GradPressureTerm;
        array_1d<double,TDim> LiquidFlux;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            PoroElementUtilities::InterpolateVariableWithComponents(BodyAcceleration,NContainer,VolumeAcceleration,GPoint);

            noalias(GradPressureTerm) = prod(trans(GradNpT),PressureVector);
            noalias(GradPressureTerm) += -LiquidDensity*BodyAcceleration;

            noalias(LiquidFlux) = -DynamicViscosityInverse*prod(mIntrinsicPermeability,GradPressureTerm);

            PoroElementUtilities::FillArray1dOutput(rOutput[GPoint],LiquidFlux);
        }
    } else if(rVariable == LIQUID_PRESSURE_GRADIENT) {

        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Defining necessary variables
        array_1d<double,TNumNodes> PressureVector;
        for(unsigned int i=0; i<TNumNodes; i++)
            PressureVector[i] = Geom[i].FastGetSolutionStepValue(LIQUID_PRESSURE);
        BoundedMatrix<double,TNumNodes, TDim> GradNpT;
        array_1d<double,TDim> GradPressure;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
        {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            noalias(GradPressure) = prod(trans(GradNpT),PressureVector);

            PoroElementUtilities::FillArray1dOutput(rOutput[GPoint],GradPressure);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput,
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
        const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

        Vector VoigtVector(strain_size);
        noalias(VoigtVector) = ZeroVector(strain_size);
        if(cl_dimension == 3) {
            VoigtVector[0] = 1.0;
            VoigtVector[1] = 1.0;
            VoigtVector[2] = 1.0;
        } else {
            VoigtVector[0] = 1.0;
            VoigtVector[1] = 1.0;
        }
        Matrix B(strain_size,TNumNodes*TDim);
        noalias(B) = ZeroMatrix(strain_size,TNumNodes*TDim);
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
        double Pressure;
        array_1d<double,TNumNodes> PressureVector;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            PressureVector[i] = Geom[i].FastGetSolutionStepValue(LIQUID_PRESSURE);
        }
        const double BiotCoefficient = Prop[BIOT_COEFFICIENT];

        //Create constitutive law parameters:
        Vector StrainVector(strain_size);
        Vector StressVector(strain_size);
        Matrix ConstitutiveMatrix(strain_size,strain_size);
        Vector Np(TNumNodes);
        Matrix GradNpT(TNumNodes,TDim);
        Matrix F = identity_matrix<double>(TDim);
        double detF = 1.0;
        ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
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

            Pressure = 0.0;
            for(unsigned int i = 0; i < TNumNodes; i++)
            {
                Pressure += NContainer(GPoint,i)*PressureVector[i];
            }

            noalias(StressVector) += -BiotCoefficient*Pressure*VoigtVector;

            rOutput[GPoint].resize(cl_dimension,cl_dimension,false );
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector);
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
    } else {
        UPlElement<TDim,TNumNodes>::CalculateOnIntegrationPoints(rVariable,rOutput,rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

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

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,rCurrentProcessInfo);

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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

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
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

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
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B),rVariables.ConstitutiveMatrix);
    noalias(rVariables.UMatrix) = prod(rVariables.UVoigtMatrix,rVariables.B)*rVariables.IntegrationCoefficient;

    //Distribute stiffness block matrix into the elemental matrix
    PoroElementUtilities::AssembleUBlockMatrix(rLeftHandSideMatrix,rVariables.UMatrix);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    noalias(rVariables.UVector) = prod(trans(rVariables.B),rVariables.VoigtVector);

    noalias(rVariables.UPMatrix) = -rVariables.BiotCoefficient*outer_prod(rVariables.UVector,rVariables.Np)*rVariables.IntegrationCoefficient;

    //Distribute coupling block matrix into the elemental matrix
    PoroElementUtilities::AssembleUPBlockMatrix(rLeftHandSideMatrix,rVariables.UPMatrix);

    noalias(rVariables.PUMatrix) = -rVariables.VelocityCoefficient*trans(rVariables.UPMatrix);

    //Distribute transposed coupling block matrix into the elemental matrix
    PoroElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix,rVariables.PUMatrix);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    noalias(rVariables.PMatrix) = rVariables.DtPressureCoefficient*rVariables.BiotModulusInverse*outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;

    //Distribute compressibility block matrix into the elemental matrix
    PoroElementUtilities::AssemblePBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT,mIntrinsicPermeability);

    noalias(rVariables.PMatrix) = rVariables.DynamicViscosityInverse*prod(rVariables.PDimMatrix,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;

    //Distribute permeability block matrix into the elemental matrix
    PoroElementUtilities::AssemblePBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rVariables.UVector) = -1.0*prod(trans(rVariables.B),rVariables.StressVector)*rVariables.IntegrationCoefficient;

    //Distribute stiffness block vector into elemental vector
    PoroElementUtilities::AssembleUBlockVector(rRightHandSideVector,rVariables.UVector);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rVariables.UVector) = rVariables.Density*prod(trans(rVariables.Nu),rVariables.BodyAcceleration)*rVariables.IntegrationCoefficient;

    //Distribute body force block vector into elemental vector
    PoroElementUtilities::AssembleUBlockVector(rRightHandSideVector,rVariables.UVector);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rVariables.UVector) = prod(trans(rVariables.B),rVariables.VoigtVector);

    noalias(rVariables.UPMatrix) = rVariables.BiotCoefficient*outer_prod(rVariables.UVector,rVariables.Np)*rVariables.IntegrationCoefficient;

    noalias(rVariables.UVector) = prod(rVariables.UPMatrix,rVariables.PressureVector);

    //Distribute coupling block vector 1 into elemental vector
    PoroElementUtilities::AssembleUBlockVector(rRightHandSideVector,rVariables.UVector);

    noalias(rVariables.PVector) = -1.0*prod(trans(rVariables.UPMatrix),rVariables.VelocityVector);

    //Distribute coupling block vector 2 into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rVariables.PMatrix) = rVariables.BiotModulusInverse*outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = -1.0*prod(rVariables.PMatrix,rVariables.DtPressureVector);

    //Distribute compressibility block vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT,mIntrinsicPermeability);

    noalias(rVariables.PMatrix) = rVariables.DynamicViscosityInverse*prod(rVariables.PDimMatrix,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = -1.0*prod(rVariables.PMatrix,rVariables.PressureVector);

    //Distribute permeability block vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT,mIntrinsicPermeability)*rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = rVariables.DynamicViscosityInverse*rVariables.LiquidDensity*
                                    prod(rVariables.PDimMatrix,rVariables.BodyAcceleration);

    //Distribute liquid body flow block vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateFluxResidual( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

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

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

        //Compute stresses
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Contributions to the right hand side
        this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, Variables);
        this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, Variables);
        this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateMixBodyForce( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

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

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

        //Compute stresses
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Contributions to the right hand side
        this->CalculateAndAddMixBodyForce(rRightHandSideVector, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateNegInternalForce( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

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

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

        //Compute stresses
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Contributions to the right hand side
        this->CalculateAndAddStiffnessForce(rRightHandSideVector, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateExplicitContributions (VectorType& rFluxResidual, VectorType& rBodyForce, VectorType& rNegInternalForces, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

    //Resetting the RHS
    if ( rFluxResidual.size() != element_size )
        rFluxResidual.resize( element_size, false );
    noalias( rFluxResidual ) = ZeroVector( element_size );
    if ( rBodyForce.size() != element_size )
        rBodyForce.resize( element_size, false );
    noalias( rBodyForce ) = ZeroVector( element_size );
    if ( rNegInternalForces.size() != element_size )
        rNegInternalForces.resize( element_size, false );
    noalias( rNegInternalForces ) = ZeroVector( element_size );

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

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,rCurrentProcessInfo);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables.GradNpT,Variables.B,Variables.StrainVector,DN_DXContainer,Variables.DisplacementVector,GPoint);

        //Compute Np, Nu and BodyAcceleration
        noalias(Variables.Np) = row(NContainer,GPoint);
        PoroElementUtilities::CalculateNuMatrix(Variables.Nu,NContainer,GPoint);
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);

        //Compute stresses
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Contributions to the right hand side
        this->CalculateAndAddCompressibilityFlow(rFluxResidual, Variables);
        this->CalculateAndAddPermeabilityFlow(rFluxResidual, Variables);
        this->CalculateAndAddFluidBodyFlow(rFluxResidual, Variables);

        this->CalculateAndAddMixBodyForce(rBodyForce, Variables);

        this->CalculateAndAddStiffnessForce(rNegInternalForces, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::SaveGPStress(Matrix& rStressContainer, const Vector& StressVector, const unsigned int& VoigtSize, const unsigned int& GPoint)
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
void UPlSmallStrainElement<2,3>::ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
{
    KRATOS_TRY

    // Triangle_2d_3 with GI_GAUSS_2

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area();
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
void UPlSmallStrainElement<2,4>::ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
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
void UPlSmallStrainElement<3,4>::ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
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
void UPlSmallStrainElement<3,8>::ExtrapolateGPValues(const Matrix& StressContainer, const unsigned int& VoigtSize)
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::InitializeElementVariables(ElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                                                                  const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //Properties variables
    const double& BulkModulusSolid = Prop[BULK_MODULUS_SOLID];
    const double& Porosity = Prop[POROSITY];
    rVariables.DynamicViscosityInverse = 1.0/Prop[DYNAMIC_VISCOSITY_LIQUID];
    rVariables.LiquidDensity = Prop[DENSITY_LIQUID];
    rVariables.Density = Porosity*rVariables.LiquidDensity + (1.0-Porosity)*Prop[DENSITY_SOLID];
    rVariables.BiotCoefficient = Prop[BIOT_COEFFICIENT];
    rVariables.BiotModulusInverse = (rVariables.BiotCoefficient-Porosity)/BulkModulusSolid + Porosity/Prop[BULK_MODULUS_LIQUID];

    //ProcessInfo variables
    rVariables.VelocityCoefficient = rCurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_LIQUID_PRESSURE_COEFFICIENT];

    //Nodal Variables
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        rVariables.PressureVector[i] = Geom[i].FastGetSolutionStepValue(LIQUID_PRESSURE);
        rVariables.DtPressureVector[i] = Geom[i].FastGetSolutionStepValue(DT_LIQUID_PRESSURE);
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
    //Auxiliary variables
    rVariables.UVoigtMatrix.resize(TNumNodes*TDim,strain_size,false);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlSmallStrainElement<TDim,TNumNodes>::CalculateKinematics(Matrix& rGradNpT,
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
void UPlSmallStrainElement<2,3>::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT)
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
void UPlSmallStrainElement<2,4>::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT)
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
void UPlSmallStrainElement<3,4>::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT)
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
void UPlSmallStrainElement<3,8>::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT)
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlSmallStrainElement<2,3>;
template class UPlSmallStrainElement<2,4>;
template class UPlSmallStrainElement<3,4>;
template class UPlSmallStrainElement<3,8>;

} // Namespace Kratos
