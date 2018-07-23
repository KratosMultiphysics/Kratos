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
#include "custom_elements/U_Pw_small_strain_interface_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwSmallStrainInterfaceElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UPwSmallStrainInterfaceElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwSmallStrainInterfaceElement<TDim,TNumNodes>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwSmallStrainInterfaceElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
int UPwSmallStrainInterfaceElement<TDim,TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const PropertiesType& Prop = this->GetProperties();

    if (this->Id() < 1)
        KRATOS_THROW_ERROR(std::logic_error, "Element found with Id 0 or negative","")

    // Verify generic variables
    int ierr = UPwElement<TDim,TNumNodes>::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Verify specific properties
    if ( MINIMUM_JOINT_WIDTH.Key() == 0 || Prop.Has( MINIMUM_JOINT_WIDTH ) == false || Prop[MINIMUM_JOINT_WIDTH] <= 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"MINIMUM_JOINT_WIDTH has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( TRANSVERSAL_PERMEABILITY.Key() == 0 || Prop.Has( TRANSVERSAL_PERMEABILITY ) == false || Prop[TRANSVERSAL_PERMEABILITY] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"TRANSVERSAL_PERMEABILITY has Key zero, is not defined or has an invalid value at element", this->Id() )

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
        ierr = Prop[CONSTITUTIVE_LAW]->Check( Prop, this->GetGeometry(), rCurrentProcessInfo );
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "A constitutive law needs to be specified for the element ", this->Id() )

    return ierr;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::Initialize()
{
    KRATOS_TRY

    UPwElement<TDim,TNumNodes>::Initialize();

    //Compute initial gap of the joint
    this->CalculateInitialGap(this->GetGeometry());

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);

    //Resizing mass matrix
    if ( rMassMatrix.size1() != element_size )
        rMassMatrix.resize( element_size, element_size, false );
    noalias( rMassMatrix ) = ZeroMatrix( element_size, element_size );

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();

    //Defining shape functions and the determinant of the jacobian at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    Vector detJContainer(NumGPoints);
    Geom.DeterminantOfJacobian(detJContainer,mThisIntegrationMethod);

    //Defining necessary variables
    double IntegrationCoefficient;
    const double& Porosity = Prop[POROSITY];
    const double Density = Porosity*Prop[DENSITY_WATER] + (1.0-Porosity)*Prop[DENSITY_SOLID];
    BoundedMatrix<double,TDim+1, TNumNodes*(TDim+1)> Nut = ZeroMatrix(TDim+1, TNumNodes*(TDim+1));
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
    BoundedMatrix<double,TDim, TDim> RotationMatrix;
    this->CalculateRotationMatrix(RotationMatrix,Geom);
    BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
    array_1d<double,TDim> LocalRelDispVector;
    array_1d<double,TDim> RelDispVector;
    const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
    double JointWidth;

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
    {
        InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

        noalias(RelDispVector) = prod(Nu,DisplacementVector);

        noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);

        this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], MinimumJointWidth,GPoint);

        InterfaceElementUtilities::CalculateNuElementMatrix(Nut,NContainer,GPoint);

        //calculating weighting coefficient for integration
        this->CalculateIntegrationCoefficient( IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Adding contribution to Mass matrix
        noalias(rMassMatrix) += Density*prod(trans(Nut),Nut)*JointWidth*IntegrationCoefficient;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //Defining necessary variables
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
    BoundedMatrix<double,TDim, TDim> RotationMatrix;
    this->CalculateRotationMatrix(RotationMatrix,Geom);
    BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
    array_1d<double,TDim> RelDispVector;
    const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
    double JointWidth;

    //Create constitutive law parameters:
    Vector StrainVector(TDim);
    Vector StressVector(TDim);
    Matrix ConstitutiveMatrix(TDim,TDim);
    Vector Np(TNumNodes);
    Matrix GradNpT(TNumNodes,TDim);
    Matrix F = identity_matrix<double>(TDim);
    double detF = 1.0;
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,rCurrentProcessInfo);
    ConstitutiveParameters.SetConstitutiveMatrix(ConstitutiveMatrix);
    ConstitutiveParameters.SetStressVector(StressVector);
    ConstitutiveParameters.SetStrainVector(StrainVector);
    ConstitutiveParameters.SetShapeFunctionsValues(Np);
    ConstitutiveParameters.SetShapeFunctionsDerivatives(GradNpT);
    ConstitutiveParameters.SetDeterminantF(detF);
    ConstitutiveParameters.SetDeformationGradientF(F);

    // Auxiliar output variables
    unsigned int NumGPoints = mConstitutiveLawVector.size();
    std::vector<double> JointWidthContainer(NumGPoints);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
    {
        InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

        noalias(RelDispVector) = prod(Nu,DisplacementVector);

        noalias(StrainVector) = prod(RotationMatrix,RelDispVector);

        JointWidthContainer[GPoint] = mInitialGap[GPoint] + StrainVector[TDim-1];

        this->CheckAndCalculateJointWidth(JointWidth, ConstitutiveParameters, StrainVector[TDim-1], MinimumJointWidth, GPoint);

        noalias(Np) = row(NContainer,GPoint);

        //compute constitutive tensor and/or stresses
        mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
    }

    if(rCurrentProcessInfo[NODAL_SMOOTHING] == true)
    {
        this->ExtrapolateGPValues(JointWidthContainer);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainInterfaceElement<2,4>::ExtrapolateGPValues (const std::vector<double>& JointWidthContainer)
{
    array_1d<double,2> DamageContainer; // 2 LobattoPoints

    for ( unsigned int i = 0;  i < 2; i++ ) // NumLobattoPoints
    {
        DamageContainer[i] = 0.0;
        DamageContainer[i] = mConstitutiveLawVector[i]->GetValue( DAMAGE_VARIABLE, DamageContainer[i] );
    }

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area();

    array_1d<double,4> NodalJointWidth;
    NodalJointWidth[0] = JointWidthContainer[0]*Area;
    NodalJointWidth[1] = JointWidthContainer[1]*Area;
    NodalJointWidth[2] = JointWidthContainer[1]*Area;
    NodalJointWidth[3] = JointWidthContainer[0]*Area;

    array_1d<double,4> NodalDamage;
    NodalDamage[0] = DamageContainer[0]*Area;
    NodalDamage[1] = DamageContainer[1]*Area;
    NodalDamage[2] = DamageContainer[1]*Area;
    NodalDamage[3] = DamageContainer[0]*Area;

    for(unsigned int i = 0; i < 4; i++) //NumNodes
    {
        rGeom[i].SetLock();
        rGeom[i].FastGetSolutionStepValue(NODAL_JOINT_WIDTH) += NodalJointWidth[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) += NodalDamage[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_JOINT_AREA) += Area;
        rGeom[i].UnSetLock();
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainInterfaceElement<3,6>::ExtrapolateGPValues (const std::vector<double>& JointWidthContainer)
{
    array_1d<double,3> DamageContainer; // 3 LobattoPoints

    for ( unsigned int i = 0;  i < 3; i++ ) // NumLobattoPoints
    {
        DamageContainer[i] = 0.0;
        DamageContainer[i] = mConstitutiveLawVector[i]->GetValue( DAMAGE_VARIABLE, DamageContainer[i] );
    }

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area();

    array_1d<double,6> NodalJointWidth;
    NodalJointWidth[0] = JointWidthContainer[0]*Area;
    NodalJointWidth[1] = JointWidthContainer[1]*Area;
    NodalJointWidth[2] = JointWidthContainer[2]*Area;
    NodalJointWidth[3] = JointWidthContainer[0]*Area;
    NodalJointWidth[4] = JointWidthContainer[1]*Area;
    NodalJointWidth[5] = JointWidthContainer[2]*Area;

    array_1d<double,6> NodalDamage;
    NodalDamage[0] = DamageContainer[0]*Area;
    NodalDamage[1] = DamageContainer[1]*Area;
    NodalDamage[2] = DamageContainer[2]*Area;
    NodalDamage[3] = DamageContainer[0]*Area;
    NodalDamage[4] = DamageContainer[1]*Area;
    NodalDamage[5] = DamageContainer[2]*Area;

    for(unsigned int i = 0; i < 6; i++) //NumNodes
    {
        rGeom[i].SetLock();
        rGeom[i].FastGetSolutionStepValue(NODAL_JOINT_WIDTH) += NodalJointWidth[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) += NodalDamage[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_JOINT_AREA) += Area;
        rGeom[i].UnSetLock();
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainInterfaceElement<3,8>::ExtrapolateGPValues (const std::vector<double>& JointWidthContainer)
{
    array_1d<double,4> DamageContainer; // 4 LobattoPoints

    for ( unsigned int i = 0;  i < 4; i++ ) // NumLobattoPoints
    {
        DamageContainer[i] = 0.0;
        DamageContainer[i] = mConstitutiveLawVector[i]->GetValue( DAMAGE_VARIABLE, DamageContainer[i] );
    }

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area();

    array_1d<double,8> NodalJointWidth;
    NodalJointWidth[0] = JointWidthContainer[0]*Area;
    NodalJointWidth[1] = JointWidthContainer[1]*Area;
    NodalJointWidth[2] = JointWidthContainer[2]*Area;
    NodalJointWidth[3] = JointWidthContainer[3]*Area;
    NodalJointWidth[4] = JointWidthContainer[0]*Area;
    NodalJointWidth[5] = JointWidthContainer[1]*Area;
    NodalJointWidth[6] = JointWidthContainer[2]*Area;
    NodalJointWidth[7] = JointWidthContainer[3]*Area;

    array_1d<double,8> NodalDamage;
    NodalDamage[0] = DamageContainer[0]*Area;
    NodalDamage[1] = DamageContainer[1]*Area;
    NodalDamage[2] = DamageContainer[2]*Area;
    NodalDamage[3] = DamageContainer[3]*Area;
    NodalDamage[4] = DamageContainer[0]*Area;
    NodalDamage[5] = DamageContainer[1]*Area;
    NodalDamage[6] = DamageContainer[2]*Area;
    NodalDamage[7] = DamageContainer[3]*Area;

    for(unsigned int i = 0; i < 8; i++) //NumNodes
    {
        rGeom[i].SetLock();
        rGeom[i].FastGetSolutionStepValue(NODAL_JOINT_WIDTH) += NodalJointWidth[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) += NodalDamage[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_JOINT_AREA) += Area;
        rGeom[i].UnSetLock();
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                                                                    std::vector<double>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if(rVariable == DAMAGE_VARIABLE)
    {
        //Variables computed on Lobatto points
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );
        std::vector<double> GPValues(NumGPoints);

        for ( unsigned int i = 0;  i < NumGPoints; i++ )
            GPValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, GPValues[i] );

        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber( GeometryData::GI_GAUSS_2 );
        if ( rValues.size() != OutputGPoints )
            rValues.resize( OutputGPoints );

        this->InterpolateOutputDoubles(rValues,GPValues);
    }
    else if(rVariable == STATE_VARIABLE)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        for ( unsigned int i = 0;  i < mConstitutiveLawVector.size(); i++ )
            rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
    }
    else if(rVariable == JOINT_WIDTH)
    {
        //Variables computed on Lobatto points
        const GeometryType& Geom = this->GetGeometry();

        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );
        std::vector<array_1d<double,3>> GPAuxValues(NumGPoints);
        this->CalculateOnIntegrationPoints(LOCAL_RELATIVE_DISPLACEMENT_VECTOR, GPAuxValues, rCurrentProcessInfo);

        std::vector<double> GPValues(NumGPoints);

        for(unsigned int i=0; i < NumGPoints; i++)
        {
            GPValues[i] = mInitialGap[i] + GPAuxValues[i][TDim-1];
        }

        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber( GeometryData::GI_GAUSS_2 );
        if ( rValues.size() != OutputGPoints )
            rValues.resize( OutputGPoints );

        this->InterpolateOutputDoubles(rValues,GPValues);
    }
    else
    {
        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = this->GetGeometry().IntegrationPointsNumber( GeometryData::GI_GAUSS_2 );
        if ( rValues.size() != OutputGPoints )
            rValues.resize( OutputGPoints );

        for(unsigned int i=0; i < OutputGPoints; i++)
        {
            rValues[i] = 0.0;
        }
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,
                                                                                    std::vector<array_1d<double,3>>& rValues,const ProcessInfo& rCurrentProcessInfo)
{
    if(rVariable == FLUID_FLUX_VECTOR || rVariable == LOCAL_STRESS_VECTOR || rVariable == LOCAL_RELATIVE_DISPLACEMENT_VECTOR || rVariable == LOCAL_FLUID_FLUX_VECTOR)
    {
        //Variables computed on Lobatto points
        const GeometryType& Geom = this->GetGeometry();
        std::vector<array_1d<double,3>> GPValues(Geom.IntegrationPointsNumber( mThisIntegrationMethod ));

        this->CalculateOnIntegrationPoints(rVariable, GPValues, rCurrentProcessInfo);

        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber( GeometryData::GI_GAUSS_2 );
        if ( rValues.size() != OutputGPoints )
            rValues.resize( OutputGPoints );

        this->InterpolateOutputValues< array_1d<double,3> >(rValues,GPValues);
    }
    else
    {
        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = this->GetGeometry().IntegrationPointsNumber( GeometryData::GI_GAUSS_2 );
        if ( rValues.size() != OutputGPoints )
            rValues.resize( OutputGPoints );

        for(unsigned int i=0; i < OutputGPoints; i++)
        {
            noalias(rValues[i]) = ZeroVector(3);
        }
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,std::vector<Matrix>& rValues,
                                                                                    const ProcessInfo& rCurrentProcessInfo)
{
    if(rVariable == PERMEABILITY_MATRIX || rVariable == LOCAL_PERMEABILITY_MATRIX)
    {
        //Variables computed on Lobatto points
        const GeometryType& Geom = this->GetGeometry();
        std::vector<Matrix> GPValues(Geom.IntegrationPointsNumber( mThisIntegrationMethod ));

        this->CalculateOnIntegrationPoints(rVariable, GPValues, rCurrentProcessInfo);

        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber( GeometryData::GI_GAUSS_2 );
        if ( rValues.size() != OutputGPoints )
            rValues.resize( OutputGPoints );

        for(unsigned int GPoint=0; GPoint<OutputGPoints; GPoint++)
            rValues[GPoint].resize(TDim,TDim,false);

        this->InterpolateOutputValues< Matrix >(rValues,GPValues);
    }
    else
    {
        //Printed on standard GiD Gauss points
        const unsigned int OutputGPoints = this->GetGeometry().IntegrationPointsNumber( GeometryData::GI_GAUSS_2 );
        if ( rValues.size() != OutputGPoints )
            rValues.resize( OutputGPoints );

        for(unsigned int i=0; i < OutputGPoints; i++)
        {
            rValues[i].resize(TDim,TDim,false);
            noalias(rValues[i]) = ZeroMatrix(TDim,TDim);
        }
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<array_1d<double,3>>& rVariable,
                                                                                std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if(rVariable == FLUID_FLUX_VECTOR)
    {
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

        //Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );
        GeometryType::JacobiansType JContainer(NumGPoints);
        Geom.Jacobian( JContainer, mThisIntegrationMethod );

        //Defining necessary variables
        array_1d<double,TNumNodes> PressureVector;
        for(unsigned int i=0; i<TNumNodes; i++)
            PressureVector[i] = Geom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        array_1d<double,TNumNodes*TDim> VolumeAcceleration;
        PoroElementUtilities::GetNodalVariableVector(VolumeAcceleration,Geom,VOLUME_ACCELERATION);
        array_1d<double,TDim> BodyAcceleration;
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
        BoundedMatrix<double,TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix,Geom);
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
        array_1d<double,TDim> LocalRelDispVector;
        array_1d<double,TDim> RelDispVector;
        const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
        double JointWidth;
        BoundedMatrix<double,TNumNodes, TDim> GradNpT;
        const double& Transversal_Permeability = Prop[TRANSVERSAL_PERMEABILITY];
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

            this->CalculateShapeFunctionsGradients< BoundedMatrix<double,TNumNodes,TDim> >(GradNpT,SFGradAuxVars,JContainer[GPoint],RotationMatrix,
                                                                                                                DN_DeContainer[GPoint],NContainer,JointWidth,GPoint);

            PoroElementUtilities::InterpolateVariableWithComponents(BodyAcceleration,NContainer,VolumeAcceleration,GPoint);

            InterfaceElementUtilities::CalculatePermeabilityMatrix(LocalPermeabilityMatrix,JointWidth,Transversal_Permeability);

            noalias(GradPressureTerm) = prod(trans(GradNpT),PressureVector);
            noalias(GradPressureTerm) += -FluidDensity*BodyAcceleration;

            noalias(LocalFluidFlux) = -DynamicViscosityInverse*prod(LocalPermeabilityMatrix,GradPressureTerm);

            noalias(FluidFlux) = prod(trans(RotationMatrix),LocalFluidFlux);

            PoroElementUtilities::FillArray1dOutput(rOutput[GPoint],FluidFlux);
        }
    }
    else if(rVariable == LOCAL_STRESS_VECTOR)
    {
        //Defining necessary variables
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
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
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
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
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            noalias(LocalStressVector) = StressVectorDynamic;

            PoroElementUtilities::FillArray1dOutput(rOutput[GPoint],LocalStressVector);
        }
    }
    else if(rVariable == LOCAL_RELATIVE_DISPLACEMENT_VECTOR)
    {
        //Defining necessary variables
        const GeometryType& Geom = this->GetGeometry();
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

            PoroElementUtilities::FillArray1dOutput(rOutput[GPoint],LocalRelDispVector);
        }
    }
    else if(rVariable == LOCAL_FLUID_FLUX_VECTOR)
    {
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );

        //Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );
        GeometryType::JacobiansType JContainer(NumGPoints);
        Geom.Jacobian( JContainer, mThisIntegrationMethod );

        //Defining necessary variables
        array_1d<double,TNumNodes> PressureVector;
        for(unsigned int i=0; i<TNumNodes; i++)
            PressureVector[i] = Geom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        array_1d<double,TNumNodes*TDim> VolumeAcceleration;
        PoroElementUtilities::GetNodalVariableVector(VolumeAcceleration,Geom,VOLUME_ACCELERATION);
        array_1d<double,TDim> BodyAcceleration;
        array_1d<double,TNumNodes*TDim> DisplacementVector;
        PoroElementUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
        BoundedMatrix<double,TDim, TDim> RotationMatrix;
        this->CalculateRotationMatrix(RotationMatrix,Geom);
        BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
        array_1d<double,TDim> LocalRelDispVector;
        array_1d<double,TDim> RelDispVector;
        const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
        double JointWidth;
        BoundedMatrix<double,TNumNodes, TDim> GradNpT;
        const double& Transversal_Permeability = Prop[TRANSVERSAL_PERMEABILITY];
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

            this->CalculateShapeFunctionsGradients< BoundedMatrix<double,TNumNodes,TDim> >(GradNpT,SFGradAuxVars,JContainer[GPoint],RotationMatrix,
                                                                                                                DN_DeContainer[GPoint],NContainer,JointWidth,GPoint);

            PoroElementUtilities::InterpolateVariableWithComponents(BodyAcceleration,NContainer,VolumeAcceleration,GPoint);

            InterfaceElementUtilities::CalculatePermeabilityMatrix(LocalPermeabilityMatrix,JointWidth,Transversal_Permeability);

            noalias(GradPressureTerm) = prod(trans(GradNpT),PressureVector);
            noalias(GradPressureTerm) += -FluidDensity*BodyAcceleration;

            noalias(LocalFluidFlux) = -DynamicViscosityInverse*prod(LocalPermeabilityMatrix,GradPressureTerm);

            PoroElementUtilities::FillArray1dOutput(rOutput[GPoint],LocalFluidFlux);
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput,
                                                                                    const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if(rVariable == PERMEABILITY_MATRIX)
    {
        const GeometryType& Geom = this->GetGeometry();
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
        const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
        double JointWidth;
        const double& Transversal_Permeability = Prop[TRANSVERSAL_PERMEABILITY];
        BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim,TDim);
        BoundedMatrix<double,TDim, TDim> PermeabilityMatrix;

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); GPoint++ )
        {
            InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

            noalias(RelDispVector) = prod(Nu,DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);

            this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], MinimumJointWidth,GPoint);

            InterfaceElementUtilities::CalculatePermeabilityMatrix(LocalPermeabilityMatrix,JointWidth,Transversal_Permeability);

            noalias(PermeabilityMatrix) = prod(trans(RotationMatrix),BoundedMatrix<double,TDim, TDim>(prod(LocalPermeabilityMatrix,RotationMatrix)));

            rOutput[GPoint].resize(TDim,TDim,false);
            noalias(rOutput[GPoint]) = PermeabilityMatrix;
        }
    }
    else if(rVariable == LOCAL_PERMEABILITY_MATRIX)
    {
        const GeometryType& Geom = this->GetGeometry();
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
        const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
        double JointWidth;
        const double& Transversal_Permeability = Prop[TRANSVERSAL_PERMEABILITY];
        BoundedMatrix<double,TDim, TDim> LocalPermeabilityMatrix = ZeroMatrix(TDim,TDim);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); GPoint++ )
        {
            InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);

            noalias(RelDispVector) = prod(Nu,DisplacementVector);

            noalias(LocalRelDispVector) = prod(RotationMatrix,RelDispVector);

            this->CalculateJointWidth(JointWidth, LocalRelDispVector[TDim-1], MinimumJointWidth,GPoint);

            InterfaceElementUtilities::CalculatePermeabilityMatrix(LocalPermeabilityMatrix,JointWidth,Transversal_Permeability);

            rOutput[GPoint].resize(TDim,TDim,false);
            noalias(rOutput[GPoint]) = LocalPermeabilityMatrix;
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< >
void UPwSmallStrainInterfaceElement<2,4>::CalculateInitialGap(const GeometryType& Geom)
{
    const double& MinimumJointWidth = this->GetProperties()[MINIMUM_JOINT_WIDTH];

    mInitialGap.resize(2);
    mIsOpen.resize(2);

    array_1d<double,3> Vx;
    noalias(Vx) = Geom.GetPoint( 3 ) - Geom.GetPoint( 0 );
    mInitialGap[0] = norm_2(Vx);
    if(mInitialGap[0] < MinimumJointWidth)
        mIsOpen[0] = false;
    else
        mIsOpen[0] = true;

    noalias(Vx) = Geom.GetPoint( 2 ) - Geom.GetPoint( 1 );
    mInitialGap[1] = norm_2(Vx);
    if(mInitialGap[1] < MinimumJointWidth)
        mIsOpen[1] = false;
    else
        mIsOpen[1] = true;
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainInterfaceElement<3,6>::CalculateInitialGap(const GeometryType& Geom)
{
    const double& MinimumJointWidth = this->GetProperties()[MINIMUM_JOINT_WIDTH];

    mInitialGap.resize(3);
    mIsOpen.resize(3);

    array_1d<double,3> Vx;
    noalias(Vx) = Geom.GetPoint( 3 ) - Geom.GetPoint( 0 );
    mInitialGap[0] = norm_2(Vx);
    if(mInitialGap[0] < MinimumJointWidth)
        mIsOpen[0] = false;
    else
        mIsOpen[0] = true;

    noalias(Vx) = Geom.GetPoint( 4 ) - Geom.GetPoint( 1 );
    mInitialGap[1] = norm_2(Vx);
    if(mInitialGap[1] < MinimumJointWidth)
        mIsOpen[1] = false;
    else
        mIsOpen[1] = true;

    noalias(Vx) = Geom.GetPoint( 5 ) - Geom.GetPoint( 2 );
    mInitialGap[2] = norm_2(Vx);
    if(mInitialGap[2] < MinimumJointWidth)
        mIsOpen[2] = false;
    else
        mIsOpen[2] = true;
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainInterfaceElement<3,8>::CalculateInitialGap(const GeometryType& Geom)
{
    const double& MinimumJointWidth = this->GetProperties()[MINIMUM_JOINT_WIDTH];

    mInitialGap.resize(4);
    mIsOpen.resize(4);

    array_1d<double,3> Vx;
    noalias(Vx) = Geom.GetPoint( 4 ) - Geom.GetPoint( 0 );
    mInitialGap[0] = norm_2(Vx);
    if(mInitialGap[0] < MinimumJointWidth)
        mIsOpen[0] = false;
    else
        mIsOpen[0] = true;

    noalias(Vx) = Geom.GetPoint( 5 ) - Geom.GetPoint( 1 );
    mInitialGap[1] = norm_2(Vx);
    if(mInitialGap[1] < MinimumJointWidth)
        mIsOpen[1] = false;
    else
        mIsOpen[1] = true;

    noalias(Vx) = Geom.GetPoint( 6 ) - Geom.GetPoint( 2 );
    mInitialGap[2] = norm_2(Vx);
    if(mInitialGap[2] < MinimumJointWidth)
        mIsOpen[2] = false;
    else
        mIsOpen[2] = true;

    noalias(Vx) = Geom.GetPoint( 7 ) - Geom.GetPoint( 3 );
    mInitialGap[3] = norm_2(Vx);
    if(mInitialGap[3] < MinimumJointWidth)
        mIsOpen[3] = false;
    else
        mIsOpen[3] = true;
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo )
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
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer(NumGPoints);
    Geom.Jacobian( JContainer, mThisIntegrationMethod );
    Vector detJContainer(NumGPoints);
    Geom.DeterminantOfJacobian(detJContainer,mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,CurrentProcessInfo);

    //Auxiliary variables
    const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
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
        this->CheckAndCalculateJointWidth(Variables.JointWidth,ConstitutiveParameters,Variables.StrainVector[TDim-1], MinimumJointWidth, GPoint);
        this->CalculateShapeFunctionsGradients< Matrix >(Variables.GradNpT,SFGradAuxVars,JContainer[GPoint],Variables.RotationMatrix,
                                                        DN_DeContainer[GPoint],NContainer,Variables.JointWidth,GPoint);

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
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
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
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);

    //Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,CurrentProcessInfo);

    //Auxiliary variables
    const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
    const double& Transversal_Permeability = Prop[TRANSVERSAL_PERMEABILITY];
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
        this->CheckAndCalculateJointWidth(Variables.JointWidth,ConstitutiveParameters,Variables.StrainVector[TDim-1], MinimumJointWidth, GPoint);
        this->CalculateShapeFunctionsGradients< Matrix >(Variables.GradNpT,SFGradAuxVars,JContainer[GPoint],Variables.RotationMatrix,
                                                        DN_DeContainer[GPoint],NContainer,Variables.JointWidth,GPoint);

        //Compute BodyAcceleration and Permeability Matrix
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);
        InterfaceElementUtilities::CalculatePermeabilityMatrix(Variables.LocalPermeabilityMatrix,Variables.JointWidth,Transversal_Permeability);

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
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
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
    ConstitutiveLaw::Parameters ConstitutiveParameters(Geom,Prop,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);

    //Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables,ConstitutiveParameters,Geom,Prop,CurrentProcessInfo);

    //Auxiliary variables
    const double& MinimumJointWidth = Prop[MINIMUM_JOINT_WIDTH];
    const double& Transversal_Permeability = Prop[TRANSVERSAL_PERMEABILITY];
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
        this->CheckAndCalculateJointWidth(Variables.JointWidth,ConstitutiveParameters,Variables.StrainVector[TDim-1], MinimumJointWidth, GPoint);
        this->CalculateShapeFunctionsGradients< Matrix >(Variables.GradNpT,SFGradAuxVars,JContainer[GPoint],Variables.RotationMatrix,
                                                        DN_DeContainer[GPoint],NContainer,Variables.JointWidth,GPoint);

        //Compute BodyAcceleration and Permeability Matrix
        PoroElementUtilities::InterpolateVariableWithComponents(Variables.BodyAcceleration,NContainer,Variables.VolumeAcceleration,GPoint);
        InterfaceElementUtilities::CalculatePermeabilityMatrix(Variables.LocalPermeabilityMatrix,Variables.JointWidth,Transversal_Permeability);

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
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::InitializeElementVariables(InterfaceElementVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                                                                  const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    //Properties variables
    const double& BulkModulusSolid = Prop[BULK_MODULUS_SOLID];
    const double& Porosity = Prop[POROSITY];
    const double BulkModulus = Prop[YOUNG_MODULUS]/(3.0*(1.0-2.0*Prop[POISSON_RATIO]));
    rVariables.DynamicViscosityInverse = 1.0/Prop[DYNAMIC_VISCOSITY];
    rVariables.FluidDensity = Prop[DENSITY_WATER];
    rVariables.Density = Porosity*rVariables.FluidDensity + (1.0-Porosity)*Prop[DENSITY_SOLID];
    rVariables.BiotCoefficient = 1.0-BulkModulus/BulkModulusSolid;
    rVariables.BiotModulusInverse = (rVariables.BiotCoefficient-Porosity)/BulkModulusSolid + Porosity/Prop[BULK_MODULUS_FLUID];

    //ProcessInfo variables
    rVariables.VelocityCoefficient = CurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtPressureCoefficient = CurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    //Nodal Variables
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        rVariables.PressureVector[i] = Geom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.DtPressureVector[i] = Geom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }
    PoroElementUtilities::GetNodalVariableVector(rVariables.DisplacementVector,Geom,DISPLACEMENT);
    PoroElementUtilities::GetNodalVariableVector(rVariables.VelocityVector,Geom,VELOCITY);
    PoroElementUtilities::GetNodalVariableVector(rVariables.VolumeAcceleration,Geom,VOLUME_ACCELERATION);

    //General Variables
    this->CalculateRotationMatrix(rVariables.RotationMatrix,Geom);
    InterfaceElementUtilities::CalculateVoigtVector(rVariables.VoigtVector);

    //Variables computed at each GP
    //Constitutive Law parameters
    rVariables.StrainVector.resize(TDim,false);
    rVariables.StressVector.resize(TDim,false);
    rVariables.ConstitutiveMatrix.resize(TDim,TDim,false);
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
    noalias(rVariables.Nu) = ZeroMatrix(TDim, TNumNodes*TDim);
    noalias(rVariables.LocalPermeabilityMatrix) = ZeroMatrix(TDim,TDim);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainInterfaceElement<2,4>::CalculateRotationMatrix(BoundedMatrix<double,2,2>& rRotationMatrix, const GeometryType& Geom)
{
    KRATOS_TRY

    //Define mid-plane points for quadrilateral_interface_2d_4
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    noalias(pmid0) = 0.5 * (Geom.GetPoint( 0 ) + Geom.GetPoint( 3 ));
    noalias(pmid1) = 0.5 * (Geom.GetPoint( 1 ) + Geom.GetPoint( 2 ));

    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = pmid1 - pmid0;
    double inv_norm_x = 1.0/norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;

    //Rotation Matrix
    rRotationMatrix(0,0) = Vx[0];
    rRotationMatrix(0,1) = Vx[1];

    // We need to determine the unitary vector in local y direction pointing towards the TOP face of the joint

    // Unitary vector in local x direction (3D)
    array_1d<double, 3> Vx3D;
    Vx3D[0] = Vx[0];
    Vx3D[1] = Vx[1];
    Vx3D[2] = 0.0;

    // Unitary vector in local y direction (first option)
    array_1d<double, 3> Vy3D;
    Vy3D[0] = -Vx[1];
    Vy3D[1] = Vx[0];
    Vy3D[2] = 0.0;

    // Vector in global z direction (first option)
    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx3D, Vy3D);

    // Vz must have the same sign as vector (0,0,1)
    if(Vz[2] > 0.0)
    {
        rRotationMatrix(1,0) = -Vx[1];
        rRotationMatrix(1,1) = Vx[0];
    }
    else
    {
        rRotationMatrix(1,0) = Vx[1];
        rRotationMatrix(1,1) = -Vx[0];
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainInterfaceElement<3,6>::CalculateRotationMatrix(BoundedMatrix<double,3,3>& rRotationMatrix, const GeometryType& Geom)
{
    KRATOS_TRY

    //Define mid-plane points for prism_interface_3d_6
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    array_1d<double, 3> pmid2;
    noalias(pmid0) = 0.5 * (Geom.GetPoint( 0 ) + Geom.GetPoint( 3 ));
    noalias(pmid1) = 0.5 * (Geom.GetPoint( 1 ) + Geom.GetPoint( 4 ));
    noalias(pmid2) = 0.5 * (Geom.GetPoint( 2 ) + Geom.GetPoint( 5 ));

    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = pmid1 - pmid0;
    double inv_norm = 1.0/norm_2(Vx);
    Vx[0] *= inv_norm;
    Vx[1] *= inv_norm;
    Vx[2] *= inv_norm;

    //Unitary vector in local z direction
    array_1d<double, 3> Vy;
    noalias(Vy) = pmid2 - pmid0;
    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    inv_norm = 1.0/norm_2(Vz);
    Vz[0] *= inv_norm;
    Vz[1] *= inv_norm;
    Vz[2] *= inv_norm;

    //Unitary vector in local y direction
    MathUtils<double>::CrossProduct( Vy, Vz, Vx);

    //Rotation Matrix
    rRotationMatrix(0,0) = Vx[0];
    rRotationMatrix(0,1) = Vx[1];
    rRotationMatrix(0,2) = Vx[2];

    rRotationMatrix(1,0) = Vy[0];
    rRotationMatrix(1,1) = Vy[1];
    rRotationMatrix(1,2) = Vy[2];

    rRotationMatrix(2,0) = Vz[0];
    rRotationMatrix(2,1) = Vz[1];
    rRotationMatrix(2,2) = Vz[2];

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainInterfaceElement<3,8>::CalculateRotationMatrix(BoundedMatrix<double,3,3>& rRotationMatrix, const GeometryType& Geom)
{
    KRATOS_TRY

    //Define mid-plane points for hexahedra_interface_3d_8
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    array_1d<double, 3> pmid2;
    noalias(pmid0) = 0.5 * (Geom.GetPoint( 0 ) + Geom.GetPoint( 4 ));
    noalias(pmid1) = 0.5 * (Geom.GetPoint( 1 ) + Geom.GetPoint( 5 ));
    noalias(pmid2) = 0.5 * (Geom.GetPoint( 2 ) + Geom.GetPoint( 6 ));

    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = pmid1 - pmid0;
    double inv_norm = 1.0/norm_2(Vx);
    Vx[0] *= inv_norm;
    Vx[1] *= inv_norm;
    Vx[2] *= inv_norm;

    //Unitary vector in local z direction
    array_1d<double, 3> Vy;
    noalias(Vy) = pmid2 - pmid0;
    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    inv_norm = 1.0/norm_2(Vz);
    Vz[0] *= inv_norm;
    Vz[1] *= inv_norm;
    Vz[2] *= inv_norm;

    //Unitary vector in local y direction
    MathUtils<double>::CrossProduct( Vy, Vz, Vx);

    //Rotation Matrix
    rRotationMatrix(0,0) = Vx[0];
    rRotationMatrix(0,1) = Vx[1];
    rRotationMatrix(0,2) = Vx[2];

    rRotationMatrix(1,0) = Vy[0];
    rRotationMatrix(1,1) = Vy[1];
    rRotationMatrix(1,2) = Vy[2];

    rRotationMatrix(2,0) = Vz[0];
    rRotationMatrix(2,1) = Vz[1];
    rRotationMatrix(2,2) = Vz[2];

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateJointWidth(double& rJointWidth, const double& NormalRelDisp,
                                                                        const double& MinimumJointWidth,const unsigned int& GPoint)
{
    rJointWidth = mInitialGap[GPoint] + NormalRelDisp;

    if(rJointWidth < MinimumJointWidth)
    {
        rJointWidth = MinimumJointWidth;
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CheckAndCalculateJointWidth(double& rJointWidth, ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                                                                double& rNormalRelDisp,const double& MinimumJointWidth,const unsigned int& GPoint)
{
    rJointWidth = mInitialGap[GPoint] + rNormalRelDisp;

    rConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY); // No contact between interfaces

    // Initally open joint
    if(mIsOpen[GPoint]==true)
    {
        if(rJointWidth < MinimumJointWidth)
        {
            rConstitutiveParameters.Reset(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY); // Contact between interfaces
            rNormalRelDisp = rJointWidth - MinimumJointWidth;
            rJointWidth = MinimumJointWidth;
        }
    }
    // Initally closed joint
    else
    {
        if(rJointWidth < 0.0)
        {
            rConstitutiveParameters.Reset(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY); // Contact between interfaces
            rNormalRelDisp = rJointWidth;
            rJointWidth = MinimumJointWidth;
        }
        else if(rJointWidth < MinimumJointWidth)
        {
            rJointWidth = MinimumJointWidth;
        }
    }
}

//----------------------------------------------------------------------------------------

template< >
template< class TMatrixType >
void UPwSmallStrainInterfaceElement<2,4>::CalculateShapeFunctionsGradients(TMatrixType& rGradNpT,SFGradAuxVariables& rAuxVariables,
                                                    const Matrix& Jacobian,const BoundedMatrix<double,2,2>& RotationMatrix,
                                                    const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint)
{
    //Quadrilateral_interface_2d_4
    rAuxVariables.GlobalCoordinatesGradients[0] = Jacobian(0,0);
    rAuxVariables.GlobalCoordinatesGradients[1] = Jacobian(1,0);
    noalias(rAuxVariables.LocalCoordinatesGradients) = prod(RotationMatrix,rAuxVariables.GlobalCoordinatesGradients);

    rGradNpT(0,0) = DN_De(0,0)/rAuxVariables.LocalCoordinatesGradients[0]; rGradNpT(0,1) = -Ncontainer(GPoint,0)/JointWidth;
    rGradNpT(1,0) = DN_De(1,0)/rAuxVariables.LocalCoordinatesGradients[0]; rGradNpT(1,1) = -Ncontainer(GPoint,1)/JointWidth;
    rGradNpT(2,0) = DN_De(2,0)/rAuxVariables.LocalCoordinatesGradients[0]; rGradNpT(2,1) = Ncontainer(GPoint,2)/JointWidth;
    rGradNpT(3,0) = DN_De(3,0)/rAuxVariables.LocalCoordinatesGradients[0]; rGradNpT(3,1) = Ncontainer(GPoint,3)/JointWidth;
}

//----------------------------------------------------------------------------------------

template< >
template< class TMatrixType >
void UPwSmallStrainInterfaceElement<3,6>::CalculateShapeFunctionsGradients(TMatrixType& rGradNpT,SFGradAuxVariables& rAuxVariables,
                                                    const Matrix& Jacobian,const BoundedMatrix<double,3,3>& RotationMatrix,
                                                    const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint)
{
    //Prism_interface_3d_6
    for(unsigned int i = 0; i < 6; i++)
    {
        rAuxVariables.ShapeFunctionsNaturalGradientsMatrix(i,0) = DN_De(i,0); rAuxVariables.ShapeFunctionsNaturalGradientsMatrix(i,1) = DN_De(i,1);
    }

    rAuxVariables.GlobalCoordinatesGradients[0] = Jacobian(0,0);
    rAuxVariables.GlobalCoordinatesGradients[1] = Jacobian(1,0);
    rAuxVariables.GlobalCoordinatesGradients[2] = Jacobian(2,0);
    noalias(rAuxVariables.LocalCoordinatesGradients) = prod(RotationMatrix,rAuxVariables.GlobalCoordinatesGradients);

    rAuxVariables.LocalCoordinatesGradientsMatrix(0,0) = rAuxVariables.LocalCoordinatesGradients[0];
    rAuxVariables.LocalCoordinatesGradientsMatrix(1,0) = rAuxVariables.LocalCoordinatesGradients[1];

    rAuxVariables.GlobalCoordinatesGradients[0] = Jacobian(0,1);
    rAuxVariables.GlobalCoordinatesGradients[1] = Jacobian(1,1);
    rAuxVariables.GlobalCoordinatesGradients[2] = Jacobian(2,1);
    noalias(rAuxVariables.LocalCoordinatesGradients) = prod(RotationMatrix,rAuxVariables.GlobalCoordinatesGradients);

    rAuxVariables.LocalCoordinatesGradientsMatrix(0,1) = rAuxVariables.LocalCoordinatesGradients[0];
    rAuxVariables.LocalCoordinatesGradientsMatrix(1,1) = rAuxVariables.LocalCoordinatesGradients[1];

    PoroElementUtilities::InvertMatrix2( rAuxVariables.LocalCoordinatesGradientsInvMatrix, rAuxVariables.LocalCoordinatesGradientsMatrix );

    noalias(rAuxVariables.ShapeFunctionsGradientsMatrix) = prod(rAuxVariables.ShapeFunctionsNaturalGradientsMatrix,rAuxVariables.LocalCoordinatesGradientsInvMatrix);

    rGradNpT(0,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(0,0); rGradNpT(0,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(0,1); rGradNpT(0,2) = -Ncontainer(GPoint,0)/JointWidth;
    rGradNpT(1,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(1,0); rGradNpT(1,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(1,1); rGradNpT(1,2) = -Ncontainer(GPoint,1)/JointWidth;
    rGradNpT(2,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(2,0); rGradNpT(2,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(2,1); rGradNpT(2,2) = -Ncontainer(GPoint,2)/JointWidth;
    rGradNpT(3,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(3,0); rGradNpT(3,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(3,1); rGradNpT(3,2) = Ncontainer(GPoint,3)/JointWidth;
    rGradNpT(4,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(4,0); rGradNpT(4,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(4,1); rGradNpT(4,2) = Ncontainer(GPoint,4)/JointWidth;
    rGradNpT(5,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(5,0); rGradNpT(5,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(5,1); rGradNpT(5,2) = Ncontainer(GPoint,5)/JointWidth;
}

//----------------------------------------------------------------------------------------

template< >
template< class TMatrixType >
void UPwSmallStrainInterfaceElement<3,8>::CalculateShapeFunctionsGradients(TMatrixType& rGradNpT,SFGradAuxVariables& rAuxVariables,
                                                    const Matrix& Jacobian,const BoundedMatrix<double,3,3>& RotationMatrix,
                                                    const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint)
{
    //Hexahedral_interface_3d_8
    for(unsigned int i = 0; i < 8; i++)
    {
        rAuxVariables.ShapeFunctionsNaturalGradientsMatrix(i,0) = DN_De(i,0); rAuxVariables.ShapeFunctionsNaturalGradientsMatrix(i,1) = DN_De(i,1);
    }

    rAuxVariables.GlobalCoordinatesGradients[0] = Jacobian(0,0);
    rAuxVariables.GlobalCoordinatesGradients[1] = Jacobian(1,0);
    rAuxVariables.GlobalCoordinatesGradients[2] = Jacobian(2,0);
    noalias(rAuxVariables.LocalCoordinatesGradients) = prod(RotationMatrix,rAuxVariables.GlobalCoordinatesGradients);

    rAuxVariables.LocalCoordinatesGradientsMatrix(0,0) = rAuxVariables.LocalCoordinatesGradients[0];
    rAuxVariables.LocalCoordinatesGradientsMatrix(1,0) = rAuxVariables.LocalCoordinatesGradients[1];

    rAuxVariables.GlobalCoordinatesGradients[0] = Jacobian(0,1);
    rAuxVariables.GlobalCoordinatesGradients[1] = Jacobian(1,1);
    rAuxVariables.GlobalCoordinatesGradients[2] = Jacobian(2,1);
    noalias(rAuxVariables.LocalCoordinatesGradients) = prod(RotationMatrix,rAuxVariables.GlobalCoordinatesGradients);

    rAuxVariables.LocalCoordinatesGradientsMatrix(0,1) = rAuxVariables.LocalCoordinatesGradients[0];
    rAuxVariables.LocalCoordinatesGradientsMatrix(1,1) = rAuxVariables.LocalCoordinatesGradients[1];

    PoroElementUtilities::InvertMatrix2( rAuxVariables.LocalCoordinatesGradientsInvMatrix, rAuxVariables.LocalCoordinatesGradientsMatrix );

    noalias(rAuxVariables.ShapeFunctionsGradientsMatrix) = prod(rAuxVariables.ShapeFunctionsNaturalGradientsMatrix,rAuxVariables.LocalCoordinatesGradientsInvMatrix);

    rGradNpT(0,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(0,0); rGradNpT(0,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(0,1); rGradNpT(0,2) = -Ncontainer(GPoint,0)/JointWidth;
    rGradNpT(1,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(1,0); rGradNpT(1,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(1,1); rGradNpT(1,2) = -Ncontainer(GPoint,1)/JointWidth;
    rGradNpT(2,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(2,0); rGradNpT(2,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(2,1); rGradNpT(2,2) = -Ncontainer(GPoint,2)/JointWidth;
    rGradNpT(3,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(3,0); rGradNpT(3,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(3,1); rGradNpT(3,2) = -Ncontainer(GPoint,3)/JointWidth;
    rGradNpT(4,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(4,0); rGradNpT(4,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(4,1); rGradNpT(4,2) = Ncontainer(GPoint,4)/JointWidth;
    rGradNpT(5,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(5,0); rGradNpT(5,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(5,1); rGradNpT(5,2) = Ncontainer(GPoint,5)/JointWidth;
    rGradNpT(6,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(6,0); rGradNpT(6,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(6,1); rGradNpT(6,2) = Ncontainer(GPoint,6)/JointWidth;
    rGradNpT(7,0) = rAuxVariables.ShapeFunctionsGradientsMatrix(7,0); rGradNpT(7,1) = rAuxVariables.ShapeFunctionsGradientsMatrix(7,1); rGradNpT(7,2) = Ncontainer(GPoint,7)/JointWidth;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables)
{
    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables)
{
    noalias(rVariables.DimMatrix) = prod(trans(rVariables.RotationMatrix),
                                        BoundedMatrix<double,TDim,TDim>(prod(rVariables.ConstitutiveMatrix,
                                        rVariables.RotationMatrix)));
    noalias(rVariables.UDimMatrix) = prod(trans(rVariables.Nu),rVariables.DimMatrix);
    noalias(rVariables.UMatrix) = prod(rVariables.UDimMatrix,rVariables.Nu)*rVariables.IntegrationCoefficient;

    //Distribute stiffness block matrix into the elemental matrix
    PoroElementUtilities::AssembleUBlockMatrix(rLeftHandSideMatrix,rVariables.UMatrix);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables)
{
    noalias(rVariables.UDimMatrix) = prod(trans(rVariables.Nu),trans(rVariables.RotationMatrix));

    noalias(rVariables.UVector) = prod(rVariables.UDimMatrix,rVariables.VoigtVector);

    noalias(rVariables.UPMatrix) = -rVariables.BiotCoefficient*outer_prod(rVariables.UVector,rVariables.Np)*rVariables.IntegrationCoefficient;

    //Distribute coupling block matrix into the elemental matrix
    PoroElementUtilities::AssembleUPBlockMatrix(rLeftHandSideMatrix,rVariables.UPMatrix);

    noalias(rVariables.PUMatrix) = -rVariables.VelocityCoefficient*trans(rVariables.UPMatrix);

    //Distribute transposed coupling block matrix into the elemental matrix
    PoroElementUtilities::AssemblePUBlockMatrix(rLeftHandSideMatrix,rVariables.PUMatrix);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables)
{
    noalias(rVariables.PMatrix) = rVariables.DtPressureCoefficient*rVariables.BiotModulusInverse*outer_prod(rVariables.Np,rVariables.Np)*rVariables.JointWidth*rVariables.IntegrationCoefficient;

    //Distribute compressibility block matrix into the elemental matrix
    PoroElementUtilities::AssemblePBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables)
{
    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT,rVariables.LocalPermeabilityMatrix);

    noalias(rVariables.PMatrix) = rVariables.DynamicViscosityInverse*prod(rVariables.PDimMatrix,trans(rVariables.GradNpT))*rVariables.JointWidth*rVariables.IntegrationCoefficient;

    //Distribute permeability block matrix into the elemental matrix
    PoroElementUtilities::AssemblePBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rVariables.PMatrix,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables)
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
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables)
{
    noalias(rVariables.UDimMatrix) = prod(trans(rVariables.Nu),trans(rVariables.RotationMatrix));

    noalias(rVariables.UVector) = -1.0*prod(rVariables.UDimMatrix,rVariables.StressVector)*rVariables.IntegrationCoefficient;

    //Distribute stiffness block vector into elemental vector
    PoroElementUtilities::AssembleUBlockVector(rRightHandSideVector,rVariables.UVector);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables)
{
    noalias(rVariables.UVector) = rVariables.Density*prod(trans(rVariables.Nu),rVariables.BodyAcceleration)*rVariables.JointWidth*rVariables.IntegrationCoefficient;

    //Distribute body force block vector into elemental vector
    PoroElementUtilities::AssembleUBlockVector(rRightHandSideVector,rVariables.UVector);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables)
{
    noalias(rVariables.UDimMatrix) = prod(trans(rVariables.Nu),trans(rVariables.RotationMatrix));

    noalias(rVariables.UVector) = prod(rVariables.UDimMatrix,rVariables.VoigtVector);

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
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables)
{
    noalias(rVariables.PMatrix) = rVariables.BiotModulusInverse*outer_prod(rVariables.Np,rVariables.Np)*rVariables.JointWidth*rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = -1.0*prod(rVariables.PMatrix,rVariables.DtPressureVector);

    //Distribute compressibility block vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables)
{
    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT,rVariables.LocalPermeabilityMatrix);

    noalias(rVariables.PMatrix) = rVariables.DynamicViscosityInverse*prod(rVariables.PDimMatrix,trans(rVariables.GradNpT))*rVariables.JointWidth*rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = -1.0*prod(rVariables.PMatrix,rVariables.PressureVector);

    //Distribute permeability block vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainInterfaceElement<TDim,TNumNodes>::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables)
{
    noalias(rVariables.PDimMatrix) = prod(rVariables.GradNpT,rVariables.LocalPermeabilityMatrix)*rVariables.JointWidth*rVariables.IntegrationCoefficient;

    noalias(rVariables.PVector) = rVariables.DynamicViscosityInverse*rVariables.FluidDensity*
                                    prod(rVariables.PDimMatrix,rVariables.BodyAcceleration);

    //Distribute fluid body flow block vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< >
void UPwSmallStrainInterfaceElement<2,4>::InterpolateOutputDoubles( std::vector<double>& rOutput, const std::vector<double>& GPValues )
{
    //Interpolation of computed values at Lobatto GP to the standard GiD gauss points

    rOutput[0] = 0.6220084679281462 * GPValues[0] + 0.16666666666666663 * GPValues[1] + 0.044658198738520435 * GPValues[1] + 0.16666666666666663 * GPValues[0];

    rOutput[1] = 0.16666666666666663 * GPValues[0] + 0.6220084679281462 * GPValues[1] + 0.16666666666666663 * GPValues[1] + 0.044658198738520435 * GPValues[0];

    rOutput[2]= 0.044658198738520435 * GPValues[0] + 0.16666666666666663 * GPValues[1] + 0.6220084679281462 * GPValues[1] + 0.16666666666666663 * GPValues[0];

    rOutput[3] = 0.16666666666666663 * GPValues[0] + 0.044658198738520435 * GPValues[1] + 0.16666666666666663 * GPValues[1] + 0.6220084679281462 * GPValues[0];
}

//----------------------------------------------------------------------------------------

template< >
void UPwSmallStrainInterfaceElement<3,6>::InterpolateOutputDoubles( std::vector<double>& rOutput, const std::vector<double>& GPValues )
{
    //Interpolation of computed values at Lobatto GP to the standard GiD gauss points

    rOutput[0] = 0.5257834230632086 * GPValues[0] + 0.13144585576580214 * GPValues[1] + 0.13144585576580214 * GPValues[2]
                + 0.14088324360345805 * GPValues[0] + 0.03522081090086451 * GPValues[1] + 0.03522081090086451 * GPValues[2];

    rOutput[1] = 0.13144585576580214 * GPValues[0] + 0.5257834230632086 * GPValues[1] + 0.13144585576580214 * GPValues[2]
                + 0.03522081090086451 * GPValues[0] + 0.14088324360345805 * GPValues[1] + 0.03522081090086451 * GPValues[2];

    rOutput[2] = 0.13144585576580214 * GPValues[0] + 0.13144585576580214 * GPValues[1] + 0.5257834230632086 * GPValues[2]
                + 0.03522081090086451 * GPValues[0] + 0.03522081090086451 * GPValues[1] + 0.14088324360345805 * GPValues[2];

    rOutput[3] = 0.14088324360345805 * GPValues[0] + 0.03522081090086451 * GPValues[1] + 0.03522081090086451 * GPValues[2]
                + 0.5257834230632086 * GPValues[0] + 0.13144585576580214 * GPValues[1] + 0.13144585576580214 * GPValues[2];

    rOutput[4] = 0.03522081090086451 * GPValues[0] + 0.14088324360345805 * GPValues[1] + 0.03522081090086451 * GPValues[2]
                + 0.13144585576580214 * GPValues[0] + 0.5257834230632086 * GPValues[1] + 0.13144585576580214 * GPValues[2];

    rOutput[5] = 0.03522081090086451 * GPValues[0] + 0.03522081090086451 * GPValues[1] + 0.14088324360345805 * GPValues[2]
                + 0.13144585576580214 * GPValues[0] + 0.13144585576580214 * GPValues[1] + 0.5257834230632086 * GPValues[2];
}

//----------------------------------------------------------------------------------------

template<>
void UPwSmallStrainInterfaceElement<3,8>::InterpolateOutputDoubles( std::vector<double>& rOutput, const std::vector<double>& GPValues )
{
    //Interpolation of computed values at Lobatto GP to the standard GiD gauss points

    rOutput[0] = 0.4905626121623441 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3]
                + 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.009437387837655926 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    rOutput[1] = 0.13144585576580212 * GPValues[0] + 0.4905626121623441 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3]
                + 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.009437387837655926 * GPValues[3];

    rOutput[2] = 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.4905626121623441 * GPValues[2] + 0.13144585576580212 * GPValues[3]
                + 0.009437387837655926 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    rOutput[3] = 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.4905626121623441 * GPValues[3]
                + 0.035220810900864506 * GPValues[0] + 0.009437387837655926 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    rOutput[4] = 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.009437387837655926 * GPValues[2] + 0.035220810900864506 * GPValues[3]
                + 0.4905626121623441 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    rOutput[5] = 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.009437387837655926 * GPValues[3]
                + 0.13144585576580212 * GPValues[0] + 0.4905626121623441 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    rOutput[6] = 0.009437387837655926 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3]
                + 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.4905626121623441 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    rOutput[7] = 0.035220810900864506 * GPValues[0] + 0.009437387837655926 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3]
                + 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.4905626121623441 * GPValues[3];
}

//----------------------------------------------------------------------------------------

template<>
template< class TValueType >
void UPwSmallStrainInterfaceElement<2,4>::InterpolateOutputValues( std::vector<TValueType>& rOutput, const std::vector<TValueType>& GPValues )
{
    //Interpolation of computed values at Lobatto GP to the standard GiD gauss points

    noalias(rOutput[0]) = 0.6220084679281462 * GPValues[0] + 0.16666666666666663 * GPValues[1] + 0.044658198738520435 * GPValues[1] + 0.16666666666666663 * GPValues[0];

    noalias(rOutput[1]) = 0.16666666666666663 * GPValues[0] + 0.6220084679281462 * GPValues[1] + 0.16666666666666663 * GPValues[1] + 0.044658198738520435 * GPValues[0];

    noalias(rOutput[2])= 0.044658198738520435 * GPValues[0] + 0.16666666666666663 * GPValues[1] + 0.6220084679281462 * GPValues[1] + 0.16666666666666663 * GPValues[0];

    noalias(rOutput[3]) = 0.16666666666666663 * GPValues[0] + 0.044658198738520435 * GPValues[1] + 0.16666666666666663 * GPValues[1] + 0.6220084679281462 * GPValues[0];
}

//----------------------------------------------------------------------------------------

template<>
template< class TValueType >
void UPwSmallStrainInterfaceElement<3,6>::InterpolateOutputValues( std::vector<TValueType>& rOutput, const std::vector<TValueType>& GPValues )
{
    //Interpolation of computed values at Lobatto GP to the standard GiD gauss points

    noalias(rOutput[0]) = 0.5257834230632086 * GPValues[0] + 0.13144585576580214 * GPValues[1] + 0.13144585576580214 * GPValues[2]
                        + 0.14088324360345805 * GPValues[0] + 0.03522081090086451 * GPValues[1] + 0.03522081090086451 * GPValues[2];

    noalias(rOutput[1]) = 0.13144585576580214 * GPValues[0] + 0.5257834230632086 * GPValues[1] + 0.13144585576580214 * GPValues[2]
                        + 0.03522081090086451 * GPValues[0] + 0.14088324360345805 * GPValues[1] + 0.03522081090086451 * GPValues[2];

    noalias(rOutput[2]) = 0.13144585576580214 * GPValues[0] + 0.13144585576580214 * GPValues[1] + 0.5257834230632086 * GPValues[2]
                        + 0.03522081090086451 * GPValues[0] + 0.03522081090086451 * GPValues[1] + 0.14088324360345805 * GPValues[2];

    noalias(rOutput[3]) = 0.14088324360345805 * GPValues[0] + 0.03522081090086451 * GPValues[1] + 0.03522081090086451 * GPValues[2]
                        + 0.5257834230632086 * GPValues[0] + 0.13144585576580214 * GPValues[1] + 0.13144585576580214 * GPValues[2];

    noalias(rOutput[4]) = 0.03522081090086451 * GPValues[0] + 0.14088324360345805 * GPValues[1] + 0.03522081090086451 * GPValues[2]
                        + 0.13144585576580214 * GPValues[0] + 0.5257834230632086 * GPValues[1] + 0.13144585576580214 * GPValues[2];

    noalias(rOutput[5]) = 0.03522081090086451 * GPValues[0] + 0.03522081090086451 * GPValues[1] + 0.14088324360345805 * GPValues[2]
                        + 0.13144585576580214 * GPValues[0] + 0.13144585576580214 * GPValues[1] + 0.5257834230632086 * GPValues[2];
}

//----------------------------------------------------------------------------------------

template<>
template< class TValueType >
void UPwSmallStrainInterfaceElement<3,8>::InterpolateOutputValues( std::vector<TValueType>& rOutput, const std::vector<TValueType>& GPValues )
{
    //Interpolation of computed values at Lobatto GP to the standard GiD gauss points

    noalias(rOutput[0]) = 0.4905626121623441 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3]
                        + 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.009437387837655926 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    noalias(rOutput[1]) = 0.13144585576580212 * GPValues[0] + 0.4905626121623441 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3]
                        + 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.009437387837655926 * GPValues[3];

    noalias(rOutput[2]) = 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.4905626121623441 * GPValues[2] + 0.13144585576580212 * GPValues[3]
                        + 0.009437387837655926 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    noalias(rOutput[3]) = 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.4905626121623441 * GPValues[3]
                        + 0.035220810900864506 * GPValues[0] + 0.009437387837655926 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    noalias(rOutput[4]) = 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.009437387837655926 * GPValues[2] + 0.035220810900864506 * GPValues[3]
                        + 0.4905626121623441 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    noalias(rOutput[5]) = 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.009437387837655926 * GPValues[3]
                        + 0.13144585576580212 * GPValues[0] + 0.4905626121623441 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3];

    noalias(rOutput[6]) = 0.009437387837655926 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.035220810900864506 * GPValues[3]
                        + 0.035220810900864506 * GPValues[0] + 0.13144585576580212 * GPValues[1] + 0.4905626121623441 * GPValues[2] + 0.13144585576580212 * GPValues[3];

    noalias(rOutput[7]) = 0.035220810900864506 * GPValues[0] + 0.009437387837655926 * GPValues[1] + 0.035220810900864506 * GPValues[2] + 0.13144585576580212 * GPValues[3]
                        + 0.13144585576580212 * GPValues[0] + 0.035220810900864506 * GPValues[1] + 0.13144585576580212 * GPValues[2] + 0.4905626121623441 * GPValues[3];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template void UPwSmallStrainInterfaceElement<2,4>::CalculateShapeFunctionsGradients< BoundedMatrix<double,4,2> >
                                                    ( BoundedMatrix<double,4,2>& rGradNpT, SFGradAuxVariables& rAuxVariables,
                                                    const Matrix& Jacobian,const BoundedMatrix<double,2,2>& RotationMatrix,
                                                    const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint );
template void UPwSmallStrainInterfaceElement<2,4>::CalculateShapeFunctionsGradients< Matrix >
                                                    ( Matrix& rGradNpT, SFGradAuxVariables& rAuxVariables,const Matrix& Jacobian,
                                                    const BoundedMatrix<double,2,2>& RotationMatrix,
                                                    const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint );
template void UPwSmallStrainInterfaceElement<3,6>::CalculateShapeFunctionsGradients< BoundedMatrix<double,6,3> >
                                                    ( BoundedMatrix<double,6,3>& rGradNpT, SFGradAuxVariables& rAuxVariables,
                                                    const Matrix& Jacobian,const BoundedMatrix<double,3,3>& RotationMatrix,
                                                    const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint );
template void UPwSmallStrainInterfaceElement<3,6>::CalculateShapeFunctionsGradients< Matrix >
                                                    ( Matrix& rGradNpT, SFGradAuxVariables& rAuxVariables,const Matrix& Jacobian,
                                                    const BoundedMatrix<double,3,3>& RotationMatrix,
                                                    const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint );
template void UPwSmallStrainInterfaceElement<3,8>::CalculateShapeFunctionsGradients< BoundedMatrix<double,8,3> >
                                                    ( BoundedMatrix<double,8,3>& rGradNpT, SFGradAuxVariables& rAuxVariables,
                                                    const Matrix& Jacobian,const BoundedMatrix<double,3,3>& RotationMatrix,
                                                    const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint );
template void UPwSmallStrainInterfaceElement<3,8>::CalculateShapeFunctionsGradients< Matrix >
                                                    ( Matrix& rGradNpT, SFGradAuxVariables& rAuxVariables,const Matrix& Jacobian,
                                                    const BoundedMatrix<double,3,3>& RotationMatrix,
                                                    const Matrix& DN_De,const Matrix& Ncontainer, const double& JointWidth,const unsigned int& GPoint );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwSmallStrainInterfaceElement<2,4>;
template class UPwSmallStrainInterfaceElement<3,6>;
template class UPwSmallStrainInterfaceElement<3,8>;

} // Namespace Kratos
