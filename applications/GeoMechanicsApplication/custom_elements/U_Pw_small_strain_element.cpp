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
Element::Pointer UPwSmallStrainElement<TDim,TNumNodes>::
    Create(IndexType NewId,
           NodesArrayType const& ThisNodes,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new UPwSmallStrainElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer UPwSmallStrainElement<TDim,TNumNodes>::
    Create(IndexType NewId,
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
    // KRATOS_INFO("0-UPwSmallStrainElement::Check()") << this->Id() << std::endl;

    // Base class checks for positive area and Id > 0
    // Verify generic variables
    int ierr = UPwBaseElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();

    if (rGeom.DomainSize() < 1.0e-15)
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

    // Verify specific properties
    if ( rProp.Has( IGNORE_UNDRAINED ) == false)
        KRATOS_ERROR << "IGNORE_UNDRAINED does not exist in the parameter list" << this->Id() << std::endl;

    bool IgnoreUndrained = rProp[IGNORE_UNDRAINED];

    if (!IgnoreUndrained) {
        if ( rProp.Has( BULK_MODULUS_FLUID ) == false || rProp[BULK_MODULUS_FLUID] < 0.0 )
            KRATOS_ERROR << "BULK_MODULUS_FLUID has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

        if ( rProp.Has( DYNAMIC_VISCOSITY ) == false || rProp[DYNAMIC_VISCOSITY] < 0.0 )
            KRATOS_ERROR << "DYNAMIC_VISCOSITY has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

        if ( rProp.Has( PERMEABILITY_XX ) == false || rProp[PERMEABILITY_XX] < 0.0 )
            KRATOS_ERROR << "PERMEABILITY_XX has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

        if ( rProp.Has( PERMEABILITY_YY ) == false || rProp[PERMEABILITY_YY] < 0.0 )
            KRATOS_ERROR << "PERMEABILITY_YY has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

        if ( rProp.Has( PERMEABILITY_XY ) == false || rProp[PERMEABILITY_XY] < 0.0 )
            KRATOS_ERROR << "PERMEABILITY_XY has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

        if (TDim > 2) {
            if ( rProp.Has( PERMEABILITY_ZZ ) == false || rProp[PERMEABILITY_ZZ] < 0.0 )
                KRATOS_ERROR << "PERMEABILITY_ZZ has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

            if ( rProp.Has( PERMEABILITY_YZ ) == false || rProp[PERMEABILITY_YZ] < 0.0 )
                KRATOS_ERROR << "PERMEABILITY_YZ has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

            if ( rProp.Has( PERMEABILITY_ZX ) == false || rProp[PERMEABILITY_ZX] < 0.0 )
                KRATOS_ERROR << "PERMEABILITY_ZX has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;
        }
    }


    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW )) 
        << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strainSize = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    if ( TDim == 2 ) {
        KRATOS_ERROR_IF_NOT( strainSize == VOIGT_SIZE_2D_PLANE_STRAIN )
        << "Wrong constitutive law used. This is a 2D element! expected strain size is "
        << VOIGT_SIZE_2D_PLANE_STRAIN
        << " But received: "
        << strainSize
        << " in element id: "
        << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT( strainSize == VOIGT_SIZE_3D )
        << "Wrong constitutive law used. This is a 3D element! expected strain size is "
        << VOIGT_SIZE_3D
        << " But received: "
        << strainSize
        << " in element id: "
        << this->Id() << std::endl;
    }

    // Check constitutive law
    if ( mConstitutiveLawVector.size() > 0 ) {
        return mConstitutiveLawVector[0]->Check( rProp, rGeom, rCurrentProcessInfo );
    }

    // Check constitutive law
    if ( mRetentionLawVector.size() > 0 ) {
        return mRetentionLawVector[0]->Check( rProp, rCurrentProcessInfo );
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::Check()") << std::endl;

    return ierr;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeSolutionStep()") << this->Id() << std::endl;

    if (!mIsInitialised) this->Initialize(rCurrentProcessInfo);

    //Defining necessary variables
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); //Note: this is for nonlocal damage

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,
                                     rCurrentProcessInfo);

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        //Compute Np, GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Initialize constitutive law
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->InitializeMaterialResponseCauchy(ConstitutiveParameters);

        // Initialize retention law
        mRetentionLawVector[GPoint]->InitializeSolutionStep(RetentionParameters);
    }

    // reset hydraulic discharge
    this->ResetHydraulicDischarge();

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeSolutionStep()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    ResetHydraulicDischarge()
{
    KRATOS_TRY;

    // reset hydraulic discharge
    GeometryType& rGeom = this->GetGeometry();
    for (unsigned int i=0; i<TNumNodes; ++i) {
        ThreadSafeNodeWrite(rGeom[i], HYDRAULIC_DISCHARGE, 0.0 );
    }

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateHydraulicDischarge(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateHydraulicDischarge()") << std::endl;

    std::vector<array_1d<double,3>> FluidFlux;
    this->CalculateOnIntegrationPoints( FLUID_FLUX_VECTOR,
                                        FluidFlux,
                                        rCurrentProcessInfo );

    const GeometryType& rGeom   = this->GetGeometry();
    const unsigned int NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( mThisIntegrationMethod );

    ElementVariables Variables;
    // gradient of shape functions and determinant of Jacobian

    Variables.GradNpTInitialConfiguration.resize(TNumNodes, TDim, false);
    Variables.GradNpT.resize(TNumNodes, TDim, false);
    (Variables.detJContainer).resize(NumGPoints, false);
    rGeom.ShapeFunctionsIntegrationPointsGradients( Variables.DN_DXContainer,
                                                   Variables.detJContainer,
                                                   mThisIntegrationMethod );

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint ) {
        noalias(Variables.GradNpT) = Variables.DN_DXContainer[GPoint];
        Variables.detJ = Variables.detJContainer[GPoint];

        //Compute weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints,
                                                  GPoint,
                                                  Variables.detJ);

        for (unsigned int node = 0; node < TNumNodes; ++node) {
            double HydraulicDischarge = 0;
            for (unsigned int iDir = 0; iDir < TDim; ++iDir) {
                HydraulicDischarge += Variables.GradNpT(node, iDir) * FluidFlux[GPoint][iDir];
            }

            HydraulicDischarge *= Variables.IntegrationCoefficient;
            HydraulicDischarge += rGeom[node].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE);
            ThreadSafeNodeWrite(this->GetGeometry()[node], HYDRAULIC_DISCHARGE, HydraulicDischarge);
        }
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateHydraulicDischarge()") << std::endl;

    KRATOS_CATCH("");
}



//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeNonLinearIteration()") << std::endl;

    //Defining necessary variables
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,this->GetProperties(),rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); //Note: this is for nonlocal damage

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,
                                     rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        this->CalculateKinematics(Variables, GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        // Compute Stress
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeNonLinearIteration()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------

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

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::FinalizeSolutionStep()") << std::endl;

    this->CalculateHydraulicDischarge(rCurrentProcessInfo);

    //Defining necessary variables
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,this->GetProperties(),rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,
                                     rCurrentProcessInfo);

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);

    if (rCurrentProcessInfo[NODAL_SMOOTHING] == true) {

        Matrix StressContainer(NumGPoints, mStressVector[0].size());

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint ) {
            this->CalculateKinematics(Variables, GPoint);

            //Compute infinitessimal strain
            this->CalculateStrain(Variables, GPoint);

            //set gauss points variables to constitutivelaw parameters
            this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

            //compute constitutive tensor and/or stresses
            noalias(Variables.StressVector) = mStressVector[GPoint];
            ConstitutiveParameters.SetStressVector(Variables.StressVector);
            mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
            mStateVariablesFinalized[GPoint] = 
                mConstitutiveLawVector[GPoint]->GetValue( STATE_VARIABLES,
                                                          mStateVariablesFinalized[GPoint] );

            // retention law
            mRetentionLawVector[GPoint]->FinalizeSolutionStep(RetentionParameters);

            this->SaveGPStress(StressContainer, mStressVector[GPoint], GPoint);
        }
        this->ExtrapolateGPValues(StressContainer);
    } else {
        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint ) {
            //Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            //Compute infinitessimal strain
            this->CalculateStrain(Variables, GPoint);

            //set gauss points variables to constitutivelaw parameters
            this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

            //compute constitutive tensor and/or stresses
            noalias(Variables.StressVector) = mStressVector[GPoint];
            ConstitutiveParameters.SetStressVector(Variables.StressVector);
            mConstitutiveLawVector[GPoint]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
            mStateVariablesFinalized[GPoint] = 
                mConstitutiveLawVector[GPoint]->GetValue( STATE_VARIABLES,
                                                          mStateVariablesFinalized[GPoint] );

            // retention law
            mRetentionLawVector[GPoint]->FinalizeSolutionStep(RetentionParameters);
        }
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::FinalizeSolutionStep()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::SaveGPStress(Matrix& rStressContainer,
                                                         const Vector& StressVector,
                                                         const unsigned int& GPoint)
{

    KRATOS_TRY

    // KRATOS_INFO("0-UPwSmallStrainElement::SaveGPStress()") << std::endl;

    for (unsigned int i = 0; i < StressVector.size(); ++i) {
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
    ExtrapolateGPValues(const Matrix& StressContainer)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::ExtrapolateGPValues()") << std::endl;

    array_1d<double, TNumNodes> DamageContainer;

    for ( unsigned int Node = 0;  Node < TNumNodes; ++Node ) {
        DamageContainer[Node] = 0.0;
        DamageContainer[Node] = mConstitutiveLawVector[Node]->GetValue( DAMAGE_VARIABLE, DamageContainer[Node] );
    }

    GeometryType& rGeom = this->GetGeometry();
    const double& Area = rGeom.Area(); // In 3D this is volume
    array_1d<Vector, TNumNodes> NodalStressVector; //List with stresses at each node
    array_1d<Matrix, TNumNodes> NodalStressTensor;

    for (unsigned int Node = 0; Node < TNumNodes; ++Node) {
        NodalStressVector[Node].resize(VoigtSize);
        NodalStressTensor[Node].resize(StressTensorSize, StressTensorSize);
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> ExtrapolationMatrix;
    this->CalculateExtrapolationMatrix(ExtrapolationMatrix);

    Matrix AuxNodalStress;
    AuxNodalStress.resize(TNumNodes, VoigtSize);
    noalias(AuxNodalStress) = prod(ExtrapolationMatrix, StressContainer);

    /* INFO:
        *
        *                  |S0-0 S1-0 S2-0|
        * AuxNodalStress = |S0-1 S1-1 S2-1|
        *                  |S0-2 S1-2 S2-2|
        *
        * S1-0 = S[1] at node 0
    */

    array_1d<double, TNumNodes> NodalDamage;
    noalias(NodalDamage) = prod(ExtrapolationMatrix,DamageContainer);

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
        noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);

        rGeom[i].SetLock();
        noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) += NodalDamage[i]*Area;
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
        rGeom[i].UnSetLock();
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::ExtrapolateGPValues()") << std::endl;
    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateOnIntegrationPoints( const Variable<double>& rVariable,
                                  std::vector<double>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateOnIntegrationPoints()") << std::endl;
    //Defining necessary variables
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != NumGPoints )
        rOutput.resize(NumGPoints);

    if (rVariable == VON_MISES_STRESS) {
        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint ) {
            ComparisonUtilities EquivalentStress;
            rOutput[GPoint] = EquivalentStress.CalculateVonMises(mStressVector[GPoint]);
        }
    } else if (rVariable == DEGREE_OF_SATURATION ||
               rVariable == EFFECTIVE_SATURATION ||
               rVariable == BISHOP_COEFICIENT ||
               rVariable == DERIVATIVE_OF_SATURATION ||
               rVariable == RELATIVE_PERMEABILITY ) {

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,
                                         rCurrentProcessInfo);

        // create general parametes of retention law
        RetentionLaw::Parameters RetentionParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint ) {
            //Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            Variables.FluidPressure = this->CalculateFluidPressure(Variables, GPoint);
            this->SetRetentionParameters(Variables, RetentionParameters);

            if (rVariable == DEGREE_OF_SATURATION)     rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateSaturation(RetentionParameters);
            if (rVariable == EFFECTIVE_SATURATION)     rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateEffectiveSaturation(RetentionParameters);
            if (rVariable == BISHOP_COEFICIENT)        rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(RetentionParameters);
            if (rVariable == DERIVATIVE_OF_SATURATION) rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(RetentionParameters);
            if (rVariable == RELATIVE_PERMEABILITY )   rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateRelativePermeability(RetentionParameters);
        }
    } else if (rVariable == HYDRAULIC_HEAD) {
        const double NumericalLimit = std::numeric_limits<double>::epsilon();
        const PropertiesType& rProp = this->GetProperties();

        //Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );

        //Defining necessary variables
        array_1d<double,TNumNodes> NodalHydraulicHead;
        for (unsigned int node=0; node < TNumNodes; ++node) {
            array_1d<double,3> NodeVolumeAcceleration;
            noalias(NodeVolumeAcceleration) = rGeom[node].FastGetSolutionStepValue(VOLUME_ACCELERATION, 0);
            const double g = norm_2(NodeVolumeAcceleration);
            if (g > NumericalLimit) {
                const double FluidWeight = g * rProp[DENSITY_WATER];

                array_1d<double,3> NodeCoordinates;
                noalias(NodeCoordinates) = rGeom[node].Coordinates();
                array_1d<double,3> NodeVolumeAccelerationUnitVector;
                noalias(NodeVolumeAccelerationUnitVector) = NodeVolumeAcceleration / g;

                const double WaterPressure = rGeom[node].FastGetSolutionStepValue(WATER_PRESSURE);
                NodalHydraulicHead[node] =- inner_prod(NodeCoordinates, NodeVolumeAccelerationUnitVector)
                                          - PORE_PRESSURE_SIGN_FACTOR  * WaterPressure / FluidWeight;
            } else {
                NodalHydraulicHead[node] = 0.0;
            }
        }

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint ) {
            double HydraulicHead = 0.0;
            for (unsigned int node = 0; node < TNumNodes; ++node)
                HydraulicHead += NContainer(GPoint, node) * NodalHydraulicHead[node];

            rOutput[GPoint] = HydraulicHead;
        }
    } else {
        if ( rOutput.size() != mConstitutiveLawVector.size() )
            rOutput.resize(mConstitutiveLawVector.size());

        for ( unsigned int i = 0;  i < mConstitutiveLawVector.size(); ++i ) {
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rOutput[i] );
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
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateOnIntegrationPoints<double,3>()") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );
    if ( rOutput.size() != NumGPoints )
        rOutput.resize(NumGPoints);

    if (rVariable == FLUID_FLUX_VECTOR) {
        const PropertiesType& rProp = this->GetProperties();

        //Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
        rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,mThisIntegrationMethod);

        //Defining necessary variables
        array_1d<double,TNumNodes> PressureVector;
        for (unsigned int i=0; i<TNumNodes; ++i)
            PressureVector[i] = rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        array_1d<double,TNumNodes*TDim> VolumeAcceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(VolumeAcceleration,rGeom,VOLUME_ACCELERATION);
        array_1d<double,TDim> BodyAcceleration;
        BoundedMatrix<double,TNumNodes, TDim> GradNpT;
        BoundedMatrix<double,TDim, TDim> PermeabilityMatrix;
        GeoElementUtilities::FillPermeabilityMatrix(PermeabilityMatrix,rProp);
        const double& DynamicViscosityInverse = 1.0/rProp[DYNAMIC_VISCOSITY];
        const double& FluidDensity = rProp[DENSITY_WATER];
        array_1d<double,TDim> GradPressureTerm;
        array_1d<double,TDim> FluidFlux;

        // create general parametes of retention law
        RetentionLaw::Parameters RetentionParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint ) {
            noalias(GradNpT) = DN_DXContainer[GPoint];

            GeoElementUtilities::
                InterpolateVariableWithComponents<TDim, TNumNodes>( BodyAcceleration,
                                                                    NContainer,
                                                                    VolumeAcceleration,
                                                                    GPoint );

            double FluidPressure = 0.0;
            for (unsigned int node = 0; node < TNumNodes; ++node)
                FluidPressure += NContainer(GPoint, node) * PressureVector[node];

            RetentionParameters.SetFluidPressure(FluidPressure);

            double RelativePermeability = mRetentionLawVector[GPoint]->
                                            CalculateRelativePermeability(RetentionParameters);

            noalias(GradPressureTerm) =  prod(trans(GradNpT),PressureVector);

            noalias(GradPressureTerm) += PORE_PRESSURE_SIGN_FACTOR * FluidDensity * BodyAcceleration;

            noalias(FluidFlux) =   PORE_PRESSURE_SIGN_FACTOR
                                 * DynamicViscosityInverse
                                 * RelativePermeability
                                 * prod(PermeabilityMatrix,GradPressureTerm);

            GeoElementUtilities::FillArray1dOutput(rOutput[GPoint], FluidFlux);
        }
    } else {
        if ( rOutput.size() != mConstitutiveLawVector.size() )
            rOutput.resize(mConstitutiveLawVector.size());

        for ( unsigned int i = 0;  i < mConstitutiveLawVector.size(); ++i ) {
            noalias(rOutput[i]) = ZeroVector(3);
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rOutput[i] );
        }
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateOnIntegrationPoints<double,3>()") << std::endl;
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
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateOnIntegrationPoints()") << rVariable << std::endl;

    //Defining necessary variables
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != NumGPoints )
        rOutput.resize(NumGPoints);

    if (rVariable == CAUCHY_STRESS_TENSOR) {
        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint ) {
            rOutput[GPoint].resize(StressTensorSize, StressTensorSize, false );
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(mStressVector[GPoint]);
        }
    } else if (rVariable == TOTAL_STRESS_TENSOR) {
        //Defining necessary variables
        const PropertiesType& rProp = this->GetProperties();

        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,rProp,rCurrentProcessInfo);
        ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,
                                         rCurrentProcessInfo);

        // create general parametes of retention law
        RetentionLaw::Parameters RetentionParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);

        Vector VoigtVector(mStressVector[0].size());
        noalias(VoigtVector) = ZeroVector(VoigtVector.size());

        for (unsigned int i=0; i < StressTensorSize; ++i) VoigtVector[i] = 1.0;

        if ( rOutput.size() != NumGPoints )
            rOutput.resize(NumGPoints);

        const bool hasBiotCoefficient = rProp.Has(BIOT_COEFFICIENT);

        Vector TotalStressVector(mStressVector[0].size());

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint ) {
            //Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            //Compute infinitessimal strain
            this->CalculateStrain(Variables, GPoint);

            //set gauss points variables to constitutivelaw parameters
            this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

            //compute constitutive tensor and/or stresses
            noalias(Variables.StressVector) = mStressVector[GPoint];
            ConstitutiveParameters.SetStressVector(Variables.StressVector);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            Variables.BiotCoefficient = CalculateBiotCoefficient(Variables, hasBiotCoefficient);

            this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

            noalias(TotalStressVector) = mStressVector[GPoint];
            noalias(TotalStressVector) += PORE_PRESSURE_SIGN_FACTOR
                                         * Variables.BiotCoefficient
                                         * Variables.BishopCoefficient
                                         * Variables.FluidPressure
                                         * VoigtVector; 

            rOutput[GPoint].resize( StressTensorSize, StressTensorSize, false );
            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(TotalStressVector);
        }
    } else if (rVariable == ENGINEERING_STRAIN_TENSOR) {
        //Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            noalias(Variables.Np) = row(Variables.NContainer, GPoint);

            Matrix J0,InvJ0;
            this->CalculateDerivativesOnInitialConfiguration(Variables.detJInitialConfiguration,
                                                             J0,
                                                             InvJ0,
                                                             Variables.GradNpTInitialConfiguration,
                                                             GPoint);

            // Calculating operator B
            this->CalculateBMatrix( Variables.B, Variables.GradNpTInitialConfiguration, Variables.Np);

            //Compute infinitessimal strain
            this->CalculateCauchyStrain( Variables );

            if ( rOutput[GPoint].size2() != TDim )
                rOutput[GPoint].resize(TDim,TDim,false );

            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(Variables.StrainVector);
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR) {
        //Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,GPoint);

            //Compute strain
            this->CalculateStrain(Variables, GPoint);

            if ( rOutput[GPoint].size2() != TDim )
                rOutput[GPoint].resize(TDim,TDim,false );

            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(Variables.StrainVector);
        }
    } else if (rVariable == PERMEABILITY_MATRIX) {

        //If the permeability of the element is a given property
        BoundedMatrix<double,TDim,TDim> PermeabilityMatrix;
        GeoElementUtilities::FillPermeabilityMatrix(PermeabilityMatrix,this->GetProperties());

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint )
        {
            rOutput[GPoint].resize(TDim,TDim,false);
            noalias(rOutput[GPoint]) = PermeabilityMatrix;
        }
    } else {
        if ( rOutput.size() != mConstitutiveLawVector.size() )
            rOutput.resize(mConstitutiveLawVector.size());

        for ( unsigned int i = 0;  i < mConstitutiveLawVector.size(); ++i ) {
            rOutput[i].resize(TDim,TDim,false);
            noalias(rOutput[i]) = ZeroMatrix(TDim,TDim);
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rOutput[i] );
        }
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateMaterialStiffnessMatrix( MatrixType& rStiffnessMatrix,
                                      const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateMaterialStiffnessMatrix()") << std::endl;

    const unsigned int N_DOF = TNumNodes * (TDim + 1);

    //Resizing stiffness matrix
    if ( rStiffnessMatrix.size1() != N_DOF )
        rStiffnessMatrix.resize( N_DOF, N_DOF, false );
    noalias( rStiffnessMatrix ) = ZeroMatrix( N_DOF, N_DOF );

    //Previous definitions
    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = IntegrationPoints.size();

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,rProp,rCurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,
                                     rCurrentProcessInfo);

    //Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        //Compute Np, GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        //Compute constitutive tensor
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //Compute weighting coefficient for integration
        Variables.IntegrationCoefficient = 
            this->CalculateIntegrationCoefficient(IntegrationPoints,
                                                  GPoint,
                                                  Variables.detJ);

        //Compute stiffness matrix
        this->CalculateAndAddStiffnessMatrix(rStiffnessMatrix, Variables);
    }
    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateMaterialStiffnessMatrix()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddGeometricStiffnessMatrix( MatrixType& rLeftHandSideMatrix,
                                             ElementVariables& rVariables,
                                             unsigned int GPoint )
{
    KRATOS_TRY

    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( mStressVector[GPoint] );
    Matrix ReducedKgMatrix = prod( rVariables.GradNpT,
                                   rVariables.IntegrationCoefficient *
                                   Matrix( prod( StressTensor, trans(rVariables.GradNpT) ) ) ); //to be optimized

    Matrix UMatrix(TNumNodes*TDim, TNumNodes*TDim);
    noalias(UMatrix) = ZeroMatrix(TNumNodes*TDim, TNumNodes*TDim);
    MathUtils<double>::ExpandAndAddReducedMatrix( UMatrix, ReducedKgMatrix, TDim );

    //Distribute stiffness block matrix into the elemental matrix
    GeoElementUtilities::AssembleUBlockMatrix(rLeftHandSideMatrix, UMatrix, TNumNodes, TDim);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateMassMatrix( MatrixType& rMassMatrix,
                         const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int N_DOF = this->GetNumberOfDOF();

    //Resizing mass matrix
    if ( rMassMatrix.size1() != N_DOF )
        rMassMatrix.resize( N_DOF, N_DOF, false );
    noalias( rMassMatrix ) = ZeroMatrix( N_DOF, N_DOF );

    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = IntegrationPoints.size();

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables( Variables,
                                      rCurrentProcessInfo );

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);

    //Defining shape functions at all integration points
    //Defining necessary variables
    BoundedMatrix<double, TDim+1, TNumNodes*(TDim+1)> Nut              = ZeroMatrix(TDim+1, TNumNodes*(TDim+1));
    BoundedMatrix<double, TDim+1, TNumNodes*(TDim+1)> AuxDensityMatrix = ZeroMatrix(TDim+1, TNumNodes*(TDim+1));
    BoundedMatrix<double, TDim+1, TDim+1>             DensityMatrix    = ZeroMatrix(TDim+1, TDim+1);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint ) {
        GeoElementUtilities::CalculateNuElementMatrix<TDim, TNumNodes>(Nut, Variables.NContainer, GPoint);

        Matrix J0, InvJ0, DNu_DX0;
        this->CalculateDerivativesOnInitialConfiguration(Variables.detJInitialConfiguration,
                                                         J0,
                                                         InvJ0,
                                                         DNu_DX0,
                                                         GPoint);

        //calculating weighting coefficient for integration
        Variables.IntegrationCoefficientInitialConfiguration = 
            this->CalculateIntegrationCoefficient(IntegrationPoints,
                                                  GPoint,
                                                  Variables.detJInitialConfiguration);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->CalculateSoilDensity(Variables);

        GeoElementUtilities::
            AssembleDensityMatrix< TDim >(DensityMatrix, Variables.Density);
        
        noalias(AuxDensityMatrix) = prod(DensityMatrix, Nut);

        //Adding contribution to Mass matrix
        noalias(rMassMatrix) +=  prod(trans(Nut), AuxDensityMatrix)
                               * Variables.IntegrationCoefficientInitialConfiguration;
    }

    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& rCurrentProcessInfo,
                  const bool CalculateStiffnessMatrixFlag,
                  const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAll()") << std::endl;

    //Previous definitions
    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = IntegrationPoints.size();

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);

    // stiffness matrix is needed to calculate Biot coefficient
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag)  ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables( Variables,
                                      rCurrentProcessInfo );

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);

    const bool hasBiotCoefficient = rProp.Has(BIOT_COEFFICIENT);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        //Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        //Compute Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim, TNumNodes>( Variables.BodyAcceleration,
                                                                Variables.NContainer,
                                                                Variables.VolumeAcceleration,
                                                                GPoint );


        //Compute constitutive tensor and stresses
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);

        //Compute weighting coefficient for integration
        Variables.IntegrationCoefficient = 
            this->CalculateIntegrationCoefficient(IntegrationPoints,
                                                  GPoint,
                                                  Variables.detJ);

        Variables.IntegrationCoefficientInitialConfiguration =
            this->CalculateIntegrationCoefficient(IntegrationPoints,
                                                  GPoint,
                                                  Variables.detJInitialConfiguration);

        //Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if (CalculateResidualVectorFlag)  this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAll()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateBiotCoefficient( const ElementVariables& rVariables,
                              const bool &hasBiotCoefficient) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateBiotCoefficient()") << std::endl;

    const PropertiesType& rProp = this->GetProperties();

    //Properties variables
    if (hasBiotCoefficient) {
        return rProp[BIOT_COEFFICIENT];
    }
    else {
        // calculate Bulk modulus from stiffness matrix
        const double BulkModulus = CalculateBulkModulus(rVariables.ConstitutiveMatrix);
        return 1.0 - BulkModulus / rProp[BULK_MODULUS_SOLID];
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateBiotCoefficient()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    InitializeBiotCoefficients( ElementVariables& rVariables,
                                const bool &hasBiotCoefficient)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeBiotCoefficients()") << std::endl;

    const PropertiesType& rProp = this->GetProperties();

    //Properties variables
    rVariables.BiotCoefficient = CalculateBiotCoefficient(rVariables, hasBiotCoefficient);

    if (!rProp[IGNORE_UNDRAINED]) {
        rVariables.BiotModulusInverse =  (rVariables.BiotCoefficient - rProp[POROSITY])/rProp[BULK_MODULUS_SOLID]
                                       + rProp[POROSITY]/rProp[BULK_MODULUS_FLUID];
    } else {
        rVariables.BiotModulusInverse =  (rVariables.BiotCoefficient - rProp[POROSITY])/rProp[BULK_MODULUS_SOLID]
                                       + rProp[POROSITY]/TINY;
    }

    rVariables.BiotModulusInverse *= rVariables.DegreeOfSaturation;
    rVariables.BiotModulusInverse -= rVariables.DerivativeOfSaturation*rProp[POROSITY];

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeBiotCoefficients()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateBulkModulus(const Matrix &ConstitutiveMatrix) const
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
                                const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeElementVariables()") << this->Id() << std::endl;

    //Properties variables
    this->InitializeProperties( rVariables );

    //ProcessInfo variables
    rVariables.VelocityCoefficient = rCurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    //Nodal Variables
    this->InitializeNodalDisplacementVariables( rVariables );
    this->InitializeNodalPorePressureVariables( rVariables );
    this->InitializeNodalVolumeAccelerationVariables( rVariables );

    //Variables computed at each GP
    noalias(rVariables.Nu) = ZeroMatrix(TDim, TNumNodes*TDim);
    rVariables.Np.resize(TNumNodes,false);
    rVariables.GradNpT.resize(TNumNodes,TDim,false);
    rVariables.F.resize(TDim,TDim,false);

    noalias(rVariables.F) = identity_matrix<double>(TDim);

    rVariables.detF = 1.0;

    //General Variables
    rVariables.VoigtVector.resize(VoigtSize);
    noalias(rVariables.VoigtVector) = ZeroVector(VoigtSize);
    for (unsigned int i=0; i < TDim; ++i) rVariables.VoigtVector[i] = 1.0;

    rVariables.B.resize(VoigtSize,TNumNodes*TDim,false);
    noalias(rVariables.B) = ZeroMatrix(VoigtSize,TNumNodes*TDim);

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    // shape functions
    (rVariables.NContainer).resize(NumGPoints, TNumNodes, false);
    rVariables.NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );

    // gradient of shape functions and determinant of Jacobian
    (rVariables.detJContainer).resize(NumGPoints,false);

    rGeom.ShapeFunctionsIntegrationPointsGradients( rVariables.DN_DXContainer,
                                                   rVariables.detJContainer,
                                                   mThisIntegrationMethod );

    //Constitutive Law parameters
    
    rVariables.StressVector.resize(VoigtSize, false);
    rVariables.StrainVector.resize(VoigtSize, false);
    rVariables.ConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);

    //Auxiliary variables
    rVariables.UVoigtMatrix.resize(TNumNodes*TDim, VoigtSize, false);

    // Retention law
    rVariables.FluidPressure = 0.0;
    rVariables.DegreeOfSaturation = 1.0;
    rVariables.DerivativeOfSaturation = 0.0;
    rVariables.RelativePermeability = 1.0;
    rVariables.BishopCoefficient = 1.0;

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeElementVariables()") << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateBMatrix(Matrix& rB,
                     const Matrix& GradNpT,
                     const Vector &Np)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateBMatrix()") << std::endl;

    unsigned int index;

    if (TDim > 2) {
        for ( unsigned int i = 0; i < TNumNodes; ++i ) {
            index = TDim * i;

            rB( INDEX_3D_XX, index + INDEX_X ) = GradNpT( i, INDEX_X );
            rB( INDEX_3D_YY, index + INDEX_Y ) = GradNpT( i, INDEX_Y );
            rB( INDEX_3D_ZZ, index + INDEX_Z ) = GradNpT( i, INDEX_Z );
            rB( INDEX_3D_XY, index + INDEX_X ) = GradNpT( i, INDEX_Y );
            rB( INDEX_3D_XY, index + INDEX_Y ) = GradNpT( i, INDEX_X );
            rB( INDEX_3D_YZ, index + INDEX_Y ) = GradNpT( i, INDEX_Z );
            rB( INDEX_3D_YZ, index + INDEX_Z ) = GradNpT( i, INDEX_Y );
            rB( INDEX_3D_XZ, index + INDEX_X ) = GradNpT( i, INDEX_Z );
            rB( INDEX_3D_XZ, index + INDEX_Z ) = GradNpT( i, INDEX_X );
        }
    } else {
        // 2D plane strain
        for ( unsigned int i = 0; i < TNumNodes; ++i ) {
            index = TDim * i;

            rB( INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X ) = GradNpT( i, INDEX_X );
            rB( INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y ) = GradNpT( i, INDEX_Y );
            rB( INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X ) = GradNpT( i, INDEX_Y );
            rB( INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y ) = GradNpT( i, INDEX_X );
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
    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix,rVariables);
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

    noalias(rVariables.UPMatrix) =  PORE_PRESSURE_SIGN_FACTOR 
                                  * rVariables.BiotCoefficient
                                  * rVariables.BishopCoefficient
                                  * outer_prod(rVariables.UVector,rVariables.Np)
                                  * rVariables.IntegrationCoefficient;

    //Distribute coupling block matrix into the elemental matrix
    GeoElementUtilities::AssembleUPBlockMatrix<TDim, TNumNodes>(rLeftHandSideMatrix, rVariables.UPMatrix);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        noalias(rVariables.PUMatrix) =  PORE_PRESSURE_SIGN_FACTOR
                                      * SaturationCoefficient
                                      * rVariables.VelocityCoefficient
                                      * trans(rVariables.UPMatrix);

        //Distribute transposed coupling block matrix into the elemental matrix
        GeoElementUtilities::AssemblePUBlockMatrix<TDim, TNumNodes>(rLeftHandSideMatrix,rVariables.PUMatrix);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddCouplingMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateCompressibilityMatrix(BoundedMatrix<double,TNumNodes,TNumNodes> &PMatrix,
                                   const ElementVariables &rVariables) const
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddCompressibilityMatrix()") << std::endl;

    noalias(PMatrix) = - PORE_PRESSURE_SIGN_FACTOR 
                       * rVariables.DtPressureCoefficient
                       * rVariables.BiotModulusInverse
                       * outer_prod(rVariables.Np, rVariables.Np)
                       * rVariables.IntegrationCoefficient;

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddCompressibilityMatrix()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix,
                                         ElementVariables& rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddCompressibilityMatrix()") << std::endl;

    this->CalculateCompressibilityMatrix(rVariables.PMatrix, rVariables);

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
    CalculatePermeabilityMatrix(BoundedMatrix<double,TNumNodes,TDim> &PDimMatrix,
                                BoundedMatrix<double,TNumNodes,TNumNodes> &PMatrix,
                                const ElementVariables &rVariables) const
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculatePermeabilityMatrix()") << std::endl;

    noalias(PDimMatrix) = - PORE_PRESSURE_SIGN_FACTOR 
                          * prod(rVariables.GradNpT, rVariables.PermeabilityMatrix);

    noalias(PMatrix) =  rVariables.DynamicViscosityInverse
                      * rVariables.RelativePermeability
                      * prod(PDimMatrix, trans(rVariables.GradNpT))
                      * rVariables.IntegrationCoefficient;

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculatePermeabilityMatrix()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddPermeabilityMatrix(MatrixType &rLeftHandSideMatrix,
                                      ElementVariables &rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddPermeabilityMatrix()") << std::endl;

    this->CalculatePermeabilityMatrix(rVariables.PDimMatrix,
                                      rVariables.PMatrix,
                                      rVariables);

    //Distribute permeability block matrix into the elemental matrix
    GeoElementUtilities::
        AssemblePBlockMatrix< TDim, TNumNodes >(rLeftHandSideMatrix, rVariables.PMatrix);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddPermeabilityMatrix()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddRHS(VectorType& rRightHandSideVector,
                       ElementVariables& rVariables,
                       unsigned int GPoint)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddRHS()") << std::endl;

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables, GPoint);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    if (!rVariables.IgnoreUndrained) {
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
                                   ElementVariables& rVariables,
                                   unsigned int GPoint )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddStiffnessForce()") << std::endl;

    noalias(rVariables.UVector) = -1.0 * prod(trans(rVariables.B), mStressVector[GPoint])
                                       * rVariables.IntegrationCoefficient;

    //Distribute stiffness block vector into elemental vector
    GeoElementUtilities::AssembleUBlockVector<TDim, TNumNodes>(rRightHandSideVector, rVariables.UVector);

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

    this->CalculateSoilGamma(rVariables);

    noalias(rVariables.UVector) =   prod(trans(rVariables.Nu), rVariables.SoilGamma)
                                  * rVariables.IntegrationCoefficientInitialConfiguration;

    //Distribute body force block vector into elemental vector
    GeoElementUtilities::AssembleUBlockVector<TDim, TNumNodes>(rRightHandSideVector,rVariables.UVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddMixBodyForce()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateSoilDensity(ElementVariables &rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateSoilDensity()") << std::endl;

    rVariables.Density = (  rVariables.DegreeOfSaturation
                          * rVariables.Porosity
                          * rVariables.FluidDensity )
                        + (1.0 - rVariables.Porosity)*rVariables.SolidDensity;

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateSoilDensity()") << std::endl;
    KRATOS_CATCH("");

}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateSoilGamma(ElementVariables &rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateSoilGamma()") << std::endl;

    this->CalculateSoilDensity(rVariables);

    noalias(rVariables.SoilGamma) = rVariables.Density * rVariables.BodyAcceleration;

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateSoilGamma()") << std::endl;
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

    noalias(rVariables.UPMatrix) = - PORE_PRESSURE_SIGN_FACTOR
                                   * rVariables.BiotCoefficient
                                   * rVariables.BishopCoefficient
                                   * outer_prod(rVariables.UVector,rVariables.Np)
                                   * rVariables.IntegrationCoefficient;

    noalias(rVariables.UVector) = prod(rVariables.UPMatrix,rVariables.PressureVector);

    //Distribute coupling block vector 1 into elemental vector
    GeoElementUtilities::AssembleUBlockVector<TDim, TNumNodes>(rRightHandSideVector, rVariables.UVector);

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        noalias(rVariables.PVector) =  PORE_PRESSURE_SIGN_FACTOR
                                     * SaturationCoefficient
                                     * prod(trans(rVariables.UPMatrix), rVariables.VelocityVector);

        //Distribute coupling block vector 2 into elemental vector
        GeoElementUtilities::AssemblePBlockVector< TDim, TNumNodes >(rRightHandSideVector, rVariables.PVector);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddCouplingTerms()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateCompressibilityFlow(BoundedMatrix<double,TNumNodes,TNumNodes> &PMatrix,
                                 array_1d<double,TNumNodes> &PVector,
                                 const ElementVariables &rVariables) const
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateCompressibilityFlow()") << std::endl;

    noalias(PMatrix) = - PORE_PRESSURE_SIGN_FACTOR 
                       * rVariables.BiotModulusInverse
                       * outer_prod(rVariables.Np,rVariables.Np)
                       * rVariables.IntegrationCoefficient;

    noalias(PVector) = - prod(PMatrix, rVariables.DtPressureVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateCompressibilityFlow()") << std::endl;
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

    this->CalculateCompressibilityFlow(rVariables.PMatrix,
                                       rVariables.PVector,
                                       rVariables);

    //Distribute compressibility block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<TDim, TNumNodes>(rRightHandSideVector, rVariables.PVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddCompressibilityFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculatePermeabilityFlow(BoundedMatrix<double,TNumNodes,TDim> &PDimMatrix,
                              BoundedMatrix<double,TNumNodes,TNumNodes> &PMatrix,
                              array_1d<double,TNumNodes> &PVector,
                              const ElementVariables &rVariables) const
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculatePermeabilityFlow()") << std::endl;

    noalias(PDimMatrix) = prod(rVariables.GradNpT, rVariables.PermeabilityMatrix);

    noalias(PMatrix) = - PORE_PRESSURE_SIGN_FACTOR 
                       * rVariables.DynamicViscosityInverse
                       * rVariables.RelativePermeability
                       * prod(PDimMatrix,trans(rVariables.GradNpT))
                       * rVariables.IntegrationCoefficient;

    noalias(PVector) = - prod(PMatrix, rVariables.PressureVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculatePermeabilityFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateAndAddPermeabilityFlow( VectorType &rRightHandSideVector,
                                     ElementVariables &rVariables )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateAndAddPermeabilityFlow()") << std::endl;

    this->CalculatePermeabilityFlow(rVariables.PDimMatrix,
                                    rVariables.PMatrix,
                                    rVariables.PVector,
                                    rVariables);

    //Distribute permeability block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<TDim, TNumNodes>(rRightHandSideVector, rVariables.PVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddPermeabilityFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateFluidBodyFlow(BoundedMatrix<double,TNumNodes,TDim> &PDimMatrix,
                           array_1d<double,TNumNodes> &PVector,
                           const ElementVariables &rVariables) const
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateFluidBodyFlow()") << std::endl;

    noalias(PDimMatrix) =  prod(rVariables.GradNpT,rVariables.PermeabilityMatrix)
                         * rVariables.IntegrationCoefficient;

    noalias(PVector) =  rVariables.DynamicViscosityInverse
                      * rVariables.FluidDensity
                      * rVariables.RelativePermeability
                      * prod(PDimMatrix,rVariables.BodyAcceleration);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateFluidBodyFlow()") << std::endl;
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

    this->CalculateFluidBodyFlow(rVariables.PDimMatrix,
                                 rVariables.PVector,
                                 rVariables);

    //Distribute fluid body flow block vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector<TDim, TNumNodes>(rRightHandSideVector,rVariables.PVector);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateAndAddFluidBodyFlow()") << std::endl;
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateStrain( ElementVariables& rVariables, const IndexType& GPoint )
{
    if (rVariables.UseHenckyStrain) {
        this->CalculateDeformationGradient(rVariables, GPoint);
        this->CalculateHenckyStrain( rVariables );
    } else {
        this->CalculateCauchyStrain( rVariables );
    }
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateCauchyStrain( ElementVariables& rVariables )
{
    noalias(rVariables.StrainVector) = prod(rVariables.B, rVariables.DisplacementVector);
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateCauchyGreenStrain( ElementVariables& rVariables )
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateCauchyGreenStrain()") << std::endl;

    //-Compute total deformation gradient
    const Matrix& F = rVariables.F;

    Matrix ETensor;
    ETensor = prod(trans(F), F);

    for (unsigned int i=0; i<TDim; ++i)
        ETensor(i,i) -= 1.0;
    ETensor *= 0.5;

    if (TDim==2) {
        Vector StrainVector;
        StrainVector = MathUtils<double>::StrainTensorToVector(ETensor);
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_XX] = StrainVector[0];
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_YY] = StrainVector[1];
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_ZZ] = 0.0;
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_XY] = StrainVector[2];
    } else {
        noalias(rVariables.StrainVector) = MathUtils<double>::StrainTensorToVector(ETensor);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateCauchyGreenStrain()") << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateCauchyAlmansiStrain(ElementVariables& rVariables )
{
   // KRATOS_INFO("0-UPwSmallStrainElement::CalculateCauchyAlmansiStrain()") << std::endl;

    //-Compute total deformation gradient
    const Matrix& F = rVariables.F;

    Matrix LeftCauchyGreen;
    LeftCauchyGreen = prod(F, trans(F));

    Matrix ETensor;
    double det;
    MathUtils<double>::InvertMatrix(LeftCauchyGreen, ETensor, det );

    for (unsigned int i=0; i<TDim; ++i)
        ETensor(i,i) = 1.0 - ETensor(i,i);

    ETensor *= 0.5;

    if (TDim==2) {
        Vector StrainVector;
        StrainVector = MathUtils<double>::StrainTensorToVector(ETensor);
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_XX] = StrainVector[0];
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_YY] = StrainVector[1];
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_ZZ] = 0.0;
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_XY] = StrainVector[2];
    } else {
        noalias(rVariables.StrainVector) = MathUtils<double>::StrainTensorToVector(ETensor);
    }

   // KRATOS_INFO("1-UPwSmallStrainElement::CalculateCauchyAlmansiStrain()") << std::endl;
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateHenckyStrain( ElementVariables& rVariables )
{
    KRATOS_TRY

    //-Compute total deformation gradient
    const Matrix& F = rVariables.F;

    Matrix CMatrix;
    CMatrix = prod(trans(F), F);

    // Declare the different matrix
    Matrix EigenValuesMatrix = ZeroMatrix(TDim, TDim);
    Matrix EigenVectorsMatrix = ZeroMatrix(TDim, TDim);

    // Decompose matrix
    MathUtils<double>::GaussSeidelEigenSystem(CMatrix, EigenVectorsMatrix, EigenValuesMatrix, 1.0e-16, 20);

    // Calculate the eigenvalues of the E matrix
    for (IndexType i = 0; i < TDim; ++i) {
        EigenValuesMatrix(i, i) = 0.5 * std::log(EigenValuesMatrix(i, i));
    }

    // Calculate E matrix
    Matrix ETensor = ZeroMatrix(TDim, TDim);
    MathUtils<double>::BDBtProductOperation(ETensor, EigenValuesMatrix, EigenVectorsMatrix);

    // Hencky Strain Calculation
    if (TDim==2) {
        Vector StrainVector;
        StrainVector = MathUtils<double>::StrainTensorToVector(ETensor);
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_XX] = StrainVector[0];
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_YY] = StrainVector[1];
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_ZZ] = 0.0;
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_XY] = StrainVector[2];
    } else {
        noalias(rVariables.StrainVector) = MathUtils<double>::StrainTensorToVector(ETensor);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateDeformationGradient( ElementVariables& rVariables,
                                  const IndexType& GPoint)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateDeformationGradient()") << std::endl;

    // calculation of derivative of shape function with respect to reference configuration
    // derivative of shape function (displacement)
    Matrix J0, InvJ0, DNu_DX0;
    double detJ0;
    this->CalculateDerivativesOnInitialConfiguration(detJ0,
                                                     J0,
                                                     InvJ0,
                                                     DNu_DX0,
                                                     GPoint);

    //Calculating current jacobian in order to find deformation gradient
    Matrix J, InvJ, DNu_DX;
    double detJ;
    this->CalculateJacobianOnCurrentConfiguration(detJ,
                                                  J,
                                                  InvJ,
                                                  GPoint);


#ifdef KRATOS_COMPILED_IN_WINDOWS
    if (detJ < 0.0) {
        KRATOS_INFO("negative detJ")
        << "ERROR:: ELEMENT ID: "
        << this->Id()
        << " INVERTED. DETJ: "
        << detJ
        << " nodes:" << this->GetGeometry()
        << std::endl;
    }
#endif

    KRATOS_ERROR_IF(detJ < 0.0)
     << "ERROR:: ELEMENT ID: "
     << this->Id()
     << " INVERTED. DETJ: "
     << detJ
     << std::endl;


    // Deformation gradient
    noalias(rVariables.F) = prod( J, InvJ0 );
    rVariables.detF = MathUtils<double>::Det(rVariables.F);

    // KRATOS_INFO("1-UpdatedLagrangianUPwDiffOrderElement::CalculateDeformationGradient()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    InitializeNodalPorePressureVariables( ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeNodalPorePressureVariables") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();

    //Nodal Variables
    for (unsigned int i=0; i<TNumNodes; ++i) {
        rVariables.PressureVector[i] = rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.DtPressureVector[i] = rGeom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeNodalPorePressureVariables") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    InitializeNodalDisplacementVariables( ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeNodalDisplacementVariables") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();

    //Nodal Variables
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.DisplacementVector, rGeom, DISPLACEMENT);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VelocityVector,     rGeom, VELOCITY);

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeNodalDisplacementVariables") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    InitializeNodalVolumeAccelerationVariables( ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeNodalVolumeAccelerationVariables") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();

    //Nodal Variables
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rVariables.VolumeAcceleration, rGeom, VOLUME_ACCELERATION);

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeNodalVolumeAccelerationVariables") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    InitializeProperties( ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::InitializeProperties") << std::endl;

    const PropertiesType& rProp = this->GetProperties();

    rVariables.IgnoreUndrained = rProp[IGNORE_UNDRAINED];
    rVariables.UseHenckyStrain = false;
    if (rProp.Has(USE_HENCKY_STRAIN))
        rVariables.UseHenckyStrain = rProp[USE_HENCKY_STRAIN];

    rVariables.ConsiderGeometricStiffness = false;
    if (rProp.Has(CONSIDER_GEOMETRIC_STIFFNESS))
        rVariables.ConsiderGeometricStiffness = rProp[CONSIDER_GEOMETRIC_STIFFNESS];

    rVariables.DynamicViscosityInverse = 1.0/rProp[DYNAMIC_VISCOSITY];
    rVariables.FluidDensity = rProp[DENSITY_WATER];
    rVariables.SolidDensity = rProp[DENSITY_SOLID];
    rVariables.Porosity     = rProp[POROSITY];
    GeoElementUtilities::FillPermeabilityMatrix(rVariables.PermeabilityMatrix, rProp);

    // KRATOS_INFO("1-UPwSmallStrainElement::InitializeProperties") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateKinematics( ElementVariables& rVariables,
                         const unsigned int &GPoint )

{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateKinematics") << std::endl;

    //Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Np) = row(rVariables.NContainer, GPoint);
    noalias(rVariables.GradNpT) = rVariables.DN_DXContainer[GPoint];

    rVariables.detJ = rVariables.detJContainer[GPoint];

    //Compute the deformation matrix B
    this->CalculateBMatrix(rVariables.B, rVariables.GradNpT, rVariables.Np);

    Matrix J0,InvJ0;
    this->CalculateDerivativesOnInitialConfiguration(rVariables.detJInitialConfiguration,
                                                     J0,
                                                     InvJ0,
                                                     rVariables.GradNpTInitialConfiguration,
                                                     GPoint);

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateKinematics") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    SetConstitutiveParameters(ElementVariables& rVariables,
                              ConstitutiveLaw::Parameters& rConstitutiveParameters)
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::SetConstitutiveParameters") << std::endl;

    rConstitutiveParameters.SetStrainVector(rVariables.StrainVector);
    // rConstitutiveParameters.SetStressVector(rVariables.StressVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rConstitutiveParameters.SetShapeFunctionsValues(rVariables.Np);
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rVariables.GradNpT);
    rConstitutiveParameters.SetDeformationGradientF(rVariables.F);
    rConstitutiveParameters.SetDeterminantF(rVariables.detF);

    // KRATOS_INFO("1-UPwSmallStrainElement::SetConstitutiveParameters") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    SetRetentionParameters(const ElementVariables& rVariables,
                           RetentionLaw::Parameters& rRetentionParameters)
{
    KRATOS_TRY

    rRetentionParameters.SetFluidPressure(rVariables.FluidPressure);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateFluidPressure(const ElementVariables &rVariables, const unsigned int &GPoint)
{
    KRATOS_TRY

    return inner_prod(rVariables.Np, rVariables.PressureVector);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateRetentionResponse( ElementVariables& rVariables,
                                RetentionLaw::Parameters& rRetentionParameters,
                                const unsigned int &GPoint )
{
    KRATOS_TRY

    rVariables.FluidPressure = CalculateFluidPressure(rVariables, GPoint);
    SetRetentionParameters(rVariables, rRetentionParameters);

    rVariables.DegreeOfSaturation = mRetentionLawVector[GPoint]->CalculateSaturation(rRetentionParameters);
    rVariables.DerivativeOfSaturation = mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(rRetentionParameters);
    rVariables.RelativePermeability = mRetentionLawVector[GPoint]->CalculateRelativePermeability(rRetentionParameters);
    rVariables.BishopCoefficient = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(rRetentionParameters);

    KRATOS_CATCH( "" )
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwSmallStrainElement<TDim,TNumNodes>::
    CalculateExtrapolationMatrix(BoundedMatrix<double,TNumNodes,TNumNodes>& rExtrapolationMatrix)
{
    KRATOS_TRY

    KRATOS_ERROR << "undefined number of nodes in CalculateExtrapolationMatrix ... TNumNodes:" 
                 << TNumNodes
                 << " element: "
                 << this->Id()
                 << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void UPwSmallStrainElement< 2, 3 >::
    CalculateExtrapolationMatrix(BoundedMatrix<double,3,3>& rExtrapolationMatrix)
{
    /// The matrix contains the shape functions at each GP evaluated at each node.
    /// Rows: nodes
    /// Columns: GP

    //Triangle_2d_3
    //GI_GAUSS_2

    rExtrapolationMatrix(0,0) = 1.6666666666666666666; rExtrapolationMatrix(0,1) = -0.33333333333333333333; rExtrapolationMatrix(0,2) = -0.33333333333333333333;
    rExtrapolationMatrix(1,0) = -0.33333333333333333333; rExtrapolationMatrix(1,1) = 1.6666666666666666666; rExtrapolationMatrix(1,2) = -0.33333333333333333333;
    rExtrapolationMatrix(2,0) = -0.33333333333333333333; rExtrapolationMatrix(2,1) = -0.33333333333333333333; rExtrapolationMatrix(2,2) = 1.6666666666666666666;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void UPwSmallStrainElement< 2, 4 >::
    CalculateExtrapolationMatrix(BoundedMatrix<double,4,4>& rExtrapolationMatrix)
{
    //Quadrilateral_2d_4
    //GI_GAUSS_2

    rExtrapolationMatrix(0,0) = 1.8660254037844386; rExtrapolationMatrix(0,1) = -0.5; rExtrapolationMatrix(0,2) = 0.13397459621556132; rExtrapolationMatrix(0,3) = -0.5;
    rExtrapolationMatrix(1,0) = -0.5; rExtrapolationMatrix(1,1) = 1.8660254037844386; rExtrapolationMatrix(1,2) = -0.5; rExtrapolationMatrix(1,3) = 0.13397459621556132;
    rExtrapolationMatrix(2,0) = 0.13397459621556132; rExtrapolationMatrix(2,1) = -0.5; rExtrapolationMatrix(2,2) = 1.8660254037844386; rExtrapolationMatrix(2,3) = -0.5;
    rExtrapolationMatrix(3,0) = -0.5; rExtrapolationMatrix(3,1) = 0.13397459621556132; rExtrapolationMatrix(3,2) = -0.5; rExtrapolationMatrix(3,3) = 1.8660254037844386;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void UPwSmallStrainElement< 3, 4 >::
    CalculateExtrapolationMatrix(BoundedMatrix<double,4,4>& rExtrapolationMatrix)
{
    //Tetrahedra_3d_4
    //GI_GAUSS_2

    rExtrapolationMatrix(0,0) = -0.309016988749894905; rExtrapolationMatrix(0,1) = -0.3090169887498949046; rExtrapolationMatrix(0,2) = -0.309016988749894905; rExtrapolationMatrix(0,3) = 1.9270509662496847144;
    rExtrapolationMatrix(1,0) = 1.9270509662496847144; rExtrapolationMatrix(1,1) = -0.30901698874989490481; rExtrapolationMatrix(1,2) = -0.3090169887498949049; rExtrapolationMatrix(1,3) = -0.30901698874989490481;
    rExtrapolationMatrix(2,0) = -0.30901698874989490473; rExtrapolationMatrix(2,1) = 1.9270509662496847143; rExtrapolationMatrix(2,2) = -0.3090169887498949049; rExtrapolationMatrix(2,3) = -0.30901698874989490481;
    rExtrapolationMatrix(3,0) = -0.3090169887498949048; rExtrapolationMatrix(3,1) = -0.30901698874989490471; rExtrapolationMatrix(3,2) = 1.9270509662496847143; rExtrapolationMatrix(3,3) = -0.30901698874989490481;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void UPwSmallStrainElement< 3, 8 >::
    CalculateExtrapolationMatrix(BoundedMatrix<double,8,8>& rExtrapolationMatrix)
{
    //Hexahedra_3d_8
    //GI_GAUSS_2

    rExtrapolationMatrix(0,0) = 2.549038105676658; rExtrapolationMatrix(0,1) = -0.6830127018922192; rExtrapolationMatrix(0,2) = 0.18301270189221927; rExtrapolationMatrix(0,3) = -0.6830127018922192;
    rExtrapolationMatrix(0,4) = -0.6830127018922192; rExtrapolationMatrix(0,5) = 0.18301270189221927; rExtrapolationMatrix(0,6) = -0.04903810567665795; rExtrapolationMatrix(0,7) = 0.18301270189221927;

    rExtrapolationMatrix(1,0) = -0.6830127018922192; rExtrapolationMatrix(1,1) = 2.549038105676658; rExtrapolationMatrix(1,2) = -0.6830127018922192; rExtrapolationMatrix(1,3) = 0.18301270189221927;
    rExtrapolationMatrix(1,4) = 0.18301270189221927; rExtrapolationMatrix(1,5) = -0.6830127018922192; rExtrapolationMatrix(1,6) = 0.18301270189221927; rExtrapolationMatrix(1,7) = -0.04903810567665795;

    rExtrapolationMatrix(2,0) = 0.18301270189221927; rExtrapolationMatrix(2,1) = -0.6830127018922192; rExtrapolationMatrix(2,2) = 2.549038105676658; rExtrapolationMatrix(2,3) = -0.6830127018922192;
    rExtrapolationMatrix(2,4) = -0.04903810567665795; rExtrapolationMatrix(2,5) = 0.18301270189221927; rExtrapolationMatrix(2,6) = -0.6830127018922192; rExtrapolationMatrix(2,7) = 0.18301270189221927;

    rExtrapolationMatrix(3,0) = -0.6830127018922192; rExtrapolationMatrix(3,1) = 0.18301270189221927; rExtrapolationMatrix(3,2) = -0.6830127018922192; rExtrapolationMatrix(3,3) = 2.549038105676658;
    rExtrapolationMatrix(3,4) = 0.18301270189221927; rExtrapolationMatrix(3,5) = -0.04903810567665795; rExtrapolationMatrix(3,6) = 0.18301270189221927; rExtrapolationMatrix(3,7) = -0.6830127018922192;

    rExtrapolationMatrix(4,0) = -0.6830127018922192; rExtrapolationMatrix(4,1) = 0.18301270189221927; rExtrapolationMatrix(4,2) = -0.04903810567665795; rExtrapolationMatrix(4,3) = 0.18301270189221927;
    rExtrapolationMatrix(4,4) = 2.549038105676658; rExtrapolationMatrix(4,5) = -0.6830127018922192; rExtrapolationMatrix(4,6) = 0.18301270189221927; rExtrapolationMatrix(4,7) = -0.6830127018922192;

    rExtrapolationMatrix(5,0) = 0.18301270189221927; rExtrapolationMatrix(5,1) = -0.6830127018922192; rExtrapolationMatrix(5,2) = 0.18301270189221927; rExtrapolationMatrix(5,3) = -0.04903810567665795;
    rExtrapolationMatrix(5,4) = -0.6830127018922192; rExtrapolationMatrix(5,5) = 2.549038105676658; rExtrapolationMatrix(5,6) = -0.6830127018922192; rExtrapolationMatrix(5,7) = 0.18301270189221927;

    rExtrapolationMatrix(6,0) = -0.04903810567665795; rExtrapolationMatrix(6,1) = 0.18301270189221927; rExtrapolationMatrix(6,2) = -0.6830127018922192; rExtrapolationMatrix(6,3) = 0.18301270189221927;
    rExtrapolationMatrix(6,4) = 0.18301270189221927; rExtrapolationMatrix(6,5) = -0.6830127018922192; rExtrapolationMatrix(6,6) = 2.549038105676658; rExtrapolationMatrix(6,7) = -0.6830127018922192;

    rExtrapolationMatrix(7,0) = 0.18301270189221927; rExtrapolationMatrix(7,1) = -0.04903810567665795; rExtrapolationMatrix(7,2) = 0.18301270189221927; rExtrapolationMatrix(7,3) = -0.6830127018922192;
    rExtrapolationMatrix(7,4) = -0.6830127018922192; rExtrapolationMatrix(7,5) = 0.18301270189221927; rExtrapolationMatrix(7,6) = -0.6830127018922192; rExtrapolationMatrix(7,7) = 2.549038105676658;
}

//----------------------------------------------------------------------------------------------------

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
