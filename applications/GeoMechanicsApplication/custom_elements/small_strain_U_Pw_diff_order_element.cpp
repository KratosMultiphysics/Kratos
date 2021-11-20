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

// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"

// Application includes
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

// Default Constructor
SmallStrainUPwDiffOrderElement::SmallStrainUPwDiffOrderElement() : Element() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Constructor 1
SmallStrainUPwDiffOrderElement::
    SmallStrainUPwDiffOrderElement( IndexType NewId, GeometryType::Pointer pGeometry ) : Element( NewId, pGeometry ) {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Constructor 2
SmallStrainUPwDiffOrderElement::
    SmallStrainUPwDiffOrderElement( IndexType NewId,
                                    GeometryType::Pointer pGeometry,
                                    PropertiesType::Pointer pProperties ) : Element( NewId, pGeometry, pProperties ) {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Destructor
SmallStrainUPwDiffOrderElement::~SmallStrainUPwDiffOrderElement() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Element::Pointer SmallStrainUPwDiffOrderElement::Create( IndexType NewId,
                                                         NodesArrayType const& ThisNodes,
                                                         PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SmallStrainUPwDiffOrderElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Element::Pointer SmallStrainUPwDiffOrderElement::Create(IndexType NewId,
                                                        GeometryType::Pointer pGeom,
                                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new SmallStrainUPwDiffOrderElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int  SmallStrainUPwDiffOrderElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::Check()") << this->Id() << std::endl;

    const GeometryType& rGeom = GetGeometry();

    if (rGeom.DomainSize() < 1.0e-15)
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

    //verify that the variables are correctly initialized
    // Verify specific properties
    const PropertiesType& rProp = this->GetProperties();

    if ( rProp.Has( IGNORE_UNDRAINED ) == false)
        KRATOS_ERROR << "IGNORE_UNDRAINED does not exist in the parameter list" << this->Id() << std::endl;

    bool IgnoreUndrained = rProp[IGNORE_UNDRAINED];

    if (!IgnoreUndrained) {
        if ( rProp.Has( PERMEABILITY_XX ) == false || rProp[PERMEABILITY_XX] < 0.0 )
            KRATOS_ERROR << "PERMEABILITY_XX has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

        if ( rProp.Has( PERMEABILITY_YY ) == false || rProp[PERMEABILITY_YY] < 0.0 )
            KRATOS_ERROR << "PERMEABILITY_YY has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

        if ( rProp.Has( PERMEABILITY_XY ) == false || rProp[PERMEABILITY_XY] < 0.0 )
            KRATOS_ERROR << "PERMEABILITY_XY has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

        if (rGeom.WorkingSpaceDimension() > 2) {
            if ( rProp.Has( PERMEABILITY_ZZ ) == false || rProp[PERMEABILITY_ZZ] < 0.0 )
                KRATOS_ERROR << "PERMEABILITY_ZZ has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

            if ( rProp.Has( PERMEABILITY_YZ ) == false || rProp[PERMEABILITY_YZ] < 0.0 )
                KRATOS_ERROR << "PERMEABILITY_YZ has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;

            if ( rProp.Has( PERMEABILITY_ZX ) == false || rProp[PERMEABILITY_ZX] < 0.0 )
                KRATOS_ERROR << "PERMEABILITY_ZX has Key zero, is not defined or has an invalid value at element" << this->Id() << std::endl;
        }
    }

    //verify that the dofs exist
    for ( unsigned int i = 0; i < rGeom.size(); ++i ) {
        if ( rGeom[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            KRATOS_ERROR << "missing variable DISPLACEMENT on node " << rGeom[i].Id() << std::endl;

        if ( rGeom[i].HasDofFor( DISPLACEMENT_X ) == false || rGeom[i].HasDofFor( DISPLACEMENT_Y ) == false || rGeom[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << rGeom[i].Id() << std::endl;

        if ( rGeom[i].SolutionStepsDataHas( WATER_PRESSURE ) == false )
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << rGeom[i].Id() << std::endl;

        if ( rGeom[i].HasDofFor( WATER_PRESSURE ) == false )
            KRATOS_ERROR << "missing the dof for the variable WATER_PRESSURE on node " << rGeom[i].Id() << std::endl;
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(rProp.Has( CONSTITUTIVE_LAW )) << "Constitutive law not provided for property " << rProp.Id() << std::endl;

    //verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    rProp.GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for (unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); ++i) {
        if (LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal)
            correct_strain_measure = true;
    }

    if ( correct_strain_measure == false )
        KRATOS_ERROR << "constitutive law is not compatible with the element type StrainMeasure_Infinitesimal " << this->Id() << std::endl;

    rProp.GetValue( CONSTITUTIVE_LAW )->Check( rProp, rGeom, rCurrentProcessInfo );

    // Verify that the constitutive law has the correct dimension
    const SizeType strainSize = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    if ( rGeom.WorkingSpaceDimension() > 2 ) {
        KRATOS_ERROR_IF_NOT( strainSize == VOIGT_SIZE_3D )
        << "Wrong constitutive law used. This is a 3D element! expected strain size is "
        << VOIGT_SIZE_3D
        << " But received: "
        << strainSize
        << " in element id: "
        << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT( strainSize == VOIGT_SIZE_2D_PLANE_STRAIN )
        << "Wrong constitutive law used. This is a 2D element! expected strain size is "
        << VOIGT_SIZE_2D_PLANE_STRAIN
        << " But received: "
        << strainSize
        << " in element id: "
        << this->Id() << std::endl;
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::Check()") << this->Id() << std::endl;

    return 0;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::Initialize()") << this->Id() << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const GeometryType::IntegrationPointsArrayType&
        IntegrationPoints = rGeom.IntegrationPoints( this->GetIntegrationMethod() );

    if ( mConstitutiveLawVector.size() != IntegrationPoints.size() )
        mConstitutiveLawVector.resize( IntegrationPoints.size() );

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL ) {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i ) {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(),
                                                           rGeom,
                                                           row(rGeom.ShapeFunctionsValues(this->GetIntegrationMethod()),i) );
        }
    }
    else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    // Retention law
    if ( mRetentionLawVector.size() != IntegrationPoints.size() )
        mRetentionLawVector.resize( IntegrationPoints.size() );

    for ( unsigned int i = 0; i < mRetentionLawVector.size(); ++i ) {
        //RetentionLawFactory::Pointer pRetentionFactory;
        mRetentionLawVector[i] = RetentionLawFactory::Clone(GetProperties());
        mRetentionLawVector[i]-> InitializeMaterial( GetProperties(),
                                                     rGeom,
                                                     row( rGeom.ShapeFunctionsValues( this->GetIntegrationMethod() ), i ) );
    }

    const SizeType NumUNodes = rGeom.PointsNumber();

    switch(NumUNodes)
    {
        case 6: //2D T6P3
            mpPressureGeometry = GeometryType::Pointer( new Triangle2D3< Node<3> >(rGeom(0), rGeom(1), rGeom(2)) );
            break;
        case 8: //2D Q8P4
            mpPressureGeometry = GeometryType::Pointer( new Quadrilateral2D4< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 9: //2D Q9P4
            mpPressureGeometry = GeometryType::Pointer( new Quadrilateral2D4< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 10: //3D T10P4
            mpPressureGeometry = GeometryType::Pointer( new Tetrahedra3D4< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 20: //3D H20P8
            mpPressureGeometry = GeometryType::Pointer( new Hexahedra3D8< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3), rGeom(4), rGeom(5), rGeom(6), rGeom(7)) );
            break;
        case 27: //3D H27P8
            mpPressureGeometry = GeometryType::Pointer( new Hexahedra3D8< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3), rGeom(4), rGeom(5), rGeom(6), rGeom(7)) );
            break;
        default:
            KRATOS_ERROR << "Unexpected geometry type for different order interpolation element" << this->Id() << std::endl;
            break;
    }

    // resize mStressVector:
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType VoigtSize  = ( Dim == 3 ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN);
    if ( mStressVector.size() != IntegrationPoints.size() ) {
       mStressVector.resize(IntegrationPoints.size());
       for (unsigned int i=0; i < mStressVector.size(); ++i) {
          mStressVector[i].resize(VoigtSize);
          std::fill(mStressVector[i].begin(), mStressVector[i].end(), 0.0);
       }
    }

    if ( mStateVariablesFinalized.size() != IntegrationPoints.size() )
        mStateVariablesFinalized.resize(IntegrationPoints.size());

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i ) {
        int nStateVariables = 0;
        nStateVariables = mConstitutiveLawVector[i]->GetValue( NUMBER_OF_UMAT_STATE_VARIABLES,
                                                               nStateVariables );
        if (nStateVariables > 0) {
            mConstitutiveLawVector[i]->SetValue( STATE_VARIABLES,
                                                 mStateVariablesFinalized[i],
                                                 rCurrentProcessInfo );
        }
    }

    mIsInitialised = true;
    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::Initialize()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    ResetConstitutiveLaw()
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::ResetConstitutiveLaw()") << std::endl;

    // erasing stress vectors
    for (unsigned int i=0; i < mStressVector.size(); ++i) {
        mStressVector[i].clear();
    }
    mStressVector.clear();

    for (unsigned int i=0; i < mStateVariablesFinalized.size(); ++i) {
        mStateVariablesFinalized[i].clear();
    }
    mStateVariablesFinalized.clear();

    KRATOS_CATCH( "" )

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::ResetConstitutiveLaw()") << std::endl;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::InitializeSolutionStep()") << std::endl;

    if (!mIsInitialised) this->Initialize(rCurrentProcessInfo);

    //Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // Set constitutive law flags:
    //ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); //Note: this is for nonlocal damage

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
        //compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables,GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables, ConstitutiveParameters);

        //compute constitutive tensor and/or stresses
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->InitializeMaterialResponseCauchy(ConstitutiveParameters);

        // retention law
        mRetentionLawVector[GPoint]->InitializeSolutionStep(RetentionParameters);
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::InitializeSolutionStep()") << std::endl;
    KRATOS_CATCH( "" )
}

//-------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    GetDofList( DofsVectorType& rElementalDofList,
                const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::GetDofList()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    // const bool IgnoreUndrained = GetProperties()[IGNORE_UNDRAINED];
    // const SizeType ElementSize = (IgnoreUndrained ? NumUNodes * Dim : NumUNodes * Dim + NumPNodes);
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    if (rElementalDofList.size() != ElementSize)
        rElementalDofList.resize(ElementSize);

    SizeType Index = 0;

    if (Dim > 2) {
        for (SizeType i = 0; i < NumUNodes; ++i) {
            rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
            rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
            rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
        }
    } else {
        for (SizeType i = 0; i < NumUNodes; ++i) {
            rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
            rElementalDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
        }
    }

    // if (!IgnoreUndrained) {
        for (SizeType i=0; i < NumPNodes; ++i)
            rElementalDofList[Index++] = GetGeometry()[i].pGetDof( WATER_PRESSURE );
    // }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::GetDofList()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                          VectorType& rRightHandSideVector,
                          const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateLocalSystem()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    // const bool IgnoreUndrained = GetProperties()[IGNORE_UNDRAINED];
    // const SizeType ElementSize = (IgnoreUndrained ? NumUNodes * Dim : NumUNodes * Dim + NumPNodes);
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != ElementSize )
        rLeftHandSideMatrix.resize( ElementSize, ElementSize, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( ElementSize, ElementSize );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != ElementSize )
        rRightHandSideVector.resize( ElementSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( ElementSize );

    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix,
                 rRightHandSideVector,
                 rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateLocalSystem()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
                           const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateLeftHandSide()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    // const bool IgnoreUndrained = GetProperties()[IGNORE_UNDRAINED];
    // const SizeType ElementSize = (IgnoreUndrained ? NumUNodes * Dim : NumUNodes * Dim + NumPNodes);
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != ElementSize )
        rLeftHandSideMatrix.resize( ElementSize, ElementSize, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( ElementSize, ElementSize );

    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = false;
    VectorType tempRightHandSideVector;

    CalculateAll(rLeftHandSideMatrix,
                 tempRightHandSideVector,
                 rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateLeftHandSide()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateRightHandSide( VectorType& rRightHandSideVector,
                            const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateRightHandSide()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    // const bool IgnoreUndrained = GetProperties()[IGNORE_UNDRAINED];
    // const SizeType ElementSize = (IgnoreUndrained ? NumUNodes * Dim : NumUNodes * Dim + NumPNodes);
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    //Resetting the RHS
    if ( rRightHandSideVector.size() != ElementSize )
        rRightHandSideVector.resize( ElementSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( ElementSize );

    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp,
                 rRightHandSideVector,
                 rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateRightHandSide()") << std::endl;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateMassMatrix( MatrixType& rMassMatrix,
                         const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateMassMatrix()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType BlockElementSize = NumUNodes * Dim;
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( this->GetIntegrationMethod() );

    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(rGeom, this->GetProperties(), rCurrentProcessInfo);

    Matrix M = ZeroMatrix(BlockElementSize, BlockElementSize);

    //Defining shape functions and the determinant of the jacobian at all integration points

    //Loop over integration points
    Matrix Nu                = ZeroMatrix( Dim , NumUNodes * Dim );
    Matrix AuxDensityMatrix  = ZeroMatrix( Dim , NumUNodes * Dim );
    Matrix DensityMatrix     = ZeroMatrix( Dim, Dim );

    for ( SizeType GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint ) {
        //compute element kinematics (Np, gradNpT, |J|, B)
        this->CalculateKinematics(Variables,GPoint);

        //calculating weighting coefficient for integration
        Variables.IntegrationCoefficientInitialConfiguration =
            this->CalculateIntegrationCoefficient(IntegrationPoints,
                                                  GPoint,
                                                  Variables.detJInitialConfiguration);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->CalculateSoilDensity(Variables);

        //Setting the shape function matrix
        SizeType Index = 0;
        for (SizeType i = 0; i < NumUNodes; ++i) {
            for (SizeType iDim = 0; iDim < Dim; ++iDim) {
                Nu(iDim, Index++) = Variables.Nu(i);
            }
        }

        GeoElementUtilities::
            AssembleDensityMatrix(DensityMatrix, Variables.Density);

        noalias(AuxDensityMatrix) = prod(DensityMatrix, Nu);

        //Adding contribution to Mass matrix
        noalias(M) += prod(trans(Nu), AuxDensityMatrix) * Variables.IntegrationCoefficientInitialConfiguration;
    }

    //Distribute mass block matrix into the elemental matrix
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    // const bool IgnoreUndrained = GetProperties()[IGNORE_UNDRAINED];
    // const SizeType ElementSize = (IgnoreUndrained ? BlockElementSize : BlockElementSize + NumPNodes);
    const SizeType ElementSize = BlockElementSize + NumPNodes;


    if ( rMassMatrix.size1() != ElementSize || rMassMatrix.size2() != ElementSize)
        rMassMatrix.resize( ElementSize, ElementSize, false );
    noalias( rMassMatrix ) = ZeroMatrix( ElementSize, ElementSize );

    for (SizeType i = 0; i < NumUNodes; ++i) {
        SizeType Index_i = i * Dim;

        for (SizeType j = 0; j < NumUNodes; ++j) {
            SizeType Index_j = j * Dim;
            for (SizeType idim = 0; idim < Dim; ++idim) {
                for (SizeType jdim = 0; jdim < Dim; ++jdim) {
                    rMassMatrix(Index_i+idim,  Index_j+jdim) += M(Index_i+idim,  Index_j+jdim);
                }
            }
        }
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateMassMatrix()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateDampingMatrix(MatrixType& rDampingMatrix,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateDampingMatrix()") << std::endl;

    // Rayleigh Method (Damping Matrix = alpha*M + beta*K)

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    // const bool IgnoreUndrained = GetProperties()[IGNORE_UNDRAINED];
    // const SizeType ElementSize = (IgnoreUndrained ? NumUNodes * Dim : NumUNodes * Dim + NumPNodes);
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    // Compute Mass Matrix
    MatrixType MassMatrix(ElementSize, ElementSize);

    this->CalculateMassMatrix(MassMatrix,rCurrentProcessInfo);

    // Compute Stiffness matrix
    MatrixType StiffnessMatrix(ElementSize, ElementSize);

    this->CalculateMaterialStiffnessMatrix(StiffnessMatrix, rCurrentProcessInfo);

    // Compute Damping Matrix
    if ( rDampingMatrix.size1() != ElementSize )
        rDampingMatrix.resize( ElementSize, ElementSize, false );
    noalias( rDampingMatrix ) = ZeroMatrix( ElementSize, ElementSize );

    const PropertiesType& rProp = this->GetProperties();

    if (rProp.Has( RAYLEIGH_ALPHA ))
        noalias(rDampingMatrix) += rProp[RAYLEIGH_ALPHA] * MassMatrix;
    else
        noalias(rDampingMatrix) += rCurrentProcessInfo[RAYLEIGH_ALPHA] * MassMatrix;

    if (rProp.Has( RAYLEIGH_BETA ))
        noalias(rDampingMatrix) += rProp[RAYLEIGH_BETA] * StiffnessMatrix;
    else
        noalias(rDampingMatrix) += rCurrentProcessInfo[RAYLEIGH_BETA] * StiffnessMatrix;

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateDampingMatrix()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    EquationIdVector(EquationIdVectorType& rResult,
                     const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::EquationIdVector()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    // const bool IgnoreUndrained = GetProperties()[IGNORE_UNDRAINED];
    // const SizeType ElementSize = (IgnoreUndrained ? NumUNodes * Dim : NumUNodes * Dim + NumPNodes);
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    if ( rResult.size() != ElementSize )
        rResult.resize( ElementSize, false );

    SizeType Index = 0;

    if (Dim > 2) {
        for ( SizeType i = 0; i < NumUNodes; ++i ) {
            rResult[Index++] = rGeom[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[Index++] = rGeom[i].GetDof( DISPLACEMENT_Y ).EquationId();
            rResult[Index++] = rGeom[i].GetDof( DISPLACEMENT_Z ).EquationId();
        }
    } else {
        for ( SizeType i = 0; i < NumUNodes; ++i ) {
            rResult[Index++] = rGeom[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[Index++] = rGeom[i].GetDof( DISPLACEMENT_Y ).EquationId();
        }
    }

    // if (!IgnoreUndrained) {
        for ( SizeType i = 0; i < NumPNodes; ++i )
            rResult[Index++] = rGeom[i].GetDof( WATER_PRESSURE ).EquationId();
    // }



    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::EquationIdVector()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    GetFirstDerivativesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::GetFirstDerivativesVector()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    if ( rValues.size() != ElementSize )
        rValues.resize( ElementSize, false );

    SizeType Index = 0;

    if ( Dim > 2 ) {
        for ( SizeType i = 0; i < NumUNodes; ++i ) {
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( VELOCITY_X, Step );
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( VELOCITY_Y, Step );
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
        }
    }
    else {
        for ( SizeType i = 0; i < NumUNodes; ++i ) {
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( VELOCITY_X, Step );
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( VELOCITY_Y, Step );
        }
    }

    for ( SizeType i = 0; i < NumPNodes; ++i )
        rValues[Index++] = 0.0;

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::GetSecondDerivativesVector()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::GetSecondDerivativesVector( Vector& rValues, int Step ) const
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::GetSecondDerivativesVector()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ElementSize = NumUNodes * Dim + NumPNodes;

    if ( rValues.size() != ElementSize )
        rValues.resize( ElementSize, false );

    SizeType Index = 0;

    if ( Dim > 2 ) {
        for ( SizeType i = 0; i < NumUNodes; ++i ) {
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
        }
    }
    else {
        for ( SizeType i = 0; i < NumUNodes; ++i ) {
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );
        }

    }

    for ( SizeType i = 0; i < NumPNodes; ++i )
        rValues[Index++] = 0.0;

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::GetSecondDerivativesVector()") << std::endl;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::FinalizeSolutionStep()") << std::endl;

    //Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    //ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
        //compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables,GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables,ConstitutiveParameters);

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

    const bool IgnoreUndrained = GetProperties()[IGNORE_UNDRAINED];

    //Assign pressure values to the intermediate nodes for post-processing
    if (!IgnoreUndrained) AssignPressureToIntermediateNodes();

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::FinalizeSolutionStep()") << std::endl;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::AssignPressureToIntermediateNodes()
{
    //Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::AssignPressureToIntermediateNodes()") << std::endl;

    GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();

    switch (NumUNodes)
    {
        case 6: //2D T6P3
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            ThreadSafeNodeWrite(rGeom[3],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[4],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[5],WATER_PRESSURE, 0.5 * (p2 + p0) );
            break;
        }
        case 8: //2D Q8P4
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            ThreadSafeNodeWrite(rGeom[4],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[5],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[6],WATER_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[7],WATER_PRESSURE, 0.5 * (p3 + p0) );
            break;
        }
        case 9: //2D Q9P4
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            ThreadSafeNodeWrite(rGeom[4],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[5],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[6],WATER_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[7],WATER_PRESSURE, 0.5 * (p3 + p0) );
            ThreadSafeNodeWrite(rGeom[8],WATER_PRESSURE, 0.25 * (p0 + p1 + p2 + p3) );
            break;
        }
        case 10: //3D T10P4
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            ThreadSafeNodeWrite(rGeom[4],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[5],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[6],WATER_PRESSURE, 0.5 * (p2 + p0) );
            ThreadSafeNodeWrite(rGeom[7],WATER_PRESSURE, 0.5 * (p0 + p3) );
            ThreadSafeNodeWrite(rGeom[8],WATER_PRESSURE, 0.5 * (p1 + p3) );
            ThreadSafeNodeWrite(rGeom[9],WATER_PRESSURE, 0.5 * (p2 + p3) );
            break;
        }
        case 20: //3D H20P8
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p4 = rGeom[4].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p5 = rGeom[5].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p6 = rGeom[6].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p7 = rGeom[7].FastGetSolutionStepValue(WATER_PRESSURE);
            // edges -- bottom
            ThreadSafeNodeWrite(rGeom[8],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[9],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[10],WATER_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[11],WATER_PRESSURE, 0.5 * (p3 + p0) );
            // edges -- middle
            ThreadSafeNodeWrite(rGeom[12],WATER_PRESSURE, 0.5 * (p4 + p0) );
            ThreadSafeNodeWrite(rGeom[13],WATER_PRESSURE, 0.5 * (p5 + p1) );
            ThreadSafeNodeWrite(rGeom[14],WATER_PRESSURE, 0.5 * (p6 + p2) );
            ThreadSafeNodeWrite(rGeom[15],WATER_PRESSURE, 0.5 * (p7 + p3) );
            // edges -- top
            ThreadSafeNodeWrite(rGeom[16],WATER_PRESSURE, 0.5 * (p4 + p5) );
            ThreadSafeNodeWrite(rGeom[17],WATER_PRESSURE, 0.5 * (p5 + p6) );
            ThreadSafeNodeWrite(rGeom[18],WATER_PRESSURE, 0.5 * (p6 + p7) );
            ThreadSafeNodeWrite(rGeom[19],WATER_PRESSURE, 0.5 * (p7 + p0) );
            break;
        }
        case 27: //3D H27P8
        {
            const double p0 = rGeom[0].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p1 = rGeom[1].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p2 = rGeom[2].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p3 = rGeom[3].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p4 = rGeom[4].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p5 = rGeom[5].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p6 = rGeom[6].FastGetSolutionStepValue(WATER_PRESSURE);
            const double p7 = rGeom[7].FastGetSolutionStepValue(WATER_PRESSURE);
            // edges -- bottom
            ThreadSafeNodeWrite(rGeom[8],WATER_PRESSURE, 0.5 * (p0 + p1) );
            ThreadSafeNodeWrite(rGeom[9],WATER_PRESSURE, 0.5 * (p1 + p2) );
            ThreadSafeNodeWrite(rGeom[10],WATER_PRESSURE, 0.5 * (p2 + p3) );
            ThreadSafeNodeWrite(rGeom[11],WATER_PRESSURE, 0.5 * (p3 + p0) );
            // edges -- middle
            ThreadSafeNodeWrite(rGeom[12],WATER_PRESSURE, 0.5 * (p4 + p0) );
            ThreadSafeNodeWrite(rGeom[13],WATER_PRESSURE, 0.5 * (p5 + p1) );
            ThreadSafeNodeWrite(rGeom[14],WATER_PRESSURE, 0.5 * (p6 + p2) );
            ThreadSafeNodeWrite(rGeom[15],WATER_PRESSURE, 0.5 * (p7 + p3) );
            // edges -- top
            ThreadSafeNodeWrite(rGeom[16],WATER_PRESSURE, 0.5 * (p4 + p5) );
            ThreadSafeNodeWrite(rGeom[17],WATER_PRESSURE, 0.5 * (p5 + p6) );
            ThreadSafeNodeWrite(rGeom[18],WATER_PRESSURE, 0.5 * (p6 + p7) );
            ThreadSafeNodeWrite(rGeom[19],WATER_PRESSURE, 0.5 * (p7 + p0) );
            // face centers
            ThreadSafeNodeWrite(rGeom[20],WATER_PRESSURE, 0.25 * (p0 + p1 + p2 + p3) );
            ThreadSafeNodeWrite(rGeom[21],WATER_PRESSURE, 0.25 * (p0 + p1 + p4 + p5) );
            ThreadSafeNodeWrite(rGeom[22],WATER_PRESSURE, 0.25 * (p1 + p2 + p5 + p6) );
            ThreadSafeNodeWrite(rGeom[23],WATER_PRESSURE, 0.25 * (p2 + p3 + p6 + p7) );
            ThreadSafeNodeWrite(rGeom[24],WATER_PRESSURE, 0.25 * (p3 + p0 + p7 + p4) );
            ThreadSafeNodeWrite(rGeom[25],WATER_PRESSURE, 0.25 * (p4 + p5 + p6 + p7) );
            // element center
            ThreadSafeNodeWrite(rGeom[26],WATER_PRESSURE, 0.125 * (p0+p1+p2+p3+p4+p5+p6+p7) );
            break;
        }
        default:
            KRATOS_ERROR << "Unexpected geometry type for different order interpolation element" << this->Id() << std::endl;
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::AssignPressureToIntermediateNodes()") << std::endl;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    SetValuesOnIntegrationPoints( const Variable<double>& rVariable,
                                  const std::vector<double>& rValues,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::SetValuesOnIntegrationPoints()") << std::endl;

    for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint )
        mConstitutiveLawVector[GPoint]->SetValue( rVariable, rValues[GPoint], rCurrentProcessInfo );

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::SetValuesOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    SetValuesOnIntegrationPoints( const Variable<Vector>& rVariable,
                                  const std::vector<Vector>& rValues,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::1-SetValuesOnIntegrationPoints()") << std::endl;

    for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint )
        mConstitutiveLawVector[GPoint]->SetValue( rVariable, rValues[GPoint], rCurrentProcessInfo );

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::1-SetValuesOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    SetValuesOnIntegrationPoints( const Variable<Matrix>& rVariable,
                                  const std::vector<Matrix>& rValues,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::2-SetValuesOnIntegrationPoints()") << std::endl;

    for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint )
        mConstitutiveLawVector[GPoint]->SetValue( rVariable, rValues[GPoint], rCurrentProcessInfo );

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::2-SetValuesOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateOnIntegrationPoints( const Variable<int>& rVariable,
                                  std::vector<int>& rValues,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints()") << std::endl;

    const unsigned int& IntegrationPointsNumber = GetGeometry().IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rValues.size() != IntegrationPointsNumber )
        rValues.resize( IntegrationPointsNumber );

    for ( unsigned int i = 0; i < IntegrationPointsNumber; ++i )
        rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                  std::vector<ConstitutiveLaw::Pointer>& rValues,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::3-CalculateOnIntegrationPoints()") << std::endl;

    if (rVariable == CONSTITUTIVE_LAW) {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        for (unsigned int i=0; i<rValues.size(); ++i)
            rValues[i] = mConstitutiveLawVector[i];
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::3-CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateOnIntegrationPoints( const Variable<double>& rVariable,
                                  std::vector<double>& rOutput,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const unsigned int& IntegrationPointsNumber = rGeom.IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != IntegrationPointsNumber )
        rOutput.resize( IntegrationPointsNumber, false );

    if ( rVariable == VON_MISES_STRESS ) {
        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            ComparisonUtilities EquivalentStress;
            rOutput[GPoint] =  EquivalentStress.CalculateVonMises(mStressVector[GPoint]);
        }
    } else if (rVariable == DEGREE_OF_SATURATION ||
             rVariable == EFFECTIVE_SATURATION ||
             rVariable == BISHOP_COEFICIENT ||
             rVariable == DERIVATIVE_OF_SATURATION ||
             rVariable == RELATIVE_PERMEABILITY )
    {
        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,
                                         rCurrentProcessInfo);

        // create general parametes of retention law
        RetentionLaw::Parameters RetentionParameters(rGeom, GetProperties(), rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mRetentionLawVector.size(); ++GPoint ) {
            //Compute Np, GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);

            Variables.FluidPressure = CalculateFluidPressure(Variables, GPoint);
            SetRetentionParameters(Variables, RetentionParameters);

            if (rVariable == DEGREE_OF_SATURATION)     rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateSaturation(RetentionParameters);
            if (rVariable == EFFECTIVE_SATURATION)     rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateEffectiveSaturation(RetentionParameters);
            if (rVariable == BISHOP_COEFICIENT)        rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(RetentionParameters);
            if (rVariable == DERIVATIVE_OF_SATURATION) rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(RetentionParameters);
            if (rVariable == RELATIVE_PERMEABILITY )   rOutput[GPoint] = mRetentionLawVector[GPoint]->CalculateRelativePermeability(RetentionParameters);
        }
    } else if (rVariable == HYDRAULIC_HEAD) {
        const double NumericalLimit = std::numeric_limits<double>::epsilon();
        const PropertiesType& rProp = this->GetProperties();
        const unsigned int NumGPoints = rGeom.IntegrationPointsNumber( this->GetIntegrationMethod() );

        //Defining the shape functions, the jacobian and the shape functions local gradients Containers
        const Matrix& NContainer = rGeom.ShapeFunctionsValues( this->GetIntegrationMethod() );
        const SizeType NumUNodes = rGeom.PointsNumber();

        //Defining necessary variables
        Vector NodalHydraulicHead = ZeroVector(NumUNodes);
        for (unsigned int node=0; node < NumUNodes; ++node) {
            Vector NodeVolumeAcceleration(3);
            noalias(NodeVolumeAcceleration) = rGeom[node].FastGetSolutionStepValue(VOLUME_ACCELERATION, 0);
            const double g = norm_2(NodeVolumeAcceleration);
            if (g > NumericalLimit) {
                const double FluidWeight = g * rProp[DENSITY_WATER];

                Vector NodeCoordinates(3);
                noalias(NodeCoordinates) = rGeom[node].Coordinates();
                Vector NodeVolumeAccelerationUnitVector(3);
                noalias(NodeVolumeAccelerationUnitVector) = NodeVolumeAcceleration / g;

                const double WaterPressure = rGeom[node].FastGetSolutionStepValue(WATER_PRESSURE);
                NodalHydraulicHead[node] =- inner_prod(NodeCoordinates, NodeVolumeAccelerationUnitVector)
                                          - PORE_PRESSURE_SIGN_FACTOR  * WaterPressure / FluidWeight;
            }
            else {
                NodalHydraulicHead[node] = 0.0;
            }
        }

        if ( rOutput.size() != NumGPoints )
            rOutput.resize(NumGPoints);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ ) {
            double HydraulicHead = 0.0;
            for (unsigned int node = 0; node < NumUNodes; ++node)
                HydraulicHead += NContainer(GPoint, node) * NodalHydraulicHead[node];

            rOutput[GPoint] = HydraulicHead;
        }
    } else {
        for ( unsigned int i = 0; i < IntegrationPointsNumber; ++i )
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rOutput[i] );
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,
                                 std::vector<array_1d<double,3>>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::1-CalculateOnIntegrationPoints(array_1d)") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const unsigned int& IntegrationPointsNumber = rGeom.IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != IntegrationPointsNumber )
        rOutput.resize( IntegrationPointsNumber );

    if ( rVariable == FLUID_FLUX_VECTOR ) {
        //Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,rCurrentProcessInfo);

        // create general parametes of retention law
        RetentionLaw::Parameters RetentionParameters(rGeom, GetProperties(), rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,GPoint);

            //Compute FluidFlux vector q [m/s]
            const SizeType Dim = rGeom.WorkingSpaceDimension();
            const SizeType NumUNodes = rGeom.PointsNumber();

            Vector BodyAcceleration = ZeroVector(Dim);
            SizeType Index = 0;
            for (SizeType i = 0; i < NumUNodes; ++i) {
                for (unsigned int idim = 0; idim < Dim; ++idim)
                    BodyAcceleration[idim] += Variables.Nu[i]*Variables.BodyAcceleration[Index++];
            }

            CalculateFluidPressure(Variables, GPoint);
            SetRetentionParameters(Variables, RetentionParameters);

            const double RelativePermeability = 
                mRetentionLawVector[GPoint]->CalculateRelativePermeability(RetentionParameters);

            Vector GradPressureTerm(Dim);
            noalias(GradPressureTerm)  =  prod(trans(Variables.DNp_DX), Variables.PressureVector);
            noalias(GradPressureTerm) +=  PORE_PRESSURE_SIGN_FACTOR
                                        * GetProperties()[DENSITY_WATER]
                                        * BodyAcceleration;

            Vector AuxFluidFlux = ZeroVector(Dim);
            AuxFluidFlux =   PORE_PRESSURE_SIGN_FACTOR
                           * Variables.DynamicViscosityInverse
                           * RelativePermeability
                           * prod(Variables.IntrinsicPermeability, GradPressureTerm );

            Vector FluidFlux = ZeroVector(3);
            for (unsigned int idim = 0; idim < Dim; ++idim)
                FluidFlux[idim] = AuxFluidFlux[idim];

            if ( rOutput[GPoint].size() != 3 )
                rOutput[GPoint].resize( 3, false );

            rOutput[GPoint] = FluidFlux;
        }
    }

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::1-CalculateOnIntegrationPoints(array_1d)") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                 std::vector<Vector>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::1-CalculateOnIntegrationPoints()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const unsigned int& IntegrationPointsNumber = rGeom.IntegrationPointsNumber( this->GetIntegrationMethod() );

    if ( rOutput.size() != IntegrationPointsNumber )
        rOutput.resize( IntegrationPointsNumber );

    if ( rVariable == CAUCHY_STRESS_VECTOR ) {
        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            if ( rOutput[GPoint].size() != mStressVector[GPoint].size() )
                rOutput[GPoint].resize( mStressVector[GPoint].size(), false );

            rOutput[GPoint] = mStressVector[GPoint];
        }
    } else if ( rVariable == ENGINEERING_STRAIN_VECTOR ) {
        //Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            noalias(Variables.Nu) = row(Variables.NuContainer, GPoint);

            Matrix J0,InvJ0;
            this->CalculateDerivativesOnInitialConfiguration(Variables.detJInitialConfiguration,
                                                             J0,
                                                             InvJ0,
                                                             Variables.DNu_DXInitialConfiguration,
                                                             GPoint);

            // Calculating operator B
            this->CalculateBMatrix( Variables.B, Variables.DNu_DXInitialConfiguration, Variables.Nu);

            //Compute infinitessimal strain
            this->CalculateCauchyStrain( Variables );

            if ( rOutput[GPoint].size() != Variables.StrainVector.size() )
                rOutput[GPoint].resize( Variables.StrainVector.size(), false );

            rOutput[GPoint] = Variables.StrainVector;
        }
    } else if ( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ) {
        //Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,GPoint);

            //Compute infinitessimal strain
            this->CalculateStrain(Variables, GPoint);

            if ( rOutput[GPoint].size() != Variables.StrainVector.size() )
                rOutput[GPoint].resize( Variables.StrainVector.size(), false );

            rOutput[GPoint] = Variables.StrainVector;
        }
    } else if ( rVariable == TOTAL_STRESS_VECTOR ) {
        //Definition of variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);

        //Create constitutive law parameters:
        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,GetProperties(),rCurrentProcessInfo);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        const PropertiesType& rProp = this->GetProperties();

        const SizeType VoigtSize = mStressVector[0].size();
        Vector VoigtVector = ZeroVector(VoigtSize);

        for (unsigned int i=0; i < rGeom.WorkingSpaceDimension(); ++i) VoigtVector[i] = 1.0;

        // create general parametes of retention law
        RetentionLaw::Parameters RetentionParameters(rGeom, GetProperties(), rCurrentProcessInfo);

        const bool hasBiotCoefficient = rProp.Has(BIOT_COEFFICIENT);

        Vector TotalStressVector(mStressVector[0].size());

        //Loop over integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,GPoint);

            //Compute infinitessimal strain
            this->CalculateStrain(Variables, GPoint);

            //set gauss points variables to constitutivelaw parameters
            this->SetConstitutiveParameters(Variables,ConstitutiveParameters);

            //compute constitutive tensor and/or stresses
            noalias(Variables.StressVector) = mStressVector[GPoint];
            ConstitutiveParameters.SetStressVector(Variables.StressVector);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            Variables.BiotCoefficient = CalculateBiotCoefficient(Variables, hasBiotCoefficient);

            this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

            noalias(TotalStressVector) = mStressVector[GPoint];
            noalias(TotalStressVector) +=  PORE_PRESSURE_SIGN_FACTOR
                                         * Variables.BiotCoefficient
                                         * Variables.BishopCoefficient
                                         * Variables.FluidPressure
                                         * VoigtVector;

            if ( rOutput[GPoint].size() != TotalStressVector.size() )
                rOutput[GPoint].resize( TotalStressVector.size(), false );

            rOutput[GPoint] = TotalStressVector;
        }
    } else {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable , rOutput[i] );
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::1-CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable,
                                 std::vector< Matrix >& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::-CalculateOnIntegrationPoints()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const unsigned int& IntegrationPointsNumber = rGeom.IntegrationPointsNumber( this->GetIntegrationMethod() );
    const unsigned int dimension = rGeom.WorkingSpaceDimension();

    if ( rOutput.size() != IntegrationPointsNumber )
        rOutput.resize( IntegrationPointsNumber );

    if ( rVariable == CAUCHY_STRESS_TENSOR ) {
        std::vector<Vector> StressVector;

        this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            if ( rOutput[GPoint].size2() != dimension )
                rOutput[GPoint].resize( dimension, dimension, false );

            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector[GPoint]);
        }
    } else if ( rVariable == ENGINEERING_STRAIN_TENSOR ) {
        std::vector<Vector> StrainVector;

        CalculateOnIntegrationPoints( ENGINEERING_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            if ( rOutput[GPoint].size2() != dimension )
                rOutput[GPoint].resize( dimension, dimension, false );

            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector[GPoint]);
        }

    } else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ) {
        std::vector<Vector> StrainVector;

        CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint ) {
            if ( rOutput[GPoint].size2() != dimension )
                rOutput[GPoint].resize( dimension, dimension, false );

            rOutput[GPoint] = MathUtils<double>::StrainVectorToTensor(StrainVector[GPoint]);
        }
    } else if (rVariable == TOTAL_STRESS_TENSOR) {
        std::vector<Vector> StressVector;

        this->CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

        //loop integration points
        for (unsigned int GPoint = 0; GPoint < mConstitutiveLawVector.size(); ++GPoint) {
            if (rOutput[GPoint].size2() != dimension)
                rOutput[GPoint].resize(dimension, dimension, false);

            rOutput[GPoint] = MathUtils<double>::StressVectorToTensor(StressVector[GPoint]);
        }

    } else {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable , rOutput[i] );
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::-CalculateOnIntegrationPoints()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& rCurrentProcessInfo,
                  bool CalculateStiffnessMatrixFlag,
                  bool CalculateResidualVectorFlag )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAll") << this->Id() << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const PropertiesType& rProp = this->GetProperties();

    //Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,rCurrentProcessInfo);

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    // Stiffness matrix is always needed to calculate Biot coefficient
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (CalculateResidualVectorFlag)  ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(rGeom, rProp, rCurrentProcessInfo);

    //Loop over integration points
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( this->GetIntegrationMethod() );

    const bool hasBiotCoefficient = rProp.Has(BIOT_COEFFICIENT);

    for ( unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint ) {
        //compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables, GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables,ConstitutiveParameters);

        //compute constitutive tensor and/or stresses
        ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        this->InitializeBiotCoefficients(Variables, hasBiotCoefficient);

        //calculating weighting coefficient for integration
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
        if (CalculateResidualVectorFlag) this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAll") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateMaterialStiffnessMatrix( MatrixType& rStiffnessMatrix,
                                      const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateMaterialStiffnessMatrix") << std::endl;

    const GeometryType& rGeom = GetGeometry();

    //Definition of variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,rCurrentProcessInfo);

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,GetProperties(),rCurrentProcessInfo);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //Loop over integration points
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( this->GetIntegrationMethod() );

    for ( unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint ) {
        //compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables,GPoint);

        //Compute infinitessimal strain
        this->CalculateStrain(Variables, GPoint);

        //set gauss points variables to constitutivelaw parameters
        this->SetConstitutiveParameters(Variables,ConstitutiveParameters);

        //compute constitutive tensor and/or stresses
        noalias(Variables.StressVector) = mStressVector[GPoint];
        ConstitutiveParameters.SetStressVector(Variables.StressVector);
        mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //calculating weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints,
                                                  GPoint,
                                                  Variables.detJ);

        //Contributions of material stiffness to the left hand side
        this->CalculateAndAddStiffnessMatrix(rStiffnessMatrix, Variables);

    }

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateMaterialStiffnessMatrix") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double SmallStrainUPwDiffOrderElement::
    CalculateBulkModulus(const Matrix &ConstitutiveMatrix) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateBulkModulus") << std::endl;

    const int IndexG = ConstitutiveMatrix.size1() - 1;
    const double M = ConstitutiveMatrix(0, 0);
    const double G = ConstitutiveMatrix(IndexG, IndexG);
    const double BulkModulus = M - (4.0/3.0)*G;

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateBulkModulus") << BulkModulus << std::endl;

    return BulkModulus;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    InitializeElementVariables( ElementVariables& rVariables,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::InitializeElementVariables") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType NumGPoints = rGeom.IntegrationPointsNumber( this->GetIntegrationMethod() );
    const SizeType Dim = rGeom.WorkingSpaceDimension();

    //Variables at all integration points
    (rVariables.NuContainer).resize(NumGPoints,NumUNodes, false);
    rVariables.NuContainer = rGeom.ShapeFunctionsValues( this->GetIntegrationMethod() );

    (rVariables.NpContainer).resize(NumGPoints,NumPNodes, false);
    rVariables.NpContainer = mpPressureGeometry->ShapeFunctionsValues( this->GetIntegrationMethod() );

    (rVariables.Nu).resize(NumUNodes, false);
    (rVariables.Np).resize(NumPNodes, false);

    (rVariables.DNu_DXContainer).resize(NumGPoints, false);
    for (SizeType i = 0; i<NumGPoints; ++i)
        ((rVariables.DNu_DXContainer)[i]).resize(NumUNodes,Dim, false);
    (rVariables.DNu_DX).resize(NumUNodes,Dim,false);
    (rVariables.DNu_DXInitialConfiguration).resize(NumUNodes,Dim,false);
    (rVariables.detJuContainer).resize(NumGPoints,false);

    try {
        rGeom.ShapeFunctionsIntegrationPointsGradients( rVariables.DNu_DXContainer,
                                                        rVariables.detJuContainer,
                                                        this->GetIntegrationMethod() );
    } catch (Kratos::Exception& e) {
        KRATOS_INFO("Original error message") << e.what() << std::endl;
#ifdef KRATOS_COMPILED_IN_WINDOWS
        KRATOS_INFO("Error in calculation of dNu/dx. Most probably the element is distorted. Element ID: ") << this->Id() << std::endl;
#endif
        KRATOS_ERROR << "In calculation of dNu/dx. Most probably the element is distorted. Element ID: " << this->Id() << std::endl;
    }

    (rVariables.DNp_DXContainer).resize(NumGPoints,false);
    for (SizeType i = 0; i<NumGPoints; ++i)
        ((rVariables.DNp_DXContainer)[i]).resize(NumPNodes,Dim,false);
    (rVariables.DNp_DX).resize(NumPNodes,Dim,false);
    Vector detJpContainer = ZeroVector(NumGPoints);

    try {
        mpPressureGeometry->ShapeFunctionsIntegrationPointsGradients( rVariables.DNp_DXContainer,
                                                                    detJpContainer,
                                                                    this->GetIntegrationMethod());
    } catch (Kratos::Exception& e) {
        KRATOS_INFO("Original error message") << e.what() << std::endl;
#ifdef KRATOS_COMPILED_IN_WINDOWS
        KRATOS_INFO("Error in calculation of dNp/dx. Most probably the element is distorted. Element ID: ") << this->Id() << std::endl;
#endif
        KRATOS_ERROR << "In calculation of dNp/dx. Most probably the element is distorted. Element ID: " << this->Id() << std::endl;
    }

    //Variables computed at each integration point
    const SizeType VoigtSize  = ( Dim == 3 ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN);

    (rVariables.B).resize(VoigtSize, NumUNodes * Dim, false);
    noalias(rVariables.B) = ZeroMatrix( VoigtSize, NumUNodes * Dim );

    (rVariables.StrainVector).resize(VoigtSize, false);
    (rVariables.ConstitutiveMatrix).resize(VoigtSize, VoigtSize, false);

    (rVariables.StressVector).resize(VoigtSize, false);

    //Needed parameters for consistency with the general constitutive law
    rVariables.detF = 1.0;
    (rVariables.F).resize(Dim, Dim, false);
    noalias(rVariables.F) = identity_matrix<double>(Dim);

    //Nodal variables
    this->InitializeNodalVariables(rVariables);

    //Properties variables
    this->InitializeProperties(rVariables);

    //ProcessInfo variables
    rVariables.VelocityCoefficient = rCurrentProcessInfo[VELOCITY_COEFFICIENT];
    rVariables.DtPressureCoefficient = rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT];

    // Retention law
    rVariables.FluidPressure = 0.0;
    rVariables.DegreeOfSaturation = 1.0;
    rVariables.DerivativeOfSaturation = 0.0;
    rVariables.RelativePermeability = 1.0;
    rVariables.BishopCoefficient = 1.0;

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::InitializeElementVariables") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::InitializeNodalVariables( ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::InitializeNodalVariables") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    Vector BodyAccelerationAux    = ZeroVector(3);
    (rVariables.BodyAcceleration).resize(NumUNodes * Dim,false);
    (rVariables.DisplacementVector).resize(NumUNodes * Dim,false);
    (rVariables.VelocityVector).resize(NumUNodes * Dim,false);

    if (Dim > 2) {
        for (SizeType i=0; i<NumUNodes; ++i) {
            SizeType Local_i = i * Dim;
            BodyAccelerationAux = rGeom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

            rVariables.BodyAcceleration[Local_i]   = BodyAccelerationAux[0];
            rVariables.DisplacementVector[Local_i] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_X);
            rVariables.VelocityVector[Local_i]     = rGeom[i].FastGetSolutionStepValue(VELOCITY_X);

            rVariables.BodyAcceleration[Local_i+1]   = BodyAccelerationAux[1];
            rVariables.DisplacementVector[Local_i+1] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
            rVariables.VelocityVector[Local_i+1]     = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);

            rVariables.BodyAcceleration[Local_i+2]   = BodyAccelerationAux[2];
            rVariables.DisplacementVector[Local_i+2] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
            rVariables.VelocityVector[Local_i+2]     = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z);
        }
    } else {
        for (SizeType i=0; i<NumUNodes; ++i) {
            SizeType Local_i = i * Dim;
            BodyAccelerationAux = rGeom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

            rVariables.BodyAcceleration[Local_i]   = BodyAccelerationAux[0];
            rVariables.DisplacementVector[Local_i] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_X);
            rVariables.VelocityVector[Local_i]     = rGeom[i].FastGetSolutionStepValue(VELOCITY_X);

            rVariables.BodyAcceleration[Local_i+1]   = BodyAccelerationAux[1];
            rVariables.DisplacementVector[Local_i+1] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
            rVariables.VelocityVector[Local_i+1]     = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);
        }
    }

    (rVariables.PressureVector).resize(NumPNodes,false);
    (rVariables.PressureDtVector).resize(NumPNodes,false);
    (rVariables.DeltaPressureVector).resize(NumPNodes,false);
    for (SizeType i=0; i<NumPNodes; ++i) {
        rVariables.PressureVector[i]   = rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.PressureDtVector[i] = rGeom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
        rVariables.DeltaPressureVector[i] =   rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE)
                                            - rGeom[i].FastGetSolutionStepValue(WATER_PRESSURE, 1);
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::InitializeNodalVariables") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
double SmallStrainUPwDiffOrderElement::
    CalculateBiotCoefficient( const ElementVariables& rVariables,
                              const bool &hasBiotCoefficient) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateBiotCoefficient()") << std::endl;

    const PropertiesType& rProp = this->GetProperties();

    //Properties variables
    if (hasBiotCoefficient) {
        return rProp[BIOT_COEFFICIENT];
    } else {
        // calculate Bulk modulus from stiffness matrix
        const double BulkModulus = CalculateBulkModulus(rVariables.ConstitutiveMatrix);
        return 1.0 - BulkModulus / rProp[BULK_MODULUS_SOLID];
    }

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateBiotCoefficient()") << std::endl;
    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    InitializeBiotCoefficients( ElementVariables &rVariables,
                                const bool &hasBiotCoefficient )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::InitializeBiotCoefficients") << std::endl;

    const PropertiesType& rProp = this->GetProperties();

    rVariables.BiotCoefficient = CalculateBiotCoefficient(rVariables, hasBiotCoefficient);

    const bool IgnoreUndrained = rProp[IGNORE_UNDRAINED];

    if (!IgnoreUndrained) {
        rVariables.BiotModulusInverse = (rVariables.BiotCoefficient - rProp[POROSITY]) / rProp[BULK_MODULUS_SOLID] 
                                    + rProp[POROSITY]/rProp[BULK_MODULUS_FLUID];
    } else {
        rVariables.BiotModulusInverse = (rVariables.BiotCoefficient - rProp[POROSITY]) / rProp[BULK_MODULUS_SOLID] 
                                       + rProp[POROSITY]/TINY;
    }

    rVariables.BiotModulusInverse *= rVariables.DegreeOfSaturation;
    rVariables.BiotModulusInverse -= rVariables.DerivativeOfSaturation*rProp[POROSITY];

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::InitializeBiotCoefficients") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::InitializeProperties( ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::InitializeProperties") << std::endl;

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const PropertiesType& rProp = this->GetProperties();

    rVariables.IgnoreUndrained = rProp[IGNORE_UNDRAINED];
    rVariables.UseHenckyStrain = false;
    if (rProp.Has(USE_HENCKY_STRAIN))
        rVariables.UseHenckyStrain = rProp[USE_HENCKY_STRAIN];

    rVariables.ConsiderGeometricStiffness = false;
    if (rProp.Has(CONSIDER_GEOMETRIC_STIFFNESS))
        rVariables.ConsiderGeometricStiffness = rProp[CONSIDER_GEOMETRIC_STIFFNESS];

    rVariables.DynamicViscosityInverse = 1.0 / rProp[DYNAMIC_VISCOSITY];
    //Setting the intrinsic permeability matrix
    (rVariables.IntrinsicPermeability).resize(dimension,dimension,false);
    rVariables.IntrinsicPermeability(0,0) = rProp[PERMEABILITY_XX];
    rVariables.IntrinsicPermeability(1,1) = rProp[PERMEABILITY_YY];
    rVariables.IntrinsicPermeability(0,1) = rProp[PERMEABILITY_XY];
    rVariables.IntrinsicPermeability(1,0) = rVariables.IntrinsicPermeability(0,1);

    if (dimension==3) {
        rVariables.IntrinsicPermeability(2,2) = rProp[PERMEABILITY_ZZ];
        rVariables.IntrinsicPermeability(2,0) = rProp[PERMEABILITY_ZX];
        rVariables.IntrinsicPermeability(1,2) = rProp[PERMEABILITY_YZ];
        rVariables.IntrinsicPermeability(0,2) = rVariables.IntrinsicPermeability(2,0);
        rVariables.IntrinsicPermeability(2,1) = rVariables.IntrinsicPermeability(1,2);
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::InitializeProperties") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateKinematics( ElementVariables& rVariables,
                         const unsigned int &GPoint )

{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateKinematics") << std::endl;

    //Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Nu) = row(rVariables.NuContainer, GPoint);
    noalias(rVariables.Np) = row(rVariables.NpContainer, GPoint);

    noalias(rVariables.DNu_DX) = rVariables.DNu_DXContainer[GPoint];
    noalias(rVariables.DNp_DX) = rVariables.DNp_DXContainer[GPoint];

    //Compute the deformation matrix B
    this->CalculateBMatrix(rVariables.B, rVariables.DNu_DX, rVariables.Nu);

    rVariables.detJ = rVariables.detJuContainer[GPoint];

    Matrix J0,InvJ0;
    this->CalculateDerivativesOnInitialConfiguration(rVariables.detJInitialConfiguration,
                                                     J0,
                                                     InvJ0,
                                                     rVariables.DNu_DXInitialConfiguration,
                                                     GPoint);


    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateKinematics") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateDerivativesOnInitialConfiguration(double& detJ,
                                               Matrix& J0,
                                               Matrix& InvJ0,
                                               Matrix& DNu_DX0,
                                               const IndexType& GPoint) const
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateDerivativesOnInitialConfiguration()") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( this->GetIntegrationMethod() );

    GeometryUtils::JacobianOnInitialConfiguration(rGeom, IntegrationPoints[GPoint], J0);
    const Matrix& DN_De = rGeom.ShapeFunctionsLocalGradients(this->GetIntegrationMethod())[GPoint];
    MathUtils<double>::InvertMatrix( J0, InvJ0, detJ );
    GeometryUtils::ShapeFunctionsGradients(DN_De, InvJ0, DNu_DX0);

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateDerivativesOnInitialConfiguration()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateBMatrix(Matrix& rB,
                     const Matrix& DNp_DX,
                     const Vector& Np)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateBMatrix()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();

    unsigned int index;

    if (Dim > 2) {
        for ( unsigned int i = 0; i < NumUNodes; ++i ) {
            index = Dim * i;

            rB( INDEX_3D_XX, index + INDEX_X ) = DNp_DX( i, INDEX_X );
            rB( INDEX_3D_YY, index + INDEX_Y ) = DNp_DX( i, INDEX_Y );
            rB( INDEX_3D_ZZ, index + INDEX_Z ) = DNp_DX( i, INDEX_Z );
            rB( INDEX_3D_XY, index + INDEX_X ) = DNp_DX( i, INDEX_Y );
            rB( INDEX_3D_XY, index + INDEX_Y ) = DNp_DX( i, INDEX_X );
            rB( INDEX_3D_YZ, index + INDEX_Y ) = DNp_DX( i, INDEX_Z );
            rB( INDEX_3D_YZ, index + INDEX_Z ) = DNp_DX( i, INDEX_Y );
            rB( INDEX_3D_XZ, index + INDEX_X ) = DNp_DX( i, INDEX_Z );
            rB( INDEX_3D_XZ, index + INDEX_Z ) = DNp_DX( i, INDEX_X );
        }
    } else {
        // 2D plane strain
        for ( unsigned int i = 0; i < NumUNodes; ++i ) {
            index = Dim * i;

            rB( INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X ) = DNp_DX( i, INDEX_X );
            rB( INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y ) = DNp_DX( i, INDEX_Y );
            rB( INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X ) = DNp_DX( i, INDEX_Y );
            rB( INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y ) = DNp_DX( i, INDEX_X );
        }
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateBMatrix()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    SetConstitutiveParameters(ElementVariables& rVariables,
                          ConstitutiveLaw::Parameters& rConstitutiveParameters)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::SetConstitutiveParameters") << std::endl;

    rConstitutiveParameters.SetStrainVector(rVariables.StrainVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);

    //Needed parameters for consistency with the general constitutive law
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rVariables.DNu_DX);
    rConstitutiveParameters.SetShapeFunctionsValues(rVariables.Nu);

    rConstitutiveParameters.SetDeterminantF(rVariables.detF);
    rConstitutiveParameters.SetDeformationGradientF(rVariables.F);

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::SetConstitutiveParameters") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
double SmallStrainUPwDiffOrderElement::
    CalculateIntegrationCoefficient(const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
                                    const IndexType& GPoint,
                                    const double& detJ)
{
    return IntegrationPoints[GPoint].Weight() * detJ;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                       ElementVariables& rVariables)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddLHS") << std::endl;

    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix,rVariables);
    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix,rVariables);
        this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddLHS") << std::endl;

    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddStiffnessMatrix( MatrixType& rLeftHandSideMatrix,
                                    ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessMatrix") << std::endl;

    Matrix StiffnessMatrix =  prod( trans(rVariables.B), Matrix(prod(rVariables.ConstitutiveMatrix, rVariables.B)) )
                            * rVariables.IntegrationCoefficient;

    //Distribute stiffness block matrix into the elemental matrix
    this->AssembleUBlockMatrix(rLeftHandSideMatrix,StiffnessMatrix);

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessMatrix") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    AssembleUBlockMatrix(Matrix &rLeftHandSideMatrix,
                        const Matrix &StiffnessMatrix) const
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::AssembleUBlockMatrix") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    SizeType Index_i, Index_j;

    for (SizeType i = 0; i < NumUNodes; ++i) {
        Index_i = i * Dim;

        for (SizeType j = 0; j < NumUNodes; ++j) {
            Index_j = j * Dim;

            for (unsigned int idim = 0; idim < Dim; ++idim) {
                for (unsigned int jdim = 0; jdim < Dim; ++jdim) {
                    rLeftHandSideMatrix(Index_i+idim, Index_j+jdim) += StiffnessMatrix(Index_i+idim, Index_j+jdim);
                }
            }
        }
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::AssembleUBlockMatrix") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddCouplingMatrix( MatrixType& rLeftHandSideMatrix,
                                   ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingMatrix") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType VoigtSize  = ( Dim == 3 ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN);

    Vector VoigtVector = ZeroVector(VoigtSize);

    for (unsigned int i=0; i < Dim; ++i) VoigtVector[i] = 1.0;

    Matrix CouplingMatrix =   PORE_PRESSURE_SIGN_FACTOR
                            * rVariables.BiotCoefficient
                            * rVariables.BishopCoefficient
                            * prod( trans(rVariables.B), Matrix(outer_prod(VoigtVector,rVariables.Np)) )
                            * rVariables.IntegrationCoefficient;

    //Distribute coupling block matrix into the elemental matrix
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for (SizeType i = 0; i<NumUNodes; ++i) {
        SizeType Index_i = i * Dim;
        for (SizeType j = 0; j<NumPNodes; ++j) {
            for (unsigned int idim = 0; idim < Dim; ++idim)
                rLeftHandSideMatrix(Index_i+idim,NumUNodes*Dim+j) += CouplingMatrix(Index_i+idim,j);
        }
    }

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        Matrix CouplingMatrixT =   PORE_PRESSURE_SIGN_FACTOR 
                                * SaturationCoefficient
                                * rVariables.VelocityCoefficient
                                * trans(CouplingMatrix);

        //Distribute transposed coupling block matrix into the elemental matrix

        for (SizeType i = 0; i<NumPNodes; ++i) {
            for (SizeType j = 0; j<NumUNodes; ++j) {
                SizeType Index_j = j * Dim;
                for (unsigned int idim = 0; idim < Dim; ++idim)
                    rLeftHandSideMatrix(NumUNodes*Dim+i, Index_j+idim) += CouplingMatrixT(i, Index_j+idim);
            }
        }
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingMatrix") << std::endl;
    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddCompressibilityMatrix( MatrixType& rLeftHandSideMatrix,
                                          ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityMatrix") << std::endl;

    Matrix CompressibilityMatrix = - PORE_PRESSURE_SIGN_FACTOR 
                                   * rVariables.DtPressureCoefficient
                                   * rVariables.BiotModulusInverse 
                                   * outer_prod(rVariables.Np,rVariables.Np)
                                   * rVariables.IntegrationCoefficient;

    //Distribute compressibility block matrix into the elemental matrix
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for (SizeType i = 0; i < NumPNodes; ++i) {
        for (SizeType j=0; j < NumPNodes; ++j) {
            rLeftHandSideMatrix(NumUNodes*Dim+i,NumUNodes*Dim+j) += CompressibilityMatrix(i,j);
        }
    }
    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityMatrix") << std::endl;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddPermeabilityMatrix( MatrixType& rLeftHandSideMatrix,
                                       ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddPermeabilityMatrix") << std::endl;

    Matrix PermeabilityMatrix = - PORE_PRESSURE_SIGN_FACTOR
                                * rVariables.DynamicViscosityInverse
                                * rVariables.RelativePermeability
                                * prod( rVariables.DNp_DX, Matrix( prod(rVariables.IntrinsicPermeability, trans(rVariables.DNp_DX)) ) )
                                * rVariables.IntegrationCoefficient ;

    //Distribute permeability block matrix into the elemental matrix
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for (SizeType i = 0; i < NumPNodes; ++i) {
        for (SizeType j=0; j < NumPNodes; ++j) {
            rLeftHandSideMatrix(NumUNodes*Dim+i,NumUNodes*Dim+j) += PermeabilityMatrix(i,j);
        }
    }
    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddPermeabilityMatrix") << std::endl;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddRHS( VectorType& rRightHandSideVector,
                        ElementVariables& rVariables,
                        unsigned int GPoint )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddRHS") << std::endl;

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables, GPoint);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    if (!rVariables.IgnoreUndrained) {
        this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

        this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

        this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddRHS") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddStiffnessForce( VectorType& rRightHandSideVector,
                                   ElementVariables& rVariables,
                                   unsigned int GPoint )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessForce") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();

    Vector StiffnessForce =  prod(trans(rVariables.B), mStressVector[GPoint]) 
                           * rVariables.IntegrationCoefficient;

    //Distribute stiffness block vector into the elemental vector
    const SizeType NumUNodes = rGeom.PointsNumber();

    for (SizeType i = 0; i < NumUNodes; ++i) {
        SizeType Index = i * Dim;
        for (SizeType idim=0; idim < Dim; ++idim) {
            rRightHandSideVector[Index+idim] -= StiffnessForce[Index+idim];
        }
    }
    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddStiffnessForce") << std::endl;

    KRATOS_CATCH( "" )

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddMixBodyForce( VectorType& rRightHandSideVector,
                                 ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddMixBodyForce") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();

    this->CalculateSoilDensity(rVariables);

    Vector BodyAcceleration = ZeroVector(Dim);
    SizeType Index = 0;
    for (SizeType i = 0; i < NumUNodes; ++i) {
        for (SizeType idim=0; idim < Dim; ++idim) {
            BodyAcceleration[idim] += rVariables.Nu[i]*rVariables.BodyAcceleration[Index++];
        }
    }

    for (SizeType i=0; i < NumUNodes; ++i) {
        Index = i * Dim;
        for (SizeType idim=0; idim < Dim; ++idim) {
            rRightHandSideVector[Index+idim] +=  rVariables.Nu[i] 
                                                * rVariables.Density
                                                * BodyAcceleration[idim]
                                                * rVariables.IntegrationCoefficientInitialConfiguration;
        }
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddMixBodyForce") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateSoilDensity(ElementVariables &rVariables)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainElement::CalculateSoilDensity()") << std::endl;
    const PropertiesType& rProp = this->GetProperties();

    rVariables.Density = (  rVariables.DegreeOfSaturation
                          * rProp[POROSITY]
                          * rProp[DENSITY_WATER] )
                         + (1.0 - rProp[POROSITY] )*rProp[DENSITY_SOLID];

    // KRATOS_INFO("1-UPwSmallStrainElement::CalculateSoilDensity()") << std::endl;
    KRATOS_CATCH("");

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddCouplingTerms( VectorType& rRightHandSideVector,
                                  ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingTerms") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType VoigtSize  = ( Dim == 3 ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN);

    Vector VoigtVector = ZeroVector(VoigtSize);
    for (SizeType idim=0; idim < Dim; ++idim)  VoigtVector[idim] = 1.0;

    Matrix CouplingMatrix = - PORE_PRESSURE_SIGN_FACTOR
                            * rVariables.BiotCoefficient
                            * rVariables.BishopCoefficient
                            * prod( trans(rVariables.B), Matrix( outer_prod(VoigtVector, rVariables.Np) ) )
                            * rVariables.IntegrationCoefficient;

    Vector CouplingForce = prod(CouplingMatrix, rVariables.PressureVector);

    //Distribute coupling block vector 1 into the elemental vector
    const SizeType NumUNodes = rGeom.PointsNumber();

    for (SizeType i = 0; i<NumUNodes; ++i) {
        SizeType Index = i * Dim;
        for (SizeType idim=0; idim < Dim; ++idim) {
            rRightHandSideVector[Index + idim] += CouplingForce[Index + idim];
        }
    }

    if (!rVariables.IgnoreUndrained) {
        const double SaturationCoefficient = rVariables.DegreeOfSaturation / rVariables.BishopCoefficient;
        Vector CouplingFlow =   PORE_PRESSURE_SIGN_FACTOR
                              * SaturationCoefficient
                              * prod(trans(CouplingMatrix),rVariables.VelocityVector);

        //Distribute coupling block vector 2 into the elemental vector
        const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
        for (SizeType i = 0; i<NumPNodes; ++i) {
            rRightHandSideVector[NumUNodes*Dim+i] += CouplingFlow[i];
        }
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddCouplingTerms") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddCompressibilityFlow( VectorType& rRightHandSideVector,
                                        ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityFlow") << std::endl;

    Matrix CompressibilityMatrix = - PORE_PRESSURE_SIGN_FACTOR
                                   * rVariables.BiotModulusInverse
                                   * outer_prod(rVariables.Np, rVariables.Np)
                                   * rVariables.IntegrationCoefficient;

    Vector CompressibilityFlow = - prod(CompressibilityMatrix, rVariables.PressureDtVector);

    //Distribute compressibility block vector into the elemental vector
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for (SizeType i = 0; i < NumPNodes; ++i) {
        rRightHandSideVector[NumUNodes*Dim+i] += CompressibilityFlow[i];
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddCompressibilityFlow") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddPermeabilityFlow( VectorType& rRightHandSideVector,
                                     ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddPermeabilityFlow") << std::endl;

    Matrix PermeabilityMatrix = - PORE_PRESSURE_SIGN_FACTOR
                                * rVariables.DynamicViscosityInverse
                                * rVariables.RelativePermeability
                                * prod( rVariables.DNp_DX, Matrix( prod(rVariables.IntrinsicPermeability, trans(rVariables.DNp_DX)) ) )
                                * rVariables.IntegrationCoefficient;

    Vector PermeabilityFlow = - prod(PermeabilityMatrix, rVariables.PressureVector);

    //Distribute permeability block vector into the elemental vector
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    for (SizeType i = 0; i < NumPNodes; ++i) {
        rRightHandSideVector[NumUNodes*Dim+i] += PermeabilityFlow[i];
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddPermeabilityFlow") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateAndAddFluidBodyFlow( VectorType& rRightHandSideVector,
                                  ElementVariables& rVariables )
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateAndAddFluidBodyFlow") << std::endl;

    Matrix GradNpTPerm =  rVariables.DynamicViscosityInverse
                        * GetProperties()[DENSITY_WATER]
                        * rVariables.RelativePermeability
                        * prod( rVariables.DNp_DX, rVariables.IntrinsicPermeability )
                        * rVariables.IntegrationCoefficient;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    Vector BodyAcceleration = ZeroVector(Dim);

    SizeType Index = 0;
    for (SizeType i = 0; i < NumUNodes; ++i) {
        for (SizeType idim=0; idim < Dim; ++idim) {
            BodyAcceleration[idim] += rVariables.Nu[i]*rVariables.BodyAcceleration[Index++];
        }
    }

    for (SizeType i = 0; i < NumPNodes; ++i) {
        rRightHandSideVector[NumUNodes*Dim+i] += inner_prod(row(GradNpTPerm,i),BodyAcceleration);
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateAndAddFluidBodyFlow") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
GeometryData::IntegrationMethod
    SmallStrainUPwDiffOrderElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
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
void SmallStrainUPwDiffOrderElement::
    CalculateCauchyStrain( ElementVariables& rVariables )
{
    noalias(rVariables.StrainVector) = prod(rVariables.B, rVariables.DisplacementVector);
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateCauchyGreenStrain( ElementVariables& rVariables )
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateCauchyGreenStrain()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();

    //-Compute total deformation gradient
    const Matrix& F = rVariables.F;

    Matrix ETensor;
    ETensor = prod(trans(F), F);

    for (unsigned int i=0; i<Dim; ++i)
        ETensor(i,i) -= 1.0;
    ETensor *= 0.5;

    if (Dim==2) {
        Vector StrainVector;
        StrainVector = MathUtils<double>::StrainTensorToVector(ETensor);
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_XX] = StrainVector[0];
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_YY] = StrainVector[1];
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_ZZ] = 0.0;
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_XY] = StrainVector[2];
    } else {
        noalias(rVariables.StrainVector) = MathUtils<double>::StrainTensorToVector(ETensor);
    }

    // KRATOS_INFO("1-UpdatedLagrangianUPwDiffOrderElement::CalculateCauchyGreenStrain()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateDeformationGradient( ElementVariables& rVariables,
                                  const IndexType& GPoint)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateDeformationGradient()") << std::endl;

    // calculation of derivative of shape function with respect to reference configuration
    // derivative of shape function (displacement)
    Matrix J0, InvJ0, DNu_DX0;
    double detJ0;
    this-> CalculateDerivativesOnInitialConfiguration(detJ0,
                                                      J0,
                                                      InvJ0,
                                                      DNu_DX0,
                                                      GPoint);

    //Calculating current jacobian in order to find deformation gradient
    Matrix J, InvJ;
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

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateCauchyAlmansiStrain(ElementVariables& rVariables )
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateCauchyAlmansiStrain()") << std::endl;

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();

    //-Compute total deformation gradient
    const Matrix& F = rVariables.F;

    Matrix LeftCauchyGreen;
    LeftCauchyGreen = prod(F, trans(F));

    Matrix ETensor;
    double det;
    MathUtils<double>::InvertMatrix(LeftCauchyGreen, ETensor, det );

    for (unsigned int i=0; i<Dim; ++i)
        ETensor(i,i) = 1.0 - ETensor(i,i);

    ETensor *= 0.5;

    if (Dim==2) {
        Vector StrainVector;
        StrainVector = MathUtils<double>::StrainTensorToVector(ETensor);
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_XX] = StrainVector[0];
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_YY] = StrainVector[1];
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_ZZ] = 0.0;
        rVariables.StrainVector[INDEX_2D_PLANE_STRAIN_XY] = StrainVector[2];
    } else {
        noalias(rVariables.StrainVector) = MathUtils<double>::StrainTensorToVector(ETensor);
    }

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateCauchyAlmansiStrain()") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    CalculateHenckyStrain( ElementVariables& rVariables )
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();

    //-Compute total deformation gradient
    const Matrix& F = rVariables.F;

    Matrix CMatrix;
    CMatrix = prod(trans(F), F);

    // Declare the different matrix
    Matrix EigenValuesMatrix = ZeroMatrix(Dim, Dim);
    Matrix EigenVectorsMatrix = ZeroMatrix(Dim, Dim);

    // Decompose matrix
    MathUtils<double>::GaussSeidelEigenSystem(CMatrix, EigenVectorsMatrix, EigenValuesMatrix, 1.0e-16, 20);

    // Calculate the eigenvalues of the E matrix
    for (IndexType i = 0; i < Dim; ++i) {
        EigenValuesMatrix(i, i) = 0.5 * std::log(EigenValuesMatrix(i, i));
    }

    // Calculate E matrix
    Matrix ETensor = ZeroMatrix(Dim, Dim);
    MathUtils<double>::BDBtProductOperation(ETensor, EigenValuesMatrix, EigenVectorsMatrix);

    // Hencky Strain Calculation
    if (Dim==2) {
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
void SmallStrainUPwDiffOrderElement::
    CalculateJacobianOnCurrentConfiguration(double& detJ,
                                            Matrix& rJ,
                                            Matrix& rInvJ,
                                            const IndexType& GPoint) const
{
    KRATOS_TRY

    // KRATOS_INFO("0-SmallStrainUPwDiffOrderElement::CalculateJacobianOnCurrentConfiguration()") << std::endl;
    const GeometryType& rGeom = this->GetGeometry();

    rJ = rGeom.Jacobian( rJ, GPoint, this->GetIntegrationMethod() );
    MathUtils<double>::InvertMatrix( rJ, rInvJ, detJ );

    // KRATOS_INFO("1-SmallStrainUPwDiffOrderElement::CalculateJacobianOnCurrentConfiguration()") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double SmallStrainUPwDiffOrderElement::
    CalculateFluidPressure(const ElementVariables &rVariables, const unsigned int &GPoint)
{
    KRATOS_TRY

    double FluidPressure = inner_prod(rVariables.Np, rVariables.PressureVector);

    return FluidPressure;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
    SetRetentionParameters(const ElementVariables& rVariables,
                           RetentionLaw::Parameters& rRetentionParameters)
{
    KRATOS_TRY

    rRetentionParameters.SetFluidPressure(rVariables.FluidPressure);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SmallStrainUPwDiffOrderElement::
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

} // Namespace Kratos
