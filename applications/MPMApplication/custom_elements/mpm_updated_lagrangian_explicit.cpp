//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "utilities/atomic_utilities.h"
#include "includes/constitutive_law.h"
#include "includes/checks.h"

// Application includes
#include "mpm_application_variables.h"
#include "custom_utilities/mpm_energy_calculation_utility.h"
#include "custom_utilities/mpm_explicit_utilities.h"
#include "custom_elements/mpm_updated_lagrangian_explicit.hpp"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianExplicit::MPMUpdatedLagrangianExplicit( )
    : Element( )
    , mMP()
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianExplicit::MPMUpdatedLagrangianExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
    : Element( NewId, pGeometry )
    , mMP()
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianExplicit::MPMUpdatedLagrangianExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
    : Element( NewId, pGeometry, pProperties )
    , mMP()
{
}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

MPMUpdatedLagrangianExplicit::MPMUpdatedLagrangianExplicit(MPMUpdatedLagrangianExplicit const& rOther)
    :Element(rOther)
    ,mMP(rOther.mMP)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
{
}

//******************************ASSIGNMENT OPERATOR***********************************
//************************************************************************************
MPMUpdatedLagrangianExplicit&  MPMUpdatedLagrangianExplicit::operator=(MPMUpdatedLagrangianExplicit const& rOther)
{
    Element::operator=(rOther);

    mMP = rOther.mMP;

    mDeformationGradientF0.clear();
    mDeformationGradientF0 = rOther.mDeformationGradientF0;

    mDeterminantF0 = rOther.mDeterminantF0;
    mConstitutiveLawVector = rOther.mConstitutiveLawVector;

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************
Element::Pointer MPMUpdatedLagrangianExplicit::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new MPMUpdatedLagrangianExplicit( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer MPMUpdatedLagrangianExplicit::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< MPMUpdatedLagrangianExplicit >(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************
Element::Pointer MPMUpdatedLagrangianExplicit::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    MPMUpdatedLagrangianExplicit NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.mMP = mMP;

    NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();

    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mDeterminantF0 = mDeterminantF0;

    return Element::Pointer( new MPMUpdatedLagrangianExplicit(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianExplicit::~MPMUpdatedLagrangianExplicit()
{
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        // Initialize parameters
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        mDeterminantF0 = 1;
        mDeformationGradientF0 = IdentityMatrix(dimension);

        // Initialize constitutive law and materials
        InitializeMaterial(rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::InitializeGeneralVariables(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const bool is_axisymmetric = (rCurrentProcessInfo.Has(IS_AXISYMMETRIC))
        ? rCurrentProcessInfo.GetValue(IS_AXISYMMETRIC)
        : false;
    const SizeType def_grad_dim = (is_axisymmetric)
        ? 3
        : dimension;

    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.detFT = 1;

    rVariables.F.resize(def_grad_dim, def_grad_dim, false );

    rVariables.F0.resize(def_grad_dim, def_grad_dim, false );

    rVariables.FT.resize(def_grad_dim, def_grad_dim, false );

    rVariables.ConstitutiveMatrix.resize(strain_size, strain_size, false );

    rVariables.StrainVector.resize(strain_size, false );

    rVariables.StressVector.resize(strain_size, false );

    rVariables.DN_DX.resize( number_of_nodes, dimension, false );

    // CurrentDisp is the unknown variable. It represents the nodal delta displacement. When it is predicted is equal to zero.
    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::SetGeneralVariables(
    GeneralVariables& rVariables,
    ConstitutiveLaw::Parameters& rValues,
    const Vector& rN)
{
    GeometryType& r_geometry = GetGeometry();

    // Variables.detF is the determinant of the incremental total deformation gradient
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    // Check if detF is negative (element is inverted)
    if(rVariables.detF<0)
    {
        KRATOS_INFO("MPMUpdatedLagrangianExplicit")<<" Element: "<<this->Id()<<std::endl;
        KRATOS_INFO("MPMUpdatedLagrangianExplicit")<<" Element position: "<< mMP.xg <<std::endl;
        KRATOS_INFO("MPMUpdatedLagrangianExplicit")<<" Element velocity: "<< mMP.velocity <<std::endl;
        const unsigned int number_of_nodes = r_geometry.PointsNumber();
        KRATOS_INFO("MPMUpdatedLagrangianExplicit") << " Shape functions: " << r_geometry.ShapeFunctionsValues() << std::endl;
        KRATOS_INFO("MPMUpdatedLagrangianExplicit") << " Quadrature points: " << r_geometry.IntegrationPointsNumber() << std::endl;
        KRATOS_INFO("MPMUpdatedLagrangianExplicit") << " Parent geometry ID: " << r_geometry.GetGeometryParent(0).Id() << std::endl;
        KRATOS_INFO("MPMUpdatedLagrangianExplicit") << " Parent geometry number of points: " << r_geometry.GetGeometryParent(0).PointsNumber() << std::endl;

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            const array_1d<double, 3> & current_position      = r_geometry[i].Coordinates();
            const array_1d<double, 3> & current_displacement  = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3> & previous_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT,1);

            KRATOS_INFO("MPMUpdatedLagrangianExplicit")<<" NODE ["<<r_geometry[i].Id()<<"]: (Current position: "<<current_position<<") "<<std::endl;
            KRATOS_INFO("MPMUpdatedLagrangianExplicit")<<" ---Current Disp: "<<current_displacement<<" (Previour Disp: "<<previous_displacement<<")"<<std::endl;
        }

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if( r_geometry[i].SolutionStepsDataHas(CONTACT_FORCE) )
            {
                const array_1d<double, 3 > & PreContactForce = r_geometry[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
                const array_1d<double, 3 > & ContactForce = r_geometry[i].FastGetSolutionStepValue(CONTACT_FORCE);
                KRATOS_INFO("MPMUpdatedLagrangianExplicit")<<" ---Contact_Force: (Pre:"<<PreContactForce<<", Current:"<<ContactForce<<") "<<std::endl;
            }
            else
            {
                KRATOS_INFO("MPMUpdatedLagrangianExplicit")<<" ---Contact_Force: NULL "<<std::endl;
            }
        }

        KRATOS_ERROR << "MPM UPDATED LAGRANGIAN DISPLACEMENT ELEMENT INVERTED: |F|<0  detF = " << rVariables.detF << std::endl;
    }

    rVariables.detFT = rVariables.detF * rVariables.detF0;
    rVariables.FT    = prod( rVariables.F, rVariables.F0 );

    rValues.SetDeterminantF(rVariables.detFT);
    rValues.SetDeformationGradientF(rVariables.FT);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rN);
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::CalculateElementalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);
    const Vector& r_N = row(GetGeometry().ShapeFunctionsValues(), 0);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // The MP_Volume (integration weight) is evaluated
    mMP.volume = mMP.mass / mMP.density;

    if (CalculateResidualVectorFlag) // if calculation of the vector is required
    {
        // Contribution to forces (in residual term) are calculated
        Vector volume_force = mMP.volume_acceleration * mMP.mass;
        this->CalculateAndAddRHS(
            rRightHandSideVector,
            Variables,
            volume_force,
            mMP.volume,
            rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::CalculateAndAddRHS(
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    Vector& rVolumeForce,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
    this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    MPMExplicitUtilities::CalculateAndAddExplicitInternalForce(rCurrentProcessInfo ,
        *this, mMP.cauchy_stress_vector, mMP.volume,
        mConstitutiveLawVector->GetStrainSize(), rRightHandSideVector);
}

//************************************************************************************
//*********************Calculate the contribution of external force*******************
void MPMUpdatedLagrangianExplicit::CalculateAndAddExternalForces(
    VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = dimension * i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[index + j] += r_N(0, i) * rVolumeForce[j];
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::CalculateExplicitStresses(
    const ProcessInfo& rCurrentProcessInfo,
    GeneralVariables& rVariables)
{
    KRATOS_TRY

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // Define the stress measure
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions = Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    // use element provided strain incremented from velocity gradient
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);

    // Compute explicit element kinematics, strain is incremented here.
    Matrix Jacobian;
    GetGeometry().Jacobian(Jacobian, 0);
    Matrix InvJ;
    double detJ;
    MathUtils<double>::InvertMatrix(Jacobian, InvJ, detJ);
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    Matrix r_DN_De = GetGeometry().ShapeFunctionLocalGradient(0);
    rVariables.DN_DX = prod(r_DN_De, InvJ); // cartesian gradients

    MPMExplicitUtilities::CalculateExplicitKinematics(rCurrentProcessInfo, *this,
        mMP.almansi_strain_vector, rVariables.F, mConstitutiveLawVector->GetStrainSize());

    rVariables.StressVector = mMP.cauchy_stress_vector;
    rVariables.StrainVector = mMP.almansi_strain_vector;

    // Update gradient deformation
    rVariables.F0 = mDeformationGradientF0; // total member def grad NOT including this increment
    rVariables.FT = prod(rVariables.F, rVariables.F0); // total def grad including this increment
    rVariables.detF = MathUtils<double>::Det(rVariables.F); // det of current increment
    rVariables.detF0 = MathUtils<double>::Det(rVariables.F0); // det of def grad NOT including this increment
    rVariables.detFT = MathUtils<double>::Det(rVariables.FT); // det of total def grad including this increment
    mDeformationGradientF0 = rVariables.FT; // update member internal total grad def
    mDeterminantF0 = rVariables.detFT; // update member internal total grad def det

    // Update MP volume
    if (rCurrentProcessInfo.GetValue(IS_COMPRESSIBLE))
    {
        mMP.density = (GetProperties()[DENSITY]) / rVariables.detFT;
        mMP.volume = mMP.mass / mMP.density;
    }

    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);

    // Set general variables to constitutivelaw parameters
    const Vector& r_N_vec = row(r_N, 0);
    this->SetGeneralVariables(rVariables, Values, r_N_vec);

    // Calculate Material Response
    /* NOTE:
    The function below will call CalculateMaterialResponseCauchy() by default and then (may)
    call CalculateMaterialResponseKirchhoff() in the constitutive_law.*/
    mConstitutiveLawVector->CalculateMaterialResponse(Values, rVariables.StressMeasure);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix = Matrix(0, 0);

    const SizeType mat_size = GetNumberOfDofs() * GetGeometry().size();
    if (rRightHandSideVector.size() != mat_size) {
        rRightHandSideVector.resize(mat_size, false);
    }
    rRightHandSideVector = ZeroVector(mat_size);

    CalculateElementalSystem(left_hand_side_matrix, rRightHandSideVector,
        rCurrentProcessInfo, false, true);
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType right_hand_side_vector = Vector(0);

    const SizeType mat_size = GetNumberOfDofs() * GetGeometry().size();
    if (rLeftHandSideMatrix.size1() != mat_size && rLeftHandSideMatrix.size2() != mat_size) {
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    CalculateElementalSystem(
        rLeftHandSideMatrix, right_hand_side_vector,
        rCurrentProcessInfo, true, false);
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType mat_size = GetNumberOfDofs() * GetGeometry().size();
    if (rLeftHandSideMatrix.size1() != mat_size && rLeftHandSideMatrix.size2() != mat_size) {
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    if (rRightHandSideVector.size() != mat_size) {
        rRightHandSideVector.resize(mat_size, false);
    }
    rRightHandSideVector = ZeroVector(mat_size);

    CalculateElementalSystem(
        rLeftHandSideMatrix, rRightHandSideVector,
        rCurrentProcessInfo, true, true);
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::AddExplicitContribution(
    const ProcessInfo& rCurrentProcessInfo)
{
    /* NOTE:
    This is moved from initially InitializeSolutionStep due to #13432 (as per @RiccardoRossi suggestion).
    This is temporary, and will be moved to an utility in the future, after restructuring of
    MPM's internal variables data structure. This function is called during predict in the schemes.

    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.*/

    GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    const bool is_explicit_central_difference = (rCurrentProcessInfo.Has(IS_EXPLICIT_CENTRAL_DIFFERENCE))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT_CENTRAL_DIFFERENCE)
        : false;

    // Calculating shape functions
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    array_1d<double,3> nodal_momentum = ZeroVector(3);
    array_1d<double,3> nodal_inertia  = ZeroVector(3);

    // Here MP contribution in terms of momentum, inertia and mass are added
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for (unsigned int j = 0; j < dimension; j++)
        {
            nodal_momentum[j] = r_N(0, i) * mMP.velocity[j] * mMP.mass;
            nodal_inertia[j] = r_N(0, i) * mMP.acceleration[j] * mMP.mass;
        }

        // Add in the predictor velocity increment for central difference explicit
        // This is the 'previous grid acceleration', which is actually
        // be the initial material point acceleration mapped to the grid.
        if (is_explicit_central_difference) {
            const double& delta_time = rCurrentProcessInfo[DELTA_TIME];
            for (unsigned int j = 0; j < dimension; j++) {
                nodal_momentum[j] += 0.5 * delta_time * (r_N(0, i) * mMP.acceleration[j]) * mMP.mass;
            }
        }

        AtomicAdd(r_geometry[i].FastGetSolutionStepValue(NODAL_MOMENTUM, 0), nodal_momentum);
        AtomicAdd(r_geometry[i].FastGetSolutionStepValue(NODAL_INERTIA, 0), nodal_inertia);
        AtomicAdd(r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0), r_N(0, i)*mMP.mass);
    }
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::FinalizeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR
    << "FinalizeSolutionStep for explicit time integration is done in the scheme";

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::FinalizeStepVariables(
    GeneralVariables & rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Update internal (historical) variables
    mDeterminantF0         = rVariables.detF* rVariables.detF0;
    mDeformationGradientF0 = prod(rVariables.F, rVariables.F0);

    mMP.cauchy_stress_vector = rVariables.StressVector;
    mMP.almansi_strain_vector = rVariables.StrainVector;

    // Delta Plastic Strains
    if (mConstitutiveLawVector->Has(MP_DELTA_PLASTIC_STRAIN))
        mConstitutiveLawVector->GetValue(MP_DELTA_PLASTIC_STRAIN, mMP.delta_plastic_strain );
    if (mConstitutiveLawVector->Has(MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN))
        mConstitutiveLawVector->GetValue(MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN, mMP.delta_plastic_volumetric_strain);
    if (mConstitutiveLawVector->Has(MP_DELTA_PLASTIC_DEVIATORIC_STRAIN))
        mConstitutiveLawVector->GetValue(MP_DELTA_PLASTIC_DEVIATORIC_STRAIN, mMP.delta_plastic_deviatoric_strain);

    // Total Plastic Strain
    if (mConstitutiveLawVector->Has(MP_EQUIVALENT_PLASTIC_STRAIN))
        mConstitutiveLawVector->GetValue(MP_EQUIVALENT_PLASTIC_STRAIN, mMP.equivalent_plastic_strain );
    if (mConstitutiveLawVector->Has(MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN))
        mConstitutiveLawVector->GetValue(MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN, mMP.accumulated_plastic_volumetric_strain);
    if (mConstitutiveLawVector->Has(MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN))
        mConstitutiveLawVector->GetValue(MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN, mMP.accumulated_plastic_deviatoric_strain);

}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::InitializeMaterial(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    GeneralVariables Variables;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        Vector N = row(GetGeometry().ShapeFunctionsValues(), 0);
        mConstitutiveLawVector->InitializeMaterial(
            GetProperties(), GetGeometry(), N);

        mMP.almansi_strain_vector = ZeroVector(mConstitutiveLawVector->GetStrainSize());
        mMP.cauchy_stress_vector = ZeroVector(mConstitutiveLawVector->GetStrainSize());

        // Resize the deformation gradient if we are axisymmetric
        if (mConstitutiveLawVector->GetStrainSize() == 4) mDeformationGradientF0 = IdentityMatrix(3);
    }
    else
        KRATOS_ERROR <<  "A constitutive law needs to be specified for the element with ID: " << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::ResetConstitutiveLaw()
{
    KRATOS_TRY
    GeneralVariables Variables;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        mConstitutiveLawVector->ResetMaterial(
            GetProperties(),
            GetGeometry(),
            row(GetGeometry().ShapeFunctionsValues(), 0));
    }

    KRATOS_CATCH( "" )
}

//*************************COMPUTE CURRENT DISPLACEMENT*******************************
//************************************************************************************
/*
This function convert the computed nodal displacement into matrix of (number_of_nodes, dimension)
*/
Matrix& MPMUpdatedLagrangianExplicit::CalculateCurrentDisp(
    Matrix & rCurrentDisp,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    rCurrentDisp = ZeroMatrix(number_of_nodes, dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3 > & current_displacement  = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rCurrentDisp(i,j) = current_displacement[j];
        }
    }

    return rCurrentDisp;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& CurrentProcessInfo) const
{
    const GeometryType& r_geometry = GetGeometry();
    int number_of_nodes = r_geometry.size();
    int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( rResult.size() != matrix_size )
        rResult.resize( matrix_size, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index] = r_geometry[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if ( dimension == 3 )
            rResult[index + 2] = r_geometry[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& CurrentProcessInfo) const
{
    const GeometryType& r_geometry = GetGeometry();
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < r_geometry.size(); i++ )
    {
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Y ) );

        if ( r_geometry.WorkingSpaceDimension() == 3 )
        {
            rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rRHSVariable == RESIDUAL_VECTOR &&
        rDestinationVariable == FORCE_RESIDUAL) {
        GeometryType& r_geometry = GetGeometry();
        const unsigned int dimension = r_geometry.WorkingSpaceDimension();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();

        for (size_t i = 0; i < number_of_nodes; ++i) {
            size_t index = dimension * i;
            array_1d<double, 3>& r_force_residual = r_geometry[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (size_t j = 0; j < dimension; ++j) {
                r_force_residual[j] += rRHSVector[index + j];
            }
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//****************MASS MATRIX*********************************************************
void MPMUpdatedLagrangianExplicit::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Call the values of the shape function for the single element
    Vector N = row(GetGeometry().ShapeFunctionsValues(), 0);

    const bool is_lumped_mass_matrix = (rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX))
        ? rCurrentProcessInfo.GetValue(COMPUTE_LUMPED_MASS_MATRIX)
        : true;

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType matrix_size = dimension * number_of_nodes;

    if ( rMassMatrix.size1() != matrix_size || rMassMatrix.size2() != matrix_size)
        rMassMatrix.resize( matrix_size, matrix_size, false );
    rMassMatrix = ZeroMatrix(matrix_size, matrix_size);

    if (!is_lumped_mass_matrix) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType j = 0; j < number_of_nodes; ++j) {
                for (IndexType k = 0; k < dimension; ++k)
                {
                    const IndexType index_i = i * dimension + k;
                    const IndexType index_j = j * dimension + k;
                    rMassMatrix(index_i, index_j) = N[i] * N[j] * mMP.mass;
                }
            }
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType k = 0; k < dimension; ++k)
            {
                const IndexType index = i * dimension + k;
                rMassMatrix(index, index) = N[i] * mMP.mass;
            }
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::GetValuesVector(
    Vector& values,
    int Step) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::GetFirstDerivativesVector(
    Vector& values,
    int Step) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::GetSecondDerivativesVector(
    Vector& values,
    int Step) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianExplicit::GetHistoricalVariables(GeneralVariables& rVariables)
{
    //Deformation Gradient F ( set to identity )
    unsigned int size =  rVariables.F.size1();
    rVariables.detF  = 1;
    rVariables.F     = IdentityMatrix(size);

    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;
}

///@}
///@name Access Get Values
///@{

void MPMUpdatedLagrangianExplicit::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == CALCULATE_EXPLICIT_MP_STRESS)
    {
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);
        this->CalculateExplicitStresses(rCurrentProcessInfo, Variables);
        this->FinalizeStepVariables(Variables, rCurrentProcessInfo);
        rValues[0] = true;
    }
    else if (rVariable == EXPLICIT_MAP_GRID_TO_MP)
    {
        MPMExplicitUtilities::UpdateGaussPointExplicit(rCurrentProcessInfo, *this);
        rValues[0] = true;
    }
    else if (rVariable == CALCULATE_MUSL_VELOCITY_FIELD)
    {
        MPMExplicitUtilities::CalculateMUSLGridVelocity(rCurrentProcessInfo, *this);
        rValues[0] = true;
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMUpdatedLagrangianExplicit::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_MATERIAL_ID) {
        rValues[0] = GetProperties().Id();
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMUpdatedLagrangianExplicit::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_DENSITY) {
        rValues[0] = mMP.density;
    }
    else if (rVariable == MP_MASS) {
        rValues[0] = mMP.mass;
    }
    else if (rVariable == MP_VOLUME) {
        rValues[0] = mMP.volume;
    }
    else if (rVariable == MP_POTENTIAL_ENERGY) {
        rValues[0] = MPMEnergyCalculationUtility::CalculatePotentialEnergy(*this);
    }
    else if (rVariable == MP_KINETIC_ENERGY) {
        rValues[0] = MPMEnergyCalculationUtility::CalculateKineticEnergy(*this);
    }
    else if (rVariable == MP_STRAIN_ENERGY) {
        rValues[0] = MPMEnergyCalculationUtility::CalculateStrainEnergy(*this);
    }
    else if (rVariable == MP_TOTAL_ENERGY) {
        rValues[0] = MPMEnergyCalculationUtility::CalculateTotalEnergy(*this);
    }
    else if (rVariable == MP_HARDENING_RATIO || rVariable == MP_EQUIVALENT_STRESS ||
        rVariable == MP_EQUIVALENT_PLASTIC_STRAIN || rVariable == MP_EQUIVALENT_PLASTIC_STRAIN_RATE ||
        rVariable == MP_TEMPERATURE) {
        rValues[0] = mConstitutiveLawVector->GetValue(rVariable, rValues[0]);
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMUpdatedLagrangianExplicit::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_COORD || rVariable == MPC_COORD) {
        rValues[0] = mMP.xg;
    }
    else if (rVariable == MP_DISPLACEMENT) {
        rValues[0] = mMP.displacement;
    }
    else if (rVariable == MP_VELOCITY) {
        rValues[0] = mMP.velocity;
    }
    else if (rVariable == MP_ACCELERATION) {
        rValues[0] = mMP.acceleration;
    }
    else if (rVariable == MP_VOLUME_ACCELERATION) {
        rValues[0] = mMP.volume_acceleration;
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMUpdatedLagrangianExplicit::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_CAUCHY_STRESS_VECTOR) {
        rValues[0] = mMP.cauchy_stress_vector;
    }
    else if (rVariable == MP_ALMANSI_STRAIN_VECTOR) {
        rValues[0] = mMP.almansi_strain_vector;
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

///@}
///@name Access Set Values
///@{

void MPMUpdatedLagrangianExplicit::SetValuesOnIntegrationPoints(
    const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_MASS) {
        mMP.mass = rValues[0];
    }
    else if (rVariable == MP_DENSITY) {
        mMP.density = rValues[0];
    }
    else if (rVariable == MP_VOLUME) {
        mMP.volume = rValues[0];
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMUpdatedLagrangianExplicit::SetValuesOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    const std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_COORD || rVariable == MPC_COORD) {
        mMP.xg = rValues[0];
    }
    else if (rVariable == MP_DISPLACEMENT) {
        mMP.displacement = rValues[0];
    }
    else if (rVariable == MP_VELOCITY) {
        mMP.velocity = rValues[0];
    }
    else if (rVariable == MP_ACCELERATION) {
        mMP.acceleration = rValues[0];
    }
    else if (rVariable == MP_VOLUME_ACCELERATION) {
        mMP.volume_acceleration = rValues[0];
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMUpdatedLagrangianExplicit::SetValuesOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    const std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_CAUCHY_STRESS_VECTOR) {
        mMP.cauchy_stress_vector = rValues[0];
    }
    else if (rVariable == MP_ALMANSI_STRAIN_VECTOR) {
        mMP.almansi_strain_vector = rValues[0];
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

///@}

/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  MPMUpdatedLagrangianExplicit::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    Element::Check(rCurrentProcessInfo);

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    // Verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;

    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
        if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient) correct_strain_measure = true;
        if (LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Velocity_Gradient) correct_strain_measure = true;
    }

    KRATOS_ERROR_IF(correct_strain_measure == false ) << "Constitutive law is not compatible with the element type: Large Displacements " << std::endl;

    // Verify that the dofs exist
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    if( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false)
    {
        KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    }
    else
    {
        // Verify that the constitutive law has the correct dimension
        if ( dimension == 2 )
        {
            KRATOS_ERROR_IF_NOT(this->GetProperties().Has( THICKNESS )) << "THICKNESS not provided for element " << this->Id() << std::endl;
        }
        else
        {
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) " << this->Id() << std::endl;
        }

        // Check constitutive law
        this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Check( this->GetProperties(), r_geometry, rCurrentProcessInfo );
    }

    return 0;

    KRATOS_CATCH( "" );
}

void MPMUpdatedLagrangianExplicit::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )

    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
    rSerializer.save("MP",mMP);
}

void MPMUpdatedLagrangianExplicit::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
    rSerializer.load("MP",mMP);
}

} // Namespace Kratos
