//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Laura Moreno Martínez
//


// System includes
#include <boost/token_functions.hpp>
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/mpm_updated_lagrangian_VP_VMS.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "mpm_application_variables.h"
#include "includes/checks.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/material_point_generator_utility.h"
#include "includes/variables.h"
#include "includes/mat_variables.h"
#include "custom_utilities/mpm_energy_calculation_utility.h"
#include "custom_utilities/mpm_explicit_utilities.h"
#include "custom_utilities/mpm_math_utilities.h"




namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianVPVMS::MPMUpdatedLagrangianVPVMS()
    : Element( )
    , mMP()
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianVPVMS::MPMUpdatedLagrangianVPVMS( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
    , mMP()
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianVPVMS::MPMUpdatedLagrangianVPVMS( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
    , mMP()
{
    mFinalizedStep = true;


}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

MPMUpdatedLagrangianVPVMS::MPMUpdatedLagrangianVPVMS( MPMUpdatedLagrangianVPVMS const& rOther)
    :Element(rOther)
    ,mMP(rOther.mMP)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    ,mFinalizedStep(rOther.mFinalizedStep)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

MPMUpdatedLagrangianVPVMS&  MPMUpdatedLagrangianVPVMS::operator=(MPMUpdatedLagrangianVPVMS const& rOther)
{
    Element::operator=(rOther);

    mMP = rOther.mMP;
    m_mp_pressure = rOther.m_mp_pressure;

    mDeformationGradientF0.clear();
    mDeformationGradientF0 = rOther.mDeformationGradientF0;

    mDeterminantF0 = rOther.mDeterminantF0;
    mConstitutiveLawVector = rOther.mConstitutiveLawVector;

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer MPMUpdatedLagrangianVPVMS::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new MPMUpdatedLagrangianVPVMS( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer MPMUpdatedLagrangianVPVMS::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< MPMUpdatedLagrangianVPVMS >(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer MPMUpdatedLagrangianVPVMS::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    MPMUpdatedLagrangianVPVMS NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.m_mp_pressure = m_mp_pressure;

    //NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();

    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mDeterminantF0 = mDeterminantF0;

    return Element::Pointer( new MPMUpdatedLagrangianVPVMS(NewElement) );
}
//*******************************DESTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianVPVMS::~MPMUpdatedLagrangianVPVMS()
{
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    //printf("    > Initialize\n");
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

void MPMUpdatedLagrangianVPVMS::GetDofList( DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo ) const
{
    //printf("    > GetDofList\n");
    rElementalDofList.resize( 0 );

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rElementalDofList.push_back( r_geometry[i].pGetDof( VELOCITY_X ) );
        rElementalDofList.push_back( r_geometry[i].pGetDof( VELOCITY_Y ) );

        if ( dimension == 3 )
        {
            rElementalDofList.push_back( r_geometry[i].pGetDof( VELOCITY_Z ) );
        }
        rElementalDofList.push_back( r_geometry[i].pGetDof( PRESSURE ));

    }
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo ) const
{
    //printf("    > EquationIdVector\n");
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension       = r_geometry.WorkingSpaceDimension();
    unsigned int element_size          = number_of_nodes * dimension + number_of_nodes;

    if ( rResult.size() != element_size )
        rResult.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension + i;
        rResult[index]     = r_geometry[i].GetDof( VELOCITY_X ).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof( VELOCITY_Y ).EquationId();

        if( dimension == 3)
        {
            rResult[index + 2] = r_geometry[i].GetDof( VELOCITY_Z ).EquationId();
            rResult[index + 3] = r_geometry[i].GetDof( PRESSURE ).EquationId();
        }
        else
        {
            rResult[index + 2] = r_geometry[i].GetDof( PRESSURE ).EquationId();
        }
    }
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::GetValuesVector( Vector& values, int Step ) const
{
    //printf("    > GetValuesVector\n");
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension       = r_geometry.WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = r_geometry[i].FastGetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
        {
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
            values[index + 3] = r_geometry[i].FastGetSolutionStepValue( PRESSURE, Step );
        }
        else
        {
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( PRESSURE, Step );
        }

    }
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::GetFirstDerivativesVector( Vector& values, int Step ) const
{
    //printf("    > GetFirstDerivativesVector\n");
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension       = r_geometry.WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );
        if ( dimension == 3 )
        {
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
            values[index + 3] = 0;
        }
        else
        {
            values[index + 2] = 0;
        }
    }
}

////************************************************************************************
////************************************************************************************

void MPMUpdatedLagrangianVPVMS::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/
    //printf("    > InitializeSolutionStep\n");
    GeometryType& r_geometry = GetGeometry();
    unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    GeneralVariables Variables;

    // Calculating shape function
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    mFinalizedStep = false;

    array_1d<double,3> aux_MP_velocity = ZeroVector(3);
    array_1d<double,3> aux_MP_acceleration = ZeroVector(3);
    array_1d<double,3> nodal_momentum = ZeroVector(3);
    array_1d<double,3> nodal_inertia = ZeroVector(3);
    double aux_MP_pressure = 0.0;

    for (unsigned int j=0; j<number_of_nodes; j++)
    {
        // These are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
        array_1d<double, 3 > nodal_acceleration = ZeroVector(3);
        if (r_geometry[j].SolutionStepsDataHas(ACCELERATION))
            nodal_acceleration = r_geometry[j].FastGetSolutionStepValue(ACCELERATION,1);

        array_1d<double, 3 > nodal_velocity = ZeroVector(3);
        if (r_geometry[j].SolutionStepsDataHas(VELOCITY))
            nodal_velocity = r_geometry[j].FastGetSolutionStepValue(VELOCITY,1);
        //array_1d<double, 3 > nodal_displacement = ZeroVector(3);
        //if (r_geometry[j].SolutionStepsDataHas(DISPLACEMENT))
        //    nodal_displacement = r_geometry[j].FastGetSolutionStepValue(DISPLACEMENT,1);

        // These are the values of nodal pressure evaluated in the initialize solution step
        const double& nodal_pressure = r_geometry[j].FastGetSolutionStepValue(PRESSURE,1);

        aux_MP_pressure += r_N(0, j) * nodal_pressure;

        for (unsigned int k = 0; k < dimension; k++)
        {
            //aux_MP_displacement[k] += r_N(0, j) * nodal_displacement[k];
            aux_MP_velocity[k] += r_N(0, j) * nodal_velocity[k];
            aux_MP_acceleration[k] += r_N(0, j) * nodal_acceleration[k];
        }
    }

    // Here MP contribution in terms of momentum, inertia, mass-pressure and mass are added
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double nodal_mpressure = r_N(0, i) * (m_mp_pressure - aux_MP_pressure) * mMP.mass;

        for (unsigned int j = 0; j < dimension; j++)
        {
            nodal_momentum[j] = r_N(0, i) * (mMP.velocity[j] - aux_MP_velocity[j]) * mMP.mass;
            nodal_inertia[j]  = r_N(0, i) * (mMP.acceleration[j] - aux_MP_acceleration[j]) * mMP.mass;
        }

        r_geometry[i].SetLock();
        r_geometry[i].FastGetSolutionStepValue(NODAL_MOMENTUM, 0)  += nodal_momentum;
        r_geometry[i].FastGetSolutionStepValue(NODAL_INERTIA, 0)   += nodal_inertia;
        r_geometry[i].FastGetSolutionStepValue(NODAL_MPRESSURE, 0) += nodal_mpressure;

        r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0) += r_N(0, i) * mMP.mass;
        r_geometry[i].UnSetLock();
    }
}

////************************************************************************************
////************************************************************************************

void MPMUpdatedLagrangianVPVMS::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    //printf("    > FinalizeSolutionStep\n");
    const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
        : false;

    KRATOS_ERROR_IF(is_explicit)
    << "FinalizeSolutionStep for explicit time integration is done in the scheme";

    // Create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);
    const Vector& r_N = row(GetGeometry().ShapeFunctionsValues(), 0);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    //Flags &ConstitutiveLawOptions=Values.GetOptions();

    //ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

    // Compute element kinematics B, F, DN_DX ...
    this->CalculateKinematics(Variables, rCurrentProcessInfo);

    // Set general variables to constitutivelaw parameters
    this->SetGeneralVariables(Variables,Values, r_N);

    // Call the constitutive law to update material variables
    //mConstitutiveLawVector->FinalizeMaterialResponse(Values, Variables.StressMeasure);

    // Call the element internal variables update
    this->FinalizeStepVariables(Variables, rCurrentProcessInfo);

    mFinalizedStep = true;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > CalculateLocalSystem\n");
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

void MPMUpdatedLagrangianVPVMS::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > CalculateRightHandSide\n");
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


void MPMUpdatedLagrangianVPVMS::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > CalculateLeftHandSide\n");
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
//****************MASS MATRIX*********************************************************

void MPMUpdatedLagrangianVPVMS::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    //printf("    > CalculateMassMatrix\n");
    // Call the values of the shape function for the single element
    Vector N = row(GetGeometry().ShapeFunctionsValues(), 0);

    const bool is_lumped_mass_matrix = (rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX))
        ? rCurrentProcessInfo.GetValue(COMPUTE_LUMPED_MASS_MATRIX)
        : true;

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType matrix_size = (dimension + 1) * number_of_nodes;

    if (rMassMatrix.size1() != matrix_size || rMassMatrix.size2() != matrix_size)
        rMassMatrix.resize(matrix_size, matrix_size, false);
    rMassMatrix = ZeroMatrix(matrix_size, matrix_size);

    if (!is_lumped_mass_matrix) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType j = 0; j < number_of_nodes; ++j) {
                for (IndexType k = 0; k < dimension; ++k)
                {
                    const IndexType index_i = i * (dimension + 1) + k;
                    const IndexType index_j = j * (dimension + 1) + k;
                    rMassMatrix(index_i, index_j) = N[i] * N[j] * mMP.mass;
                }
            }
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            for (IndexType k = 0; k < dimension; ++k)
            {
                const IndexType index = i * (dimension + 1) + k;
                rMassMatrix(index, index) = N[i] * mMP.mass;
            }
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::AddExplicitContribution(const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable, const Variable<array_1d<double, 3>>& rDestinationVariable, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    //printf("    > AddExplicitContribution\n");
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

/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  MPMUpdatedLagrangianVPVMS::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    //printf("    > Check\n");
   const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
    ? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
    : false;
    KRATOS_ERROR_IF(is_explicit)
    << "Explicit time integration not implemented for Updated Lagrangian UP MPM Element";

    Element::Check(rCurrentProcessInfo);

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    
    // Verify that the constitutive law exists
    //if( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false)
    //{
    //    KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    //}
    //else
    //{
    //    // Verify that the constitutive law has the correct dimension
    //    if ( dimension == 2 )
    //    {
    //        KRATOS_ERROR_IF_NOT(this->GetProperties().Has( THICKNESS )) << "THICKNESS not provided for element " << this->Id() << std::endl;
    //    }
    //    else
    //    {
    //        KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) " << this->Id() << std::endl;
    //    }
//
    //    // Check constitutive law
    //    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Check( this->GetProperties(), r_geometry, rCurrentProcessInfo );
    //}

    // Verify compatibility with the constitutive law
    //ConstitutiveLaw::Features LawFeatures;
    //this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    //bool correct_strain_measure = false;
    //for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    //{
    //    if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient) correct_strain_measure = true;
    //    if (is_explicit && LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Velocity_Gradient) correct_strain_measure = true;
//
    //}
    //if (true)
    //{
//
    //}
//
    //KRATOS_ERROR_IF(correct_strain_measure == false ) << "Constitutive law is not compatible with the element type: Large Displacements " << std::endl;
    //KRATOS_ERROR_IF(LawFeatures.mOptions.IsNot(ConstitutiveLaw::U_P_LAW)) << "Constitutive law is not compatible with the U-P element type: Large Displacements U_P" << std::endl;

    // Verify that the dofs exist
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)

        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z, rnode)
    }

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

//void MPMUpdatedLagrangianVPVMS::CalculateExplicitStresses(const ProcessInfo& rCurrentProcessInfo,
//    GeneralVariables& rVariables)
//{
//    KRATOS_TRY
//
//    // Create constitutive law parameters:
//    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
//
//    // Define the stress measure
//    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;
//
//    // Set constitutive law flags:
//    Flags& ConstitutiveLawOptions = Values.GetOptions();
//    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
//
//    // use element provided strain incremented from velocity gradient
//    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
//    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
//
//    // Compute explicit element kinematics, strain is incremented here.
//    Matrix Jacobian;
//    GetGeometry().Jacobian(Jacobian, 0);
//    Matrix InvJ;
//    double detJ;
//    MathUtils<double>::InvertMatrix(Jacobian, InvJ, detJ);
//    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
//    Matrix r_DN_De = GetGeometry().ShapeFunctionLocalGradient(0);
//    rVariables.DN_DX = prod(r_DN_De, InvJ); // cartesian gradients
//
//    MPMExplicitUtilities::CalculateExplicitKinematics(rCurrentProcessInfo, *this,
//        mMP.almansi_strain_vector, rVariables.F, mConstitutiveLawVector->GetStrainSize());
//
//    rVariables.StressVector = mMP.cauchy_stress_vector;
//    rVariables.StrainVector = mMP.almansi_strain_vector;
//
//    // Update gradient deformation
//    rVariables.F0 = mDeformationGradientF0; // total member def grad NOT including this increment
//    rVariables.FT = prod(rVariables.F, rVariables.F0); // total def grad including this increment
//    rVariables.detF = MathUtils<double>::Det(rVariables.F); // det of current increment
//    rVariables.detF0 = MathUtils<double>::Det(rVariables.F0); // det of def grad NOT including this increment
//    rVariables.detFT = MathUtils<double>::Det(rVariables.FT); // det of total def grad including this increment
//    mDeformationGradientF0 = rVariables.FT; // update member internal total grad def
//    mDeterminantF0 = rVariables.detFT; // update member internal total grad def det
//
//    // Update MP volume
//    if (rCurrentProcessInfo.GetValue(IS_COMPRESSIBLE))
//    {
//        mMP.density = (GetProperties()[DENSITY]) / rVariables.detFT;
//        mMP.volume = mMP.mass / mMP.density;
//    }
//
//    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
//
//    // Set general variables to constitutivelaw parameters
//    const Vector& r_N_vec = row(r_N, 0);
//    this->SetGeneralVariables(rVariables, Values, r_N_vec);
//
//    // Calculate Material Response
//    /* NOTE:
//    The function below will call CalculateMaterialResponseCauchy() by default and then (may)
//    call CalculateMaterialResponseKirchhoff() in the constitutive_law.*/
//    mConstitutiveLawVector->CalculateMaterialResponse(Values, rVariables.StressMeasure);
//
//    KRATOS_CATCH("")
//}
//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::CalculateOnIntegrationPoints(const Variable<bool>& rVariable,
    std::vector<bool>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > CalculateOnIntegrationPoints\n");
    //if (rValues.size() != 1)
    //    rValues.resize(1);
//
    //if (rVariable == CALCULATE_EXPLICIT_MP_STRESS)
    //{
    //    GeneralVariables Variables;
    //    this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);
    //    this->CalculateExplicitStresses(rCurrentProcessInfo, Variables);
    //    this->FinalizeStepVariables(Variables, rCurrentProcessInfo);
    //    rValues[0] = true;
    //}
    //else if (rVariable == EXPLICIT_MAP_GRID_TO_MP)
    //if (rVariable == EXPLICIT_MAP_GRID_TO_MP)
    //{
    //    MPMExplicitUtilities::UpdateGaussPointExplicit(rCurrentProcessInfo, *this);
    //    rValues[0] = true;
    //}
    //else if (rVariable == CALCULATE_MUSL_VELOCITY_FIELD)
    //{
    //    MPMExplicitUtilities::CalculateMUSLGridVelocity(rCurrentProcessInfo, *this);
    //    rValues[0] = true;
    //}
    //else
    //{
    //    KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    //}
}

void MPMUpdatedLagrangianVPVMS::CalculateOnIntegrationPoints(const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > CalculateOnIntegrationPoints\n");
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

void MPMUpdatedLagrangianVPVMS::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > CalculateOnIntegrationPoints\n");
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_DENSITY) {
        rValues[0] = mMP.density;
    }
    else if (rVariable == MP_PRESSURE) {
        rValues[0] = m_mp_pressure;
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

void MPMUpdatedLagrangianVPVMS::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > CalculateOnIntegrationPoints\n");
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

void MPMUpdatedLagrangianVPVMS::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    x CalculateOnIntegrationPoints\n");
//    if (rValues.size() != 1)
//        rValues.resize(1);
// 
//    if (rVariable == MP_CAUCHY_STRESS_VECTOR) {
//        rValues[0] = mMP.cauchy_stress_vector;
//    }
//    else if (rVariable == MP_ALMANSI_STRAIN_VECTOR) {
//        rValues[0] = mMP.almansi_strain_vector;
//    }
//    else
//    {
//        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
//    }
}

void MPMUpdatedLagrangianVPVMS::SetValuesOnIntegrationPoints(const Variable<int>& rVariable,
    const std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    x SetValuesOnIntegrationPoints(int)\n");
}

void MPMUpdatedLagrangianVPVMS::SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > SetValuesOnIntegrationPoints\n");
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_MASS) {
        mMP.mass = rValues[0];
    }
    else if (rVariable == MP_PRESSURE) {
        m_mp_pressure = rValues[0];
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

void MPMUpdatedLagrangianVPVMS::SetValuesOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
    const std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > SetValuesOnIntegrationPoints\n");
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

void MPMUpdatedLagrangianVPVMS::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
    const std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    x SetValuesOnIntegrationPoints(stress/strain vectors)\n");
//    KRATOS_ERROR_IF(rValues.size() > 1)
//        << "Only 1 value per integration point allowed! Passed values vector size: "
//        << rValues.size() << std::endl;
// 
//    if (rVariable == MP_CAUCHY_STRESS_VECTOR) {
//        mMP.cauchy_stress_vector = rValues[0];
//    }
//    else if (rVariable == MP_ALMANSI_STRAIN_VECTOR) {
//        mMP.almansi_strain_vector = rValues[0];
//    }
//    else
//    {
//        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
//    }
}

//************************************************************************************
//************************************************************************************

// template <class TElementData>
void MPMUpdatedLagrangianVPVMS::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > Calculate\n");
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}


void MPMUpdatedLagrangianVPVMS::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    /// Lumped projection terms
    //if (rVariable == RESPROJ_DISPL) {
    //    this->CalculateProjections(rCurrentProcessInfo);
    //} else {
    //printf("    > Calculate\n");
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    //}
}

// template <class TElementData>
void MPMUpdatedLagrangianVPVMS::Calculate(
    const Variable<Vector>& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > Calculate\n");
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

// template <class TElementData>
void MPMUpdatedLagrangianVPVMS::Calculate(
    const Variable<Matrix>& rVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > Calculate\n");
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}



void MPMUpdatedLagrangianVPVMS::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    //printf("    > InitializeGeneralVariable\n");
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    //const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const bool is_axisymmetric = (rCurrentProcessInfo.Has(IS_AXISYMMETRIC)) 
        ? rCurrentProcessInfo.GetValue(IS_AXISYMMETRIC)
        : false;
    const SizeType def_grad_dim = (is_axisymmetric)
        ? 3
        : dimension;
  
    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.detFT = 1;

    rVariables.B.resize(3, number_of_nodes * dimension, false ); //3?

    rVariables.F.resize(def_grad_dim, def_grad_dim, false );

    rVariables.F0.resize(def_grad_dim, def_grad_dim, false );

    rVariables.FT.resize(def_grad_dim, def_grad_dim, false );

    //rVariables.StrainVector.resize(strain_size, false );

    //rVariables.StressVector.resize(strain_size, false );

    rVariables.DN_DX.resize( number_of_nodes, dimension, false );

    // CurrentVel is the unknown variable. It represents the nodal delta velocity. When it is predicted is equal to zero.
    rVariables.CurrentVel = CalculateCurrentVel(rVariables.CurrentVel, rCurrentProcessInfo);
    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);

    // Initialize stabilization parameters
    rVariables.tau1  = 0;
    rVariables.tau2  = 0;

    // Set Pressure and Pressure Gradient in gauss points
    rVariables.PressureGP = 0;
    rVariables.PressureGradient = ZeroVector(dimension);

    //Set Velocity and Velocity Gradient in gauss points
    rVariables.VelocityGP = ZeroVector(dimension);
    rVariables.VelocityGradient = ZeroMatrix(dimension, dimension);

    // Set dynamic coefficients for stabilization
    rVariables.DynamicCoefficient = 0;
    rVariables.DynamicRHS = ZeroVector(dimension);

    // Set Identity matrices
    rVariables.Identity = IdentityMatrix(dimension);

    // Set Body forces
    rVariables.BodyForceMP = ZeroVector(dimension);
    
    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************


void MPMUpdatedLagrangianVPVMS::CalculateElementalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    //printf("    > CalculateElementalSystem\n");
    CalculateResidualVectorFlag = true; //temporary solution
    // Create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    const Vector& r_N = row(GetGeometry().ShapeFunctionsValues(), 0);
    const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
        : false;

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    //Flags &ConstitutiveLawOptions=Values.GetOptions();
    //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    if (!is_explicit)
    {
        //ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables, rCurrentProcessInfo);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables, Values, r_N);
        
        // Calculate Material Response
        /* NOTE:
        The function below will call CalculateMaterialResponseCauchy() by default and then (may)
        call CalculateMaterialResponseKirchhoff() in the constitutive_law.*/
        //mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

        /* NOTE:SetSpecificVariables
        The material points will have constant mass as defined at the beginning.
        However, the density and volume (integration weight) are changing every time step.*/
        //mMP.density = (GetProperties()[DENSITY]) / Variables.detFT;
        mMP.density = (GetProperties()[DENSITY])/ Variables.detFT;

        // Compute other variables needed for stabilization
        SetSpecificVariables(Variables,rCurrentProcessInfo);

        // Compute stabilization parameters
        CalculateTaus(rCurrentProcessInfo.GetValue(STABILIZATION_TYPE),Variables,rCurrentProcessInfo);
    }

    // The MP_Volume (integration weight) is evaluated
    mMP.volume = mMP.mass / mMP.density;

    if (CalculateStiffnessMatrixFlag && !is_explicit) // if calculation of the matrix is required
    {
        // Contributions to stiffness matrix calculated on the reference configuration
        this->CalculateAndAddLHS(
            rLeftHandSideMatrix,
            Variables,
            mMP.volume,
            rCurrentProcessInfo);
    }

    if (CalculateResidualVectorFlag) // if calculation of the vector is required
    {
        // Contribution to forces (in residual term) are calculated
        //#BODYFORCE
        //KRATOS_WATCH(mMP.volume_acceleration);
        //KRATOS_WATCH(mMP.mass);
        //KRATOS_WATCH(Variables.BodyForceMP);
        Vector volume_force = (mMP.volume_acceleration * mMP.mass ) + (Variables.BodyForceMP * mMP.mass); // caso statico: il primo termine deve essere zero.
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

void MPMUpdatedLagrangianVPVMS::SetSpecificVariables(GeneralVariables& rVariables,const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > SetSpecificVariables\n");
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    // Set Pressure and Pressure Gradient in gauss points
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rVariables.PressureGP += r_N(0,j) * r_geometry[j].FastGetSolutionStepValue(PRESSURE,1);
        for ( unsigned int i = 0; i < dimension; i++ )
        {
            rVariables.PressureGradient[i] += rVariables.DN_DX(j,i) * r_geometry[j].FastGetSolutionStepValue(PRESSURE,1);
        }
    }
    // ??
    // Set Velocity and Velocity Gradient in gauss points
    for ( unsigned int i = 0; i < number_of_nodes; i++)
    {
        rVariables.VelocityGP[0] += r_N(0,i) * r_geometry[i].FastGetSolutionStepValue(VELOCITY_X,1);
        rVariables.VelocityGP[1] += r_N(0,i) * r_geometry[i].FastGetSolutionStepValue(VELOCITY_Y,1);
        if (dimension == 3)
        {
            rVariables.VelocityGP[2] += r_N(0,i) * r_geometry[i].FastGetSolutionStepValue(VELOCITY_Z,1);
        }
        for (unsigned int j = 0; j < dimension; j++)
        {
            rVariables.VelocityGradient(j,0) = rVariables.DN_DX(i,j) * r_geometry[i].FastGetSolutionStepValue(VELOCITY_X,1);
            rVariables.VelocityGradient(j,1) = rVariables.DN_DX(i,j) * r_geometry[i].FastGetSolutionStepValue(VELOCITY_Y,1);
            if (dimension == 3)
            {
                rVariables.VelocityGradient(j,2) = rVariables.DN_DX(i,j) * r_geometry[i].FastGetSolutionStepValue(VELOCITY_Z,1);
            }
        }
    }

    // Set if the model is dynamic
    const bool is_dynamic = rCurrentProcessInfo.Has(IS_DYNAMIC)
        ? rCurrentProcessInfo.GetValue(IS_DYNAMIC)
        : false;
    const int current_step = rCurrentProcessInfo.GetValue(STEP);
    if (is_dynamic && current_step > 1) ComputeDynamicTerms(rVariables,rCurrentProcessInfo);

    // Compute Shear modulus and Bulk Modulus
    if (GetProperties().Has(YOUNG_MODULUS) && GetProperties().Has(POISSON_RATIO))
    {
        const double& young_modulus = GetProperties()[YOUNG_MODULUS];
        const double& poisson_ratio = GetProperties()[POISSON_RATIO];
        rVariables.ShearModulus = young_modulus / (2.0 * (1.0 + poisson_ratio));
        rVariables.BulkModulus  = young_modulus / (3.0 * (1.0 - 2.0 * poisson_ratio));

        const double tolerance = 10.e-7;
        const bool check = bool( (poisson_ratio > 0.5-tolerance ) || (poisson_ratio < (-1.0 + tolerance)) );
        if (POISSON_RATIO.Key() == 0 || check==true)  rVariables.BulkModulus= 1e16;

        if (rVariables.BulkModulus!=rVariables.BulkModulus)
            rVariables.BulkModulus= 1e16;
    }
    else if (GetProperties().Has(DYNAMIC_VISCOSITY))
    {
        rVariables.ShearModulus = GetProperties()[DYNAMIC_VISCOSITY];
        rVariables.BulkModulus = GetProperties()[BULK_MODULUS];
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

// Calulate element size depending on the dimension of worspace

void MPMUpdatedLagrangianVPVMS::ComputeElementSize(double& ElemSize){

    KRATOS_TRY
    //printf("    > ComputeElementSize\n");
    GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    double Area = r_geometry.GetGeometryParent(0).Area();

    if (dimension == 2)
    {
        ElemSize = 1.128379167 * sqrt(Area);
    }
    else if (dimension== 3)
    {
        ElemSize = 0.60046878 * pow(Area,0.333333333333333333333);
    }

        KRATOS_CATCH("")
}

// Calculate stabilization parameters

void MPMUpdatedLagrangianVPVMS::CalculateTaus(const int& stabilization_type,
    GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    //printf("    > CalculateTaus\n");
    // Add computations for the tau stabilization

    //const double constant1= 10.0;
    const double constant1 = rCurrentProcessInfo.GetValue(CONSTANT1);
    const double constant2 = rCurrentProcessInfo.GetValue(CONSTANT2);
    double characteristic_element_size;
    ComputeElementSize(characteristic_element_size);

    rVariables.tau1 = constant1 * pow(characteristic_element_size,2) / (rVariables.ShearModulus);
    rVariables.tau2 = constant2 * rVariables.ShearModulus;
    //rVariables.tau1 = constant1 * mMP.density * pow(characteristic_element_size,2) / (rVariables.ShearModulus);
    //rVariables.tau2 = constant2 * pow(characteristic_element_size,2) / rVariables.tau1;
    KRATOS_CATCH( "" )
}

void MPMUpdatedLagrangianVPVMS::ComputeDynamicTerms(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)//°
{
    // ONLY FOR NEWMARK APPROACH, to be revisited
    //printf("    > ComputeDynamicTerms\n");
    KRATOS_TRY
    GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];
    const double beta = 0.25; // Should be read from the general scheme mpm_residual_based_bossak_scheme

    array_1d<double,3> aux_MP_velocity = ZeroVector(3);
    array_1d<double,3> aux_MP_acceleration = ZeroVector(3);
    array_1d<double,3> aux_MP_displacement = ZeroVector(3);
    array_1d<double,3> current_MP_velocity = ZeroVector(3);


    for (unsigned int j=0; j<number_of_nodes; j++)
    {
        // These are the values of nodal displacement and nodal acceleration evaluated in the initialize solution step
        array_1d<double, 3 > previous_acceleration = ZeroVector(3);
        //if (r_geometry[j].SolutionStepsDataHas(ACCELERATION))
            previous_acceleration = r_geometry[j].FastGetSolutionStepValue(ACCELERATION,1);

        array_1d<double, 3 > previous_velocity = ZeroVector(3);
        //if (r_geometry[j].SolutionStepsDataHas(VELOCITY))
            previous_velocity = r_geometry[j].FastGetSolutionStepValue(VELOCITY,1);

        array_1d<double, 3 > previous_displacement = ZeroVector(3);
        previous_displacement = r_geometry[j].FastGetSolutionStepValue(DISPLACEMENT,1);

        for (unsigned int k = 0; k < dimension; k++)
        {
            aux_MP_velocity[k] += r_N(0, j) * previous_velocity[k];
            aux_MP_acceleration[k] += r_N(0, j) * previous_acceleration[k];
            aux_MP_displacement[k] += r_N(0, j) * previous_displacement[k];
        }
    }

    rVariables.DynamicCoefficient = 1 / (beta * delta_time * delta_time);
    const double coeff1 = 1 / (delta_time * 0.25);
    const double coeff2 = (0.5 - beta) / beta;

    rVariables.DynamicRHS =ZeroVector(3);

    for (unsigned int idime = 0; idime < dimension; idime++)
    {
           rVariables.DynamicRHS[idime] -= rVariables.DynamicCoefficient * aux_MP_displacement[idime];
           rVariables.DynamicRHS[idime] -= coeff1 * aux_MP_displacement[idime];
           rVariables.DynamicRHS[idime] += coeff1 * aux_MP_velocity[idime] + coeff2 * aux_MP_acceleration[idime];
/*
        rVariables.DynamicRHS[idime] -= rVariables.DynamicCoefficient * mMP.displacement[idime];
        rVariables.DynamicRHS[idime] -= coeff1 * mMP.displacement[idime];
        rVariables.DynamicRHS[idime] += coeff1 * mMP.velocity[idime] + coeff2 * mMP.acceleration[idime];*/
    }

//     rVariables.DynamicCoefficient =0;

 KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::CalculateAndAddRHS(
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    Vector& rVolumeForce,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > CalculateAndAddRHS\n");
    //KRATOS_WATCH(rRightHandSideVector);
    // Operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
    CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // To be removed
    //CalculateAndAddStabilizedVelocity( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight);

    //if (rCurrentProcessInfo.GetValue(STABILIZATION_TYPE)==1)
    // Operation performed: rRightHandSideVector -= Stabilized Pressure Forces
    CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight);
    //KRATOS_WATCH(rRightHandSideVector);
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY
    //printf("    > CalculateAndAddExternalForces\n");
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index_i = dimension * i + i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {           
            rRightHandSideVector[index_i + j] += r_N(0, i) * rVolumeForce[j]; // * rIntegrationWeight;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
//° mancano i termini dinamici
void MPMUpdatedLagrangianVPVMS::CalculateAndAddStabilizedVelocity(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY
    //printf("    > CalculateAndAddStabilizedVelocity\n");
    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void MPMUpdatedLagrangianVPVMS::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY
    //printf("    > CalculateAndAddStabilizedPressure\n");
    GeometryType& r_geometry = GetGeometry();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int index_i = dimension;
    
    Vector resizedVolumeForce;
    if (dimension==2) {
        resizedVolumeForce.resize(3, false);
        resizedVolumeForce = rVolumeForce; // * rIntegrationWeight;
        resizedVolumeForce.resize(2, true);
    } else {
        KRATOS_ERROR << "3d problem not yet implemented " << std::endl;
    }
    
    Vector StabilizedPressure;
    StabilizedPressure.resize(number_of_nodes, false);
    StabilizedPressure = prod(rVariables.DN_DX, resizedVolumeForce);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {   
        rRightHandSideVector[index_i] -= rVariables.tau1 * StabilizedPressure[i] / mMP.density;
        index_i += (dimension + 1);
    }
    //KRATOS_WATCH(rVariables.StressVector);
    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************
//° mancano i termini dinamici
void MPMUpdatedLagrangianVPVMS::CalculateAndAddLHS(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > CalculateAndAddLHS\n");
    /*
        K   Bp
        Bv  0
    */
    // Operation performed: add K to the rLefsHandSideMatrix
    CalculateAndAddK( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    //KRATOS_WATCH(rLeftHandSideMatrix);

    // Operation performed: add Bp to the rLefsHandSideMatrix
    CalculateAndAddBp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    //KRATOS_WATCH(rLeftHandSideMatrix);
    
    // Operation performed: add Bv to the rLefsHandSideMatrix
    CalculateAndAddBv( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    //KRATOS_WATCH(rLeftHandSideMatrix);

    //CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kuu stabilization to the rLefsHandSideMatrix
    CalculateAndAddSv( rLeftHandSideMatrix, rVariables, rIntegrationWeight, rCurrentProcessInfo ); //°i termini dinamici o vanno omessi, o si aggiungono altrove
    
    // Operations performed: add Kpp stabilization to the rLefsHandSideMatrix
    CalculateAndAddSp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    //KRATOS_WATCH(rLeftHandSideMatrix);

}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    double bulk_modulus  = 1.0;


   if (GetProperties().Has(YOUNG_MODULUS) && GetProperties().Has(POISSON_RATIO))
    {
        const double& young_modulus = GetProperties()[YOUNG_MODULUS];
        const double& poisson_ratio    = GetProperties()[POISSON_RATIO];
        bulk_modulus  = young_modulus/(3.0*(1.0-2.0*poisson_ratio));
        // Check if Bulk Modulus is not NaN
        if (bulk_modulus != bulk_modulus)
            bulk_modulus = 1e16;
    }
   else if (GetProperties().Has(DYNAMIC_VISCOSITY)) {
        bulk_modulus  = GetProperties()[BULK_MODULUS];
    }


    //double delta_coefficient = rVariables.detF0 - 1; //FLUID

    unsigned int indexpi = dimension;

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int indexpj = dimension;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
             rLeftHandSideMatrix(indexpi,indexpj)  -= ((1.0)/(bulk_modulus)) * r_N(0, i) * r_N(0, j) * rIntegrationWeight;
	    // FLUID-UP	
           //rLeftHandSideMatrix(indexpi, indexpj) -= ((1.0) / (bulk_modulus)) * r_N(0, i) * r_N(0, j) * rIntegrationWeight / (delta_coefficient * (rVariables.detF0 / rVariables.detF));

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}

void MPMUpdatedLagrangianVPVMS::CalculateAndAddK(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight
                                             )
{
    KRATOS_TRY
    //printf("    > CalculateAndAddK\n");
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    const int size = number_of_nodes * dimension;
    Matrix reduced_K = prod(rVariables.DN_DX, rIntegrationWeight * trans(rVariables.DN_DX));
    Matrix K = ZeroMatrix(size,size);
    MathUtils<double>::ExpandAndAddReducedMatrix( K, reduced_K, dimension );
    unsigned int indexi = 0;
    unsigned int indexj = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(indexi+i,indexj+j) += rVariables.ShearModulus * K(indexi,indexj);                    
                    indexj++;
                }
            }
            indexi++;
        }
    }

    KRATOS_CATCH( "" )
}

void MPMUpdatedLagrangianVPVMS::CalculateAndAddBp (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY
    //printf("    > CalculateAndAddBp\n");
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    unsigned int index_j;
    unsigned int index_i; 
    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        index_j  = dimension;
        index_i = dimension * i + i;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(index_i+k, index_j) += rVariables.DN_DX ( i, k ) * r_N(0, j) * rIntegrationWeight;
            }
            index_j += (dimension + 1);
        }
    }
    

    KRATOS_CATCH( "" )
}

void MPMUpdatedLagrangianVPVMS::CalculateAndAddBv (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY
    //printf("    > CalculateAndAddBv\n");
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    // Assemble components considering added DOF matrix system
    unsigned int indexi = dimension;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            unsigned int indexj = dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(indexi,indexj+k) += r_N(0, i) * rVariables.DN_DX ( j, k ) * rIntegrationWeight;
            }
        }
        indexi += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}

void MPMUpdatedLagrangianVPVMS::CalculateAndAddSv (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight,
        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    //printf("    > CalculateAndAddSv\n");
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(i*(dimension+1) + idim, j*(dimension+1) + jdim) -= rVariables.tau2 * rVariables.DN_DX(i,idim) * rVariables.DN_DX(j,jdim) * rIntegrationWeight;
                }
            }
        }
    }

    KRATOS_CATCH( "" )
}

void MPMUpdatedLagrangianVPVMS::CalculateAndAddSp (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)

{
    KRATOS_TRY
    //printf("    > CalculateAndAddSp\n");
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix K = prod(rVariables.DN_DX, rIntegrationWeight * trans(rVariables.DN_DX));

    unsigned int indexi = 0;
    unsigned int indexj = 0;
    unsigned int lhs1 = 0;
    unsigned int lhs2 = 0;
    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        indexj=0;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            rLeftHandSideMatrix(dimension + i*(dimension+1), dimension + j*(dimension+1)) -= rVariables.tau1 * K(indexi,indexj);
            indexj++;     
        }
        indexi++;
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

//void MPMUpdatedLagrangianVPVMS::SetGeneralVariables(GeneralVariables& rVariables, const Vector& rN)
void MPMUpdatedLagrangianVPVMS::SetGeneralVariables(GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues, const Vector& rN)
{
    //printf("    > SetGeneralVariables\n");
    GeometryType& r_geometry = GetGeometry();
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    // Variables.detF is the determinant of the incremental total deformation gradient
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    // Check if detF is negative (element is inverted)
    if(rVariables.detF<0)
    {
        KRATOS_INFO("MPMUpdatedLagrangianVPVMS")<<" Element: "<<this->Id()<<std::endl;
        KRATOS_INFO("MPMUpdatedLagrangianVPVMS")<<" Element position: "<< mMP.xg <<std::endl;
        KRATOS_INFO("MPMUpdatedLagrangianVPVMS")<<" Element velocity: "<< mMP.velocity <<std::endl;

        KRATOS_INFO("MPMUpdatedLagrangianVPVMS") << " Shape functions: " << r_geometry.ShapeFunctionsValues() << std::endl;
        KRATOS_INFO("MPMUpdatedLagrangianVPVMS") << " Quadrature points: " << r_geometry.IntegrationPointsNumber() << std::endl;
        KRATOS_INFO("MPMUpdatedLagrangianVPVMS") << " Parent geometry ID: " << r_geometry.GetGeometryParent(0).Id() << std::endl;
        KRATOS_INFO("MPMUpdatedLagrangianVPVMS") << " Parent geometry number of points: " << r_geometry.GetGeometryParent(0).PointsNumber() << std::endl;

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            const array_1d<double, 3> & current_position      = r_geometry[i].Coordinates();
            const array_1d<double, 3> & current_displacement  = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3> & previous_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            const array_1d<double, 3> & current_velocity      = r_geometry[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3> & previous_velocity     = r_geometry[i].FastGetSolutionStepValue(VELOCITY,1);

            
            KRATOS_INFO("MPMUpdatedLagrangianVMS")<<" NODE ["<<r_geometry[i].Id()<<"]: (Current position: "<<current_position<<") "<<std::endl;
            KRATOS_INFO("MPMUpdatedLagrangianVMS")<<" ---Current Disp: "<<current_displacement<<" (Previous Disp: "<<previous_displacement<<")"<<std::endl;
            KRATOS_INFO("MPMUpdatedLagrangianVMS")<<" ---Current Vel: "<<current_velocity<<" (Previous Vel: "<<previous_velocity<<")"<<std::endl;

        }

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if( r_geometry[i].SolutionStepsDataHas(CONTACT_FORCE) )
            {
                const array_1d<double, 3 > & PreContactForce = r_geometry[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
                const array_1d<double, 3 > & ContactForce = r_geometry[i].FastGetSolutionStepValue(CONTACT_FORCE);
                KRATOS_INFO("MPMUpdatedLagrangianVMS")<<" ---Contact_Force: (Pre:"<<PreContactForce<<", Current:"<<ContactForce<<") "<<std::endl;
            }
            else
            {
                KRATOS_INFO("MPMUpdatedLagrangianVMS")<<" ---Contact_Force: NULL "<<std::endl;
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
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rN);

    // Body forces
    rVariables.BodyForceMP = ZeroVector(3);
    array_1d<double, 3 > nodal_body_force = ZeroVector(3);
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        if (r_geometry[j].SolutionStepsDataHas(BODY_FORCE))
        {

            nodal_body_force = r_geometry[j].FastGetSolutionStepValue(BODY_FORCE,0);
            for (unsigned int k = 0; k < dimension; k++)
            {
                rVariables.BodyForceMP[k] += r_N(0, j) * nodal_body_force[k];
            }
        }
    }

}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void MPMUpdatedLagrangianVPVMS::CalculateKinematics(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)//°
{
    KRATOS_TRY
    //printf("    > CalculateKinematics\n");
    // Define the stress measure
    //rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    // Calculating the reference jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n/d£]
    Matrix Jacobian;
    GetGeometry().Jacobian(Jacobian, 0);

    // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    double detJ;
    MathUtils<double>::InvertMatrix( Jacobian, InvJ, detJ);

    // Calculating the current jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n+1/d£]
    Matrix jacobian;
    GetGeometry().Jacobian(jacobian, 0, GetGeometry().GetDefaultIntegrationMethod(), -1.0 * rVariables.CurrentDisp);

    // Calculating the inverse of the jacobian and the parameters needed [d£/(dx_n+1)]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( jacobian, Invj, detJ); //overwrites detJ

    // Compute cartesian derivatives [dN/dx_n+1]
    Matrix r_DN_De = GetGeometry().ShapeFunctionLocalGradient(0);
    rVariables.DN_DX = prod(r_DN_De, Invj); //overwrites DX now is the current position dx

    /* NOTE::
    Deformation Gradient F [(dx_n+1 - dx_n)/dx_n] is to be updated in constitutive law parameter as total deformation gradient.
    The increment of total deformation gradient can be evaluated in 2 ways, which are:
    1. By: noalias( rVariables.F ) = prod( jacobian, InvJ);
    2. By means of the gradient of nodal displacement: using this second expression quadratic convergence is not guarantee

    (NOTICE: Here, we are using method no. 1)
    */

    // Update Deformation gradient
    noalias( rVariables.F ) = prod( jacobian, InvJ);

    // Determinant of the previous Deformation Gradient F_n
    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;

    // Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);

    KRATOS_CATCH( "" )
}
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rF,
        Matrix& rDN_DX)
{
    KRATOS_TRY
    //printf("    x CalculateDeformationMatrix\n");
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rB.clear(); // Set all components to zero

    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;
            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );
        }
    }
    else if( dimension == 3 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 3 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 2 ) = rDN_DX( i, 2 );

            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

            rB( 4, index + 1 ) = rDN_DX( i, 2 );
            rB( 4, index + 2 ) = rDN_DX( i, 1 );

            rB( 5, index + 0 ) = rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rDN_DX( i, 0 );
        }
    }
    else
    {
        KRATOS_ERROR << "Dimension given is wrong!" << std::endl;
    }

    KRATOS_CATCH( "" )
}

//*************************COMPUTE CURRENT DISPLACEMENT*******************************
//************************************************************************************
/*
This function convert the computed nodal displacement into matrix of (number_of_nodes, dimension)
*/
Matrix& MPMUpdatedLagrangianVPVMS::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    //printf("    > CalculateCurrentDisp\n");
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

//*************************COMPUTE CURRENT VELOCITY*******************************
//************************************************************************************
/*
This function convert the computed nodal veloity into matrix of (number_of_nodes, dimension)
*/
Matrix& MPMUpdatedLagrangianVPVMS::CalculateCurrentVel(Matrix & rCurrentVel, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    //printf("    > CalculateCurrentVel\n");
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    rCurrentVel = ZeroMatrix(number_of_nodes, dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3 > & current_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rCurrentVel(i,j) = current_velocity[j];
        }
    }

    return rCurrentVel;

    KRATOS_CATCH( "" )
}
//*************************COMPUTE ALMANSI STRAIN*************************************
//************************************************************************************
// Almansi Strain: E = 0.5 (I - U^(-2))
//void MPMUpdatedLagrangianVPVMS::CalculateAlmansiStrain(const Matrix& rF,
//        Vector& rStrainVector )
//{
//    KRATOS_TRY
//
//    KRATOS_CATCH( "" )
//}

//*************************COMPUTE GREEN-LAGRANGE STRAIN*************************************
//************************************************************************************
// Green-Lagrange Strain: E = 0.5 * (U^2 - I) = 0.5 * (C - I)
//void MPMUpdatedLagrangianVPVMS::CalculateGreenLagrangeStrain(
//    const Matrix& rF,
//    Vector& rStrainVector)
//{
//    KRATOS_TRY
//
//    KRATOS_CATCH( "" )
//}

//*************************DECIMAL CORRECTION OF STRAINS******************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::DecimalCorrection(Vector& rVector)
{
    KRATOS_TRY
    //printf("    > DecimalCorrection\n");
    for ( unsigned int i = 0; i < rVector.size(); i++ )
    {
        if( rVector[i]*rVector[i]<1e-24 )
        {
            rVector[i]=0;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    > FinalizeStepVariables\n");
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    // Update internal (historical) variables
    mDeterminantF0         = rVariables.detF* rVariables.detF0;
    mDeformationGradientF0 = prod(rVariables.F, rVariables.F0);

   // mMP.cauchy_stress_vector = rVariables.StressVector;
  // mMP.almansi_strain_vector = rVariables.StrainVector;

    // Delta Plastic Strains
/*     if (mConstitutiveLawVector->Has(MP_DELTA_PLASTIC_STRAIN))
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
 */
    const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
        : false;
    if (!is_explicit) this->UpdateGaussPoint(rVariables, rCurrentProcessInfo);

    //// Evaluation of the pressure on the material point
    //double nodal_mean_stress = 0.0;
    //for (unsigned int i = 0; i < number_of_nodes; i++)
    //    nodal_mean_stress += r_geometry[i].FastGetSolutionStepValue( PRESSURE ) * r_N(0, i);
//
    //// Evaluation of the mean stress on the material point
    //double mean_stress = 0.0;
    //for (unsigned int i = 0; i < dimension; i++)
    //    mean_stress += rVariables.StressVector[i];
    //mean_stress /= dimension;
//
    //Vector stress_vector = ZeroVector(voigtsize);
    //stress_vector = rVariables.StressVector;
    //for (unsigned int i = 0; i < dimension; i++)
    //    stress_vector[i] += (nodal_mean_stress - mean_stress);
//
    //mMP.cauchy_stress_vector = stress_vector;

}

//************************************************************************************
//************************************************************************************
/**
 * The position of the Gauss points/Material points is updated
 */

void MPMUpdatedLagrangianVPVMS::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    //printf("    > UpdateGaussPoint\n");
    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo); //°?
    rVariables.CurrentVel = CalculateCurrentVel(rVariables.CurrentVel, rCurrentProcessInfo);    
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    unsigned int dimension = r_geometry.WorkingSpaceDimension();

    //array_1d<double,3> delta_xg = ZeroVector(3);
    array_1d<double,3> MP_xg = ZeroVector(3);
    array_1d<double,3> MP_acceleration = ZeroVector(3);
    //array_1d<double,3> MP_velocity = ZeroVector(3);
    array_1d<double,3> delta_velocity = ZeroVector(3);
    double MP_pressure = 0.0;
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];

    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (r_N(0, i) > std::numeric_limits<double>::epsilon())
        {
            auto r_geometry = GetGeometry();
            array_1d<double, 3 > nodal_acceleration = ZeroVector(3);
            if (r_geometry[i].SolutionStepsDataHas(ACCELERATION))
                nodal_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION);

            const double& nodal_pressure = r_geometry[i].FastGetSolutionStepValue(PRESSURE, 0);
            MP_pressure += r_N(0, i) * nodal_pressure;

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                delta_velocity[j] += r_N(0, i) * rVariables.CurrentVel(i,j);
                MP_acceleration[j] += r_N(0, i) * nodal_acceleration[j];

                /* NOTE: The following interpolation techniques have been tried:
                    MP_velocity[j]      += rVariables.N[i] * nodal_velocity[j];
                    MP_acceleration[j]  += nodal_inertia[j]/(rVariables.N[i] * MP_mass * MP_number);
                    MP_velocity[j]      += nodal_momentum[j]/(rVariables.N[i] * MP_mass * MP_number);
                    MP_velocity[j]      += delta_time * rVariables.N[i] * nodal_acceleration[j];
                */
            }
        }
    }

    /* NOTE:
    Another way to update the MP velocity (see paper Guilkey and Weiss, 2003).
    This assume newmark (or trapezoidal, since n.gamma=0.5) rule of integration*/
    //mMP.velocity = mMP.velocity + 0.5 * delta_time * (MP_acceleration + mMP.acceleration);
    mMP.xg += delta_time * mMP.velocity + pow(delta_time, 2)/4 * (mMP.acceleration + MP_acceleration);
    mMP.displacement += delta_time * mMP.velocity + pow(delta_time, 2)/4 * (mMP.acceleration + MP_acceleration);
    /* NOTE: The following interpolation techniques have been tried:
        MP_acceleration = 4/(delta_time * delta_time) * delta_xg - 4/delta_time * MP_PreviousVelocity;
        MP_velocity = 2.0/delta_time * delta_xg - MP_PreviousVelocity;
    */

    // Update the MP Pressure
    m_mp_pressure = MP_pressure;

    // Update the MP Position
    //mMP.xg += delta_xg ;
    mMP.velocity += delta_velocity;

    //Update the MP Acceleration
    mMP.acceleration = MP_acceleration;

    // Update the MP total displacement
    //mMP.displacement += delta_xg;

    KRATOS_CATCH( "" )
}

void MPMUpdatedLagrangianVPVMS::InitializeMaterial(const ProcessInfo& rCurrentProcessInfo)
{
    //printf("    x InitializeMaterial\n");
//    KRATOS_TRY
//    GeneralVariables Variables;
//
//    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
//    {
//        mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();
//        Vector N = row(GetGeometry().ShapeFunctionsValues(), 0);
//        mConstitutiveLawVector->InitializeMaterial(
//            GetProperties(), GetGeometry(), N);
//
//        mMP.almansi_strain_vector = ZeroVector(mConstitutiveLawVector->GetStrainSize());
//        mMP.cauchy_stress_vector = ZeroVector(mConstitutiveLawVector->GetStrainSize());
//        
//        // Resize the deformation gradient if we are axisymmetric
//        if (mConstitutiveLawVector->GetStrainSize() == 4) mDeformationGradientF0 = IdentityMatrix(3);
//    }
//    else
//        KRATOS_ERROR <<  "A constitutive law needs to be specified for the element with ID: " << this->Id() << std::endl;
//
//    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::ResetConstitutiveLaw()
{
    //printf("    x ResetConstitutiveLaw\n");
//    KRATOS_TRY
//    GeneralVariables Variables;//

//    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
//    {
//        mConstitutiveLawVector->ResetMaterial(
//            GetProperties(),
//            GetGeometry(),
//            row(GetGeometry().ShapeFunctionsValues(), 0));
//    }//

//    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianVPVMS::GetHistoricalVariables( GeneralVariables& rVariables )
{
    //printf("    > GetHistoricalVariables\n");
    //Deformation Gradient F ( set to identity )
    unsigned int size =  rVariables.F.size1();
    rVariables.detF  = 1;
    rVariables.F     = IdentityMatrix(size);

    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;
}

//************************************************************************************
//************************************************************************************

double& MPMUpdatedLagrangianVPVMS::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    //printf("    > CalculateIntegrationWeight\n");
    if( dimension == 2 )
        rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
}

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& MPMUpdatedLagrangianVPVMS::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY
    //printf("    > CalculateVolumeChange\n");
    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

void MPMUpdatedLagrangianVPVMS::CalculateDeformationGradient(const Matrix& rDN_DX, Matrix& rF, Matrix& rDisplacement,
    const bool IsAxisymmetric)
{
    KRATOS_TRY
    //printf("    > CalculateDeformationGradient\n");
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (IsAxisymmetric)
    {
        // Compute radius
        const double current_radius = MPMMathUtilities<double>::CalculateRadius(GetGeometry().ShapeFunctionsValues(), GetGeometry());
        const double initial_radius = MPMMathUtilities<double>::CalculateRadius(GetGeometry().ShapeFunctionsValues(), GetGeometry(), Initial);

        rF = IdentityMatrix(3);

        if (dimension == 2)
        {
            for (IndexType i = 0; i < GetGeometry().PointsNumber(); ++i)
            {
                rF(0, 0) += rDisplacement(i, 0) * rDN_DX(i, 0);
                rF(0, 1) += rDisplacement(i, 0) * rDN_DX(i, 1);
                rF(1, 0) += rDisplacement(i, 1) * rDN_DX(i, 0);
                rF(1, 1) += rDisplacement(i, 1) * rDN_DX(i, 1);
            }

            rF(2, 2) = current_radius / initial_radius;
        }
        else KRATOS_ERROR << "Dimension given is wrong!" << std::endl;
    }
    else
    {
        /* NOTE::
    Deformation Gradient F [(dx_n+1 - dx_n)/dx_n] is to be updated in constitutive law parameter as total deformation gradient.
    The increment of total deformation gradient can be evaluated in 2 ways, which are:
    1. By: noalias( rVariables.F ) = prod( jacobian, InvJ);
    2. By means of the gradient of nodal displacement: using this second expression quadratic convergence is not guarantee

    (NOTICE: Here, we are using method no. 2)*/ //°°

    // METHOD 1: Update Deformation gradient: F [dx_n+1/dx_n] = [dx_n+1/d£] [d£/dx_n]
    // noalias( rVariables.F ) = prod( jacobian, InvJ);

    // METHOD 2: Update Deformation gradient: F_ij = δ_ij + u_i,j

        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        Matrix I = IdentityMatrix(dimension);
        Matrix gradient_displacement = ZeroMatrix(dimension, dimension);
        gradient_displacement = prod(trans(rDisplacement), rDN_DX);
        noalias(rF) = (I + gradient_displacement);
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

//void MPMUpdatedLagrangianVPVMS::ComputeResidual(GeneralVariables& rVariables, Vector& rVolumeForce, Vector& rResidualV, double& rResidualP)
//{
//
//    GeometryType& r_geometry = this->GetGeometry();
//    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
//
//    for (unsigned int k = 0; k < dimension; ++k) rResidualV[k] =  rVolumeForce[k] + rVariables.PressureGradient[k];
//    rResidualP =  -rVariables.PressureGP/rVariables.BulkModulus + (rVariables.detFT*rVariables.detF0-1); //°
//
//}
//
void MPMUpdatedLagrangianVPVMS::save( Serializer& rSerializer ) const
{
    //printf("    > save\n");
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMUpdatedLagrangianVPVMS )
    rSerializer.save("Pressure",m_mp_pressure);
}

void MPMUpdatedLagrangianVPVMS::load( Serializer& rSerializer )
{
    //printf("    > load\n");
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMUpdatedLagrangianVPVMS )
    rSerializer.load("Pressure",m_mp_pressure);
}

} // Namespace Kratos

