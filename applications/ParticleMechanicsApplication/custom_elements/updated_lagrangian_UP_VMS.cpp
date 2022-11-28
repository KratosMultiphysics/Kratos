//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Laura Moreno
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_UP_VMS.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application_variables.h"
#include "includes/checks.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{


    const unsigned int UpdatedLagrangianUPVMS::msIndexVoigt3D6C [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };
    const unsigned int UpdatedLagrangianUPVMS::msIndexVoigt2D4C [4][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1} };
    const unsigned int UpdatedLagrangianUPVMS::msIndexVoigt2D3C [3][2] = { {0, 0}, {1, 1}, {0, 1} };

///**
//* Flags related to the element computation
//*/
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_RHS_VECTOR,                 0 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_LHS_MATRIX,                 1 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUPVMS::UpdatedLagrangianUPVMS()
    : UpdatedLagrangian()
    , m_mp_pressure(1.0)
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianUPVMS::UpdatedLagrangianUPVMS( IndexType NewId, GeometryType::Pointer pGeometry )
    : UpdatedLagrangian( NewId, pGeometry )
    , m_mp_pressure(1.0)
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUPVMS::UpdatedLagrangianUPVMS( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : UpdatedLagrangian( NewId, pGeometry, pProperties )
    , m_mp_pressure(1.0)
{
    mFinalizedStep = true;


}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianUPVMS::UpdatedLagrangianUPVMS( UpdatedLagrangianUPVMS const& rOther)
    :UpdatedLagrangian(rOther)
    , m_mp_pressure(rOther.m_mp_pressure)
     //,mDeformationGradientF0(rOther.mDeformationGradientF0)
     //,mDeterminantF0(rOther.mDeterminantF0)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianUPVMS&  UpdatedLagrangianUPVMS::operator=(UpdatedLagrangianUPVMS const& rOther)
{
    UpdatedLagrangian::operator=(rOther);

    m_mp_pressure = rOther.m_mp_pressure;

    return *this;
}
//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianUPVMS::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangianUPVMS( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer UpdatedLagrangianUPVMS::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< UpdatedLagrangianUPVMS >(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianUPVMS::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    UpdatedLagrangianUPVMS NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.m_mp_pressure = m_mp_pressure;

    NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();

    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mDeterminantF0 = mDeterminantF0;

    return Element::Pointer( new UpdatedLagrangianUPVMS(NewElement) );
}
//*******************************DESTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianUPVMS::~UpdatedLagrangianUPVMS()
{
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialize parameters
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    mDeterminantF0 = 1;
    mDeformationGradientF0 = IdentityMatrix(dimension);

    // Initialize constitutive law and materials
    InitializeMaterial();

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    UpdatedLagrangian::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);
    GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int voigt_dimension = (dimension-1)*(dimension-1)+2;


    // Initialize stabilization parameters
    rVariables.tau1  = 0;
    rVariables.tau2  = 0;

    // Set Pressure and Pressure Gradient in gauss points
    rVariables.PressureGP = 0;
    rVariables.PressureGradient = ZeroVector(dimension);

    // Set dynamic coefficients for stabilization
    rVariables.DynamicCoefficient = 0;
    rVariables.DynamicRHS = ZeroVector(dimension);

    // Set Identity matrices
    rVariables.Identity = IdentityMatrix(dimension);
    rVariables.TensorIdentityMatrix = ZeroMatrix(voigt_dimension,voigt_dimension);

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************


void UpdatedLagrangianUPVMS::CalculateElementalSystem(
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
    const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
        : false;

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    if (!is_explicit)
    {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables, rCurrentProcessInfo);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables, Values, r_N);

        // Calculate Material Response
        /* NOTE:
        The function below will call CalculateMaterialResponseCauchy() by default and then (may)
        call CalculateMaterialResponseKirchhoff() in the constitutive_law.*/
        mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

        /* NOTE:SetSpecificVariables
        The material points will have constant mass as defined at the beginning.
        However, the density and volume (integration weight) are changing every time step.*/
        // Update MP_Density
        mMP.density = (GetProperties()[DENSITY]) / Variables.detFT;


        // Compute other variables needed for stabilization
        SetSpecificVariables(Variables,rCurrentProcessInfo);

        // Compute stabilization parameters
        CalculateTaus(rCurrentProcessInfo.GetValue(STABILIZATION_TYPE),Variables);

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
/**
 * The position of the Gauss points/Material points is updated
 */

void UpdatedLagrangianUPVMS::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    unsigned int dimension = r_geometry.WorkingSpaceDimension();

    array_1d<double,3> delta_xg = ZeroVector(3);
    array_1d<double,3> MP_acceleration = ZeroVector(3);
    array_1d<double,3> MP_velocity = ZeroVector(3);
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
                delta_xg[j] += r_N(0, i) * rVariables.CurrentDisp(i,j);
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
    mMP.velocity = mMP.velocity + 0.5 * delta_time * (MP_acceleration + mMP.acceleration);

    /* NOTE: The following interpolation techniques have been tried:
        MP_acceleration = 4/(delta_time * delta_time) * delta_xg - 4/delta_time * MP_PreviousVelocity;
        MP_velocity = 2.0/delta_time * delta_xg - MP_PreviousVelocity;
    */

    // Update the MP Pressure
    m_mp_pressure = MP_pressure;

    // Update the MP Position
    mMP.xg += delta_xg ;

    //Update the MP Acceleration
    mMP.acceleration = MP_acceleration;

    // Update the MP total displacement
    mMP.displacement += delta_xg;

    KRATOS_CATCH( "" )
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateKinematics(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Define the stress measure
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

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

void UpdatedLagrangianUPVMS::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rF,
        Matrix& rDN_DX)
{
    KRATOS_TRY

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

////************************************************************************************
////************************************************************************************

void UpdatedLagrangianUPVMS::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/

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

        // These are the values of nodal pressure evaluated in the initialize solution step
        const double& nodal_pressure = r_geometry[j].FastGetSolutionStepValue(PRESSURE,1);

        aux_MP_pressure += r_N(0, j) * nodal_pressure;

        for (unsigned int k = 0; k < dimension; k++)
        {
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
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::SetSpecificVariables(GeneralVariables& rVariables,const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    // Set Pressure and Pressure Gradient in gauss points
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rVariables.PressureGP += r_N(0,j) * r_geometry[j].FastGetSolutionStepValue(PRESSURE);
        for ( unsigned int i = 0; i < dimension; i++ )
        {
            rVariables.PressureGradient[i] += rVariables.DN_DX(j,i) * r_geometry[j].FastGetSolutionStepValue(PRESSURE);
        }
    }

    ConvertPressureGradientInVoigt(rVariables.PressureGradient,rVariables.PressureGradientVoigt);

    CalculateTensorIdentityMatrix(rVariables,rVariables.TensorIdentityMatrix);

    const bool is_dynamic = rCurrentProcessInfo.Has(IS_DYNAMIC)
        ? rCurrentProcessInfo.GetValue(IS_DYNAMIC)
        : false;

    //if (is_dynamic) ComputeDynamicTerms(rVariables,rCurrentProcessInfo);


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
        rVariables.BulkModulus = 1e16;
    }



    // Compute Residual if the stabilization is OSGS type
    if (rCurrentProcessInfo.GetValue(STABILIZATION_TYPE)==3)
    {
     Vector volume_force = mMP.volume_acceleration * mMP.mass;
     ComputeResidual(rVariables,volume_force,rVariables.GPResidualU,rVariables.GPResidualP);
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

// Calulate element size depending on the dimension of worspace

void UpdatedLagrangianUPVMS::ComputeElementSize(double& ElemSize){

    KRATOS_TRY

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

void UpdatedLagrangianUPVMS::CalculateTaus(const int& stabilization_type,
    GeneralVariables& rVariables)
{
    KRATOS_TRY

    // Add computations for the tau stabilization

    const double constant1=1.0;
    const double constant2=1.0;
    double characteristic_element_size;
    ComputeElementSize(characteristic_element_size);

    rVariables.tau1 = constant1 *  pow(characteristic_element_size,2) / (2 * rVariables.ShearModulus);
    rVariables.tau2 = 2 * constant2 * rVariables.ShearModulus;

    if ( stabilization_type == 1)
     {
         rVariables.tau1 = 0;
         rVariables.tau2 = 0;
     }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************


void UpdatedLagrangianUPVMS::CalculateTensorIdentityMatrix (GeneralVariables& rVariables, Matrix& rTensorIdentityMatrix)
{

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double voigt_dimension = (dimension - 1)*(dimension - 1)+2;
    rTensorIdentityMatrix=ZeroMatrix(voigt_dimension,voigt_dimension);

    if (dimension == 2)
    {
        for(unsigned int i=0; i<voigt_dimension; i++)
        {
            for(unsigned int j=0; j<voigt_dimension; j++)
            {
            rTensorIdentityMatrix( i, j ) = TensorIdentityComponent(rTensorIdentityMatrix( i, j ), rVariables,
                                          this->msIndexVoigt2D3C[i][0],this-> msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
            }
        }

    }
    else
    {
        for(unsigned int i=0; i<voigt_dimension; i++)
        {
            for(unsigned int j=0; j<voigt_dimension; j++)
            {
            rTensorIdentityMatrix( i, j ) = TensorIdentityComponent(rTensorIdentityMatrix( i, j ), rVariables,
                                          msIndexVoigt3D6C[i][0], msIndexVoigt3D6C[i][1], msIndexVoigt3D6C[j][0], msIndexVoigt3D6C[j][1]);
            }
        }

    }

}

//************************************************************************************
//************************************************************************************

double& UpdatedLagrangianUPVMS::TensorIdentityComponent (double& rCabcd, GeneralVariables& rVariables,
    const unsigned int& a, const unsigned int& b, const unsigned int& c, const unsigned int& d)
{

    double IdotI = rVariables.Identity(a, b) * rVariables.Identity(c, d);
    double Isym = (rVariables.Identity(a, c) * rVariables.Identity(b, d) +
        rVariables.Identity(a, d) * rVariables.Identity(b, c)) / 2.0;

    rCabcd += IdotI;
    rCabcd -= 2 * Isym;

    return rCabcd;
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::ConvertPressureGradientInVoigt(Vector& PressureGradient,Vector& PressureGradientVoigt)
{
    KRATOS_TRY
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int voigt_dimension=(dimension-1)*(dimension-1)+2;
    PressureGradientVoigt=ZeroVector(voigt_dimension);

   if( dimension == 2 )
    {
        PressureGradientVoigt(0) = PressureGradient(0);
        PressureGradientVoigt(1) = PressureGradient(1);
        PressureGradientVoigt(2) = 0.5*(PressureGradient(0) + PressureGradient(1));

    }
    else if( dimension == 3 )
    {

        PressureGradientVoigt(0) = PressureGradient(0);
        PressureGradientVoigt(1) = PressureGradient(1);
        PressureGradientVoigt(2) = PressureGradient(2);
        PressureGradientVoigt(3) = 0.5*(PressureGradient(0) + PressureGradient(1));
        PressureGradientVoigt(4) = 0.5*(PressureGradient(2) + PressureGradient(1));
        PressureGradientVoigt(5) = 0.5*(PressureGradient(2) + PressureGradient(0));

    }

    //return PressureGradientVoigt;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::ComputeDynamicTerms(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    // ONLY FOR NEWMARK APPROACH

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
    array_1d<double,3> current_MP_displacement = ZeroVector(3);


    for (unsigned int j=0; j<number_of_nodes; j++)
    {
        // These are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
        array_1d<double, 3 > previous_acceleration = ZeroVector(3);
        if (r_geometry[j].SolutionStepsDataHas(ACCELERATION))
            previous_acceleration = r_geometry[j].GetSolutionStepValue(ACCELERATION,1);

        array_1d<double, 3 > previous_velocity = ZeroVector(3);
        if (r_geometry[j].SolutionStepsDataHas(VELOCITY))
            previous_velocity = r_geometry[j].GetSolutionStepValue(VELOCITY,1);

        array_1d<double, 3 > previous_displacement = ZeroVector(3);
        previous_displacement = r_geometry[j].GetSolutionStepValue(DISPLACEMENT,1);

        array_1d<double, 3 > nodal_displacement = ZeroVector(3);
        nodal_displacement = r_geometry[j].FastGetSolutionStepValue(DISPLACEMENT,0);

        for (unsigned int k = 0; k < dimension; k++)
        {
            aux_MP_velocity[k] += r_N(0, j) * previous_velocity[k];
            aux_MP_acceleration[k] += r_N(0, j) * previous_acceleration[k];
            aux_MP_displacement[k] += r_N(0, j) * previous_displacement[k];
            current_MP_displacement[k] += r_N(0, j) * nodal_displacement[k];
        }
    }

    rVariables.DynamicCoefficient = 1 / (beta * delta_time * delta_time);

    const double coeff1 = 1 / (delta_time * 0.25);
    const double coeff2 = (0.5 - beta) / beta;

    rVariables.DynamicRHS =ZeroVector(3);

    for (unsigned int idime = 0; idime < dimension; idime++)
    {
        //rVariables.DynamicRHS[idime] -= rVariables.DynamicCoefficient * current_MP_displacement[idime];
        rVariables.DynamicRHS[idime] += rVariables.DynamicCoefficient * aux_MP_displacement[idime];
        rVariables.DynamicRHS[idime] += coeff1 * aux_MP_velocity[idime] + coeff2 * aux_MP_acceleration[idime];
    }


 KRATOS_CATCH( "" )
}


void UpdatedLagrangianUPVMS::ComputeResidual(GeneralVariables& rVariables, Vector& rVolumeForce, Vector& rResidualU, double& rResidualP)
{
    // Only for OSGS stabilization
    KRATOS_TRY

    rResidualU = rVolumeForce - rVariables.PressureGradient + mMP.density *mMP.acceleration;
    rResidualP = rVariables.PressureGP/rVariables.BulkModulus -(1.0 - 1.0 / rVariables.detFT);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddRHS(
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    Vector& rVolumeForce,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    //rVariables.detF0   *= rVariables.detF;
    //double determinant_F = rVariables.detF;
    //rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    // Operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
    CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // Operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
    CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // Operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
    CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // Operation performed: rRightHandSideVector -= Stabilized terms of the momentum equation
    CalculateAndAddStabilizedDisplacement( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight);

    // Operation performed: rRightHandSideVector -= Stabilized Pressure Forces
    CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight);

    //rVariables.detF     = determinant_F;
    //rVariables.detF0   /= rVariables.detF;

}
//************************************************************************************
//*********************Calculate the contribution of external force*******************

void UpdatedLagrangianUPVMS::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
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
        int index_up = dimension * i + i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[index_up + j] += r_N(0, i) * rVolumeForce[j];
        }
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    VectorType internal_forces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index_up = dimension * i + i;
        unsigned int index_u  = dimension * i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[index_up + j] -= internal_forces[index_u + j];
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int index_p = dimension;
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();



    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {

        rRightHandSideVector[index_p] += ((rVariables.PressureGP/rVariables.BulkModulus)-(1.0 - 1.0 / rVariables.detFT)) * r_N(0, i) * rIntegrationWeight;

        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddStabilizedDisplacement(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    //unsigned int index_p = dimension;
    Vector Testf1 = prod(rVariables.DN_DX, rVariables.PressureGradient);
    // Vector Stab1 = rVolumeForce + rVariables.PressureGradient;
    Vector Testf2   = prod(prod(trans(rVariables.B),rVariables.TensorIdentityMatrix),rVariables.PressureGradientVoigt);
    unsigned int indexi  = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index_up = dimension * i + i;
        unsigned int index_u  = dimension * i;

        for ( unsigned int jdim = 0; jdim < dimension; jdim ++ )
        {
            rRightHandSideVector[index_up + jdim] += rVariables.tau1 * (rVolumeForce[jdim] - rVariables.PressureGradient[jdim] + rVariables.DynamicRHS[jdim]) * Testf1(i) * rIntegrationWeight;
            rRightHandSideVector[index_up + jdim] += rVariables.tau1 * (rVolumeForce[jdim] - rVariables.PressureGradient[jdim] + rVariables.DynamicRHS[jdim]) * Testf2(indexi)  * rIntegrationWeight;

            rRightHandSideVector[index_up + jdim] += rVariables.tau2  * (-(1.0 - 1.0 / rVariables.detFT)) * rVariables.DN_DX(i,jdim) *rIntegrationWeight;
            rRightHandSideVector[index_up + jdim] += rVariables.tau2  * ((rVariables.PressureGP/rVariables.BulkModulus)) * rVariables.DN_DX(i,jdim) *rIntegrationWeight;
            indexi++;
        }
    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int index_p = dimension;
    Vector aux_vector;

    aux_vector = - rVariables.PressureGradient + rVolumeForce + rVariables.DynamicRHS;


    Vector Stab1 = prod(rVariables.DN_DX,aux_vector);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rRightHandSideVector[index_p] -= rVariables.tau1  *  Stab1(i) * rIntegrationWeight;
        //for ( unsigned int idime = 0; idime < dimension; idime++ ) {
        //    rRightHandSideVector[index_p] += rVariables.tau1  *  rVariables.DN_DX(i,idime)*(rVariables.PressureGradient[idime] + rVolumeForce[idime] - mMP.density* rVariables.DynamicRHS[idime]) * rIntegrationWeight;
        //}

        rRightHandSideVector[index_p] += rVariables.tau2  * ((rVariables.PressureGP/rVariables.BulkModulus)-(1.0 - 1.0 / rVariables.detFT)) * r_N(0, i) * (1/rVariables.BulkModulus) * rIntegrationWeight;

        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddLHS(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    //rVariables.detF0   *= rVariables.detF;
    //double determinant_F = rVariables.detF;
    //rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    // Operation performed: add Km to the rLefsHandSideMatrix
    CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kg to the rLefsHandSideMatrix
    if (!rCurrentProcessInfo.Has(IGNORE_GEOMETRIC_STIFFNESS))
    {
        CalculateAndAddKuugUP(rLeftHandSideMatrix, rVariables, rIntegrationWeight);
    }

    // Operation performed: add Kup to the rLefsHandSideMatrix
    CalculateAndAddKup( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operations performed: add Kpu to the rLefsHandSideMatrix
    CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kpp to the rLefsHandSideMatrix
    CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kuu stabilization to the rLefsHandSideMatrix
    CalculateAndAddKuuStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kup stabilization to the rLefsHandSideMatrix
    CalculateAndAddKupStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operations performed: add Kpu stabilization to the rLefsHandSideMatrix
    CalculateAndAddKpuStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kpp stabilization to the rLefsHandSideMatrix
    CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //rVariables.detF     = determinant_F;
    //rVariables.detF0   /= rVariables.detF;

}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight
                                             )
{
    KRATOS_TRY

    Matrix Kuum = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) );

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexi = 0;
    unsigned int indexj = 0;

    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuum(indexi,indexj);
                    indexj++;
                }
            }
            indexi++;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddKuugUP(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const int size = number_of_nodes * dimension;

    Matrix stress_tensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
    Matrix reduced_Kg = prod( rVariables.DN_DX, rIntegrationWeight * Matrix( prod( stress_tensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized
    Matrix Kuug = ZeroMatrix(size, size);
    MathUtils<double>::ExpandAndAddReducedMatrix( Kuug, reduced_Kg, dimension );

    // Assemble components considering added DOF matrix system
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
                    rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuug(indexi,indexj);
                    indexj++;
                }
            }
            indexi++;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index_p  = dimension;
        unsigned int index_up = dimension * i + i;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(index_up + k, index_p) += rVariables.DN_DX(i, k) * r_N(0, j) * rIntegrationWeight;
            }
            index_p += (dimension + 1);
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    // Assemble components considering added DOF matrix system
    unsigned int index_p = dimension;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            unsigned int index_up = dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(index_p, index_up + k) += r_N(0, i) * rVariables.DN_DX(j, k) * rIntegrationWeight;
            }
        }
        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

        // ATTENTION: class not used in the current implementation!!

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    //double delta_coefficient = rVariables.detF0 - 1;

    unsigned int indexpi = dimension;

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int indexpj = dimension;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            rLeftHandSideMatrix(indexpi, indexpj) -= (1.0 / rVariables.BulkModulus) * r_N(0, i) * r_N(0, j) * rIntegrationWeight; // / (delta_coefficient * (rVariables.detF0 / rVariables.detF))

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }


    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddKuuStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)

{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    Vector Kuustab1 = prod(rVariables.DN_DX, rVariables.PressureGradient);
    Vector Kuustab2 = prod( Matrix(trans(prod(rVariables.TensorIdentityMatrix,rVariables.B))),rVariables.PressureGradientVoigt);
    Vector testf1   = Kuustab1;
    Vector testf2   = prod( Matrix(prod(trans(rVariables.B),rVariables.TensorIdentityMatrix)),rVariables.PressureGradientVoigt);


    unsigned int indexi = 0;
    unsigned int indexj = 0;

    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(indexi+i,indexj+j)+= rVariables.tau1 *  rVariables.DynamicCoefficient * r_N( 0 , i ) * rIntegrationWeight * testf1(j);
                    rLeftHandSideMatrix(indexi+i,indexj+j)+= rVariables.tau1 *  rVariables.DynamicCoefficient * r_N( 0 , i ) * rIntegrationWeight * testf2(indexj);
                    rLeftHandSideMatrix(indexi+i,indexj+j)-= rVariables.tau1 * Kuustab1(i) * rIntegrationWeight * testf1(j);
                    rLeftHandSideMatrix(indexi+i,indexj+j)-= rVariables.tau1 * Kuustab2(indexi) * rIntegrationWeight * testf1(j);
                    rLeftHandSideMatrix(indexi+i,indexj+j)-= rVariables.tau1 * Kuustab1(i) * rIntegrationWeight * testf2(indexj);
                    rLeftHandSideMatrix(indexi+i,indexj+j)-= rVariables.tau1 * Kuustab2(indexi) * rIntegrationWeight * testf2(indexj);
                    rLeftHandSideMatrix(indexi+i,indexj+j)+= rVariables.tau2 * rVariables.DN_DX(i,idim) * rIntegrationWeight * rVariables.DN_DX(j,jdim);
                    indexj++;
                }
            }
            indexi++;
        }
    }

    KRATOS_CATCH( "" )
}



//***********************************************************************************
//***********************************************************************************
void UpdatedLagrangianUPVMS::CalculateAndAddKupStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    Vector Stab1 = prod(rVariables.DN_DX, rVariables.PressureGradient);
    Vector Stab2 = prod(Matrix(trans(prod(rVariables.TensorIdentityMatrix,rVariables.B))),rVariables.PressureGradientVoigt);

    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index_p  = dimension;
        unsigned int index_up = dimension * i + i;

        unsigned int indexj = 0;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(index_up + k, index_p) += rVariables.tau1 *  rVariables.DynamicCoefficient * r_N(0 , j) * rVariables.DN_DX(i, k) * rIntegrationWeight;
                rLeftHandSideMatrix(index_up + k, index_p) -= rVariables.tau1 * Stab1(j) * rVariables.DN_DX(i, k) * rIntegrationWeight;
                rLeftHandSideMatrix(index_up + k, index_p) -= rVariables.tau1 * Stab2(indexj) * rVariables.DN_DX(i, k) * rIntegrationWeight;

                rLeftHandSideMatrix(index_up + k, index_p) += rVariables.tau2 * (1/rVariables.BulkModulus)* rVariables.DN_DX(i, k) * r_N(0, j) * rIntegrationWeight;
                indexj++;
            }
            index_p += (dimension + 1);
        }
    }


    KRATOS_CATCH( "" )
}



//***********************************************************************************
//***********************************************************************************
void UpdatedLagrangianUPVMS::CalculateAndAddKpuStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    Vector Testf1 = prod(rVariables.DN_DX, rVariables.PressureGradient);
    Vector Testf2   = prod(Matrix(prod(trans(rVariables.B),rVariables.TensorIdentityMatrix)),rVariables.PressureGradientVoigt);

    // Assemble components considering added DOF matrix system
    unsigned int index_p = dimension;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexj = 0;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            unsigned int index_up = dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(index_p, index_up + k) -= rVariables.tau1 * rVariables.DN_DX(i, k) * Testf1(j) * rIntegrationWeight;
                rLeftHandSideMatrix(index_p, index_up + k) -= rVariables.tau1 * rVariables.DN_DX(i, k) * Testf2(indexj) * rIntegrationWeight;

                indexj++;
                rLeftHandSideMatrix(index_p, index_up + k) -= rVariables.tau2 * (1.0 / rVariables.BulkModulus) * r_N(0, i) * rVariables.DN_DX(j, k) * rIntegrationWeight;
            }
        }
        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void UpdatedLagrangianUPVMS::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY;
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    unsigned int indexpi = dimension;
    Matrix Stab1 = prod(rVariables.DN_DX, trans(rVariables.DN_DX));


    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int indexpj = dimension;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            rLeftHandSideMatrix(indexpi, indexpj) -= rVariables.tau1 * Stab1(i,j)  * rIntegrationWeight;
            rLeftHandSideMatrix(indexpi, indexpj) -= rVariables.tau2 * (1.0 / pow(rVariables.BulkModulus,2)) *r_N(0, i) * r_N(0, j) * rIntegrationWeight;

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }



    KRATOS_CATCH( "" )
}



//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangianUPVMS::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension       = r_geometry.WorkingSpaceDimension();
    unsigned int element_size          = number_of_nodes * dimension + number_of_nodes;

    if ( rResult.size() != element_size )
        rResult.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension + i;
        rResult[index]     = r_geometry[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if( dimension == 3)
        {
            rResult[index + 2] = r_geometry[i].GetDof( DISPLACEMENT_Z ).EquationId();
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

void UpdatedLagrangianUPVMS::GetDofList( DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo ) const
{
    rElementalDofList.resize( 0 );

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    for ( unsigned int i = 0; i < r_geometry.size(); i++ )
    {
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Y ) );

        if( dimension == 3 )
            rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Z ) );

        rElementalDofList.push_back( r_geometry[i].pGetDof( PRESSURE ));
    }
}


//************************************************************************************
//****************MASS MATRIX*********************************************************

void UpdatedLagrangianUPVMS::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

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

void UpdatedLagrangianUPVMS::GetValuesVector( Vector& values, int Step ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension       = r_geometry.WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
        {
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
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

void UpdatedLagrangianUPVMS::GetFirstDerivativesVector( Vector& values, int Step ) const
{
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
            values[index + 3] = 0;
        }
        else
        {
            values[index + 2] = 0;
        }
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::GetSecondDerivativesVector( Vector& values, int Step ) const
{
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

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::GetHistoricalVariables( GeneralVariables& rVariables )
{
    //Deformation Gradient F ( set to identity )
    UpdatedLagrangian::GetHistoricalVariables(rVariables);
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUPVMS::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    double voigtsize = 3;
    if ( dimension == 3)
        voigtsize = 6;

    UpdatedLagrangian::FinalizeStepVariables( rVariables, rCurrentProcessInfo);

    // Evaluation of the pressure on the material point
    double nodal_mean_stress = 0.0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
        nodal_mean_stress += r_geometry[i].FastGetSolutionStepValue( PRESSURE ) * r_N(0, i);

    // Evaluation of the mean stress on the material point
    double mean_stress = 0.0;
    for (unsigned int i = 0; i < dimension; i++)
        mean_stress += rVariables.StressVector[i];
    mean_stress /= dimension;

    Vector stress_vector = ZeroVector(voigtsize);
    stress_vector = rVariables.StressVector;
    for (unsigned int i = 0; i < dimension; i++)
        stress_vector[i] += (nodal_mean_stress - mean_stress);

    mMP.cauchy_stress_vector = stress_vector;

}

void UpdatedLagrangianUPVMS::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_PRESSURE) {
        rValues[0] = m_mp_pressure;
    }
    else {
        UpdatedLagrangian::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void UpdatedLagrangianUPVMS::SetValuesOnIntegrationPoints(
    const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_PRESSURE) {
        m_mp_pressure = rValues[0];
    }
    else {
        UpdatedLagrangian::SetValuesOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

//************************************************************************************
//************************************************************************************
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int UpdatedLagrangianUPVMS::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
    ? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
    : false;
    KRATOS_ERROR_IF(is_explicit)
    << "Explicit time integration not implemented for Updated Lagrangian UP MPM Element";

    int correct = 0;
    correct = UpdatedLagrangian::Check(rCurrentProcessInfo);

    // Verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetLawFeatures(LawFeatures);

    KRATOS_ERROR_IF(LawFeatures.mOptions.IsNot(ConstitutiveLaw::U_P_LAW)) << "Constitutive law is not compatible with the U-P element type: Large Displacements U_P" << std::endl;

    return correct;

    KRATOS_CATCH( "" );
}

void UpdatedLagrangianUPVMS::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangian )
    rSerializer.save("Pressure",m_mp_pressure);
}

void UpdatedLagrangianUPVMS::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangian )
    rSerializer.load("Pressure",m_mp_pressure);
}

} // Namespace Kratos

