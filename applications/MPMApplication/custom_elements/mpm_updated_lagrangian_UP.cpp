//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta and Alessandro Contri
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/mpm_updated_lagrangian_UP.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "mpm_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{

///**
//* Flags related to the element computation
//*/
//KRATOS_CREATE_LOCAL_FLAG( MPMUpdatedLagrangian, COMPUTE_RHS_VECTOR,                 0 );
//KRATOS_CREATE_LOCAL_FLAG( MPMUpdatedLagrangian, COMPUTE_LHS_MATRIX,                 1 );
//KRATOS_CREATE_LOCAL_FLAG( MPMUpdatedLagrangian, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
//KRATOS_CREATE_LOCAL_FLAG( MPMUpdatedLagrangian, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianUP::MPMUpdatedLagrangianUP()
    : MPMUpdatedLagrangian()
    , m_mp_pressure(1.0)
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianUP::MPMUpdatedLagrangianUP( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMUpdatedLagrangian( NewId, pGeometry )
    , m_mp_pressure(1.0)
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianUP::MPMUpdatedLagrangianUP( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : MPMUpdatedLagrangian( NewId, pGeometry, pProperties )
    , m_mp_pressure(1.0)
{
    mFinalizedStep = true;


}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

MPMUpdatedLagrangianUP::MPMUpdatedLagrangianUP( MPMUpdatedLagrangianUP const& rOther)
    :MPMUpdatedLagrangian(rOther)
    , m_mp_pressure(rOther.m_mp_pressure)
     //,mDeformationGradientF0(rOther.mDeformationGradientF0)
     //,mDeterminantF0(rOther.mDeterminantF0)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

MPMUpdatedLagrangianUP&  MPMUpdatedLagrangianUP::operator=(MPMUpdatedLagrangianUP const& rOther)
{
    MPMUpdatedLagrangian::operator=(rOther);

    m_mp_pressure = rOther.m_mp_pressure;

    return *this;
}
//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer MPMUpdatedLagrangianUP::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new MPMUpdatedLagrangianUP( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer MPMUpdatedLagrangianUP::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< MPMUpdatedLagrangianUP >(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer MPMUpdatedLagrangianUP::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    MPMUpdatedLagrangianUP NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.m_mp_pressure = m_mp_pressure;

    NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();

    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mDeterminantF0 = mDeterminantF0;

    return Element::Pointer( new MPMUpdatedLagrangianUP(NewElement) );
}
//*******************************DESTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianUP::~MPMUpdatedLagrangianUP()
{
}


//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUP::Initialize(const ProcessInfo& rCurrentProcessInfo)
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

void MPMUpdatedLagrangianUP::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    MPMUpdatedLagrangian::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************
/**
 * The position of the Gauss points/Material points is updated
 */

void MPMUpdatedLagrangianUP::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
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


void MPMUpdatedLagrangianUP::CalculateKinematics(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
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

void MPMUpdatedLagrangianUP::CalculateDeformationMatrix(Matrix& rB,
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

void MPMUpdatedLagrangianUP::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
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

void MPMUpdatedLagrangianUP::CalculateAndAddRHS(
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    Vector& rVolumeForce,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    rVariables.detF0   *= rVariables.detF;
    double determinant_F = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    // Operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
    CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // Operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
    CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // Operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
    CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    if (rCurrentProcessInfo.GetValue(STABILIZATION_TYPE)==1)
        // Operation performed: rRightHandSideVector -= Stabilized Pressure Forces
        CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, rIntegrationWeight);

    rVariables.detF     = determinant_F;
    rVariables.detF0   /= rVariables.detF;

}
//************************************************************************************
//*********************Calculate the contribution of external force*******************

void MPMUpdatedLagrangianUP::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
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

void MPMUpdatedLagrangianUP::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
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

//******************************************************************************************************************
//******************************************************************************************************************
double& MPMUpdatedLagrangianUP::CalculatePUCoefficient(double& rCoefficient, GeneralVariables & rVariables)
{
    KRATOS_TRY

    // TODO: Check what is the meaning of this function
    rCoefficient = rVariables.detF0 - 1;

    return rCoefficient;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

double& MPMUpdatedLagrangianUP::CalculatePUDeltaCoefficient(double &rDeltaCoefficient, GeneralVariables & rVariables)
{

    KRATOS_TRY

    // TODO: Check what is the meaning of this function
    rDeltaCoefficient = 1.0;

    return rDeltaCoefficient;

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUP::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int index_p = dimension;
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    const double& young_modulus = GetProperties()[YOUNG_MODULUS];
    const double& poisson_ratio    = GetProperties()[POISSON_RATIO];
    double bulk_modulus  = young_modulus/(3.0*(1.0-2.0*poisson_ratio));

    // Check if Bulk Modulus is not NaN
    if (bulk_modulus != bulk_modulus)
        bulk_modulus = 1.e16;

    double delta_coefficient = 0;
    delta_coefficient = this->CalculatePUDeltaCoefficient( delta_coefficient, rVariables );

    double coefficient = 0;
    coefficient = this->CalculatePUCoefficient( coefficient, rVariables );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            const double& pressure = r_geometry[j].FastGetSolutionStepValue(PRESSURE);

            rRightHandSideVector[index_p] += (1.0/(delta_coefficient * bulk_modulus)) * r_N(0, i) * r_N(0, j) * pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //2D-3D
        }

        rRightHandSideVector[index_p] -=  coefficient/delta_coefficient * r_N(0, i) * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUP::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int index_p = dimension;

    // Stabilization alpha parameters
    double alpha_stabilization  = 1.0;
    if (GetProperties().Has(YOUNG_MODULUS) && GetProperties().Has(POISSON_RATIO))
    {
        const double& young_modulus = GetProperties()[YOUNG_MODULUS];
        const double& poisson_ratio = GetProperties()[POISSON_RATIO];
        const double& shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio));

        double factor_value = 8.0; //JMR deffault value
        if (dimension == 3)
            factor_value = 10.0; //JMC deffault value

    	alpha_stabilization = alpha_stabilization * factor_value / shear_modulus;
    }
    else
    {
        KRATOS_ERROR << "No pressure stabilization existing if YOUNG_MODULUS and POISSON_RATIO or DYNAMIC_VISCOSITY is specified" << std::endl;
    }


    double consistent = 1.0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            const double& pressure = r_geometry[j].FastGetSolutionStepValue(PRESSURE);

            if( dimension == 2 )
            {
                consistent = (-1) * alpha_stabilization / 36.0;
                if (i == j)
                    consistent = 2 * alpha_stabilization / 36.0;

                rRightHandSideVector[index_p] += consistent * pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D
            }
            else
            {
                consistent = (-1) * alpha_stabilization / 80.0;
                if (i == j)
                    consistent = 3 * alpha_stabilization / 80.0;

                rRightHandSideVector[index_p] += consistent * pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //3D
            }
        }
        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUP::CalculateAndAddLHS(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    rVariables.detF0   *= rVariables.detF;
    double determinant_F = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    // Operation performed: add Km to the rLefsHandSideMatrix
    CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kg to the rLefsHandSideMatrix
    if (!rCurrentProcessInfo.Has(IGNORE_GEOMETRIC_STIFFNESS))
    {
        CalculateAndAddKuugUP(rLeftHandSideMatrix, rVariables, rIntegrationWeight);
    }

    // Operation performed: add Kup to the rLefsHandSideMatrix
    CalculateAndAddKup( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kpu to the rLefsHandSideMatrix
    CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kpp to the rLefsHandSideMatrix
    CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    if (rCurrentProcessInfo.GetValue(STABILIZATION_TYPE)==1)
        // Operation performed: add Kpp_Stab to the rLefsHandSideMatrix
        CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    rVariables.detF     = determinant_F;
    rVariables.detF0   /= rVariables.detF;

}
//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUP::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
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

void MPMUpdatedLagrangianUP::CalculateAndAddKuugUP(MatrixType& rLeftHandSideMatrix,
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

void MPMUpdatedLagrangianUP::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
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
                rLeftHandSideMatrix(index_up+k,index_p) +=  rVariables.DN_DX ( i, k ) * r_N(0, j) * rIntegrationWeight * rVariables.detF;
            }
            index_p += (dimension + 1);
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUP::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
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
                rLeftHandSideMatrix(index_p,index_up+k) += r_N(0, i) * rVariables.DN_DX ( j, k ) * rIntegrationWeight * rVariables.detF;
            }
        }
        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUP::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    const double& young_modulus = GetProperties()[YOUNG_MODULUS];
    const double& poisson_ratio = GetProperties()[POISSON_RATIO];
    double bulk_modulus = young_modulus / (3.0 * (1.0 - 2.0 * poisson_ratio));

    // Check if Bulk Modulus is not NaN
    if (bulk_modulus != bulk_modulus)
        bulk_modulus = 1.e16;

    unsigned int indexpi = dimension;

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int indexpj = dimension;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            rLeftHandSideMatrix(indexpi,indexpj)  -= ((1.0)/(bulk_modulus)) * r_N(0, i) * r_N(0, j) * rIntegrationWeight /((rVariables.detF0/rVariables.detF));

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************
// I changed the constant matrix in the stabilized term:
// as in MPM the position of the integration points does not coincide with the
// position of the Gauss points the first matrix is substitute with the product of the
// shape function values of each integration point
//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUP::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables & rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    // Stabilization alpha parameters
    double alpha_stabilization = 1.0;
    if (GetProperties().Has(YOUNG_MODULUS) && GetProperties().Has(POISSON_RATIO))
    {
        const double& young_modulus = GetProperties()[YOUNG_MODULUS];
        const double& poisson_ratio = GetProperties()[POISSON_RATIO];
        const double& shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio));

        double factor_value = 8.0; //JMR deffault value
        if (dimension == 3)
            factor_value = 10.0; //JMC deffault value

        alpha_stabilization = alpha_stabilization * factor_value / shear_modulus;
    }
        else
    {
        KRATOS_ERROR << "No pressure stabilization existing if YOUNG_MODULUS and POISSON_RATIO or DYNAMIC_VISCOSITY is specified" << std::endl;
    }

    double consistent = 1.0;
    unsigned int indexpi = dimension;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            if( dimension == 2 )  //consistent 2D
            {
                consistent = (-1) * alpha_stabilization / 36.0;
                if (indexpi == indexpj)
                    consistent = 2 * alpha_stabilization / 36.0;

                rLeftHandSideMatrix(indexpi, indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //2D
            }
            else
            {
                consistent = (-1) * alpha_stabilization / 80.0;
                if (indexpi == indexpj)
                    consistent = 3 * alpha_stabilization / 80.0;

                rLeftHandSideMatrix(indexpi, indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //3D
            }
            indexpj += (dimension + 1);
        }
        indexpi += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}


//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& MPMUpdatedLagrangianUP::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUP::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo ) const
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

void MPMUpdatedLagrangianUP::GetDofList( DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo ) const
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

void MPMUpdatedLagrangianUP::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
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

void MPMUpdatedLagrangianUP::GetValuesVector( Vector& values, int Step ) const
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

void MPMUpdatedLagrangianUP::GetFirstDerivativesVector( Vector& values, int Step ) const
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

void MPMUpdatedLagrangianUP::GetSecondDerivativesVector( Vector& values, int Step ) const
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

void MPMUpdatedLagrangianUP::GetHistoricalVariables( GeneralVariables& rVariables )
{
    //Deformation Gradient F ( set to identity )
    MPMUpdatedLagrangian::GetHistoricalVariables(rVariables);
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUP::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    double voigtsize = 3;
    if ( dimension == 3)
        voigtsize = 6;

    MPMUpdatedLagrangian::FinalizeStepVariables( rVariables, rCurrentProcessInfo);

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

void MPMUpdatedLagrangianUP::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_PRESSURE) {
        rValues[0] = m_mp_pressure;
    }
    else {
        MPMUpdatedLagrangian::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMUpdatedLagrangianUP::SetValuesOnIntegrationPoints(
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
        MPMUpdatedLagrangian::SetValuesOnIntegrationPoints(
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
int MPMUpdatedLagrangianUP::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
    ? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
    : false;
    KRATOS_ERROR_IF(is_explicit)
    << "Explicit time integration not implemented for Updated Lagrangian UP MPM Element";

    int correct = 0;
    correct = MPMUpdatedLagrangian::Check(rCurrentProcessInfo);

    // Verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetLawFeatures(LawFeatures);

    KRATOS_ERROR_IF(LawFeatures.mOptions.IsNot(ConstitutiveLaw::U_P_LAW)) << "Constitutive law is not compatible with the U-P element type: Large Displacements U_P" << std::endl;

    return correct;

    KRATOS_CATCH( "" );
}

void MPMUpdatedLagrangianUP::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMUpdatedLagrangian )
    rSerializer.save("Pressure",m_mp_pressure);
}

void MPMUpdatedLagrangianUP::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMUpdatedLagrangian )
    rSerializer.load("Pressure",m_mp_pressure);
}

} // Namespace Kratos

