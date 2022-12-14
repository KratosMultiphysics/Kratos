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


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUPVMS::UpdatedLagrangianUPVMS()
    : UpdatedLagrangianUP()
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianUPVMS::UpdatedLagrangianUPVMS( IndexType NewId, GeometryType::Pointer pGeometry )
    : UpdatedLagrangianUP( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUPVMS::UpdatedLagrangianUPVMS( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : UpdatedLagrangianUP( NewId, pGeometry, pProperties )
{
    mFinalizedStep = true;


}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianUPVMS::UpdatedLagrangianUPVMS( UpdatedLagrangianUPVMS const& rOther)
    :UpdatedLagrangianUP(rOther)

{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianUPVMS&  UpdatedLagrangianUPVMS::operator=(UpdatedLagrangianUPVMS const& rOther)
{
    UpdatedLagrangianUP::operator=(rOther);


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

//     if ( stabilization_type == 1)
//      {
//           rVariables.tau1 = 0;
//          rVariables.tau2 = 0;
//      }

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
//         unsigned int index_u  = dimension * i;

        for ( unsigned int jdim = 0; jdim < dimension; jdim ++ )
        {
            rRightHandSideVector[index_up + jdim] -= rVariables.tau1 * (-rVolumeForce[jdim] - rVariables.PressureGradient[jdim] + rVariables.DynamicRHS[jdim]) * Testf1(i) * rIntegrationWeight;
            rRightHandSideVector[index_up + jdim] -= rVariables.tau1 * (-rVolumeForce[jdim] - rVariables.PressureGradient[jdim] + rVariables.DynamicRHS[jdim]) * Testf2(indexi)  * rIntegrationWeight;

            rRightHandSideVector[index_up + jdim] += rVariables.tau2  * ((rVariables.PressureGP/rVariables.BulkModulus)-(1.0 - 1.0 / rVariables.detFT) ) * rVariables.DN_DX(i,jdim) *rIntegrationWeight; //(rVariables.PressureGP/rVariables.BulkModulus)+
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

//     aux_vector = - rVariables.PressureGradient + rVolumeForce + rVariables.DynamicRHS;
//     Vector Stab1 = prod(rVariables.DN_DX,aux_vector);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
//         rRightHandSideVector[index_p] -= rVariables.tau1  *  Stab1(i) * rIntegrationWeight;
        for ( unsigned int idime = 0; idime < dimension; idime++ ) {
           rRightHandSideVector[index_p] -= rVariables.tau1  *  rVariables.DN_DX(i,idime)*(-rVariables.PressureGradient[idime] - rVolumeForce[idime] + rVariables.DynamicRHS[idime]) * rIntegrationWeight;
        }

        rRightHandSideVector[index_p] -= rVariables.tau2  * ((rVariables.PressureGP/rVariables.BulkModulus)-(1.0 - 1.0 / rVariables.detFT)) * r_N(0, i) * (1/rVariables.BulkModulus) * rIntegrationWeight;

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

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    unsigned int indexpi = dimension;

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int indexpj = dimension;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            rLeftHandSideMatrix(indexpi, indexpj) -= (1.0 / rVariables.BulkModulus) * r_N(0, i) * r_N(0, j) * rIntegrationWeight;

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
            indexj = 0;
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

                rLeftHandSideMatrix(index_up + k, index_p) -= rVariables.tau2 * (1/rVariables.BulkModulus)* rVariables.DN_DX(i, k) * r_N(0, j) * rIntegrationWeight;
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
            rLeftHandSideMatrix(indexpi, indexpj) += rVariables.tau2 * (1.0 / pow(rVariables.BulkModulus,2)) *r_N(0, i) * r_N(0, j) * rIntegrationWeight;

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }

    KRATOS_CATCH( "" )
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

