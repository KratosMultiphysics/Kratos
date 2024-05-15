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
#include "custom_elements/mpm_updated_lagrangian.hpp"
#include "custom_elements/mpm_updated_lagrangian_UP_VMS.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "mpm_application_variables.h"
#include "includes/checks.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/material_point_generator_utility.h"

namespace Kratos
{


    const unsigned int MPMUpdatedLagrangianUPVMS::msIndexVoigt3D6C [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };
    const unsigned int MPMUpdatedLagrangianUPVMS::msIndexVoigt2D4C [4][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1} };
    const unsigned int MPMUpdatedLagrangianUPVMS::msIndexVoigt2D3C [3][2] = { {0, 0}, {1, 1}, {0, 1} };


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianUPVMS::MPMUpdatedLagrangianUPVMS()
    : MPMUpdatedLagrangianUP()
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianUPVMS::MPMUpdatedLagrangianUPVMS( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMUpdatedLagrangianUP( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianUPVMS::MPMUpdatedLagrangianUPVMS( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : MPMUpdatedLagrangianUP( NewId, pGeometry, pProperties )
{
    mFinalizedStep = true;


}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

MPMUpdatedLagrangianUPVMS::MPMUpdatedLagrangianUPVMS( MPMUpdatedLagrangianUPVMS const& rOther)
    : MPMUpdatedLagrangianUP(rOther)

{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

MPMUpdatedLagrangianUPVMS&  MPMUpdatedLagrangianUPVMS::operator=(MPMUpdatedLagrangianUPVMS const& rOther)
{
    MPMUpdatedLagrangianUP::operator=(rOther);


    return *this;
}
//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer MPMUpdatedLagrangianUPVMS::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new MPMUpdatedLagrangianUPVMS( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer MPMUpdatedLagrangianUPVMS::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< MPMUpdatedLagrangianUPVMS >(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer MPMUpdatedLagrangianUPVMS::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    MPMUpdatedLagrangianUPVMS NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.m_mp_pressure = m_mp_pressure;

    NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();

    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mDeterminantF0 = mDeterminantF0;

    return Element::Pointer( new MPMUpdatedLagrangianUPVMS(NewElement) );
}
//*******************************DESTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianUPVMS::~MPMUpdatedLagrangianUPVMS()
{
}


//************************************************************************************
//************************************************************************************

// template <class TElementData>
void MPMUpdatedLagrangianUPVMS::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}


void MPMUpdatedLagrangianUPVMS::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput, const ProcessInfo& rCurrentProcessInfo) {
    // Lumped projection terms
    if (rVariable == RESPROJ_DISPL) {
        this->CalculateProjections(rCurrentProcessInfo);
    } else {
        BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

// template <class TElementData>
void MPMUpdatedLagrangianUPVMS::Calculate(
    const Variable<Vector>& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

// template <class TElementData>
void MPMUpdatedLagrangianUPVMS::Calculate(
    const Variable<Matrix>& rVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}



void MPMUpdatedLagrangianUPVMS::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    MPMUpdatedLagrangian::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);
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
    rVariables.DiscreteAcceleration=ZeroVector(dimension);

    // Set Identity matrices
    rVariables.Identity = IdentityMatrix(dimension);
    rVariables.TensorIdentityMatrix = ZeroMatrix(voigt_dimension,voigt_dimension);

    // Set variables for OSGS stabilization
    rVariables.ResProjDisplGP = ZeroVector(dimension);
    rVariables.ResProjPressGP = 0;

    // Set Body forces
    rVariables.BodyForceMP = ZeroVector(dimension);

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************


void MPMUpdatedLagrangianUPVMS::CalculateElementalSystem(
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
        //#BODYFORCE
        Vector volume_force = (mMP.volume_acceleration * mMP.mass ) + (Variables.BodyForceMP * mMP.mass);
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

void MPMUpdatedLagrangianUPVMS::SetSpecificVariables(GeneralVariables& rVariables,const ProcessInfo& rCurrentProcessInfo)
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

    // Set the identity matrix tensor
    CalculateTensorIdentityMatrix(rVariables,rVariables.TensorIdentityMatrix);

    // Set if the model is dynamic
    const bool is_dynamic = rCurrentProcessInfo.Has(IS_DYNAMIC)
        ? rCurrentProcessInfo.GetValue(IS_DYNAMIC)
        : false;
    const int current_step = rCurrentProcessInfo.GetValue(STEP);
    rVariables.DiscreteAcceleration =ZeroVector(3);
    if (rCurrentProcessInfo.GetValue(STABILIZATION_TYPE) ==2) ComputeDynamicTerms(rVariables,rCurrentProcessInfo);


    // Compute Residual Projection in integration points
    if (rCurrentProcessInfo.GetValue(STABILIZATION_TYPE) == 3)
    {
        array_1d<double, 3 > nodal_resprojdispl = ZeroVector(3);
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            rVariables.ResProjPressGP += r_N(0,j) * r_geometry[j].FastGetSolutionStepValue(RESPROJ_PRESS);

            nodal_resprojdispl = r_geometry[j].FastGetSolutionStepValue(RESPROJ_DISPL,0);
            for (unsigned int k = 0; k < dimension; k++)
            {
               rVariables.ResProjDisplGP[k] += r_N(0, j) * nodal_resprojdispl[k];
            }
        }

    }

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

void MPMUpdatedLagrangianUPVMS::ComputeElementSize(double& ElemSize){

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

void MPMUpdatedLagrangianUPVMS::CalculateTaus(const int& stabilization_type,
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

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************


void MPMUpdatedLagrangianUPVMS::CalculateTensorIdentityMatrix (GeneralVariables& rVariables, Matrix& rTensorIdentityMatrix)
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

double& MPMUpdatedLagrangianUPVMS::TensorIdentityComponent (double& rCabcd, GeneralVariables& rVariables,
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

void MPMUpdatedLagrangianUPVMS::ConvertPressureGradientInVoigt(Vector& PressureGradient,Vector& PressureGradientVoigt)
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

void MPMUpdatedLagrangianUPVMS::ComputeDynamicTerms(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
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
    
    array_1d<double, 3 > previous_displacement = ZeroVector(3);
    array_1d<double, 3 > previous_velocity = ZeroVector(3);
    array_1d<double, 3 > previous_displacement_it = ZeroVector(3);
    array_1d<double, 3 > previous_acceleration = ZeroVector(3);


    for (unsigned int j=0; j<number_of_nodes; j++)
    {
        // These are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
        previous_acceleration = r_geometry[j].FastGetSolutionStepValue(ACCELERATION,1);
        previous_velocity = r_geometry[j].FastGetSolutionStepValue(VELOCITY,1);
        previous_displacement_it = r_geometry[j].FastGetSolutionStepValue(DISPLACEMENT,0);
        previous_displacement = r_geometry[j].FastGetSolutionStepValue(DISPLACEMENT,1);

        for (unsigned int k = 0; k < dimension; k++)
        {
            aux_MP_velocity[k] += r_N(0, j) * previous_velocity[k];
            aux_MP_acceleration[k] += r_N(0, j) * previous_acceleration[k];
            aux_MP_displacement[k] += r_N(0, j) * previous_displacement[k];
            current_MP_displacement[k] += r_N(0, j) * previous_displacement_it[k];
        }
    }

    const double density_mp=(GetProperties()[DENSITY]) / rVariables.detFT;
    rVariables.DynamicCoefficient = density_mp / (beta * delta_time * delta_time);
    const double coeff1 = 1 / (delta_time * 0.25);
    const double coeff2 = (0.5 - beta) / beta;

    rVariables.DynamicRHS =ZeroVector(3);
    rVariables.DiscreteAcceleration =ZeroVector(3);

    for (unsigned int idime = 0; idime < dimension; idime++)
    {
           rVariables.DynamicRHS[idime] -= rVariables.DynamicCoefficient * current_MP_displacement[idime];
           rVariables.DynamicRHS[idime] += rVariables.DynamicCoefficient * aux_MP_displacement[idime];
           rVariables.DiscreteAcceleration[idime]=rVariables.DynamicCoefficient * aux_MP_displacement[idime];
           rVariables.DynamicRHS[idime] += coeff1 * density_mp * aux_MP_velocity[idime] + coeff2 * density_mp * aux_MP_acceleration[idime];

    }

    //std::cout<<"dynamicRHS "<<rVariables.DynamicRHS[1]<<"\n";
    //std::cout<<"dynamicCoefficinet "<<rVariables.DynamicCoefficient <<"\n";
   

    //rVariables.DynamicRHS=mMP.density*rVariables.DynamicRHS;
    //rVariables.DynamicCoefficient=mMP.density*rVariables.DynamicRHS;

    //rVariables.DynamicRHS =ZeroVector(3);
   // rVariables.DynamicCoefficient=0;

    

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateAndAddRHS(
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
//************************************************************************************




void MPMUpdatedLagrangianUPVMS::CalculateAndAddStabilizedDisplacement(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    Vector Testf1 = prod(rVariables.DN_DX, rVariables.PressureGradient);
    Vector Testf2   = prod(prod(trans(rVariables.B),rVariables.TensorIdentityMatrix),rVariables.PressureGradientVoigt);
    double VolumetricStrainFunction = this->CalculateVolumetricStrainFunction( VolumetricStrainFunction, rVariables );

    unsigned int indexi  = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index_up = dimension * i + i;
        for ( unsigned int jdim = 0; jdim < dimension; jdim ++ )
        {
            rRightHandSideVector[index_up + jdim] -= rVariables.tau1 * (-rVolumeForce[jdim] - rVariables.PressureGradient[jdim] - rVariables.DynamicRHS[jdim] + rVariables.ResProjDisplGP[jdim]) * Testf1(i) * rIntegrationWeight;
            rRightHandSideVector[index_up + jdim] -= rVariables.tau1 * (-rVolumeForce[jdim] - rVariables.PressureGradient[jdim] - rVariables.DynamicRHS[jdim] + rVariables.ResProjDisplGP[jdim]) * Testf2(indexi)  * rIntegrationWeight;

            rRightHandSideVector[index_up + jdim] += rVariables.tau2  * ((rVariables.PressureGP/rVariables.BulkModulus)- VolumetricStrainFunction  + rVariables.ResProjPressGP) * rVariables.DN_DX(i,jdim) *rIntegrationWeight;
            indexi++; //
        }
    }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
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
    double VolumetricStrainFunction = this->CalculateVolumetricStrainFunction( VolumetricStrainFunction, rVariables );
    double functionJ = this->CalculateFunctionFromLinearizarionOfVolumetricStrain( functionJ, rVariables );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
         for ( unsigned int idime = 0; idime < dimension; idime++ ) {
            rRightHandSideVector[index_p] -= rVariables.tau1  * functionJ * rVariables.DN_DX(i,idime)*(-rVariables.PressureGradient[idime] - rVolumeForce[idime] - rVariables.DynamicRHS[idime] + rVariables.ResProjDisplGP[idime]) * rIntegrationWeight;
         }

        rRightHandSideVector[index_p] -= rVariables.tau2  * ((rVariables.PressureGP/rVariables.BulkModulus) - VolumetricStrainFunction + rVariables.ResProjPressGP) * r_N(0, i) * (1/rVariables.BulkModulus) * rIntegrationWeight;

        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateAndAddLHS(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{

    // Operation performed: add Km to the rLefsHandSideMatrix
    MPMUpdatedLagrangianUP::CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kg to the rLefsHandSideMatrix
    if (!rCurrentProcessInfo.Has(IGNORE_GEOMETRIC_STIFFNESS))
    {
        CalculateAndAddKuugUP(rLeftHandSideMatrix, rVariables, rIntegrationWeight);
    }

    // Operation performed: add Kup to the rLefsHandSideMatrix
    MPMUpdatedLagrangianUP::CalculateAndAddKup( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operations performed: add Kpu to the rLefsHandSideMatrix
    MPMUpdatedLagrangianUP::CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kpp to the rLefsHandSideMatrix
    MPMUpdatedLagrangianUP::CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

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

void MPMUpdatedLagrangianUPVMS::CalculateAndAddKuuStab (MatrixType& rLeftHandSideMatrix,
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
    double functionJ = this->CalculateFunctionFromLinearizarionOfVolumetricStrain( functionJ, rVariables );


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
                    rLeftHandSideMatrix(indexi+i,indexj+j)+= rVariables.tau2 * functionJ *rVariables.DN_DX(i,idim) * rIntegrationWeight * rVariables.DN_DX(j,jdim);
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
void MPMUpdatedLagrangianUPVMS::CalculateAndAddKupStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    Vector Stab1 = prod(rVariables.DN_DX, rVariables.PressureGradient);
    Vector Stab2 = prod(Matrix(trans(prod(rVariables.TensorIdentityMatrix,rVariables.B))),rVariables.PressureGradientVoigt);
    double functionJ = this->CalculateFunctionFromLinearizarionOfVolumetricStrain( functionJ, rVariables );

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
                rLeftHandSideMatrix(index_up + k, index_p) -= rVariables.tau1 * Stab1(j) * functionJ * rVariables.DN_DX(i, k) * rIntegrationWeight;
                rLeftHandSideMatrix(index_up + k, index_p) -= rVariables.tau1 * Stab2(indexj) * functionJ * rVariables.DN_DX(i, k) * rIntegrationWeight;
                rLeftHandSideMatrix(index_up + k, index_p) -= rVariables.tau2 * functionJ * (1 / rVariables.BulkModulus)* rVariables.DN_DX(i, k) * r_N(0, j) * rIntegrationWeight;
                indexj++;
            }
            index_p += (dimension + 1);
        }
    }


    KRATOS_CATCH( "" )
}



//***********************************************************************************
//***********************************************************************************
void MPMUpdatedLagrangianUPVMS::CalculateAndAddKpuStab (MatrixType& rLeftHandSideMatrix,
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

void MPMUpdatedLagrangianUPVMS::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY;
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    unsigned int indexpi = dimension;
    Matrix Stab1 = prod(rVariables.DN_DX, trans(rVariables.DN_DX));
    double functionJ = this->CalculateFunctionFromLinearizarionOfVolumetricStrain( functionJ, rVariables );

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int indexpj = dimension;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            rLeftHandSideMatrix(indexpi, indexpj) -= rVariables.tau1 * functionJ * Stab1(i,j)  * rIntegrationWeight;
            rLeftHandSideMatrix(indexpi, indexpj) += rVariables.tau2 * (1.0 / pow(rVariables.BulkModulus,2)) *r_N(0, i) * r_N(0, j) * rIntegrationWeight;

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
// Calculate Projections for OSGS stabilization
//***********************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateProjections(const ProcessInfo &rCurrentProcessInfo)
{

    GeometryType& r_geometry = this->GetGeometry();
    const unsigned int number_of_nodes  = r_geometry.size();
    const unsigned int dimension        = r_geometry.WorkingSpaceDimension();
    const Vector& r_N       = row(GetGeometry().ShapeFunctionsValues(), 0);
    Vector MomentumRes      = ZeroVector(3);
    double ConservRes       = 0.0;
    VectorType momentum_rhs = ZeroVector(number_of_nodes*dimension);
    VectorType conserv_rhs  = ZeroVector(number_of_nodes);
    VectorType NodalArea    = ZeroVector(number_of_nodes);

     // Set the Initialize General Variables
    GeneralVariables Variables;
    this-> InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // Set Variables
    this-> CalculateKinematics(Variables, rCurrentProcessInfo);
    this-> SetGeneralVariables(Variables, Values, r_N);
    this-> SetSpecificVariables(Variables,rCurrentProcessInfo);


    mMP.volume = mMP.mass / mMP.density;
    Vector volume_force = (mMP.volume_acceleration * mMP.mass ) + (Variables.BodyForceMP * mMP.mass);


    //std::cout<< "Nodes"<< number_of_nodes <<"\n";

    // Compute the Residual
    this->ComputeResidual(Variables,volume_force,MomentumRes,ConservRes);

    // Loop over the nodes for projecting
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        double W = mMP.volume*r_N(i);
        unsigned int row = i*dimension;
        for (unsigned int d = 0; d < dimension; d++)
            momentum_rhs[row+d] += W*MomentumRes[d];
        NodalArea[i] += W;
        conserv_rhs[i] += W*ConservRes;
    }

     // Add to Global nodal variables
    for (SizeType i = 0; i < number_of_nodes; ++i)
      {
          r_geometry[i].SetLock(); // So it is safe to write in the node in OpenMP
          array_1d<double,3>& rMomValue = r_geometry[i].FastGetSolutionStepValue(RESPROJ_DISPL);
          unsigned int row = i*dimension;
          for (unsigned int d = 0; d < dimension; d++) rMomValue[d] += momentum_rhs[row + d];

          r_geometry[i].FastGetSolutionStepValue(RESPROJ_PRESS) += conserv_rhs[i];
          r_geometry[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
          r_geometry[i].UnSetLock(); // Free the node for other threads
      }

  
   
}


void MPMUpdatedLagrangianUPVMS::ComputeResidual(GeneralVariables& rVariables, Vector& rVolumeForce, Vector& rResidualU, double& rResidualP)
{

    GeometryType& r_geometry = this->GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    for (unsigned int k = 0; k < dimension; ++k) rResidualU[k] =  rVolumeForce[k] + rVariables.PressureGradient[k] +rVariables.DiscreteAcceleration[k]; //
    rResidualP =  -rVariables.PressureGP/rVariables.BulkModulus + (rVariables.detFT*rVariables.detF0-1);

    //std::cout<<mMP.acceleration[1]<<"\n";

}


void MPMUpdatedLagrangianUPVMS::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMUpdatedLagrangian )
    rSerializer.save("Pressure",m_mp_pressure);
}

void MPMUpdatedLagrangianUPVMS::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMUpdatedLagrangian )
    rSerializer.load("Pressure",m_mp_pressure);
}
} // Namespace Kratos

