//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pavel Ryzhakov and Julio Marti
//

//#define GRADPN_FORM
//#define STOKES

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/hypo_elastic_solid3d.h"
#include "utilities/math_utils.h"
#include "ULF_application.h"
#include "utilities/geometry_utilities.h"
#include "includes/kratos_flags.h"

namespace Kratos
{
//static variables



//************************************************************************************
//************************************************************************************
HypoElasticSolid3D::HypoElasticSolid3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
HypoElasticSolid3D::HypoElasticSolid3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer HypoElasticSolid3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Element::Pointer(new HypoElasticSolid3D(NewId, GetGeometry().Create(ThisNodes), pProperties));

    KRATOS_CATCH("");
}

HypoElasticSolid3D::~HypoElasticSolid3D()
{
}

//************************************************************************************
//************************************************************************************
void HypoElasticSolid3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    if (rLeftHandSideMatrix.size1() != 12)
        rLeftHandSideMatrix.resize(12, 12, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(12, 12);
    //fill in the RHS
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void HypoElasticSolid3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if(rRightHandSideVector.size() != 12)
        rRightHandSideVector.resize(12,false);

    noalias(rRightHandSideVector) = ZeroVector(12);

    BoundedMatrix<double,4,3> msDN_Dx;
    array_1d<double,4> msN; //dimension = number of nodes
    double current_vol;

    const double& density = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
                                         GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
                                         GetGeometry()[2].FastGetSolutionStepValue(DENSITY)+ 
                                         GetGeometry()[3].FastGetSolutionStepValue(DENSITY));

    unsigned int number_of_nodes = GetGeometry().size();

    //calculate current area
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, current_vol);

    //writing the body force
    const array_1d<double,3>& body_force = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+
                                           GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +
                                           GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE)+
                                           GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE));

    for(unsigned int i = 0; i<number_of_nodes; i++)
    {
        //rRightHandSideVector[i*3] = body_force[0]* density * mV0 * 0.25;
        //rRightHandSideVector[i*3+1] = body_force[1] * density * mV0 * 0.25;
        //rRightHandSideVector[i*3+2] = body_force[2] * density * mV0 * 0.25;
	rRightHandSideVector[i*3] = body_force[0]* density * 0.25;
        rRightHandSideVector[i*3+1] = body_force[1] * density  * 0.25;
        rRightHandSideVector[i*3+2] = body_force[2] * density  * 0.25;
    }

    //get the value of Cauchy stress at the Gauss point. It is given by: 
    const BoundedMatrix<double,3,3> & CauchyStress=this->GetValue(CAUCHY_STRESS_TENSOR);

//TO DO WITH JULIO!!!!!

    //dN1/dx*SigmaX + 0 + 0 + dN1/dy*TauXY + dN1/dz*TauXZ + 0
    rRightHandSideVector[0] -= msDN_Dx(0,0)*CauchyStress(0,0) + msDN_Dx(0,1)*CauchyStress(0,1) + msDN_Dx(0,2)*CauchyStress(0,2);
    //0  +   dN1/dy*SigmaY + 0 + dN1/dx*TauXY + 0 + dN1/dz*TauYZ
    rRightHandSideVector[1] -= msDN_Dx(0,1)*CauchyStress(1,1) + msDN_Dx(0,0)*CauchyStress(0,1) + msDN_Dx(0,2)*CauchyStress(1,2);
    //0  +     0    + dN1/dz*SigmaZ + 0 + dN1/dx*TauXZ + dN1/dy*TauYZ
    rRightHandSideVector[2] -= msDN_Dx(0,2)*CauchyStress(2,2) + msDN_Dx(0,0)*CauchyStress(0,2) + msDN_Dx(0,1)*CauchyStress(1,2);


   //dN2/dx*SigmaX + 0 + 0 + dN2/dy*TauXY + dN2/dz*TauXZ + 0
    rRightHandSideVector[3] -= msDN_Dx(1,0)*CauchyStress(0,0) + msDN_Dx(1,1)*CauchyStress(0,1) + msDN_Dx(1,2)*CauchyStress(0,2);
    //0  +   dN2/dy*SigmaY + 0 + dN2/dx*TauXY + 0 + dN2/dz*TauYZ
    rRightHandSideVector[4] -= msDN_Dx(1,1)*CauchyStress(1,1) + msDN_Dx(1,0)*CauchyStress(0,1) + msDN_Dx(1,2)*CauchyStress(1,2);
    //0  +     0    + dN2/dz*SigmaZ + 0 + dN2/dx*TauXZ + dN2/dy*TauYZ
    rRightHandSideVector[5] -= msDN_Dx(1,2)*CauchyStress(2,2) + msDN_Dx(1,0)*CauchyStress(0,2) + msDN_Dx(1,1)*CauchyStress(1,2);


    //dN3/dx*SigmaX + 0 + 0 + dN3/dy*TauXY + dN3/dz*TauXZ + 0
    rRightHandSideVector[6] -= msDN_Dx(2,0)*CauchyStress(0,0) + msDN_Dx(2,1)*CauchyStress(0,1) + msDN_Dx(2,2)*CauchyStress(0,2);
    //0  +   dN3/dy*SigmaY + 0 + dN2/dx*TauXY + 0 + dN3/dz*TauYZ
    rRightHandSideVector[7] -= msDN_Dx(2,1)*CauchyStress(1,1) + msDN_Dx(2,0)*CauchyStress(0,1) + msDN_Dx(2,2)*CauchyStress(1,2);
    //0  +     0    + dN3/dz*SigmaZ + 0 + dN3/dx*TauXZ + dN3/dy*TauYZ
    rRightHandSideVector[8] -= msDN_Dx(2,2)*CauchyStress(2,2) + msDN_Dx(2,0)*CauchyStress(0,2) + msDN_Dx(2,1)*CauchyStress(1,2);

    //dN4/dx*SigmaX + 0 + 0 + dN4/dy*TauXY + dN4/dz*TauXZ + 0
    rRightHandSideVector[9] -= msDN_Dx(3,0)*CauchyStress(0,0) + msDN_Dx(3,1)*CauchyStress(0,1) + msDN_Dx(3,2)*CauchyStress(0,2);
    //0  +   dN4/dy*SigmaY + 0 + dN4/dx*TauXY + 0 + dN4/dz*TauYZ
    rRightHandSideVector[10] -= msDN_Dx(3,1)*CauchyStress(1,1) + msDN_Dx(3,0)*CauchyStress(0,1) + msDN_Dx(3,2)*CauchyStress(1,2);
    //0  +     0    + dN4/dz*SigmaZ + 0 + dN4/dx*TauXZ + dN4/dy*TauYZ
    rRightHandSideVector[11] -= msDN_Dx(3,2)*CauchyStress(2,2) + msDN_Dx(3,0)*CauchyStress(0,2) + msDN_Dx(3,1)*CauchyStress(1,2);

    rRightHandSideVector*=current_vol;
    //KRATOS_WATCH(rRightHandSideVector)

}

//************************************************************************************
//************************************************************************************
void HypoElasticSolid3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const double& density = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
                                         GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
                                         GetGeometry()[2].FastGetSolutionStepValue(DENSITY)+ 
                                         GetGeometry()[3].FastGetSolutionStepValue(DENSITY));
    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int NumberOfNodes = GetGeometry().size();

    if(rMassMatrix.size1() != 12)
        rMassMatrix.resize(12,12,false);

    noalias(rMassMatrix) = ZeroMatrix(12,12);

    double nodal_mass = mV0 * density * 0.25;

    for(unsigned int i=0; i<NumberOfNodes; i++)
    {
        for(unsigned int j=0; j<dimension; j++)
        {
            unsigned int index = i*dimension + j;
            rMassMatrix(index,index) = nodal_mass;
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************



//void HypoElasticSolid3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
void HypoElasticSolid3D::CalculateLocalVelocityContribution(MatrixType& rDampingMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    if (rDampingMatrix.size1() != 12)
        rDampingMatrix.resize(12, 12, false);

    noalias(rDampingMatrix) = ZeroMatrix(12, 12);

    //fill in the damping matrix
    BoundedMatrix<double,6,12> msB = ZeroMatrix(6,12);
    BoundedMatrix<double,6,6> ms_constitutive_matrix;
    BoundedMatrix<double,6,12> ms_temp;

    BoundedMatrix<double,4,3> msDN_Dx;
    array_1d<double,4> msN; //dimension = number of nodes

    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //getting data for the given geometry
    double current_vol;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, current_vol);

    double NU = GetProperties()[POISSON_RATIO];
    double E = GetProperties()[YOUNG_MODULUS];

    double dt = rCurrentProcessInfo[DELTA_TIME];
    //Lame constants. note that for the hypoelastic solid the Lame constants must be multiplied by the timestep
    //const 
    double MU=0.5*dt*E/(1.0+NU);
    //const 
    double LAMBDA=NU*E*dt/((1.0+NU)*(1.0-2.0*NU));
    //const 
    double KAPPA=LAMBDA+0.6666666*MU;

    unsigned int number_of_nodes = GetGeometry().size();


    //SHEAR CONTRIBUTION TO THE "DAMPING" MATRIX
    for(unsigned int i = 0; i<number_of_nodes; i++)
    {
        unsigned int start = dim*i;

        msB(0,start) =	msDN_Dx(i,0);
        msB(1,start+1)=	msDN_Dx(i,1);
        msB(2,start+2)= msDN_Dx(i,2);
        msB(3,start) =	msDN_Dx(i,1);
        msB(3,start+1) = msDN_Dx(i,0);
        msB(4,start) =	msDN_Dx(i,2);
        msB(4,start+2) = msDN_Dx(i,0);
        msB(5,start+1)= msDN_Dx(i,2);
        msB(5,start+2) = msDN_Dx(i,1);
    }

    //constitutive tensor
    ms_constitutive_matrix(0,0) = (4.0/3.0)*MU;
    ms_constitutive_matrix(0,1) = -2.0/3.0*MU;
    ms_constitutive_matrix(0,2) = -2.0/3.0*MU;
    ms_constitutive_matrix(0,3) = 0.0;
    ms_constitutive_matrix(0,4) = 0.0;
    ms_constitutive_matrix(0,5) = 0.0;

    ms_constitutive_matrix(1,0) = -2.0/3.0*MU;
    ms_constitutive_matrix(1,1) = 4.0/3.0*MU;
    ms_constitutive_matrix(1,2) = -2.0/3.0*MU;
    ms_constitutive_matrix(1,3) = 0.0;
    ms_constitutive_matrix(1,4) = 0.0;
    ms_constitutive_matrix(1,5) = 0.0;

    ms_constitutive_matrix(2,0) = -2.0/3.0*MU;
    ms_constitutive_matrix(2,1) = -2.0/3.0*MU;
    ms_constitutive_matrix(2,2) = 4.0/3.0*MU;
    ms_constitutive_matrix(2,3) = 0.0;
    ms_constitutive_matrix(2,4) = 0.0;
    ms_constitutive_matrix(2,5) = 0.0;

    ms_constitutive_matrix(3,0) = 0.0;
    ms_constitutive_matrix(3,1) = 0.0;
    ms_constitutive_matrix(3,2) = 0.0;
    ms_constitutive_matrix(3,3) = MU;
    ms_constitutive_matrix(3,4) = 0.0;
    ms_constitutive_matrix(3,5) = 0.0;

    ms_constitutive_matrix(4,0) = 0.0;
    ms_constitutive_matrix(4,1) = 0.0;
    ms_constitutive_matrix(4,2) = 0.0;
    ms_constitutive_matrix(4,3) = 0.0;
    ms_constitutive_matrix(4,4) = MU;
    ms_constitutive_matrix(4,5) = 0.0;

    ms_constitutive_matrix(5,0) = 0.0;
    ms_constitutive_matrix(5,1) = 0.0;
    ms_constitutive_matrix(5,2) = 0.0;
    ms_constitutive_matrix(5,3) = 0.0;
    ms_constitutive_matrix(5,4) = 0.0;
    ms_constitutive_matrix(5,5) = MU;

    //calculating viscous contributions
    ms_temp = prod( ms_constitutive_matrix , msB);
    noalias(rDampingMatrix) = prod( trans(msB) , ms_temp);

    
    //now we reuse the constitutive tensor
    //only the 3x3 left corner block is filled, the rest of the valyes are zero
    ms_constitutive_matrix=ZeroMatrix(6,6);
    ms_constitutive_matrix(0,0) = KAPPA;
    ms_constitutive_matrix(0,1) = KAPPA ;
    ms_constitutive_matrix(0,2) = KAPPA ;
    
    ms_constitutive_matrix(1,0) = KAPPA;
    ms_constitutive_matrix(1,1) = KAPPA ;
    ms_constitutive_matrix(1,2) = KAPPA ;

    ms_constitutive_matrix(2,0) = KAPPA;
    ms_constitutive_matrix(2,1) = KAPPA ;
    ms_constitutive_matrix(2,2) = KAPPA ;

    //calculating volumetric contribution
    ms_temp = prod( ms_constitutive_matrix , msB);
    //ms_temp*=dt;
    rDampingMatrix+= prod( trans(msB) , ms_temp);
    
    rDampingMatrix *= current_vol;

    //Now calculate an additional contribution to the residual: r -= rDampingMatrix * (v)
    array_1d< double, 12 > Vel;
    Vel[0] = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X);
    Vel[1] = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y);
    Vel[2] = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z);

    Vel[3] = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X);
    Vel[4] = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Y);
    Vel[5] = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Z);
    
    Vel[6] = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_X);
    Vel[7] = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_Y);
    Vel[8] = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_Z);

    Vel[9] = this->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY_X);
    Vel[10] = this->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY_Y);
    Vel[11] = this->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY_Z);

    noalias(rRightHandSideVector) -= prod(rDampingMatrix, Vel);
    //KRATOS_WATCH(rDampingMatrix)

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void HypoElasticSolid3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    //save original Area
    mV0 = GeometryUtils::CalculateVolume3D(GetGeometry());
    GetGeometry()[0].Set(STRUCTURE, true);
    GetGeometry()[1].Set(STRUCTURE, true);
    GetGeometry()[2].Set(STRUCTURE, true);
    GetGeometry()[3].Set(STRUCTURE, true);

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
void HypoElasticSolid3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 3;

    if(rResult.size() != number_of_nodes*dim)
        rResult.resize(number_of_nodes*dim,false);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i*dim] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
        rResult[i*dim+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        rResult[i*dim+2] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
    }

}

//************************************************************************************
//************************************************************************************
void HypoElasticSolid3D::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 3;

    if(rElementalDofList.size() != number_of_nodes*dim)
        rElementalDofList.resize(number_of_nodes*dim);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rElementalDofList[i*dim] = GetGeometry()[i].pGetDof(VELOCITY_X);
        rElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
        rElementalDofList[i*dim+2] = GetGeometry()[i].pGetDof(VELOCITY_Z);
    }
}

//************************************************************************************
//************************************************************************************
void HypoElasticSolid3D::GetValuesVector(Vector& rValues, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = 3;
    //const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes*dim;

    if(rValues.size() != MatSize)	
          rValues.resize(MatSize,false);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = i*dim;
        const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,Step);
        rValues[index] = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];
    }
}


//************************************************************************************
//************************************************************************************
void HypoElasticSolid3D::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes*dim;
    if(rValues.size() != MatSize)   rValues.resize(MatSize,false);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = i*dim;
        const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,Step);
        rValues[index] = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];
    }
}
//************************************************************************************
//************************************************************************************
void HypoElasticSolid3D::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes*dim;
    if(rValues.size() != MatSize) rValues.resize(MatSize,false);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = i*dim;
        const array_1d<double,3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION,Step);
        rValues[index] = acc[0];
        rValues[index + 1] = acc[1];
        rValues[index + 2] = acc[2];
    }
}

//************************************************************************************
//************************************************************************************



void  HypoElasticSolid3D::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo)
{
    if(rVariable == CAUCHY_STRESS_TENSOR)
    {
    BoundedMatrix<double,3,3> CauchyStress=ZeroMatrix(3,3);
    BoundedMatrix<double,3,3> HistoricalCauchyStress=ZeroMatrix(3,3);

    const array_1d<double,3>& v0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& v1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& v2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& v3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);
    
    double current_area;
    BoundedMatrix<double,4,3> msDN_Dx;
    array_1d<double,4> msN; //dimension = number of nodes
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, current_area);

    double NU = GetProperties()[POISSON_RATIO];
    double E = GetProperties()[YOUNG_MODULUS];

    double dt = rCurrentProcessInfo[DELTA_TIME];
    //Lame constants. note that for the hypoelastic solid the Lame constants must be multiplied by the timestep
    //const 
    double MU=0.5*dt*E/(1.0+NU);
    //const 
    double LAMBDA=NU*E*dt/((1.0+NU)*(1.0-2.0*NU));
    //const 
    double KAPPA=LAMBDA+0.6666666*MU;

    //Symmetric stress tensor
    CauchyStress(0,0)=2.0*(msDN_Dx(0,0)*v0[0]+msDN_Dx(1,0)*v1[0]+msDN_Dx(2,0)*v2[0] + msDN_Dx(3,0)*v3[0]);
    CauchyStress(0,1)=msDN_Dx(0,1)*v0[0] + msDN_Dx(1,1)*v1[0] +  msDN_Dx(2,1)*v2[0] + msDN_Dx(3,1)*v3[0]    +       msDN_Dx(0,0)*v0[1]+msDN_Dx(1,0)*v1[1]+msDN_Dx(2,0)*v2[1] + msDN_Dx(3,0)*v3[1] ;
    CauchyStress(0,2)=msDN_Dx(0,2)*v0[0] + msDN_Dx(1,2)*v1[0] + msDN_Dx(2,2)*v2[0] + msDN_Dx(3,2)*v3[0]     +       msDN_Dx(0,0)*v0[2]+msDN_Dx(1,0)*v1[2]+msDN_Dx(2,0)*v2[2] + msDN_Dx(3,0)*v3[2] ;

    CauchyStress(1,0)=msDN_Dx(0,1)*v0[0] + msDN_Dx(1,1)*v1[0] +  msDN_Dx(2,1)*v2[0] + msDN_Dx(3,1)*v3[0]    +       msDN_Dx(0,0)*v0[1]+msDN_Dx(1,0)*v1[1]+msDN_Dx(2,0)*v2[1] + msDN_Dx(3,0)*v3[1] ;
    CauchyStress(1,1)=2.0*(msDN_Dx(0,1)*v0[1]+msDN_Dx(1,1)*v1[1]+msDN_Dx(2,1)*v2[1] + msDN_Dx(3,1)*v3[1]);
    CauchyStress(1,2)=msDN_Dx(0,2)*v0[1] + msDN_Dx(1,2)*v1[1] + msDN_Dx(2,2)*v2[1] + msDN_Dx(3,2)*v3[1]     +       msDN_Dx(0,1)*v0[2] + msDN_Dx(1,1)*v1[2]+msDN_Dx(2,1)*v2[2] + msDN_Dx(3,1)*v3[2];

    CauchyStress(2,0)=msDN_Dx(0,2)*v0[0] + msDN_Dx(1,2)*v1[0] + msDN_Dx(2,2)*v2[0] + msDN_Dx(3,2)*v3[0]     +       msDN_Dx(0,0)*v0[2]+msDN_Dx(1,0)*v1[2]+msDN_Dx(2,0)*v2[2] + msDN_Dx(3,0)*v3[2] ;
    CauchyStress(2,1)=msDN_Dx(0,2)*v0[1] + msDN_Dx(1,2)*v1[1] + msDN_Dx(2,2)*v2[1] + msDN_Dx(3,2)*v3[1]     +       msDN_Dx(0,1)*v0[2] + msDN_Dx(1,1)*v1[2]+msDN_Dx(2,1)*v2[2] + msDN_Dx(3,1)*v3[2];
    CauchyStress(2,2)=2.0*(msDN_Dx(0,2)*v0[2]+msDN_Dx(1,2)*v1[2]+msDN_Dx(2,2)*v2[2] + msDN_Dx(3,2)*v3[2]);

    CauchyStress*=MU;

    //KRATOS_WATCH(CauchyStress)

    //adding the volumetric part
    double div_v = msDN_Dx(0,0)*v0[0] + msDN_Dx(0,1)*v0[1] + msDN_Dx(0,2)*v0[2] ;
    div_v+=       msDN_Dx(1,0)*v1[0] + msDN_Dx(1,1)*v1[1] + msDN_Dx(1,2)*v1[2];
    div_v+=	  msDN_Dx(2,0)*v2[0] + msDN_Dx(2,1)*v2[1] + msDN_Dx(2,2)*v2[2];
    div_v+=	  msDN_Dx(3,0)*v3[0] + msDN_Dx(3,1)*v3[1] + msDN_Dx(3,2)*v3[2];

    CauchyStress(0,0)+=KAPPA*div_v;
    CauchyStress(1,1)+=KAPPA*div_v;
    CauchyStress(2,2)+=KAPPA*div_v;

    HistoricalCauchyStress=this->GetValue(CAUCHY_STRESS_TENSOR);
    CauchyStress+=HistoricalCauchyStress;
    //KRATOS_WATCH("After adding divergence")
    //KRATOS_WATCH(CauchyStress)

    this->SetValue(CAUCHY_STRESS_TENSOR, CauchyStress);

    }
    else 
       KRATOS_ERROR << "Wrong variable. Calculate function of hypoelastic element is meant to compute Cauchy stress only." << std::endl;


}

} // Namespace Kratos


