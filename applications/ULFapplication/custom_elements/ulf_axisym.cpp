/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-05-16 13:59:01 $
//   Revision:            $Revision: 1.6 $
//
//

//#define GRADPN_FORM
//#define STOKES

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/ulf_axisym.h"
#include "utilities/math_utils.h"
#include "ULF_application.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
//static variables



//************************************************************************************
//************************************************************************************
UlfAxisym::UlfAxisym(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
UlfAxisym::UlfAxisym(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer UlfAxisym::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Element::Pointer(new UlfAxisym(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

UlfAxisym::~UlfAxisym()
{
}

//************************************************************************************
//************************************************************************************
void UlfAxisym::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BoundedMatrix<double,4,6> msB;
    BoundedMatrix<double,4,4> ms_constitutive_matrix;
    BoundedMatrix<double,4,6> ms_temp;
    BoundedMatrix<double,3,2> msDN_Dx;
    array_1d<double,3> msN; //dimension = number of nodes

    const double& density = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
                                         GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
                                         GetGeometry()[2].FastGetSolutionStepValue(DENSITY));
    double K = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS)+
                            GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS) +
                            GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS));
//		double K = GetProperties()[BULK_MODULUS];
    K *= density;

    int is_int = GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
    is_int += GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
    is_int += GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
    

  int is_wall = GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE);
    is_wall += GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);
    is_wall += GetGeometry()[2].FastGetSolutionStepValue(IS_STRUCTURE);
   

    unsigned int dim = 2;
    unsigned int number_of_nodes = 3;

    if(rLeftHandSideMatrix.size1() != 6)
        rLeftHandSideMatrix.resize(6,6,false);

    if(rRightHandSideVector.size() != 6)
        rRightHandSideVector.resize(6,false);

    //calculate current area
    double current_area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, current_area);
    double r=0.33333333333*(GetGeometry()[0].X()+GetGeometry()[1].X()+GetGeometry()[2].X());


    if (is_int!=0 && is_wall!=0)
	{
	//current_area=0.0;
	//mA0=0.0;

	current_area*=0.0000000001;
	mA0*=0.0000000001;
	}


    //writing the body force
    const array_1d<double,3>& body_force = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+
                                           GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +
                                           GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE));

    //axisym coef: 2*pi*r
    for(unsigned int i = 0; i<number_of_nodes; i++)
    {
        //rRightHandSideVector[i*2] = 6.28*r*body_force[0]* density * mA0 * 0.3333333333333;
        //rRightHandSideVector[i*2+1] = 6.28*r*body_force[1] * density * mA0 * 0.3333333333333;
	
	rRightHandSideVector[i*2] = r*body_force[0]* density * mA0 * 0.3333333333333;
        rRightHandSideVector[i*2+1] = r*body_force[1] * density * mA0 * 0.3333333333333;
    }

    //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
    // rLeftHandSideMatrix += Laplacian * nu;
    //filling matrix B
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = dim*i;
        msB(0,index+0)=msDN_Dx(i,0);   msB(0,index+1)= 0.0;
        msB(1,index+0)=0.0;            msB(1,index+1)= msDN_Dx(i,1);
	//the third raw are the "new" axisymmetric entries corresponding to the hoop strains.. so we directly proceed to the fourth row.. third raw is filled separetely
        msB(3,index+0)= msDN_Dx(i,1);  msB(3,index+1)= msDN_Dx(i,0);
    }
    //and now we fill the 3rd raw of the axisymmetric Strain-Displ matrix B: N/r (N/x_center in our case)    
    double norm_r = r*r;
    norm_r=sqrt(norm_r);
    //to avoid the singularity at the central axis (where r=0), one of the tricks is to assume that there epsilon_theta=epsilon_r, or msB(2,i)=msB(0,i)
    
    
    if (norm_r>0.0000000001)
	{
	//r=norm_r;
	msB(2,0)=msN[0]/r;  msB(2,1)=0.0;   	msB(2,2)=msN[1]/r; msB(2,3)=0.0;   	msB(2,4)=msN[2]/r; msB(2,5)=0.0;
	}
    else
	{
	msB(2,0)=msB(0,0);  msB(2,1)=msB(0,1); msB(2,2)=msB(0,2); msB(2,3)=msB(0,3);   msB(2,4)=msB(0,4); msB(2,5)=msB(0,5);
	}    
   
    //constitutive tensor
    /*
    ms_constitutive_matrix(0,0) = K;
    ms_constitutive_matrix(0,1) = K ;
    ms_constitutive_matrix(0,2) = 0.0;
    ms_constitutive_matrix(1,0) = K;
    ms_constitutive_matrix(1,1) = K;
    ms_constitutive_matrix(1,2) = 0.0;
    ms_constitutive_matrix(2,0) = 0.0;
    ms_constitutive_matrix(2,1) = 0.0;
    ms_constitutive_matrix(2,2) = 0.0;
    */
    //TO CHECK!!! IS THIS CONST MATRIX CORRECT??

    ms_constitutive_matrix(0,0) = K;
    ms_constitutive_matrix(0,1) = K;
    ms_constitutive_matrix(0,2) = K;
    ms_constitutive_matrix(0,3) = 0;

    ms_constitutive_matrix(1,0) = K;
    ms_constitutive_matrix(1,1) = K;
    ms_constitutive_matrix(1,2) = K;
    ms_constitutive_matrix(1,3) = 0;

    ms_constitutive_matrix(2,0) = K;
    ms_constitutive_matrix(2,1) = K;
    ms_constitutive_matrix(2,2) = K;
    ms_constitutive_matrix(2,3) = 0;

    ms_constitutive_matrix(3,0) = 0.0;
    ms_constitutive_matrix(3,1) = 0.0;
    ms_constitutive_matrix(3,2) = 0.0;
    ms_constitutive_matrix(3,3) = 0.0;


    //calculating viscous contributions
    ms_temp = prod( ms_constitutive_matrix , msB);
    noalias(rLeftHandSideMatrix) = prod( trans(msB) , ms_temp);
    rLeftHandSideMatrix *= -current_area;
    //axisymetric formulation=> K=2pi r (BTDB) A
    //2*pi=6.28
    //rLeftHandSideMatrix *=6.28*r;
    rLeftHandSideMatrix *=r;


    //get the old pressure
    double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
    double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
    double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);

    //adding pressure gradient
    double pavg = p0old + p1old + p2old; //calculate old pressure over the element     
    pavg *= 0.33333333333333333333333 * current_area;

    //add to pavg the increment of pressure inside the time step
    double dp = K * (current_area - mA0);
    pavg += dp;
    //axisymetric formulation=> f=2pi r (BTDB) A
    //2*pi=6.28
    //pavg*=6.28*r;
    pavg*=r;

    //KRATOS_WATCH(msDN_Dx)
    rRightHandSideVector[0] += (msDN_Dx(0,0)+0.33333/r)*pavg;
    rRightHandSideVector[1] += msDN_Dx(0,1)*pavg;
    rRightHandSideVector[2] += (msDN_Dx(1,0)+0.33333/r)*pavg;
    rRightHandSideVector[3] += msDN_Dx(1,1)*pavg;
    rRightHandSideVector[4] += (msDN_Dx(2,0)+0.33333/r)*pavg;
    rRightHandSideVector[5] += msDN_Dx(2,1)*pavg;
	
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void UlfAxisym::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
KRATOS_THROW_ERROR(std::logic_error,"function not implemeneted for this element","");



}

//************************************************************************************
//************************************************************************************
void UlfAxisym::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    int is_int = GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
    is_int += GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
    is_int += GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
    

  int is_wall = GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE);
    is_wall += GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);
    is_wall += GetGeometry()[2].FastGetSolutionStepValue(IS_STRUCTURE);


    const double& density = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
                                         GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
                                         GetGeometry()[2].FastGetSolutionStepValue(DENSITY));
    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int NumberOfNodes = GetGeometry().size();
    unsigned int MatSize = dimension * NumberOfNodes;

    if(rMassMatrix.size1() != MatSize)
        rMassMatrix.resize(MatSize,MatSize,false);

    noalias(rMassMatrix) = ZeroMatrix(MatSize,MatSize);

    if (is_int!=0 && is_wall!=0)
	{
	//mA0=0.0;
	mA0*=0.0000000001;
	}


    double nodal_mass = mA0 * density * 0.333333333333333333;

    for(unsigned int i=0; i<NumberOfNodes; i++)
    {
        for(unsigned int j=0; j<dimension; j++)
        {
            unsigned int index = i*dimension + j;
            rMassMatrix(index,index) = nodal_mass;
        }
    }
    //axisymetric formulation
    double r=0.33333333333*(GetGeometry()[0].X()+GetGeometry()[1].X()+GetGeometry()[2].X());
    //rMassMatrix*=6.28*r;


    rMassMatrix*=r;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void UlfAxisym::CalculateDampingMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int is_int = GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
    is_int += GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
    is_int += GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
    

  int is_wall = GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE);
    is_wall += GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);
    is_wall += GetGeometry()[2].FastGetSolutionStepValue(IS_STRUCTURE);

    //KRATOS_WATCH("ULF_AXiSYM")
    BoundedMatrix<double,4,6> msB;
    BoundedMatrix<double,4,4> ms_constitutive_matrix;
    BoundedMatrix<double,4,6> ms_temp;
    BoundedMatrix<double,3,2> msDN_Dx;
    array_1d<double,3> msN; //dimension = number of nodes

    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    if(rDampMatrix.size1() != 6)
        rDampMatrix.resize(6,6,false);


    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_Dx, msN, Area);


    if (is_int!=0 && is_wall!=0)
	{
	
	Area*=0.0000000001;
	}

    //getting properties
    const double& nu = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY)+
                                    GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY) +
                                    GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY));

    const double& density = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
                                         GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
                                         GetGeometry()[2].FastGetSolutionStepValue(DENSITY));
    //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
    // rLeftHandSideMatrix += Laplacian * nu;
    //filling matrix B
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = dim*i;
        msB(0,index+0)=msDN_Dx(i,0);
        msB(0,index+1)= 0.0;
        msB(1,index+0)=0.0;
        msB(1,index+1)= msDN_Dx(i,1);
	//the thirrd raw are the "new" axisymmetric entries corresponding to the hoop strains.. so we directly proceed to the fourth row.. third raw is filled separetely
        msB(3,index+0)= msDN_Dx(i,1);
        msB(3,index+1)= msDN_Dx(i,0);
    }
    //and now we fill the 3rd raw of the axisymmetric Strain-Displ matrix B: N/r (N/x_center in our case)
    double r=0.33333333333*(GetGeometry()[0].X()+GetGeometry()[1].X()+GetGeometry()[2].X());
    double norm_r = r*r;
    norm_r=sqrt(norm_r);
	//KRATOS_WATCH(r)
	//KRATOS_WATCH(norm_r)
    //to avoid the singularity at the central axis (where r=0), one of the tricks is to assume that there epsilon_theta=epsilon_r, or msB(2,i)=msB(0,i)
    
    if (norm_r>0.0000000001)
	{
	//r=norm_r;
	msB(2,0)=msN[0]/r;  msB(2,1)=0.0;   msB(2,2)=msN[1]/r; msB(2,3)=0.0;   msB(2,4)=msN[2]/r; msB(2,5)=0.0;
	}
    else
	{
	msB(2,0)=msB(0,0);  msB(2,1)=msB(0,1); msB(2,2)=msB(0,2); msB(2,3)=msB(0,3);  msB(2,4)=msB(0,4); msB(2,5)=msB(0,5);
	}
  
	//KRATOS_WATCH(msB)
    


    //constitutive tensor
    const double& a = nu*density;
    ms_constitutive_matrix(0,0) = (4.0/3.0)*a; 	 ms_constitutive_matrix(0,1) = (-2.0/3.0)*a; 	ms_constitutive_matrix(0,2) = (-2.0/3.0)*a; ms_constitutive_matrix(0,3) = 0.0;
    ms_constitutive_matrix(1,0) = (-2.0/3.0)*a;  ms_constitutive_matrix(1,1) = (4.0/3.0)*a;  	ms_constitutive_matrix(1,2) = (-2.0/3.0)*a; ms_constitutive_matrix(1,3) = 0.0;
    ms_constitutive_matrix(2,0) = (-2.0/3.0)*a;  ms_constitutive_matrix(2,1) = (-2.0/3.0)*a;    ms_constitutive_matrix(2,2) = (4.0/3.0)*a;  ms_constitutive_matrix(2,3) = 0.0;
    ms_constitutive_matrix(3,0) = 0.0;     	 ms_constitutive_matrix(3,1) = 0.0;		ms_constitutive_matrix(3,2) = 0.0;	    ms_constitutive_matrix(3,3) = a;

    //calculating viscous contributions
    ms_temp = prod( ms_constitutive_matrix , msB);
    noalias(rDampMatrix) = prod( trans(msB) , ms_temp);

    rDampMatrix *= Area;
    //axisymetric formulation=> K=2pi r (BTDB) A
    //2*pi=6.28
    //rDampMatrix *=6.28*r;
    rDampMatrix *=r;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure
void UlfAxisym::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    //save original Area
    mA0 = GeometryUtils::CalculateVolume2D(GetGeometry());
    //save the original radius
    r0=0.33333333333*(GetGeometry()[0].X()+GetGeometry()[1].X()+GetGeometry()[2].X());


    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
void UlfAxisym::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 2;

    if(rResult.size() != number_of_nodes*dim)
        rResult.resize(number_of_nodes*dim,false);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i*dim] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[i*dim+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
    }

}

//************************************************************************************
//************************************************************************************
void UlfAxisym::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 2;

    if(ElementalDofList.size() != number_of_nodes*dim)
        ElementalDofList.resize(number_of_nodes*dim);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i*dim] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        ElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
    }
}

//************************************************************************************
//************************************************************************************
void UlfAxisym::GetValuesVector(Vector& values, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes*dim;
    if(values.size() != MatSize)	values.resize(MatSize,false);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = i*dim;
        const array_1d<double,3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,Step);
        values[index] = disp[0];
        values[index + 1] = disp[1];
    }
}


//************************************************************************************
//************************************************************************************
void UlfAxisym::GetFirstDerivativesVector(Vector& values, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes*dim;
    if(values.size() != MatSize)   values.resize(MatSize,false);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = i*dim;
        const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,Step);
        values[index] = vel[0];
        values[index + 1] = vel[1];
    }
}
//************************************************************************************
//************************************************************************************
void UlfAxisym::GetSecondDerivativesVector(Vector& values, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes*dim;
    if(values.size() != MatSize) values.resize(MatSize,false);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = i*dim;
        const array_1d<double,3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION,Step);
        values[index] = acc[0];
        values[index + 1] = acc[1];
    }
}

//************************************************************************************
//************************************************************************************
void  UlfAxisym::Calculate(const Variable<double >& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
{
    if(rVariable == PRESSURE)
    {
	    int is_int = GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
	    is_int += GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
	    is_int += GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
	    

	  int is_wall = GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE);
	    is_wall += GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);
	    is_wall += GetGeometry()[2].FastGetSolutionStepValue(IS_STRUCTURE);

        array_1d<double,3> msN; //dimension = number of nodes

        const double density = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
                                            GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
                                            GetGeometry()[2].FastGetSolutionStepValue(DENSITY));

        //const double& density = GetProperties()[DENSITY];
        double K = 0.333333333*(GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS)+
                                GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS) +
                                GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS));
        //		double K = GetProperties()[BULK_MODULUS];
       K *= density;

      double r=0.33333333333*(GetGeometry()[0].X()+GetGeometry()[1].X()+GetGeometry()[2].X());
      double norm_r = r*r;
      norm_r=sqrt(norm_r);

        double current_area = GeometryUtils::CalculateVolume2D(GetGeometry());

    if (is_int!=0 && is_wall!=0)
	{
//	current_area=0.0;
//	mA0=0.0;
	current_area*=0.0000000001;
	mA0*=0.0000000001;


	}

        

        //save in msN the old pressure
        msN[0] = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
        msN[1] = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
        msN[2] = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
        double diag_term = current_area/6.0;
        double out_term = current_area/12.0;



        //add to pavg the increment of pressure inside the time step
        //double dp_area = r * K * 0.3333333333333333333*(current_area - mA0);
	double dp_area = K * 0.3333333333333333333*(r*current_area - r0*mA0)/r;

        //P = Mconsistent*pold + dp
        		GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) += dp_area;
        		GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) += dp_area;
        		GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) += dp_area;
        //GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) += (dp_area + diag_term*msN[0] + out_term*msN[1]  + out_term*msN[2]) ;
        //GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) += (dp_area + out_term*msN[0] + diag_term*msN[1]  + out_term*msN[2]) ;
        //GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) += (dp_area + out_term*msN[0] + out_term*msN[1]  + diag_term*msN[2]) ;

    }
    else if (rVariable == IS_FLUID)
    {
        GetGeometry()[0].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;
        GetGeometry()[1].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;
        GetGeometry()[2].FastGetSolutionStepValue(IS_FLUID) = 1.0 ;

    }


}

} // Namespace Kratos


