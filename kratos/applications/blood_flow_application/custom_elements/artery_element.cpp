/* b
==============================================================================
KratosIncompressibleFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
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
//   Date:                $Date: 2009-01-13 15:39:56 $
//   Revision:            $Revision: 1.14 $
//
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/artery_element.h"
#include "utilities/math_utils.h"
#include "blood_flow_application.h"
#include "iostream"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

ArteryElement::ArteryElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

ArteryElement::ArteryElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer ArteryElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Element::Pointer(new ArteryElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

ArteryElement::~ArteryElement()
{
}

//************************************************************************************
//************************************************************************************

void ArteryElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense to computer the system matrix for an explicit element", "");
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ArteryElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //resize the vector to the correct size
    if (rRightHandSideVector.size() != 4)
        rRightHandSideVector.resize(4,false);
    
    double h_int = rCurrentProcessInfo[DELTA_TIME];
	//KRATOS_WATCH(this->Id())
	//KRATOS_WATCH(h_int)
    
        const double A0_actual = GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
	    const double A1_actual = GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA);
	     
	     if((A0_actual == 0.00) || (A1_actual== 0.00))
	     {
            //KRATOS_WATCH("------sssssssssssssssssssssssssss-----------");
            //KRATOS_WATCH("----");
            KRATOS_WATCH(A0_actual);
            KRATOS_WATCH(A1_actual);
            KRATOS_WATCH(GetProperties().Id());
                std::cout << "this element is" << this->Id() <<std::endl;
            KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, Please check your model", "");
	     }
	        
            //const int kk =GetProperties().Id();
       //KRATOS_WATCH(h_int)
            //KRATOS_WATCH(A0_actual);
        
    //get data as needed
    const double dynamic_viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    const double density = GetProperties()[DENSITY];
    //const double E = 0.5*(GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS) + GetGeometry()[1].FastGetSolutionStepValue(YOUNG_MODULUS));
    //const double nu =0.5*(GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO) + GetGeometry()[1].FastGetSolutionStepValue(POISSON_RATIO));
    //const double H0=GetGeometry()[0].FastGetSolutionStepValue(THICKNESS);

    const double pi = 3.14159265;
    const double coriolis_coefficient = 1.1;
    //const double kr_coefficient = 1.0;
    const double kr_coefficient = 8.0*pi* (dynamic_viscosity/density);

    //const double kinematic_viscosity = dynamic_viscosity/density;
    //const double H0 = 0.5*(GetGeometry()[0].FastGetSolutionStepValue(THICKNESS) + GetGeometry()[1].FastGetSolutionStepValue(THICKNESS));;
    //const double beta = E*H0*1.77245385/(1.0-nu*nu);


    std::vector< array_1d<double,2> > Fj(2);
    std::vector< array_1d<double,2> > Sj(2);
    std::vector< boost::numeric::ublas::bounded_matrix<double, 2,2 > > Hj(2);
    std::vector< boost::numeric::ublas::bounded_matrix<double, 2,2 > > Suj(2);

    boost::numeric::ublas::bounded_matrix<double, 2,2 > M1, M2;
    
    M1(0,0) = -0.5; M1(0,1) = -0.5; 
    M1(1,0) = 0.5; M1(1,1) = 0.5;
    
    //Consistent Matrix type 1
    M2(0,0) = 0.333333333333333333333333; M2(0,1) = 0.166666666666666666666667;
    M2(1,0) = 0.166666666666666666666667; M2(1,1) = 0.333333333333333333333333; 

    //Consistent Matrix type 2
    //M2(0,0) = 0.5; M2(0,1) = 0;
    //M2(1,0) = 0; M2(1,1) = 0.5;
    //std::cout << "this element is" << this->Id() <<std::endl;
    M2 *= mL; // mL*M2; //fabs(mL);
    //KRATOS_WATCH(M2);
    //sleep(2000);

   //getch();
   //loop on nodes
    for (unsigned int i=0; i<2; i++)
    {
        const double& A = GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA);
        const double beta=GetGeometry()[i].FastGetSolutionStepValue(BETA);
        const double C = (beta*sqrt(A*A*A))/(3.0*density*GetGeometry()[i].GetValue(NODAL_AREA));
        //KRATOS_WATCH(C);GetValue
        //KRATOS_WATCH(beta);
        //KRATOS_WATCH(mL);
        const double flow = GetGeometry()[i].FastGetSolutionStepValue(FLOW);
//       if(Id() == 1)
//       {
//           KRATOS_WATCH("FLOW EN EL ELEMENTO");
//           KRATOS_WATCH(flow);
//       }

//KRATOS_WATCH(A);
//KRATOS_WATCH(beta);
//KRATOS_WATCH(C);


        Fj[i][0] = flow;
        Fj[i][1] = C + ((coriolis_coefficient*flow*flow)/A);

        Sj[i][0] = 0.0;
        Sj[i][1] = (-kr_coefficient*flow)/(A);
        //KRATOS_WATCH(Fj[i]);
        //KRATOS_WATCH(Sj[i]);
        Hj[i](0,0)   = 0.0;
        Hj[i](0,1)   = 1.0;
        Hj[i](1,0)   = (-coriolis_coefficient*pow(flow/A,2)) + ((beta* sqrt(A))/(2.0*density*GetGeometry()[i].GetValue(NODAL_AREA)));
        Hj[i](1,1)   = 2.0*coriolis_coefficient*(flow/A);
        //KRATOS_WATCH(Hj[i]);
        Suj[i](0,0) = 0.0;
        Suj[i](0,1) = 0.0;
        Suj[i](1,0) = kr_coefficient*flow/(A*A);;
        Suj[i](1,1) = -kr_coefficient/A;
        //KRATOS_WATCH(Suj[i]);
    }

    
    array_1d<double,2> Fder;

    Fder[0] = (Fj[1][0] - Fj[0][0]) / mL;
    Fder[1] = (Fj[1][1] - Fj[0][1]) / mL;
//    if (this->Id() == 500)
//           KRATOS_WATCH(Fder);

    boost::numeric::ublas::bounded_matrix<double, 2,2 > Fw;
    boost::numeric::ublas::bounded_matrix<double, 2,2 > Sw;

    array_1d<double,2> aaa;
    noalias(aaa) = prod(Hj[0] ,Sj[0]);
    Fw(0,0) = aaa[0]; Fw(1,0) = aaa[1];
 
    noalias(aaa) = prod(Hj[1] ,Sj[1]);
    Fw(0,1) = aaa[0]; Fw(1,1) = aaa[1];

    noalias(aaa) = prod(Suj[0] ,Sj[0]);
    Sw(0,0) = aaa[0]; Sw(1,0) = aaa[1];
 
    noalias(aaa) = prod(Suj[1] ,Sj[1]);
    Sw(0,1) = aaa[0]; Sw(1,1) = aaa[1];

  // KRATOS_WATCH(Fw);
  // KRATOS_WATCH(Sw);
    for (unsigned int j=0; j<2; j++)
    {
        for (unsigned int i=0; i<2; i++)
        {
          //Hemos guardado esto transpuesto (arriba esta al reves)
          Fw(i,j) = (h_int*0.5*Fw(i,j)) + Fj[j][i];
          Sw(i,j) = (h_int*0.5*Sw(i,j)) + Sj[j][i];
        }
    }

//KRATOS_WATCH(Fw);
//KRATOS_WATCH(Sw);
    boost::numeric::ublas::bounded_matrix<double, 2,2 > F2ord;
    boost::numeric::ublas::bounded_matrix<double, 2,2 > S2ord;

    noalias(aaa) = prod( Hj[0],Fder);
    F2ord(0,0) = aaa[0]; F2ord(1,0) = aaa[1]; 
    noalias(aaa) = prod( Hj[1],Fder);
    F2ord(0,1) = aaa[0]; F2ord(1,1) = aaa[1]; 
    
    noalias(aaa) = prod( Suj[0],Fder);
    S2ord(0,0) = aaa[0]; S2ord(1,0) = aaa[1]; 
    noalias(aaa) = prod( Suj[1],Fder);
    S2ord(0,1) = aaa[0]; S2ord(1,1) = aaa[1]; 
 //KRATOS_WATCH(M1);
 //KRATOS_WATCH(M2);
    array_1d<double,2> tmp;
    //array_1d<double,2> aux;
    for (unsigned int i=0; i<2; i++)
    {
      aaa[0] = Fw(i,0); aaa[1] = Fw(i,1);
      noalias(tmp) = prod(M1,aaa);
      Fw(i,0) = tmp[0]; Fw(i,1) = tmp[1];
 
      aaa[0] = Sw(i,0); aaa[1] = Sw(i,1);
      noalias(tmp) = prod(M2,aaa);
      Sw(i,0) = tmp[0]; Sw(i,1) = tmp[1];
      
      aaa[0] = F2ord(i,0); aaa[1] = F2ord(i,1);
      noalias(tmp) = prod(M1,aaa);
      F2ord(i,0) = tmp[0]; F2ord(i,1) = tmp[1];
      
      aaa[0] = S2ord(i,0); aaa[1] = S2ord(i,1);
      noalias(tmp) = prod(M2,aaa);
      S2ord(i,0) = tmp[0]; S2ord(i,1) = tmp[1];
      
    }
//     Fw = prod(M1,Fw);
//     Sw = prod(M2,Sw);
//     F2ord = prod(M1,F2ord);
//     S2ord = prod(M2,S2ord);

//    if(Id() == 1)
//    {
//    KRATOS_WATCH(Fw);
//KRATOS_WATCH(Sw);
// KRATOS_WATCH(F2ord);
//KRATOS_WATCH(S2ord);
//    }
    //now let's compute the rhs
    array_1d<double,2> rhs;
    
    for (unsigned int i=0; i<2; i++)
    {
       rhs[0] = (Fw(i,0) + Sw(i,0)) - ((h_int*0.5)*((F2ord(i,0)) + S2ord(i,0)));
       rhs[1] = (Fw(i,1) + Sw(i,1)) - ((h_int*0.5)*((F2ord(i,1)) + S2ord(i,1)));
       //KRATOS_WATCH(rhs);
       //KRATOS_WATCH(GetProperties().Id());
       //std::cout << "this element is" << this->Id() <<std::endl;
       rRightHandSideVector[i] = rhs[0];
       rRightHandSideVector[i+2] = rhs[1];         
    }
//    if(Id() == 1)
//    {
// KRATOS_WATCH(rRightHandSideVector);
//    }
// KRATOS_ERROR(std::logic_error,"i want it to go out here","");




    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ArteryElement::Initialize()
{
    KRATOS_TRY
    const double pi = 3.14159265;
    double radius = GetProperties()[RADIUS];
    array_1d<double,2> A0;
    array_1d<double,2> c0;
    const double r0 =  radius; //GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
    A0[0] = pi*r0*r0;

    const double r1 =  radius; //GetGeometry()[1].FastGetSolutionStepValue(RADIUS);
    A0[1] = pi*r1*r1;


    const double H0 = GetProperties()[THICKNESS];
    const double E = GetProperties()[YOUNG_MODULUS];
    const double nu = GetProperties()[POISSON_RATIO];
    const double pressure = GetProperties()[PRESSURE]; //Initial Pressure (Dyastolic pressure)
    const double blood_density = GetProperties()[DENSITY];
    
//    KRATOS_WATCH(GetProperties().Id());
//    KRATOS_WATCH(pressure);
//    KRATOS_WATCH(this->Id());

    double beta=E*H0*1.77245385/(1.0-nu*nu);

    //compute the lenght of the element
    array_1d<double,3> lvec = GetGeometry()[1].Coordinates();
    lvec -= GetGeometry()[0].Coordinates();

    mL = norm_2(lvec);

    c0[0]=sqrt(beta/(2*blood_density*A0[0]));
    c0[1]=sqrt(beta/(2*blood_density*A0[1]));

//    if(c0[0] == 0.00 || c0[1] == 0.00)
//     {
//       KRATOS_WATCH(c0[0]);
//       KRATOS_WATCH(c0[1]);
//       KRATOS_WATCH(GetProperties().Id());
//       KRATOS_WATCH(this->Id());
//       KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, in boundary 1-2 conditions used:son", "");
//     }

    //save area to the nodes. as well as its nodal mass
    GetGeometry()[0].SetLock();
    GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += 0.5*mL;
    GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) = A0[0];
    GetGeometry()[0].FastGetSolutionStepValue(RADIUS) = radius;
    GetGeometry()[0].FastGetSolutionStepValue(THICKNESS) = H0;
    GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS) = E;
    GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO) = nu;
    GetGeometry()[0].FastGetSolutionStepValue(BETA) = beta;
    GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) = pressure;
    GetGeometry()[0].FastGetSolutionStepValue(SYSTOLIC_PRESSURE) = pressure;
    GetGeometry()[0].FastGetSolutionStepValue(DYASTOLIC_PRESSURE) = pressure;
    GetGeometry()[0].FastGetSolutionStepValue(AVERAGE_PRESSURE) = pressure;
    GetGeometry()[0].FastGetSolutionStepValue(C0) = c0[0];
    GetGeometry()[0].GetValue(NODAL_AREA) = A0[0];//here we store the initial area
    GetGeometry()[0].GetValue(PRESSURE) = pressure;//here we store the initial area
    GetGeometry()[0].UnSetLock();
    //KRATOS_WATCH(GetGeometry()[0].Id())
    //KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS))
//    if(GetGeometry()[0].FastGetSolutionStepValue(C0) == 0.00)
//    {
//       KRATOS_WATCH(GetProperties().Id());
//       KRATOS_WATCH(blood_density);
//       KRATOS_WATCH(c0[0]);
//       KRATOS_WATCH(c0[1]);
//       KRATOS_WATCH(beta);
//       KRATOS_WATCH(A0[0]);
//       KRATOS_WATCH(A0[1]);
//        KRATOS_WATCH(GetProperties().Id());
//        KRATOS_WATCH(this->Id());
//        KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, in boundary 1-2 conditions used:son", "");
//    }

    GetGeometry()[1].SetLock();
    GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += 0.5*mL;
    GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA) = A0[1];
    GetGeometry()[1].FastGetSolutionStepValue(RADIUS) = radius;
    GetGeometry()[1].FastGetSolutionStepValue(THICKNESS) = H0;
    GetGeometry()[1].FastGetSolutionStepValue(YOUNG_MODULUS) = E;
    GetGeometry()[1].FastGetSolutionStepValue(POISSON_RATIO) = nu;
    GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) = pressure;
    GetGeometry()[1].FastGetSolutionStepValue(SYSTOLIC_PRESSURE) = pressure;
    GetGeometry()[1].FastGetSolutionStepValue(DYASTOLIC_PRESSURE) = pressure;
    GetGeometry()[1].FastGetSolutionStepValue(AVERAGE_PRESSURE) = pressure;
    GetGeometry()[1].FastGetSolutionStepValue(BETA) = beta;
    GetGeometry()[1].FastGetSolutionStepValue(C0) = c0[1];
    GetGeometry()[1].GetValue(NODAL_AREA) = A0[1]; //here we store the initial area
    GetGeometry()[1].GetValue(PRESSURE) = pressure;//here we store the initial area
    GetGeometry()[1].UnSetLock();
    //KRATOS_WATCH(GetGeometry()[1].Id())
    //KRATOS_WATCH(GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS))
            //KRATOS_WATCH("ddddd32d32")
//    if(GetGeometry()[1].FastGetSolutionStepValue(C0) == 0.00)
//     {
       //KRATOS_WATCH(c0[1]);
//       KRATOS_WATCH(GetProperties().Id());
//       KRATOS_WATCH(this->Id());
//       KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, in boundary 1-2 conditions used:son", "");
//     }

     //const double kk= GetGeometry()[0].FastGetSolutionStepValue(FLOW);
     //const double kkk=GetGeometry()[1].FastGetSolutionStepValue(FLOW);

     //KRATOS_WATCH(kk);
     //KRATOS_WATCH(kkk);
     KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void ArteryElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void ArteryElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit element", "");
}

//************************************************************************************
//************************************************************************************

void ArteryElement::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit element", "");
}


//************************************************************************************
//************************************************************************************
void ArteryElement::Calculate(const Variable<double >& rVariable,
                              double& Output,
                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //the variable error_ratio is here the norm of the subscale velocity as computed at the level of the gauss point
    if (rVariable == ERROR_RATIO)
    {

    }
    KRATOS_CATCH("");

}

int ArteryElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //check the area

    //check if if is in the XY plane

    //check that no variable has zero key


    return 0;


    KRATOS_CATCH("");
}


} // Namespace Kratos



