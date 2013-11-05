/* b
==============================================================================
KratosIncompressibleFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Condition Analysis
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
#include "custom_conditions/artery_outlet_condition.h"
#include "utilities/math_utils.h"
#include "blood_flow_application.h"
#include "boost/numeric/ublas/lu.hpp"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

ArteryOutletFreeCondition::ArteryOutletFreeCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

ArteryOutletFreeCondition::ArteryOutletFreeCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer ArteryOutletFreeCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Condition::Pointer(new ArteryOutletFreeCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

ArteryOutletFreeCondition::~ArteryOutletFreeCondition()
{
}

//************************************************************************************
//************************************************************************************

void ArteryOutletFreeCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense to computer the system matrix for an explicit condition", "");
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ArteryOutletFreeCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //resize the vector to the correct size
    if (rRightHandSideVector.size() != 4)
        rRightHandSideVector.resize(4,false);

    const int Tipo_outlet = 2;
    //double h_int = rCurrentProcessInfo[DELTA_TIME];

    //get data as needed
    //const double dynamic_viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    const double density = GetProperties()[DENSITY];
    const double E = GetGeometry()[1].FastGetSolutionStepValue(YOUNG_MODULUS);
    const double nu = GetGeometry()[1].FastGetSolutionStepValue(POISSON_RATIO);
    //const double pi = 3.14159265;
    const double coriolis_coefficient = 1.1;
    //const double kr_coefficient = 1.0;

    //const double kinematic_viscosity = dynamic_viscosity/density;
    const double H0 = GetGeometry()[0].FastGetSolutionStepValue(THICKNESS);
    const double beta = E*H0*1.77245385/(1.0-nu*nu);

    const double A1 = GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA);
    const double& A = A1;//UpdateArea(beta, density);

    // Fortran::
    //An=U_n(1,jj) --> equivaldría a A1
    //Qn=U_n(2,jj) --> equivaldría a GetGeometry()[0].FastGetSolutionStepValue(FLOW);

    if(Tipo_outlet == 1) {
        //const double flow1 = GetGeometry()[0].FastGetSolutionStepValue(FLOW); // NOTE: HERE we have to put the corrected value
        const double flow1 = GetGeometry()[1].FastGetSolutionStepValue(FLOW); // Edu
        const double flow2 = GetGeometry()[0].FastGetSolutionStepValue(FLOW);
        const double flow = 2* flow2 - flow1;
        double A0 = GetGeometry()[0].GetValue(NODAL_AREA);
        const double C = beta*sqrt(A*A*A)/(3.0*density*A0);
        rRightHandSideVector[0] = 0; //-flow;
        //double temp = -(C + coriolis_coefficient*flow*flow/(A1));
        rRightHandSideVector[1] = 0; //-(C + coriolis_coefficient*flow*flow/(A1));
        rRightHandSideVector[2] = -flow;
        rRightHandSideVector[3] = -(C + coriolis_coefficient*flow*flow/(A1));
        }
    else if (Tipo_outlet == 2){
        KRATOS_WATCH("FREE_OUTLET");
        // Fortran::
        //W1=(Qn/An)+(4.d0*c0(k)*(An**0.25d0))
        //An_1=rtnewt('out',An,Qn,An,par1,c0(k),A0(k),Rt,P_OUT,TipoOutlet)
        //Qn_1=(par1*(dsqrt(An_1)-dsqrt(A0(k)))+P_init-P_OUT)*(1/Rt)
        //const double flow_n =  GetGeometry()[0].FastGetSolutionStepValue(FLOW);
        //const double A_n = GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
        //const double C = beta*sqrt(A0*A0*A0)/(3.0*density*A0);
        //const double Wave_Speed1 =(flow_n/A_n)+4*C;
        const double flow1 = GetGeometry()[1].FastGetSolutionStepValue(FLOW); // Edu
        const double flow2 = GetGeometry()[0].FastGetSolutionStepValue(FLOW);
        const double flow = 2* flow1 - flow2;

        double A0 = GetGeometry()[1].GetValue(NODAL_AREA);
        const double Area_1 = GetGeometry()[1].GetValue(NODAL_AREA);
        const double Area_2 = GetGeometry()[0].GetValue(NODAL_AREA);
        //An_1=(2*(An_1**0.25d0)-(An**0.25d0))**4
	double Area_Aux = pow((2*(pow(Area_1,0.25))-pow(Area_2,0.25)),4.0);
	KRATOS_WATCH(Area_Aux);
	KRATOS_WATCH(A1);
	KRATOS_WATCH(A);
	KRATOS_WATCH(Area_2);
	KRATOS_WATCH(Area_1);
	Area_Aux=Area_1;
        const double C = beta*sqrt(Area_Aux*Area_Aux*Area_Aux)/(3.0*density*A0);

        rRightHandSideVector[0] = 0; //-flow;
        //double temp = -(C + coriolis_coefficient*flow*flow/(Area_Aux));
        rRightHandSideVector[1] = 0; //-(C + coriolis_coefficient*flow*flow/(A1));
        rRightHandSideVector[2] = -flow;
        rRightHandSideVector[3] = -(C + coriolis_coefficient*flow*flow/(Area_Aux));
	
	
	//Fr(1)=Qb(aa)
	//C_boun=beta0(k)/(3.d0*blood_dens*A0(k))*(Ab(aa))**1.5d0
	//Fr(2)=(Coriolis*(Qb(aa)**2.d0)/Ab(aa))+C_boun
	//rhs_glob(1,jj(2))=rhs_glob(1,jj(2))-(Fr(1))/node_coord(jj(2))%mass !-h_int/2.d0*Fu_boun(1))
	//rhs_glob(2,jj(2))=rhs_glob(2,jj(2))-(Fr(2))/node_coord(jj(2))%mass !-h_int/2.d0*Fu_boun(2))
	
        }
    else if (Tipo_outlet == 3){

        }
    else {
    }



    KRATOS_CATCH("")
}

double ArteryOutletFreeCondition::UpdateArea(double Beta, double Density)
{
    KRATOS_TRY

    const int max_iteration = 10;
//    const double p_init = GetProperties()[PRESSURE];
    const double p_init = GetGeometry()[1].FastGetSolutionStepValue(DYASTOLIC_PRESSURE);
      double& A = GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA);
    const double flow =  GetGeometry()[1].FastGetSolutionStepValue(FLOW);
    double initial_area = GetGeometry()[1].GetValue(NODAL_AREA);
    const double par1 = Beta / initial_area;
    const double par2 = sqrt(Beta / (2.00*Density*initial_area));
    const double terminal_resistence = GetGeometry()[1].FastGetSolutionStepValue(TERMINAL_RESISTANCE);
    const double w1 = flow / A + 4.00*par2*pow(A,0.25);
    //const double p_init = 10640.00;
    KRATOS_WATCH("REVISAR PRESSSIIIIIIIIIIIIIIIIIONESSSSSSSSSS");
    double x = A;
    for(int i = 0 ; i < max_iteration ; i++)
    {
        double f = (-4.00 * par2 * pow(x, 1.25) * terminal_resistence) + (w1 * x * terminal_resistence) - (par1 * (sqrt(x) - sqrt(initial_area))) - p_init;
        double df= (-5.00 * par2 * pow(x, 0.25) * terminal_resistence) + (w1 * terminal_resistence) - (par1 * 0.5 / sqrt(x));

        double dx = f/df;
        x-= dx;
        if(fabs(dx) < 1e-6)
            break;
    }
    A = x;
    if(x != GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA))
    {
        std::cout << "ERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR" << std::endl;
        KRATOS_WATCH(x);
        KRATOS_WATCH(GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA));
        KRATOS_ERROR(std::logic_error,"artery_outlet_conditions", "");
    }

    return A;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ArteryOutletFreeCondition::Initialize()
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void ArteryOutletFreeCondition::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void ArteryOutletFreeCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}

//************************************************************************************
//************************************************************************************

void ArteryOutletFreeCondition::GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}


//************************************************************************************
//************************************************************************************
void ArteryOutletFreeCondition::Calculate(const Variable<double >& rVariable,
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

int ArteryOutletFreeCondition::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //check the area

    //check if if is in the XY plane

    //check that no variable has zero key


    return 0;


    KRATOS_CATCH("");
}



} // Namespace Kratos



