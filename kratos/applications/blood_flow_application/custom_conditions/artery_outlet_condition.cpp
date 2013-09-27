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
#include <iostream>

namespace Kratos
{

//************************************************************************************
//************************************************************************************

ArteryOutletCondition::ArteryOutletCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

ArteryOutletCondition::ArteryOutletCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer ArteryOutletCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Condition::Pointer(new ArteryOutletCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

ArteryOutletCondition::~ArteryOutletCondition()
{
}

//************************************************************************************
//************************************************************************************

void ArteryOutletCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense to computer the system matrix for an explicit condition", "");
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ArteryOutletCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //double total_time=rCurrentProcessInfo[TIME];
    //int nodo= GetProperties().Id();

    //KRATOS_WATCH(total_time)

    //resize the vector to the correct size
    if (rRightHandSideVector.size() != 2)
        rRightHandSideVector.resize(2,false);

    //double h_int = rCurrentProcessInfo[DELTA_TIME];

    //const int Propiedad=GetProperties().Id();

    //get data as needed
    // EDU:: TOMA LAS PROPIEDADES DE LA PROPIEDAD XXX DEL MPDA
    //const double dynamic_viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    const double density = GetProperties()[DENSITY];
    const double p_venous= GetProperties()[PRESSURE_VENOUS];
    const double p_init = GetGeometry()[0].FastGetSolutionStepValue(SYSTOLIC_PRESSURE);


    //GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
    //KRATOS_WATCH(p_init)
    //const double initial_area = GetGeometry()[0].GetValue(NODAL_AREA);
    const double A0 = GetGeometry()[0].GetValue(NODAL_AREA);
    //
    //const double E = GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS);
    //const double nu = GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO);
    //const double pi = 3.14159265;
    const double coriolis_coefficient = 1.1;
    //const double kr_coefficient = 1.0;

    const double terminal_resistence = GetProperties()[TERMINAL_RESISTANCE];
    
    //const double kinematic_viscosity = dynamic_viscosity/density;
    //const double H0 = GetGeometry()[0].FastGetSolutionStepValue(THICKNESS);
    //const double beta = E*H0*1.77245385/(1.0-nu*nu);
    const double beta = GetGeometry()[0].FastGetSolutionStepValue(BETA);
    const double par1 = beta / A0;

    //KRATOS_WATCH(Propiedad)
    //const double A1 = GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);    
	const double A1 = GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
    //KRATOS_WATCH(A1)
	const double A = UpdateArea(beta, A1);

    //const double flow1 =  GetGeometry()[0].FastGetSolutionStepValue(FLOW); // NOTE: HERE we have to put the corrected value
    //const double flow2 = GetGeometry()[0].FastGetSolutionStepValue(FLOW);
    
    const double flow1 =  (par1*(sqrt(A)-sqrt(A0))+p_init-p_venous)*(1/terminal_resistence); // NOTE: HERE we have to put the corrected value
    const double flow2 =  GetGeometry()[0].FastGetSolutionStepValue(FLOW);
    
    const double flow = 2* flow1 - flow2;
    //const double flow = flow1;

    //KRATOS_WATCH(par1)
    //KRATOS_WATCH(A)

    //KRATOS_WATCH(flow1)
    //KRATOS_WATCH(flow2)

    double C = (beta)/(3.0*density*A0)*sqrt(A*A*A);
    //KRATOS_WATCH(C)

    //C = beta/(3.0*density*GetGeometry()[0].GetValue(NODAL_AREA))*sqrt(A*A*A);
    //KRATOS_WATCH(C)¡'

    //std::cout << "outlet: " << std::endl;
    //KRATOS_WATCH(A0);
    //KRATOS_WATCH(A1);
    //KRATOS_WATCH(C);
    //KRATOS_WATCH(A);
    
    rRightHandSideVector[0] = - flow;
    //double temp = -(C + coriolis_coefficient*flow*flow/(A1));
    //EDU::CREO QUE LA EXPRESION ES A en vez de A1:
    //rRightHandSideVector[1] = -(C + (coriolis_coefficient*flow*flow/(A1)));
    rRightHandSideVector[1] = -(C + ((coriolis_coefficient*flow*flow)/A));


    KRATOS_CATCH("")
}

double ArteryOutletCondition::UpdateArea(double Beta, double A)
{

    KRATOS_TRY

    const int max_iteration = 100;
    const double p_init = GetGeometry()[0].FastGetSolutionStepValue(SYSTOLIC_PRESSURE);
    //const double A0=GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
    const double flow =  GetGeometry()[0].FastGetSolutionStepValue(FLOW);
    const double initial_area = GetGeometry()[0].GetValue(NODAL_AREA);
    const double par1 = Beta / initial_area;
    //const double par2 = sqrt(Beta / (2.00*Density*initial_area));
    const double c0 = GetGeometry()[0].FastGetSolutionStepValue(C0);
    const double terminal_resistence = GetProperties()[TERMINAL_RESISTANCE];
    const double p_venous=GetProperties()[PRESSURE_VENOUS];
    const double w1 = (flow / A) + (4.00*c0*pow(A,0.25));
//KRATOS_WATCH(flow)
//KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(FLOW,1))
//KRATOS_WATCH(A)
//KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA,1))
//KRATOS_WATCH(c0)
    //const double p_init = 10640.00;
    //std::cout << "ENTRDA" << std::endl;
    double x = A;
    for(int i = 0 ; i < max_iteration ; i++)
    {
//KRATOS_WATCH(x)
//KRATOS_WATCH(GetGeometry()[0].Id())
//KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA))
//KRATOS_WATCH(terminal_resistence)
//KRATOS_WATCH(w1)
//KRATOS_WATCH(p_init)
//KRATOS_WATCH(initial_area)
//KRATOS_WATCH(par1)
//KRATOS_WATCH((-4.00 * c0 * pow(x, 1.25) * terminal_resistence))
//KRATOS_WATCH((w1 * x * terminal_resistence))
//KRATOS_WATCH((par1 * (sqrt(x) - sqrt(initial_area) )))
        double f = (-4.00 * c0 * pow(x, 1.25) * terminal_resistence) + (w1 * x * terminal_resistence) - (par1 * (sqrt(x) - sqrt(initial_area))) - p_init + p_venous;
        double df= (-5.00 * c0 * pow(x, 0.25) * terminal_resistence) + (w1 * terminal_resistence) - (par1 * 0.5*pow(x, -0.5));
//KRATOS_WATCH(f)
//KRATOS_WATCH(df)
        double dx = f/df;
        x-= dx;
        if(fabs(dx) < 1e-6)
            A = x;
            break;
//        else
//            std::cout << "NO CONVERGEEEEEEEEEEEEEEEE::: Artery_outlet" << std::endl;
//            KRATOS_WATCH(x);
//            KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA));
//            KRATOS_ERROR(std::runtime_error, "Artery_outlet", "");
    }
    A = x;
    //KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA));
    //KRATOS_WATCH(A);

//    if(x != GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA))
//    {
//        std::cout << "ERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR" << std::endl;
//        KRATOS_WATCH(x);
//        KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA));
//        KRATOS_ERROR(std::logic_error,"artery_outlet_conditions", "");
//    }

    if(A == 0.00)
    {
        KRATOS_WATCH(A);
        KRATOS_WATCH(GetProperties().Id());
        KRATOS_WATCH(this->Id());
        KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, Please check your inlet boundary conditions used", "");
    }
    //std::cout << "SALIDA" << std::endl;
    //KRATOS_WATCH(A);

    return A;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ArteryOutletCondition::Initialize()
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void ArteryOutletCondition::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void ArteryOutletCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}

//************************************************************************************
//************************************************************************************

void ArteryOutletCondition::GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}


//************************************************************************************
//************************************************************************************
void ArteryOutletCondition::Calculate(const Variable<double >& rVariable,
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

int ArteryOutletCondition::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //check the area

    //check if if is in the XY plane

    //check that no variable has zero key


    return 0;


    KRATOS_CATCH("");
}



} // Namespace Kratos



