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
//   Last modified by:    $Author: Eduardo Soudah $
//   Date:                $Date: 2009-01-13 15:39:56 $
//   Revision:            $Revision: 1.14 $
//
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/artery_inlet_condition.h"
#include "utilities/math_utils.h"
#include "blood_flow_application.h"
#include "boost/numeric/ublas/lu.hpp"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

ArteryInletConditionPressure::ArteryInletConditionPressure(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

ArteryInletConditionPressure::ArteryInletConditionPressure(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer ArteryInletConditionPressure::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Condition::Pointer(new ArteryInletConditionPressure(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

ArteryInletConditionPressure::~ArteryInletConditionPressure()
{
}

//************************************************************************************
//************************************************************************************

void ArteryInletConditionPressure::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense to computer the system matrix for an explicit condition", "");
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ArteryInletConditionPressure::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //resize the vector to the correct size
    if (rRightHandSideVector.size() != 2)
        rRightHandSideVector.resize(2,false);


    std::cout << "PRESSURE ::::::: Artery_Inlet _ NUEVO" << std::endl;

    //double h_int = rCurrentProcessInfo[DELTA_TIME];
    // get data as needed
    // EDU:: ESTA TOMANDO LAS PROPIEDADES DEL ELEMENTO 1 (COMO ESTA DEFINIDO EN EL MPDA)
    // const double dynamic_viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    // EDU:: AQUI YA TOMA LOS VALORES CON EL UPDATE PARA LOS SIGUIENTE PASOS DE TIEMPO
    //const double E = GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS);
    //const double nu = GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO);
    //const double pi = 3.14159265;
    const double coriolis_coefficient = 1.1;
    //const double kr_coefficient = 1.0;
    const double density = GetProperties()[DENSITY];
    //const double kinematic_viscosity = dynamic_viscosity/density;
    //const double H0 = GetGeometry()[0].FastGetSolutionStepValue(THICKNESS);;
    //const double beta = E*H0*1.77245385/(1.0-nu*nu);
    const double beta =GetGeometry()[0].FastGetSolutionStepValue(BETA);
    //const double A1 = GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
    const double A0 = GetGeometry()[0].GetValue(NODAL_AREA);
    KRATOS_WATCH("A0_kratos");
    KRATOS_WATCH(A0);
    const double A = GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
    KRATOS_WATCH("A_kratos");
    KRATOS_WATCH(A);
    const double Q1 = GetGeometry()[0].FastGetSolutionStepValue(FLOW);
    KRATOS_WATCH("Q1_kratos");
    KRATOS_WATCH(Q1);
    const double flow = UpdateFlow(beta,Q1);
    KRATOS_WATCH("Flow_kratos");
    KRATOS_WATCH(flow);
    //const double A = A1; // No hago update del area
    //const double& flow = GetGeometry()[0].FastGetSolutionStepValue(FLOW);
    //flow = 2*flow -
    //const double flow3 = GetGeometry()[0].FastGetSolutionStepValue(FLOW);
    //const double& flow5 = GetGeometry()[0].GetSolutionStepValue(FLOW);
    const double C = (beta*sqrt(A*A*A)/(3.0*density*A0));
    //const double flow2 = GetGeometry()[0].GetValue(FLOW,0);
    //const double flow4 = GetGeometry()[0].GetValue(FLOW,1);
    //const int kk =GetProperties().Id();
    //KRATOS_WATCH(GetGeometry()[0].Id())

    //std::cout << "inlet: " << std::endl;
    //KRATOS_WATCH(flow);
    //KRATOS_WATCH(flow2);
    //KRATOS_WATCH(flow3);
    //KRATOS_WATCH(flow4);
    //KRATOS_WATCH(flow5);
    //KRATOS_WATCH(A);
    rRightHandSideVector[0] = flow;
    //double temp = C + coriolis_coefficient*flow*flow/(A);
    rRightHandSideVector[1] = C + (coriolis_coefficient*flow*flow/(A));
    //TODO:CHEQUEAR EL TERMINO DE LA MASA y el SIGNO
    // RHS = RHS + FR/MASA

    KRATOS_CATCH("")
}


double ArteryInletConditionPressure::UpdateFlow(double Beta, double Q)
{

    KRATOS_TRY

    const int max_iteration = 100;
    const double flow =  GetGeometry()[0].FastGetSolutionStepValue(FLOW);
    const double A=GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
    const double A0=GetGeometry()[0].GetValue(NODAL_AREA);
    const double density = GetProperties()[DENSITY];
    const double par2 = sqrt(Beta / (2.00*density*A0));
    const double w1 = (flow / A) - 4.00*par2*pow(A,0.25);

    //std::cout << "FLOWWWWWWW ::::::: Artery_Inlet" << std::endl;
    //KRATOS_WATCH(flow);

    double x = Q;
    for(int i = 0 ; i < max_iteration ; i++)
    {
        double f = pow(((w1 - (x/A0))/(4*par2)),4) - A0;
                //4.00 * par2 * pow(x, 1.25) + w1 * x - flow;
        double df=-((pow(((w1-x/A0)/4),3))/(A0*pow(par2,4)));
                //5.00 * par2 * pow(x, 0.25) + w1 ;
        double dx = f/df;
        x -= dx;
        if(fabs(dx) < 1e-6)
            Q = x;
            KRATOS_WATCH("NEWTON_R");
            KRATOS_WATCH(Q);
            break;
//        else
//            //std::cout << "NO CONVERGEEEEEEEEEEEEEEEE:: Artery_Inlet" << std::endl;
//            KRATOS_WATCH(x);
//            KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA));
//            std::cout << "NO CONVERGEEEEEEEEEEEEEEEE::: Artery_Inlet" << std::endl;
//            KRATOS_ERROR(std::runtime_error, "Artery_Inlet", "");
    }
    Q = x;
    return Q;
    KRATOS_CATCH("");
}


//************************************************************************************
//************************************************************************************
void ArteryInletConditionPressure::Initialize()
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void ArteryInletConditionPressure::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void ArteryInletConditionPressure::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}

//************************************************************************************
//************************************************************************************

void ArteryInletConditionPressure::GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}


//************************************************************************************
//************************************************************************************
void ArteryInletConditionPressure::Calculate(const Variable<double >& rVariable,
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

int ArteryInletConditionPressure::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //check the area

    //check if if is in the XY plane

    //check that no variable has zero key


    return 0;


    KRATOS_CATCH("");
}

} // Namespace Kratos



