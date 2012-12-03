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
#include "custom_conditions/artery_inlet_condition.h"
#include "utilities/math_utils.h"
#include "blood_flow_application.h"
#include "boost/numeric/ublas/lu.hpp"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

ArteryInletCondition::ArteryInletCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

ArteryInletCondition::ArteryInletCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer ArteryInletCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Condition::Pointer(new ArteryInletCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

ArteryInletCondition::~ArteryInletCondition()
{
}

//************************************************************************************
//************************************************************************************

void ArteryInletCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense to computer the system matrix for an explicit condition", "");
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ArteryInletCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //resize the vector to the correct size
    if (rRightHandSideVector.size() != 2)
        rRightHandSideVector.resize(2,false);


    //double h_int = rCurrentProcessInfo[DELTA_TIME];


    //get data as needed
    const double dynamic_viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    const double density = GetProperties()[DENSITY];
    const double E = GetProperties()[YOUNG_MODULUS];
    const double nu = GetProperties()[POISSON_RATIO];
    //const double pi = 3.14159265;
    const double coriolis_coefficient = 1.0001;
    //const double kr_coefficient = 1.0;

    //const double kinematic_viscosity = dynamic_viscosity/density;
    const double beta = E*mH0*1.77245385/(1.0-nu*nu);

    const double& A = UpdateArea(beta, density);

    const double& flow = GetGeometry()[0].FastGetSolutionStepValue(FLOW);
    const double C = beta*sqrt(A*A*A)/(3.0*density*mInitialArea[0]);
//    std::cout << "inlet: " << std::endl;
//    KRATOS_WATCH(flow);
//    KRATOS_WATCH(A);

    rRightHandSideVector[0] = flow;
    double temp = C + coriolis_coefficient*flow*flow/(A);
    rRightHandSideVector[1] = C + coriolis_coefficient*flow*flow/(A);


    KRATOS_CATCH("")
}


double ArteryInletCondition::UpdateArea(double Beta, double Density)
{
    const int max_iteration = 10;
    double& A = GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
    const double flow =  GetGeometry()[0].FastGetSolutionStepValue(FLOW);
    const double par2 = sqrt(Beta / (2.00*Density*mInitialArea[0]));
    const double w1 = flow / A - 4.00*par2*pow(A,0.25);

    double x = A;
    for(int i = 0 ; i < max_iteration ; i++)
    {
        double f = 4.00 * par2 * pow(x, 1.25) + w1 * x - flow;
        double df= 5.00 * par2 * pow(x, 0.25) + w1 ;

        double dx = f/df;
        x -= dx;
        if(fabs(dx) < 1e-10)
            break;
    }

    A = x;

    return A;
}


//************************************************************************************
//************************************************************************************
void ArteryInletCondition::Initialize()
{
    KRATOS_TRY


    const double pi = 3.14159265;
    double radius = GetProperties()[RADIUS];

    const double r0 =  radius; //GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
    mInitialArea[0] = pi*r0*r0;

    mH0 = GetProperties()[THICKNESS];

    const double E = GetProperties()[YOUNG_MODULUS];
    const double nu = GetProperties()[POISSON_RATIO];

    mBeta = E*mH0*1.77245385/(1.0-nu*nu);


    //save area to the nodes. as well as its nodal mass
    GetGeometry()[0].SetLock();
    GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) = mInitialArea[0];
    GetGeometry()[0].FastGetSolutionStepValue(RADIUS) = radius;
    GetGeometry()[0].UnSetLock();

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void ArteryInletCondition::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void ArteryInletCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}

//************************************************************************************
//************************************************************************************

void ArteryInletCondition::GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}


//************************************************************************************
//************************************************************************************
void ArteryInletCondition::Calculate(const Variable<double >& rVariable,
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

int ArteryInletCondition::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //check the area

    //check if if is in the XY plane

    //check that no variable has zero key


    return 0;


    KRATOS_CATCH("");
}

} // Namespace Kratos



