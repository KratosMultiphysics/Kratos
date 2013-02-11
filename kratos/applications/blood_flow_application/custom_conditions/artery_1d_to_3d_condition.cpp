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
#include "custom_conditions/artery_11_condition.h"
#include "utilities/math_utils.h"
#include "blood_flow_application.h"
#include "boost/numeric/ublas/lu.hpp"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

Artery1Dto3DCondition::Artery1Dto3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

Artery1Dto3DCondition::Artery1Dto3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer Artery1Dto3DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Condition::Pointer(new Artery1Dto3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

Artery1Dto3DCondition::~Artery1Dto3DCondition()
{
}

//************************************************************************************
//************************************************************************************

void Artery1Dto3DCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense to computer the system matrix for an explicit condition", "");
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void Artery1Dto3DCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

            //get data as needed
            const double dynamic_viscosity = GetProperties()[DYNAMIC_VISCOSITY];
            const double density = GetProperties()[DENSITY];
            const double E = GetProperties()[YOUNG_MODULUS];
            const double nu = GetProperties()[POISSON_RATIO];
	    //            const double pi = 3.14159265;
            const double coriolis_coefficient = 1.0001;
	    //            const double kr_coefficient = 1.0;
	    const double p3D = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);

	    //            const double kinematic_viscosity = dynamic_viscosity/density;
            const double beta = E*mH0*1.77245385/(1.0-nu*nu);


            //resize the vector to the correct size
            if (rRightHandSideVector.size() != 4)
                rRightHandSideVector.resize(4,false);

            array_1d<double,2> area;
            array_1d<double,2> flow;
            array_1d<double,2> wave_velocity;
	    array_1d<double,2> initial_wave_velocity;
            array_1d<double,2> artery_property;
            array_1d<double,2> coef;



            //loop on nodes
            for (unsigned int i=0; i<2; i++)
            {
                Node<3>& r_node = GetGeometry()[i];
                area[i] = r_node.FastGetSolutionStepValue(NODAL_AREA);
                flow[i] = r_node.FastGetSolutionStepValue(FLOW);
                artery_property[i]=mBeta/mInitialArea[i];
                coef[i] = sqrt(mBeta/(2*density * mInitialArea[i]));
            }
 	    flow[1] = flow[0]; //impose continuity of flow
            wave_velocity[0] = (flow[0] / area[0]) + 4.00 * coef[0] * pow(area[0],0.25);
            wave_velocity[1] = 1e20; //huge number, shall be infinity
            
            initial_wave_velocity[0] = (flow[0] / area[0]) + 4.00 * coef[0] * pow(mInitialArea[0],0.25);
            initial_wave_velocity[1] = 1e20; //huge number, shall be infinity

            //KRATOS_WATCH(wave_velocity);


            double convergence;
            unsigned int max_iterations = 100;
            double tolerance = 1e-6;
	    
	    Vector rhs(2,0.0);
	    Matrix lhs(2,2,0.0);

            for(unsigned int i = 0 ; i < max_iterations ; i++)
            {
		//construct rhs
		rhs[0] = initial_wave_velocity[0] - (flow[0] / area[0]) + 4.00 * coef[0] * pow(area[0],0.25);
		rhs[1] = mBeta*(sqrt(area[0]) - sqrt(mInitialArea[0]))/mInitialArea[0] - p3D -(1.0/2.0)*density*flow[0]*flow[0]*(1.0/pow(mInitialArea[1],2)-1.0/pow(area[0],2) );

		//construct tangent
		lhs(0,0) = -1.0/area[0];
		lhs(0,1) = flow[0]/pow(area[0],2) - sqrt(mBeta/(density*mInitialArea[0])) / (sqrt(2.0)*pow(area[0],0.75));
 		lhs(1,0) = -density* (1.0/pow(mInitialArea[1],2) - 1.0 /pow(area[0],2)) *flow[0];
		lhs(1,1) = mBeta/(2.0*sqrt(area[0]*mInitialArea[0])) - density*flow[0]*flow[0]/pow(area[0],3) ;
		
                permutation_matrix<double> permutation(2);
                array_1d<double,2> delta_x = rhs;
				
                lu_factorize(lhs, permutation);
                lu_substitute(lhs,permutation, delta_x);
		
                convergence = norm_2(delta_x);

		flow[0] += delta_x[0];
		area[0] += delta_x[1];

		// we have to add the relative convergence check
                if(convergence < tolerance)
                    break;
            }
            
            //next is needed as FLOW has to be fixed for 3D inlet
            GetGeometry()[1].FastGetSolutionStepValue(FLOW) = flow[0];
	    
	    
            
            double A1 = area[0];
            double A0 = mInitialArea[0];
            double C = beta*sqrt(A1*A1*A1)/(3.0*density*A0);

            rRightHandSideVector[0] = -flow[0];
            //double temp = (C + coriolis_coefficient*flow[0]*flow[0]/(A1));
            rRightHandSideVector[1] = -(C + coriolis_coefficient*flow[0]*flow[0]/(A1));

            A1 = area[1];
	    A0 = mInitialArea[1];
            C = beta*sqrt(A1*A1*A1)/(3.0*density*A0);

            rRightHandSideVector[2] = flow[1];
            //temp = (C + coriolis_coefficient*flow[1]*flow[1]/(A1));
            rRightHandSideVector[3] = (C + coriolis_coefficient*flow[1]*flow[1]/(A1));



    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void Artery1Dto3DCondition::Initialize()
{
    KRATOS_TRY


    const double pi = 3.14159265;
    double radius = GetProperties()[RADIUS];

    const double r0 =  radius; //GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
    mInitialArea[0] = pi*r0*r0;

    const double r1 =  radius; //GetGeometry()[1].FastGetSolutionStepValue(RADIUS);
    mInitialArea[1] = pi*r1*r1;


    mH0 = GetProperties()[THICKNESS];

    const double E = GetProperties()[YOUNG_MODULUS];
    const double nu = GetProperties()[POISSON_RATIO];

    mBeta = E*mH0*1.77245385/(1.0-nu*nu);


    //save area to the nodes. as well as its nodal mass
    GetGeometry()[0].SetLock();
    GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) = mInitialArea[0];
    GetGeometry()[0].FastGetSolutionStepValue(RADIUS) = radius;
    GetGeometry()[0].UnSetLock();
    GetGeometry()[1].SetLock();
    GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA) = mInitialArea[1];
    GetGeometry()[1].FastGetSolutionStepValue(RADIUS) = radius;
    GetGeometry()[1].UnSetLock();




    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void Artery1Dto3DCondition::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void Artery1Dto3DCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}

//************************************************************************************
//************************************************************************************

void Artery1Dto3DCondition::GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}


//************************************************************************************
//************************************************************************************
void Artery1Dto3DCondition::Calculate(const Variable<double >& rVariable,
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

int Artery1Dto3DCondition::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //check the area

    //check if if is in the XY plane

    //check that no variable has zero key


    return 0;


    KRATOS_CATCH("");
}

void Artery1Dto3DCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY



    KRATOS_CATCH("");
}

} // Namespace Kratos



