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
#include <iostream>

namespace Kratos
{

//************************************************************************************
//************************************************************************************

Artery11Condition::Artery11Condition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

Artery11Condition::Artery11Condition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
{
   KRATOS_TRY
   
//    const double pi = 3.14159265;
//    //double radius = GetProperties()[RADIUS];

//    const double r0 =  GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
//    //const double r0 = GetGeometry()[0].GetSolutionStepValue(RADIUS);
//    mInitialArea[0] = pi*r0*r0;

//    //KRATOS_WATCH(mInitialArea[0]);
//    //KRATOS_WATCH(r0);
//    const double r1 =  GetGeometry()[1].FastGetSolutionStepValue(RADIUS);
//    //const double r1 =  GetGeometry()[1].GetSolutionStepValue(RADIUS);
//    mInitialArea[1] = pi*r1*r1;

//    std::cout<<"hola";

//    //KRATOS_WATCH(mInitialArea[0]);
    
//    //const double mH0 = GetProperties()[THICKNESS];
//    //const double E = GetProperties()[YOUNG_MODULUS];
//    //const double nu = GetProperties()[POISSON_RATIO];

//    //save area to the nodes. as well as its nodal mass
//    GetGeometry()[0].SetLock();
//    GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) = mInitialArea[0];
//    //GetGeometry()[0].FastGetSolutionStepValue(RADIUS) = radius;
//    GetGeometry()[0].FastGetSolutionStepValue(RADIUS) = r0;
//    GetGeometry()[0].UnSetLock();
    
//    GetGeometry()[1].SetLock();
//    GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA) = mInitialArea[1];
//    //GetGeometry()[1].FastGetSolutionStepValue(RADIUS) = radius;
//    GetGeometry()[1].FastGetSolutionStepValue(RADIUS) = r1;
//    GetGeometry()[1].UnSetLock();
    
    KRATOS_CATCH("");
}

Condition::Pointer Artery11Condition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Condition::Pointer(new Artery11Condition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

Artery11Condition::~Artery11Condition()
{
}

//************************************************************************************
//************************************************************************************

void Artery11Condition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense to computer the system matrix for an explicit condition", "");
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void Artery11Condition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

            //get data as needed
            // const double dynamic_viscosity = GetProperties()[DYNAMIC_VISCOSITY];
            const double density = GetProperties()[DENSITY];
            //const double density = GetSolutionStepValue(DENSITY);
            //const double pi = 3.14159265;
            const double coriolis_coefficient = 1.1;
 
            //resize the vector to the correct size
            if (rRightHandSideVector.size() != 4)
                rRightHandSideVector.resize(4,false);

            using namespace boost::numeric::ublas;

            array_1d<double,4> f_out;
            Matrix jacobian = ZeroMatrix(4,4);

            array_1d<double,2> area;
            array_1d<double,2> flow;
            array_1d<double,2> wave_velocity;
            array_1d<double,2> artery_property;
            array_1d<double,2> coef;
            array_1d<double,2> beta;
	    
            //loop on nodes
            for (unsigned int i=0; i<2; i++)
            {
                //std::cout<<"hola";
                //KRATOS_WATCH(GetGeometry()[i].Id())
                //const double density = GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
//                const double E = GetGeometry()[i].FastGetSolutionStepValue(YOUNG_MODULUS);
//                const double nu =GetGeometry()[i].FastGetSolutionStepValue(POISSON_RATIO);
//                const double H0 =GetGeometry()[i].FastGetSolutionStepValue(THICKNESS);
//                beta[i] = E*H0*1.77245385/(1.0-nu*nu);

                Node<3>& r_node = GetGeometry()[i];
                area[i] = r_node.FastGetSolutionStepValue(NODAL_AREA);
                flow[i] = r_node.FastGetSolutionStepValue(FLOW);
                coef[i] = r_node.FastGetSolutionStepValue(C0);
                beta[i] = r_node.FastGetSolutionStepValue(BETA);
                artery_property[i]=beta[i]/area[i];
                //coef[i] = sqrt(beta[i]/(2*density * GetGeometry()[i].GetValue(NODAL_AREA)));
                //KRATOS_WATCH(r_node.GetId())
                if(coef[i] == 0.00)
                 {
                   KRATOS_WATCH(GetProperties().Id());
                   KRATOS_WATCH(this->Id());
                   KRATOS_WATCH(GetGeometry()[i].Id());
                   KRATOS_WATCH(r_node.FastGetSolutionStepValue(C0))
                   //KRATOS_WATCH(r_node.FastGetSolutionStepValue(c0))
                   KRATOS_ERROR(std::runtime_error, "Zero Coef found, in boundary 1-1 conditions used", "");
                 }
            }



            //KRATOS_WATCH(GetGeometry().Info());
            //flow[1] = flow[0];
            //wave_velocity[0] = (flow[0] / area[0]) + 4.00 * sqrt(beta[0] / (2.00 * density * GetGeometry()[0].GetValue(NODAL_AREA))) * pow(area[0],0.25);
            //wave_velocity[1] = (flow[1] / area[1]) - 4.00 * sqrt(beta[1] / (2.00 * density * GetGeometry()[1].GetValue(NODAL_AREA))) * pow(area[1],0.25);
	    wave_velocity[0] = (flow[0] / area[0]) + 4.00 * coef(0) * pow(area[0],0.25);
            wave_velocity[1] = (flow[1] / area[1]) - 4.00 * coef(1) * pow(area[1],0.25);

//KRATOS_WATCH(flow);
//KRATOS_WATCH(wave_velocity);


            double convergence;
            unsigned int max_iterations = 1000;
            double tolerance = 1e-9;

            for(unsigned int i = 0 ; i < max_iterations ; i++)
            {
                CalculateFunctional4(f_out, area, flow, artery_property, coef, wave_velocity, density);
                CalculateJacobian4(jacobian, area, flow, artery_property, coef, wave_velocity, density);
                //KRATOS_WATCH(f_out);
                //KRATOS_WATCH(jacobian);

                permutation_matrix<double> permutation(4);
                array_1d<double,4> delta_x = -f_out;		
//                 array_1d<int,4> fixity  = ZeroVector(4);
//                 if(GetGeometry()[0].IsFixed(NODAL_AREA) == true) fixity[0] = 1.0;
//                 if(GetGeometry()[1].IsFixed(NODAL_AREA) == true) fixity[1] = 1.0;
//                 if(GetGeometry()[0].IsFixed(FLOW) == true) fixity[2] = 1.0;
//                 if(GetGeometry()[1].IsFixed(FLOW) == true) fixity[3] = 1.0;
// 
//                 // fixity[0] = 1.0;p
//                 // fixity[3] = 1.0;
// 
//                 //KRATOS_WATCH(jacobian);
//                 //KRATOS_WATCH(delta_x);
//                 for(unsigned int i =0; i<4; i++)
//                     {
//                         if(fixity[i] == 1)
//                         {
//                           for(unsigned int j =0; j<4; j++)
//                           {
//                           jacobian(j,i) = 0.0;
//                           jacobian(i,j) = 0.0;
//                           }
//                           delta_x[i] = 0.0;
//                           jacobian(i,i) = 1.0;
//                         }
//                     }

                //KRATOS_WATCH(jacobian);
                //KRATOS_WATCH(f_out);
                //KRATOS_WATCH(area);
                //KRATOS_WATCH(flow);
                bool singular = lu_factorize(jacobian, permutation);
		if(singular)
		    KRATOS_ERROR(std::logic_error,"singular jacobian found in 11 condition with id",this->Id());
		
                lu_substitute(jacobian,permutation, delta_x);
                convergence = norm_2(delta_x);
                area[0] += delta_x[0];
                area[1] += delta_x[1];
                flow[0] += delta_x[2];
                flow[1] += delta_x[3];
                // KRATOS_WATCH(flow);
                // we have to add the relative convergence check
                if(convergence < tolerance)
                    break;
//                else
//                    std::cout << "NO CONVERGEEEEEEEEEEEEEEEE::: Artery11" << std::endl;
//                    KRATOS_ERROR(std::runtime_error, "Artery11", "");
//                    //KRATOS_WATCH(x);
//                    KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA));
            }
	    
            //node 0
            double A1 = area[0];
            double A0 = GetGeometry()[0].GetValue(NODAL_AREA);
            double C = (beta[0]*sqrt(A1*A1*A1)) / (3.0*density*A0);

            rRightHandSideVector[0] = -flow[0];
            rRightHandSideVector[1] = -(C + (coriolis_coefficient*flow[0]*flow[0]/(A1)));

                if(A1 == 0.00 || A0 == 0.00)
                {
                KRATOS_WATCH("FATHER")
                KRATOS_WATCH(A1);
                KRATOS_WATCH(GetProperties().Id());
                KRATOS_WATCH(this->Id());
                KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, Please check your in boundary 1-1 conditions used", "");
                }
	    
	    
            //node 1
            A1 = area[1];
            A0 = GetGeometry()[1].GetValue(NODAL_AREA);
            C = beta[1]*sqrt(A1*A1*A1)/(3.0*density*A0);
            rRightHandSideVector[2] = flow[1];
            rRightHandSideVector[3] = (C + ((coriolis_coefficient*flow[1]*flow[1]))/(A1));

	    if(A1 == 0.00 || A0 == 0.00)
	    {
		KRATOS_WATCH("SON")
		KRATOS_WATCH(A1);
		KRATOS_WATCH(GetProperties().Id());
		KRATOS_WATCH(this->Id());
		KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, Please check your in boundary 1-1 conditions used", "");
	    }
	    
	    
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void Artery11Condition::Initialize()
{
    KRATOS_TRY

    //compute the lenght of the element
//    array_1d<double,3> lvec = GetGeometry()[1].Coordinates();
//    lvec -= GetGeometry()[0].Coordinates();

//    mL = norm_2(lvec);

//    //save area to the nodes. as well as its nodal mass
//    GetGeometry()[0].SetLock();
//    GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += 0.5*mL;
//    GetGeometry()[0].UnSetLock();
//    KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS))

//    GetGeometry()[1].SetLock();
//    GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += 0.5*mL;
//    GetGeometry()[1].UnSetLock();
//    KRATOS_WATCH(GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS))
//            KRATOS_WATCH("------")
     KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void Artery11Condition::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void Artery11Condition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}

//************************************************************************************
//************************************************************************************

void Artery11Condition::GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}


//************************************************************************************
//************************************************************************************
void Artery11Condition::Calculate(const Variable<double >& rVariable,
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

int Artery11Condition::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //check the area

    //check if if is in the XY plane

    //check that no variable has zero key


    return 0;


    KRATOS_CATCH("");
}

void Artery11Condition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY



    KRATOS_CATCH("");
}



void Artery11Condition::CalculateFunctional4(array_1d<double,4>& rFunctional,
                         array_1d<double, 2> const& Area,
                         array_1d<double, 2> const& Flow, array_1d<double,2> const& ArteryProperty, array_1d<double,2> const& Coef,
                         array_1d<double, 2> const& WaveVelocity,
                         double BloodDensity)
{

  KRATOS_TRY
 // KRATOS_WATCH(Area)
  
    const double qa0 = Flow[0]/Area[0];
    const double qa1 = Flow[1]/Area[1];

    rFunctional[0]=Flow[0] + (4.00 * Coef[0]* pow(Area[0],1.25)) - WaveVelocity[0]*Area[0];
    rFunctional[1]=(ArteryProperty[0]*(sqrt(Area[0])-sqrt(GetGeometry()[0].GetValue(NODAL_AREA))))-
             (ArteryProperty[1]*(sqrt(Area[1])-sqrt(GetGeometry()[1].GetValue(NODAL_AREA))))+
             (0.50*BloodDensity*(qa0*qa0 - qa1*qa1));
    rFunctional[2]=Flow[0]-Flow[1];
    rFunctional[3]=Flow[1]-4.00*Coef[1] * pow(Area[1],1.25)-(WaveVelocity[1]*Area[1]);
    KRATOS_CATCH("");

}

void Artery11Condition::CalculateJacobian4(Matrix& rJacobian,
                         array_1d<double, 2> const& Area,
                         array_1d<double, 2> const& Flow, array_1d<double,2> const& ArteryProperty, array_1d<double,2> const& Coef,
                         array_1d<double, 2> const& WaveVelocity,
                         double BloodDensity)
{


    //const double qa0 = Flow[0]/Area[0];
    //const double qa1 = Flow[1]/Area[1];

    rJacobian(0,0)= 5.00 * Coef[0] * pow(Area[0],0.25) - WaveVelocity[0];
    rJacobian(0,1)= 0.00;
    rJacobian(0,2)= 1.00;
    rJacobian(0,3)= 0.00;

    rJacobian(1,0)=  ((0.50 * ArteryProperty[0] / sqrt(Area[0])) - ((BloodDensity* Flow[0] * Flow[0]) / (Area[0]*Area[0]*Area[0])));
    rJacobian(1,1)= -((0.50 * ArteryProperty[1] / sqrt(Area[1])) + ((BloodDensity* Flow[1] * Flow[1]) / (Area[1]*Area[1]*Area[1])));
    rJacobian(1,2)=  (BloodDensity * Flow[0]) / (Area[0] * Area[0]);
    rJacobian(1,3)= -(BloodDensity * Flow[1]) / (Area[1] * Area[1]);

    rJacobian(2,0)= 0.00;
    rJacobian(2,1)= 0.00;
    rJacobian(2,2)= 1.00;
    rJacobian(2,3)= -1.00;

    rJacobian(3,0)= 0.00;
    rJacobian(3,1)= (-5.00 * Coef[1] * pow(Area[1],0.25)) - WaveVelocity[1];
    rJacobian(3,2)= 0.00;
    rJacobian(3,3)= 1.00;
}

} // Namespace Kratos



