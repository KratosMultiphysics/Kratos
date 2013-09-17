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
#include "custom_conditions/artery_12_condition.h"
#include "utilities/math_utils.h"
#include "blood_flow_application.h"
#include "boost/numeric/ublas/lu.hpp"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

Artery12Condition::Artery12Condition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

Artery12Condition::Artery12Condition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
{


}

Condition::Pointer Artery12Condition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Condition::Pointer(new Artery12Condition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

Artery12Condition::~Artery12Condition()
{
}

//************************************************************************************
//************************************************************************************

void Artery12Condition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense to computer the system matrix for an explicit condition", "");
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void Artery12Condition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

            //get data as needed
            //const double dynamic_viscosity = GetProperties()[DYNAMIC_VISCOSITY];
            const double density = GetProperties()[DENSITY];
            //const double E = GetProperties()[YOUNG_MODULUS];
            //const double nu = GetProperties()[POISSON_RATIO];
            //const double pi = 3.14159265;
            const double coriolis_coefficient = 1.1; // TODO
            //const double kr_coefficient = 1.0;
            //const double kinematic_viscosity = dynamic_viscosity/density;
	    

            //resize the vector to the correct size
            if (rRightHandSideVector.size() != 6)
                rRightHandSideVector.resize(6,false);

            using namespace boost::numeric::ublas;

            array_1d<double,6> f_out;
            Matrix jacobian = ZeroMatrix(6,6);

            array_1d<double,3> area;
            array_1d<double,3> flow;
            array_1d<double,3> wave_velocity;
            array_1d<double,3> artery_property;
            array_1d<double,3> coef;
            array_1d<double,3> beta;
            array_1d<double,3> InitialArea;

            //loop on nodes
            for (unsigned int i=0; i<3; i++)
            {            
                //const double E = GetGeometry()[i].FastGetSolutionStepValue(YOUNG_MODULUS);
                //const double nu = GetGeometry()[i].FastGetSolutionStepValue(POISSON_RATIO);
                //const double mH0 = GetGeometry()[i].FastGetSolutionStepValue(THICKNESS);
                InitialArea[i] = GetGeometry()[i].GetValue(NODAL_AREA);

                Node<3>& r_node = GetGeometry()[i];
                area[i] = r_node.FastGetSolutionStepValue(NODAL_AREA);
                flow[i] = r_node.FastGetSolutionStepValue(FLOW);
                coef[i] = r_node.FastGetSolutionStepValue(C0);
                beta[i] = r_node.FastGetSolutionStepValue(BETA);
                //artery_property[i]=beta[i]/GetGeometry()[i].GetValue(NODAL_AREA);
                //coef[i] = sqrt(beta[i]/(2*density * GetGeometry()[i].GetValue(NODAL_AREA)));
                //const double r0 = r_node.FastGetSolutionStepValue(RADIUS);
                //beta[i] = E*mH0*1.77245385/(1.0-nu*nu);
                //artery_property[i]=mBeta/mInitialArea[i];
                artery_property[i]=beta[i]/InitialArea[i];
                //coef[i] = sqrt(mBeta/(2*density * mInitialArea[i]));
                //coef[i] = sqrt(beta[i]/(2*density * InitialArea[i]));
            }

            //wave_velocity[0] = (flow[0] / area[0]) + 4.00 * sqrt(mBeta / (2.00 * density * mInitialArea[0])) * pow(area[0],0.25);
            //wave_velocity[1] = (flow[1] / area[1]) - 4.00 * sqrt(mBeta / (2.00 * density * mInitialArea[1])) * pow(area[1],0.25);
            //wave_velocity[2] = (flow[2] / area[2]) - 4.00 * sqrt(mBeta / (2.00 * density * mInitialArea[2])) * pow(area[2],0.25);

            wave_velocity[0] = (flow[0] / area[0]) + (4.00 * (sqrt(beta[0] / (2.00 * density * InitialArea[0]))) * pow(area[0],0.25));
            wave_velocity[1] = (flow[1] / area[1]) - (4.00 * (sqrt(beta[1] / (2.00 * density * InitialArea[1]))) * pow(area[1],0.25));
            wave_velocity[2] = (flow[2] / area[2]) - (4.00 * (sqrt(beta[2] / (2.00 * density * InitialArea[2]))) * pow(area[2],0.25));


            double convergence;
            unsigned int max_iterations = 100;
            double tolerance = 1e-5;

            for(unsigned int i = 0 ; i < max_iterations ; i++)
            {
                CalculateFunctional6(f_out, area, InitialArea, flow, artery_property, coef, wave_velocity, density);
                CalculateJacobian6(jacobian, area, flow, artery_property, coef, wave_velocity, density);

                //KRATOS_WATCH(f_out);
                //KRATOS_WATCH(jacobian);

                permutation_matrix<double> permutation(6);
                array_1d<double,6> delta_x = -f_out;

                //lu_factorize(jacobian, permutation);
                bool singular = lu_factorize(jacobian, permutation);
                if(singular)
                    KRATOS_ERROR(std::logic_error,"singular jacobian found in 12 condition with id",this->Id());

                lu_substitute(jacobian,permutation, delta_x);
                convergence = norm_2(delta_x);

                area[0] += delta_x[0];
                area[1] += delta_x[1];
                area[2] += delta_x[2];
                flow[0] += delta_x[3];
                flow[1] += delta_x[4];
                flow[2] += delta_x[5];
                //KRATOS_WATCH("AAAAA::::SSSS");
                //KRATOS_WATCH(flow);

                // we have to add the relative convergence check
                if(convergence < tolerance)
                    break;
            }


             double A1 = area[0];
             double A0 = InitialArea[0];
             //double E = GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS);
             //double nu = GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO);
             //double mH0 = GetGeometry()[0].FastGetSolutionStepValue(THICKNESS);
             beta[0] = GetGeometry()[0].FastGetSolutionStepValue(BETA);
             //beta[0] = E*mH0*1.77245385/(1.0-nu*nu);
             double C = beta[0]*sqrt(A1*A1*A1)/(3.0*density*A0);                 
             rRightHandSideVector[0] = -flow[0];
             //double temp = (C + coriolis_coefficient*flow[0]*flow[0]/(A1));
             rRightHandSideVector[1] = -(C + (coriolis_coefficient*flow[0]*flow[0]/(A1)));

             if(A1 == 0.00 || A0 == 0.00)
             {
                 KRATOS_WATCH(A0);
                 KRATOS_WATCH(A1);
                 KRATOS_WATCH(GetProperties().Id());
                 KRATOS_WATCH(this->Id());
                 KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, in boundary 1-2 conditions used: father", "");
             }
	     	     
             A1 = area[1];
             A0 = InitialArea[1];
             //E = GetGeometry()[1].FastGetSolutionStepValue(YOUNG_MODULUS);
             //nu = GetGeometry()[1].FastGetSolutionStepValue(POISSON_RATIO);
             //mH0 = GetGeometry()[1].FastGetSolutionStepValue(THICKNESS);
             //beta[1] = E*mH0*1.77245385/(1.0-nu*nu);
             beta[1]=GetGeometry()[1].FastGetSolutionStepValue(BETA);
             C = beta[1]*sqrt(A1*A1*A1)/(3.0*density*A0);
             rRightHandSideVector[2] = flow[1];
             //temp = (C + coriolis_coefficient*flow[1]*flow[1]/(A1));
             rRightHandSideVector[3] = (C + (coriolis_coefficient*flow[1]*flow[1]/(A1)));

	     if(A1 == 0.00 || A0 == 0.00)
	      {
            KRATOS_WATCH(A0);
            KRATOS_WATCH(A1);
            KRATOS_WATCH(GetProperties().Id());
            KRATOS_WATCH(this->Id());
            KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, in boundary 1-2 conditions used:son", "");
	      }	     
	     
             A1 = area[2];
             A0 = InitialArea[2];
             //E = GetGeometry()[2].FastGetSolutionStepValue(YOUNG_MODULUS);
             //nu = GetGeometry()[2].FastGetSolutionStepValue(POISSON_RATIO);
             //mH0 = GetGeometry()[2].FastGetSolutionStepValue(THICKNESS);
             beta[2]=GetGeometry()[2].FastGetSolutionStepValue(BETA);
             //beta[2] = E*mH0*1.77245385/(1.0-nu*nu);
             C = beta[2]*sqrt(A1*A1*A1)/(3.0*density*A0);
             rRightHandSideVector[4] = flow[2];
             //temp = (C + coriolis_coefficient*flow[2]*flow[2]/(A1));
             rRightHandSideVector[5] = (C + (coriolis_coefficient*flow[2]*flow[2]/(A1)));

	     if(A1 == 0.00 || A0 == 0.00)
	      {
            KRATOS_WATCH(A0);
            KRATOS_WATCH(A1);
            KRATOS_WATCH(GetProperties().Id());
            KRATOS_WATCH(this->Id());
            KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, in boundary 1-2 conditions used:son", "");
	      }	     
	     	     
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void Artery12Condition::Initialize()
{
    KRATOS_TRY

        //compute the lenght of the element
//        array_1d<double,3> lvec = GetGeometry()[1].Coordinates();
//        lvec -= GetGeometry()[0].Coordinates();

//        mL = norm_2(lvec);

//        //save area to the nodes. as well as its nodal mass
//        GetGeometry()[0].SetLock();
//        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += 0.5*mL;
//        GetGeometry()[0].UnSetLock();
//        KRATOS_WATCH(GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS))

//        GetGeometry()[1].SetLock();
//        GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += 0.5*mL;
//        GetGeometry()[1].UnSetLock();
//        KRATOS_WATCH(GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS))
//        KRATOS_WATCH("------")

    KRATOS_CATCH("");
    }

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void Artery12Condition::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void Artery12Condition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}

//************************************************************************************
//************************************************************************************

void Artery12Condition::GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}


//************************************************************************************
//************************************************************************************
void Artery12Condition::Calculate(const Variable<double >& rVariable,
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

int Artery12Condition::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //check the area

    //check if if is in the XY plane

    //check that no variable has zero key


    return 0;


    KRATOS_CATCH("");
}

void Artery12Condition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY


    KRATOS_CATCH("");
}



void Artery12Condition::CalculateFunctional6(array_1d<double,6>& rFunctional,
                         array_1d<double, 3> const& Area,
                         array_1d<double, 3> const& InitialArea,
                         array_1d<double, 3> const& Flow, array_1d<double,3> const& ArteryProperty, array_1d<double,3> const& Coef,
                         array_1d<double, 3> const& WaveVelocity,
                         double BloodDensity)
{


    const double qa0 = Flow[0]/Area[0];
    const double qa1 = Flow[1]/Area[1];
    const double qa2 = Flow[2]/Area[2];

    rFunctional[0]=Flow[0] + 4.00 * Coef[0]* pow(Area[0],1.25) - WaveVelocity[0]*Area[0];

    rFunctional[1]=ArteryProperty[0]*(sqrt(Area[0])-sqrt(InitialArea[0]))-
             ArteryProperty[1]*(sqrt(Area[1])-sqrt(InitialArea[1]))+
             0.50*BloodDensity*(qa0*qa0 - qa1*qa1);

    rFunctional[2]=ArteryProperty[0]*(sqrt(Area[0])-sqrt(InitialArea[0]))-
             ArteryProperty[2]*(sqrt(Area[2])-sqrt(InitialArea[2]))+
             0.50*BloodDensity*(qa0*qa0 - qa2*qa2);

    rFunctional[3]=Flow[0]-Flow[1]-Flow[2];
    rFunctional[4]=Flow[1]-4.00*Coef[1] * pow(Area[1],1.25)-WaveVelocity[1]*Area[1];
    rFunctional[5]=Flow[2]-4.00*Coef[2] * pow(Area[2],1.25)-WaveVelocity[2]*Area[2];


}

void Artery12Condition::CalculateJacobian6(Matrix& rJacobian,
                         array_1d<double, 3> const& Area,
                         array_1d<double, 3> const& Flow, array_1d<double,3> const& ArteryProperty, array_1d<double,3> const& Coef,
                         array_1d<double, 3> const& WaveVelocity,
                         double BloodDensity)
{


    const double qa0 = Flow[0]/Area[0];
    const double qa1 = Flow[1]/Area[1];
    const double qa2 = Flow[2]/Area[2];


    rJacobian(0,0)= 5.00 * Coef[0] * pow(Area[0],0.25) - WaveVelocity[0];
    rJacobian(0,1)= 0.00;
    rJacobian(0,2)= 0.00;
    rJacobian(0,3)= 1.00;
    rJacobian(0,4)= 0.00;
    rJacobian(0,5)= 0.00;

    rJacobian(1,0)= 0.50 * ArteryProperty[0] / sqrt(Area[0]) - BloodDensity* (qa0 * qa0 / Area[0]);
    rJacobian(1,1)= -0.50 * ArteryProperty[1] / sqrt(Area[1]) + BloodDensity* (qa1 * qa1 / Area[1]);
    rJacobian(1,2) = 0.00;
    rJacobian(1,3)= BloodDensity * Flow[0] / (Area[0] * Area[0]);
    rJacobian(1,4)= -BloodDensity * Flow[1] / (Area[1] * Area[1]);
    rJacobian(1,5) = 0.00;

    rJacobian(2,0)= 0.50 * ArteryProperty[0] / sqrt(Area[0]) - BloodDensity* (qa0 * qa0 / Area[0]);
    rJacobian(2,1) = 0.00;
    rJacobian(2,2)= -0.50 * ArteryProperty[2] / sqrt(Area[2]) + BloodDensity* (qa2 * qa2 / Area[2]);
    rJacobian(2,3)= BloodDensity * Flow[0] / (Area[0] * Area[0]);
    rJacobian(2,4) = 0.00;
    rJacobian(2,5)= -BloodDensity * Flow[2] / (Area[2] * Area[2]);

    rJacobian(3,0)= 0.00;
    rJacobian(3,1)= 0.00;
    rJacobian(3,2)= 0.00;
    rJacobian(3,3)= 1.00;
    rJacobian(3,4)= -1.00;
    rJacobian(3,5)= -1.00;

    rJacobian(4,0)= 0.00;
    rJacobian(4,1)= -5.00 * Coef[1] * pow(Area[1],0.25) - WaveVelocity[1];
    rJacobian(4,2)= 0.00;
    rJacobian(4,3)= 0.00;
    rJacobian(4,4)= 1.00;
    rJacobian(4,5)= 0.00;

    rJacobian(5,0)= 0.00;
    rJacobian(5,1)= 0.00;
    rJacobian(5,2)= -5.00 * Coef[2] * pow(Area[2],0.25) - WaveVelocity[2];
    rJacobian(5,3)= 0.00;
    rJacobian(5,4)= 0.00;
    rJacobian(5,5)= 1.00;



}

} // Namespace Kratos



