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
#include "custom_conditions/artery_3d_to_1d_condition.h"
#include "utilities/math_utils.h"
#include "blood_flow_application.h"
#include "boost/numeric/ublas/lu.hpp"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

Artery3Dto1DCondition::Artery3Dto1DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

Artery3Dto1DCondition::Artery3Dto1DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer Artery3Dto1DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Condition::Pointer(new Artery3Dto1DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

Artery3Dto1DCondition::~Artery3Dto1DCondition()
{
}

//************************************************************************************
//************************************************************************************

void Artery3Dto1DCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense to computer the system matrix for an explicit condition", "");
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void Artery3Dto1DCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //get data as needed
    //const double dynamic_viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    const double density = GetProperties()[DENSITY];
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

    //copy to node 1 from 0
    GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS) = GetGeometry()[1].FastGetSolutionStepValue(YOUNG_MODULUS);
    GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO) = GetGeometry()[1].FastGetSolutionStepValue(POISSON_RATIO);
    GetGeometry()[0].FastGetSolutionStepValue(THICKNESS) = GetGeometry()[1].FastGetSolutionStepValue(THICKNESS);
    GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS);
    GetGeometry()[0].GetValue(NODAL_AREA) = GetGeometry()[1].GetValue(NODAL_AREA);
    //here we should change the
    GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) = GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA);
    //assign flow
    GetGeometry()[1].SetLock();
    GetGeometry()[1].FastGetSolutionStepValue(FLOW) = GetGeometry()[0].FastGetSolutionStepValue(FLOW);
    GetGeometry()[1].UnSetLock();
    //loop on nodes
    for (unsigned int i=0; i<2; i++)
    {
        //const double E = GetGeometry()[i].FastGetSolutionStepValue(YOUNG_MODULUS);
        //const double nu =GetGeometry()[i].FastGetSolutionStepValue(POISSON_RATIO);
        //const double H0 =GetGeometry()[i].FastGetSolutionStepValue(THICKNESS);
        //beta[i]=GetGeometry()[i].FastGetSolutionStepValue(BETA);
        //beta[i] = E*H0*1.77245385/(1.0-nu*nu);

        Node<3>& r_node = GetGeometry()[i];
        coef[i] = r_node.FastGetSolutionStepValue(C0);
        beta[i] = r_node.FastGetSolutionStepValue(BETA);
        area[i] = r_node.FastGetSolutionStepValue(NODAL_AREA);
        flow[i] = r_node.FastGetSolutionStepValue(FLOW);
        artery_property[i]=beta[i]/GetGeometry()[i].GetValue(NODAL_AREA);
        //coef[i] = sqrt(beta[i]/(2*density * GetGeometry()[i].GetValue(NODAL_AREA)));
        KRATOS_WATCH(r_node.GetId())
    }

    flow[1] = flow[0];
    wave_velocity[0] = (flow[0] / area[0]) + 4.00 * sqrt(beta[0] / (2.00 * density * GetGeometry()[0].GetValue(NODAL_AREA))) * pow(area[0],0.25);
    wave_velocity[1] = (flow[1] / area[1]) - 4.00 * sqrt(beta[1] / (2.00 * density * GetGeometry()[1].GetValue(NODAL_AREA))) * pow(area[1],0.25);
    //KRATOS_WATCH(wave_velocity);
    double convergence;
    unsigned int max_iterations = 100;
    double tolerance = 1e-6;
    for(unsigned int i = 0 ; i < max_iterations ; i++)
    {
        CalculateFunctional4(f_out, area, flow, artery_property, coef, wave_velocity, density);
        CalculateJacobian4(jacobian, area, flow, artery_property, coef, wave_velocity, density);
        //KRATOS_WATCH(f_out);
        //KRATOS_WATCH(jacobian);
        permutation_matrix<double> permutation(4);
        array_1d<double,4> delta_x = -f_out;
        array_1d<int,4> fixity  = ZeroVector(4);
        if(GetGeometry()[0].IsFixed(NODAL_AREA) == true) fixity[0] = 1.0;
        if(GetGeometry()[1].IsFixed(NODAL_AREA) == true) fixity[1] = 1.0;
        if(GetGeometry()[0].IsFixed(FLOW) == true) fixity[2] = 1.0;
        if(GetGeometry()[1].IsFixed(FLOW) == true) fixity[3] = 1.0;
        // fixity[0] = 1.0;
        // fixity[3] = 1.0;
        //KRATOS_WATCH(jacobian);
        //KRATOS_WATCH(delta_x);
        for(unsigned int i =0; i<4; i++)
        {
            if(fixity[i] == 1)
                {
                for(unsigned int j =0; j<4; j++)
                    {
                        jacobian(j,i) = 0.0;
                        jacobian(i,j) = 0.0;
                    }
                delta_x[i] = 0.0;
                jacobian(i,i) = 1.0;
                }
        }
        //		KRATOS_WATCH(jacobian);
        //		KRATOS_WATCH(f_out);
        //		KRATOS_WATCH(area);
        //KRATOS_WATCH(flow);
        lu_factorize(jacobian, permutation);
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
    }

    //node 0
    double A1 = area[0];
    const double m = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
    rRightHandSideVector[0] = flow[0] * m;
    rRightHandSideVector[1] = area[0] * m;
    //node 1
    A1 = area[1];
    A0 = GetGeometry()[1].GetValue(NODAL_AREA);
    C = beta[1]*sqrt(A1*A1*A1)/(3.0*density*A0);
    rRightHandSideVector[2] = flow[1] * m;
    rRightHandSideVector[3] = (C + coriolis_coefficient*flow[1]*flow[1]/(A1)) * m;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void Artery3Dto1DCondition::Initialize()
{
    KRATOS_TRY
    
    GetGeometry()[1].Fix(FLOW);

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure

void Artery3Dto1DCondition::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void Artery3Dto1DCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}

//************************************************************************************
//************************************************************************************

void Artery3Dto1DCondition::GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error, "method not implemented (it does not make sense for an explicit condition", "");
}


//************************************************************************************
//************************************************************************************
void Artery3Dto1DCondition::Calculate(const Variable<double >& rVariable,
                              double& Output,
                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");

}

int Artery3Dto1DCondition::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //check the area

    //check if if is in the XY plane

    //check that no variable has zero key


    return 0;


    KRATOS_CATCH("");
}

void Artery3Dto1DCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY



    KRATOS_CATCH("");
}


void Artery3Dto1DCondition::CalculateFunctional4(array_1d<double,4>& rFunctional,
                         array_1d<double, 2> const& Area,
                         array_1d<double, 2> const& Flow, array_1d<double,2> const& ArteryProperty, array_1d<double,2> const& Coef,
                         array_1d<double, 2> const& WaveVelocity,
                         double BloodDensity)
{


    const double qa0 = Flow[0]/Area[0];
    const double qa1 = Flow[1]/Area[1];

    rFunctional[0]=Flow[0] + 4.00 * Coef[0]* pow(Area[0],1.25) - WaveVelocity[0]*Area[0];

    rFunctional[1]=ArteryProperty[0]*(sqrt(Area[0])-sqrt(GetGeometry()[0].GetValue(NODAL_AREA)))-
		     ArteryProperty[1]*(sqrt(Area[1])-sqrt(GetGeometry()[1].GetValue(NODAL_AREA)))+
		     0.50*BloodDensity*(qa0*qa0 - qa1*qa1);
    rFunctional[2]=Flow[0]-Flow[1];
    rFunctional[3]=Flow[1]-4.00*Coef[1] * pow(Area[1],1.25)-WaveVelocity[1]*Area[1];


}

void Artery3Dto1DCondition::CalculateJacobian4(Matrix& rJacobian,
                         array_1d<double, 2> const& Area,
                         array_1d<double, 2> const& Flow, array_1d<double,2> const& ArteryProperty, array_1d<double,2> const& Coef,
                         array_1d<double, 2> const& WaveVelocity,
                         double BloodDensity)
{


    const double qa0 = Flow[0]/Area[0];
    const double qa1 = Flow[1]/Area[1];


    rJacobian(0,0)= 5.00 * Coef[0] * pow(Area[0],0.25) - WaveVelocity[0];
    rJacobian(0,1)= 0.00;
    rJacobian(0,2)= 1.00;
    rJacobian(0,3)= 0.00;

    rJacobian(1,0)= 0.50 * ArteryProperty[0] / sqrt(Area[0]) - BloodDensity* (qa0 * qa0 / Area[0]);
    rJacobian(1,1)= -0.50 * ArteryProperty[1] / sqrt(Area[1]) + BloodDensity* (qa1 * qa1 / Area[1]);
    rJacobian(1,2)= BloodDensity * Flow[0] / (Area[0] * Area[0]);
    rJacobian(1,3)= -BloodDensity * Flow[1] / (Area[1] * Area[1]);

    rJacobian(2,0)= 0.00;
    rJacobian(2,1)= 0.00;
    rJacobian(2,2)= 1.00;
    rJacobian(2,3)= -1.00;

    rJacobian(3,0)= 0.00;
    rJacobian(3,1)= -5.00 * Coef[1] * pow(Area[1],0.25) - WaveVelocity[1];
    rJacobian(3,2)= 0.00;
    rJacobian(3,3)= 1.00;


}


} // Namespace Kratos



