/*
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
esoudah@cimne.upc.edu
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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-10-13 06:58:23 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_ARTERY_TIME_INTEGRATOR_INCLUDED )
#define  KRATOS_ARTERY_TIME_INTEGRATOR_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
//#include "geometries/tetrahedra_3d_4.h"
#include "blood_flow_application.h"

namespace Kratos
{
class ArteryTimeIntegrator
{
public:

    void Initialize(ModelPart& ThisModelPart)
    {
        //const ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();

        //reset the acceleration on the nodes
        //std::cout << "ini:------------------------------------------------------------------ ";
        ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            in->FastGetSolutionStepValue(NODAL_MASS) = 0.0;
        }

        //compute the projection (first step)::
        //DESDE AQUI VOY AL Initialize DEL ELEMENT(ARTERY ELEMENT)
        //CARGO LAS PROPIEDADES DEL ELEMENTO
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i!=ThisModelPart.ElementsEnd(); i++)
        {
            i->Initialize();
        }

        //compute the projection (first step)::
        for(ModelPart::ConditionIterator i = ThisModelPart.ConditionsBegin();i!=ThisModelPart.ConditionsEnd(); i++)
        {
            // DESDE AQUI VOY a::  Artery12Condition::Initialize()
            i->Initialize();
        }
    }


    void SolveStep(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        //double total_time=ThisModelPart.GetProcessInfo()[TIME];

        InitializeWorkArray( ThisModelPart );

        //step1
        //std::cout << "PASO AL SIGUIENTE PASO:::::: ";
        //std::cout << std::endl;
        //KRATOS_WATCH("1")
        //if (total_time > 0.0140)
            //KRATOS_WATCH ("");

        //CheckSolution (ThisModelPart, 1);
        //KRATOS_WATCH("2")
        ComputeRHS( ThisModelPart );
        //KRATOS_WATCH("3")
        ComputeRHS_Actualize( ThisModelPart );
        //KRATOS_WATCH("4")
        //CheckSolution (ThisModelPart, 2);
        KRATOS_CATCH("")
    }

    double EstimateDeltaTime(ModelPart& ThisModelPart, double CFL, double minLength)
    {
        KRATOS_TRY

        //KRATOS_ERROR(std::logic_error,"should not enter in EstimateDeltaTime","");
        double delta_time = 1.00e12;
        //double Umax=0.0;
        double element_dt=0.0;

// 	for(ModelPart::NodesContainerType::iterator i_node=ThisModelPart.NodesBegin(); i_node != ThisModelPart.NodesEnd(); i_node++)
// 	  {
// 	    const double& AAAA = i_node->GetValue(NODAL_AREA);
//             KRATOS_WATCH(AAAA);
// 	  }

        for(ModelPart::ElementsContainerType::iterator i_element = ThisModelPart.ElementsBegin(); i_element !=ThisModelPart.ElementsEnd(); i_element++)
        {
            const double coriolis_coefficient = 1.1;

            // EDU:PODRIAMOS SABER CUAL ES LA LONGITUD DEL ELEMENTO MINIMA
            //double h = i_element->GetGeometry().Length();

            //h=0.0008;

            const double& A0 = i_element->GetGeometry()[0].GetValue(NODAL_AREA);
            const double& A1 = i_element->GetGeometry()[1].GetValue(NODAL_AREA);

            const double& A0_actual = i_element->GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
            const double& A1_actual = i_element->GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA);

            //KRATOS_WATCH(i_element->GetGeometry()[0])

            //KRATOS_WATCH(i_element->GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA));
            //KRATOS_WATCH(i_element->GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA));
            //KRATOS_WATCH(h);
            //KRATOS_WATCH(A1);
            //KRATOS_WATCH(A0_actual);
            //KRATOS_WATCH(A1_actual);

            if((A0_actual == 0.00) || (A1_actual== 0.00))
                KRATOS_WATCH(A0_actual);
            //KRATOS_WATCH(i_element->GetGeometry()[0]);

            if((A0_actual == 0.00) || (A1_actual== 0.00))
                KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, Please check your model or boundary conditions used", "");

            const double v0 = i_element->GetGeometry()[0].FastGetSolutionStepValue(FLOW) / A0_actual;
            const double v1 = i_element->GetGeometry()[1].FastGetSolutionStepValue(FLOW) / A1_actual;
            //double v = std::max(v0,v1);

            const double density = i_element->GetProperties()[DENSITY];

            const double E0 = i_element->GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS);
            const double nu0 = i_element->GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO);
            const double thickness0 = i_element->GetGeometry()[0].FastGetSolutionStepValue(THICKNESS);
            double beta = E0*thickness0*1.77245385/(1.0-nu0*nu0);


            /* KRATOS_WATCH(v0);
              KRATOS_WATCH(v1);
            */
            const double C1_0 = sqrt(beta / (2.0 * density * A0)) * pow(A0_actual,0.25);
            //const double c_alpha = sqrt(c_max*c_max + v*v * coriolis_coefficient* (coriolis_coefficient-1));
            const double c_alpha_0 = sqrt(C1_0*C1_0 + (v0*v0 * coriolis_coefficient* (coriolis_coefficient-1)));
            const double la_0 = coriolis_coefficient * v0 + c_alpha_0;
// 	    Umax = std::max(la_0, Umax);
// 	    element_dt = fabs(CFL * minLength / Umax);
            element_dt = (CFL * minLength) / la_0;
            delta_time = std::min(delta_time, element_dt);


// 	    		c_1=dsqrt(beta0(k)/(2.d0*blood_dens*A0(k)))*(U_stored(1,j))**0.25d0
// 			c_alfa=dsqrt(c_1**2+((U_stored(2,j)/U_stored(1,j))**2)*coriolis*(coriolis-1))
// 			lambda=coriolis*U_stored(2,j)/U_stored(1,j)+c_alfa
// 			h_array=CFL*lengthmin/lambda
// 			h=min(h,h_array)

            const double E1 = i_element->GetGeometry()[1].FastGetSolutionStepValue(YOUNG_MODULUS);
            const double nu1 = i_element->GetGeometry()[1].FastGetSolutionStepValue(POISSON_RATIO);
            const double thickness1 = i_element->GetGeometry()[1].FastGetSolutionStepValue(THICKNESS);
            beta = E1*thickness1*1.77245385/(1.0-nu1*nu1);

// 	    KRATOS_WATCH(E1);
// 	    KRATOS_WATCH(nu1);
// 	    KRATOS_WATCH(thickness1);
// 	    KRATOS_WATCH(beta);
// 	    KRATOS_WATCH(density);

            const double C1_1 = sqrt(beta / (2.0 * density * A1)) * pow(A1_actual,0.25);
            //const double c_alpha = sqrt(c_max*c_max + v*v * coriolis_coefficient* (coriolis_coefficient-1));
            const double c_alpha_1 = sqrt(C1_1*C1_1 + v1*v1*coriolis_coefficient*(coriolis_coefficient-1));
            const double la_1 = coriolis_coefficient * v1 + c_alpha_1;

            /*Umax = std::max(la_1, Umax);
            element_dt = fabs(CFL * minLength / Umax);
            */
            element_dt = (CFL * minLength) /la_1;
            delta_time = std::min(delta_time, element_dt);
            //delta_time=0.001;
            //const double c_max = std::max(C1_0, C1_1);
            //const double c_max = std::max(la_0, la_1);

            //const double la = coriolis_coefficient * v + c_alpha;
            //la = coriolis_coefficient * v + c_max;

            //double element_dt = fabs(CFL * h / c_max);

//            KRATOS_WATCH(i_element->Id());
//            KRATOS_WATCH(beta);
//            KRATOS_WATCH(C1_0);
//            KRATOS_WATCH(C1_1);
//            KRATOS_WATCH(v);
//            KRATOS_WATCH(h)

            //delta_time = std::min(delta_time, element_dt);
            //KRATOS_WATCH(delta_time)
        }
//
        //KRATOS_WATCH(delta_time)
        delta_time = fabs(delta_time);
        return delta_time;
        //KRATOS_WATCH("DELTA_TIME");
        //KRATOS_WATCH(delta_time);

        KRATOS_CATCH("")
    }

    double Element_minLength (ModelPart& ThisModelPart)
    {

        KRATOS_TRY
        KRATOS_WATCH("Version Arteries Solver 1D :: 1.0")
        double minLength=1e+12;

        for(ModelPart::ElementsContainerType::iterator i_element = ThisModelPart.ElementsBegin(); i_element !=ThisModelPart.ElementsEnd(); i_element++)
        {
            double h = i_element->GetGeometry().Length();
            minLength = std::min(minLength, h);
        }
        KRATOS_WATCH(minLength);
        return minLength;

        KRATOS_CATCH("")
    }

    double CheckCardiacCovergence (ModelPart& ThisModelPart, double time_cardiac_cycle )
    {

        KRATOS_TRY
        bool Check_Convergence=false;
        double Aux_convergencia=0;
        double Aux_convergencia2=0;
        double Aux_convergencia3=0;
        double Norma_Conver= 0;
        double Norma=0;
        double dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
        double t = ThisModelPart.GetProcessInfo()[TIME];
        //double total_time = ThisModelPart.GetProcessInfo()[final_time];
        //KRATOS_WATCH (time_cardiac_cycle);
        //KRATOS_WATCH (dt);
        //KRATOS_WATCH (t);
        //KRATOS_WATCH(fabs((time_cardiac_cycle/3)-t));
        //KRATOS_WATCH (dt/2);

        if (fabs((time_cardiac_cycle/3)-t) < (dt/2))

        {
        //PRINT*, "Analisys de Convergencia para el tiempo"
            for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
            {
                //double pressure=initial_pressure+(beta/original_area)*(sqrt(nodal_area)-sqrt(original_area));
                double pressure = in->FastGetSolutionStepValue(PRESSURE);
                double pressureold = in->FastGetSolutionStepValue(PRESSURE,1);
//                double FLOWWW = in->GetSolutionStepValue(FLOW);
//                KRATOS_WATCH(FLOWWW);
//                KRATOS_WATCH(pressure);
//                KRATOS_WATCH(pressureold);
                Aux_convergencia=Aux_convergencia+((pressure-pressureold)*(pressure-pressureold));
                Aux_convergencia2=Aux_convergencia2+((pressureold)*(pressureold));
            }
            Check_Convergence = true;
        }
	
        if (Check_Convergence == true)
        {
            Aux_convergencia3=sqrt(Aux_convergencia/Aux_convergencia2);
            Norma_Conver= fabs(Norma-Aux_convergencia3);
            if ( (Aux_convergencia3 < 0.01) && (Norma_Conver < 0.01) )
            {
                  time_cardiac_cycle = 10000000; //Solution Converge (this value allow to stop the calculation cyc_num>Num_T)
                  KRATOS_ERROR(std::runtime_error, "THE MODEL HAS CONVERGED", "");
            }
            //KRATOS_WATCH (Aux_convergencia);
            //KRATOS_WATCH (Aux_convergencia2);
            KRATOS_WATCH(Norma_Conver);
            Norma = Aux_convergencia3;
            Norma_Conver= 0;
            Aux_convergencia=0;
            Aux_convergencia2=0;
            Aux_convergencia3=0;
            Check_Convergence =false;
        }
        return time_cardiac_cycle;
        KRATOS_CATCH("")
    }

    
    void ComputePressure(ModelPart& ThisModelPart, double dyastolic_pressure)
    {

        KRATOS_TRY

// 	for(ModelPart::ElementsContainerType::iterator i_element = ThisModelPart.ElementsBegin(); i_element !=ThisModelPart.ElementsEnd(); i_element++)
// 	  {
// 	    const double& A0 = i_element->GetGeometry()[0].GetValue(NODAL_AREA);
//             const double& A1 = i_element->GetGeometry()[1].GetValue(NODAL_AREA);
//
// 	    const double& A0_actual = i_element->GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
//             const double& A1_actual = i_element->GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA);
//
// 	    //const double density = i_element->GetProperties()[DENSITY];
//
// 	    double E0 = i_element->GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS);
//             double nu0 = i_element->GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO);
//             double thickness0 = i_element->GetGeometry()[0].FastGetSolutionStepValue(THICKNESS);
// 	    double beta = E0*thickness0*1.77245385/(1.0-nu0*nu0);
//
// 	    double pressure=(beta/A0)*(sqrt(A0_actual)-sqrt(A0));
// 	    //i_element->GetGeometry()[0].SetSolutionStepValue(PRESSURE,0,pressure);
//
// 	    E0 = i_element->GetGeometry()[1].FastGetSolutionStepValue(YOUNG_MODULUS);
//             nu0 = i_element->GetGeometry()[1].FastGetSolutionStepValue(POISSON_RATIO);
//             thickness0 = i_element->GetGeometry()[1].FastGetSolutionStepValue(THICKNESS);
// 	    beta = E0*thickness0*1.77245385/(1.0-nu0*nu0);
//
// 	    pressure=(beta/A1)*(sqrt(A1_actual)-sqrt(A1));
// 	    //i_element->GetGeometry()[1].SetSolutionStepValue(PRESSURE,0,pressure);
// 	  }

        for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
        {
            double original_area = in->GetValue(NODAL_AREA);
            //double initial_pressure = in->GetValue(PRESSURE); //sytolic pressure. All nodes must have the same systolic pressure. Read from the config file.
            //KRATOS_WATCH(initial_pressure);
            //double E0 = in->FastGetSolutionStepValue(YOUNG_MODULUS);
            //double nu0 = in->FastGetSolutionStepValue(POISSON_RATIO);
            //double thickness0 = in->FastGetSolutionStepValue(THICKNESS);

            double beta = in->FastGetSolutionStepValue(BETA);
            //double beta2 = in->GetValue(BETA);
            //double beta3 = in->GetSolutionStepValue(BETA);
            //double beta4 = E0*thickness0*1.77245385/(1.0-nu0*nu0);
            double nodal_area = in->FastGetSolutionStepValue(NODAL_AREA);
            //double nodal_flow = in->FastGetSolutionStepValue(FLOW);
            //KRATOS_WATCH(original_area);
            //        KRATOS_WATCH(nu0);
            //            KRATOS_WATCH(thickness0);
            // 	    KRATOS_WATCH(beta);
                        //KRATOS_WATCH("-------------------------------------------------");
                        //KRATOS_WATCH(nodal_area);
            // 	    KRATOS_WATCH("-------------------------------------------------");
            // 	    KRATOS_WATCH(kkk);

            //KRATOS_WATCH(NODAL_AREA)
            //KRATOS_WATCH(dyastolic_pressure);
            //KRATOS_WATCH(beta);
            //KRATOS_WATCH(beta2);
            //KRATOS_WATCH(beta3);
            //KRATOS_WATCH(beta4);
            //KRATOS_WATCH(nodal_area);
           // KRATOS_WATCH(original_area);
            double pressure=dyastolic_pressure+((beta/original_area)*(sqrt(nodal_area)-sqrt(original_area)));
            //KRATOS_WATCH(pressure);
            in->FastGetSolutionStepValue(PRESSURE)=pressure;
            //FastSetSolutionStepValue(PRESSURE,0,pressure);
            //KRATOS_WATCH("------------------------------------");
        }
        KRATOS_CATCH("")
    }


private:

    void ComputeRHS(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();

        //reset the acceleration on the nodes
        ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            noalias(in->FastGetSolutionStepValue(RHS)) = ZeroVector(3);
        }

        //compute the projection (first step)
        Vector rhs_el(4);
        //KRATOS_WATCH(rhs_el)

        //ELEMENTOS
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i!=ThisModelPart.ElementsEnd(); i++)
        {
            i->CalculateRightHandSide(rhs_el,CurrentProcessInfo);
            //KRATOS_WATCH(rhs_el)
            //int id = i->Id();
            //KRATOS_WATCH (id);
            i->CalculateRightHandSide(rhs_el,CurrentProcessInfo);
            //KRATOS_WATCH ('--------------------------------------------');

            //KRATOS_WATCH(rhs_el);

            Geometry< Node<3> >& geom = i->GetGeometry();
            for(unsigned int j=0; j<geom.size(); j++)
            {
                int base = 2*j;
                array_1d<double,3>& rhs_node = geom[j].FastGetSolutionStepValue(RHS);
                //nodal_mass[j] = geom[j].FastGetSolutionStepValue(NODAL_MASS);
                //double nodal_mass = geom[j].GetSolutionStepValue(NODAL_MASS);
                 //std::cout << "old: ";KRATOS_WATCH(rhs_node)
                 //rhs_el[base] = rhs_el[base]/nodal_mass[0];
                 //rhs_el[base+1] = rhs_el[base]/nodal_mass[1];
                 //KRATOS_WATCH('antes contribution')
                 //KRATOS_WATCH(rhs_node)
                 //KRATOS_WATCH(rhs_el)
                 rhs_node[0] += rhs_el[base];  ///nodal_mass[0];
                 rhs_node[1] += rhs_el[base+1]; ///nodal_mass[1];
                 //KRATOS_WATCH('despues contribution')
                 //KRATOS_WATCH(rhs_node)
            }
         }

      // KRATOS_WATCH("----------------BOUNDARIES--------------------------------------------")
        //CONDICIONES DE CONTORNO
        for(ModelPart::ConditionsContainerType::iterator i = ThisModelPart.ConditionsBegin();
                i!=ThisModelPart.ConditionsEnd(); i++)
        {
            //KRATOS_WATCH(i)
            i->CalculateRightHandSide(rhs_el,CurrentProcessInfo);
            Geometry< Node<3> >& geom = i->GetGeometry();
            //array_1d<double,2> nodal_mass;
            //int id = i->Id();
            //KRATOS_WATCH (id);
            for(unsigned int j=0; j<geom.size(); j++)
            {
                int base = 2*j;
                array_1d<double,3>& rhs_node = geom[j].FastGetSolutionStepValue(RHS);
                //KRATOS_WATCH(geom[j]);
                //nodal_mass[j] = geom[j].FastGetSolutionStepValue(NODAL_MASS);
                //double nodal_mass = geom[j].GetSolutionStepValue(NODAL_MASS);
                //std::cout << "old: ";KRATOS_WATCH(rhs_el[base]);
                //std::cout << "old: ";KRATOS_WATCH(rhs_el[base+1])
                //rhs_el[base] = rhs_el[base]/nodal_mass[0];
                //rhs_el[base+1] = rhs_el[base]/nodal_mass[1];
                //KRATOS_WATCH(rhs_el)
                rhs_node[0] += rhs_el[base];  ///nodal_mass[0];
                rhs_node[1] += rhs_el[base+1]; ///nodal_mass[1];
                //KRATOS_WATCH(rhs_node)
                    //KRATOS_WATCH ('CONDITION');
                    //KRATOS_WATCH (id);
                    //KRATOS_WATCH(geom[j].FastGetSolutionStepValue(FLOW));
                    //KRATOS_WATCH(geom[j].FastGetSolutionStepValue(NODAL_AREA));
            }

        }
		//KRATOS_WATCH("----------------FINISHED--------------------------------------------")

        KRATOS_CATCH("")
    }

    void ComputeRHS_Actualize(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;
        double dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
        for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
        {
            //int id = in->Id();
            //KRATOS_WATCH (id);
            double nodal_mass = in->FastGetSolutionStepValue(NODAL_MASS);
            //KRATOS_WATCH (nodal_mass)
            //array_1d<double,3>& pepe = in->FastGetSolutionStepValue(RHS);
            array_1d<double,3> aux = in->FastGetSolutionStepValue(RHS);
            //KRATOS_WATCH (id);
            //KRATOS_WATCH (aux);
            //aux /= nodal_mass;
            aux *= dt;
            aux /= nodal_mass;
            //aux *= dt;

            //KRATOS_WATCH (pepe)
            //int id = in->Id();
            //KRATOS_WATCH (id);
            //KRATOS_WATCH (aux);

            double& nodal_area = in->FastGetSolutionStepValue(NODAL_AREA);
            double& nodal_flow = in->FastGetSolutionStepValue(FLOW);
            //double lhs1=aux[0] / nodal_mass;
            //double lhs2=aux[1] / nodal_mass;
            //KRATOS_WATCH(nodal_area)
            //KRATOS_WATCH (nodal_area);
            //KRATOS_WATCH (nodal_flow);
            //NODAL_AREA is never prescribed
            if(in->IsFixed(NODAL_AREA) == false)
                nodal_area += aux[0] ;
                //KRATOS_WATCH(nodal_area)

            //FLOW
            //KRATOS_WATCH(nodal_flow)
            if(in->IsFixed(FLOW) == false)
                nodal_flow += aux[1];
                //KRATOS_WATCH(nodal_flow)
    //           nodal_area = in->FastGetSolutionStepValue(NODAL_AREA);
    //           nodal_flow = in->FastGetSolutionStepValue(FLOW);

            //in->FastGetSolutionStepValue(VELOCITY_X) = nodal_flow / nodal_area;
//            if (id == 1 || id == 100 || id == 200 || id == 300 || id == 400 || id == 500)
//            {
//                KRATOS_WATCH (id);
//                KRATOS_WATCH( in->FastGetSolutionStepValue(FLOW));
//                KRATOS_WATCH( in->FastGetSolutionStepValue(NODAL_AREA));
//            }
            //KRATOS_WATCH(nodal_area)
            //KRATOS_WATCH(FLOW)
        }

        KRATOS_CATCH("")                
    }


    void InitializeWorkArray( ModelPart& ThisModelPart )
    {
        //
        ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            array_1d<double,3>& work = in->FastGetSolutionStepValue(WORK);
            work[0] = in->FastGetSolutionStepValue(NODAL_AREA,1);
            work[1] = in->FastGetSolutionStepValue(FLOW,1);
        }

    }


//    void ActualizeWorkArray(ModelPart& ThisModelPart, double dt)
//    {
//        //
//        ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
//        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
//        {
//            noalias(in->FastGetSolutionStepValue(WORK)) += dt * in->FastGetSolutionStepValue(RHS);
//        }

//    }

//    //takes in account BC
//    void ActualizeSolution(ModelPart& ThisModelPart, double dt)
//    {
//        array_1d<double,3> aux;
//        ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
//        KRATOS_WATCH("-------------------------------------------------->>>>>>>>>>>>>>>>>>>><<<<");
//        //set to the value at the beginning of the step
//        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
//        {
//            if(in->IsFixed(NODAL_AREA) == false)
//                in->FastGetSolutionStepValue(NODAL_AREA) = in->FastGetSolutionStepValue(NODAL_AREA,1);

//            //FLOW
//            if(in->IsFixed(FLOW) == false)
//                in->FastGetSolutionStepValue(FLOW) = in->FastGetSolutionStepValue(FLOW,1);
//        }

//        //update the solution

//        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
//        {
//            double nodal_mass = in->GetSolutionStepValue(NODAL_MASS);
//            aux = in->FastGetSolutionStepValue(RHS);
//            aux *= dt;
//            // double nodal_area = in->FastGetSolutionStepValue(NODAL_AREA);
//            //double nodal_flow = in->FastGetSolutionStepValue(FLOW);

//            if(in->IsFixed(NODAL_AREA) == false)
//                in->FastGetSolutionStepValue(NODAL_AREA) += aux[0] / nodal_mass;

//            //FLOW
//            if(in->IsFixed(FLOW) == false)
//                in->FastGetSolutionStepValue(FLOW) += aux[1] / nodal_mass;

//            //nodal_area = in->FastGetSolutionStepValue(NODAL_AREA);
//            //nodal_flow = in->FastGetSolutionStepValue(FLOW);
//        }

//    }

    void CheckSolution(ModelPart& ThisModelPart, int Var_Aux)
    {
        KRATOS_TRY;
        for(ModelPart::ElementsContainerType::iterator i_element = ThisModelPart.ElementsBegin(); i_element !=ThisModelPart.ElementsEnd(); i_element++)
        {
            const double A0_actual = i_element->GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA);
            const double A1_actual = i_element->GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA);
            if((A0_actual == 0.00) || (A1_actual== 0.00))
             {
                KRATOS_WATCH(A0_actual);
                KRATOS_WATCH(A1_actual);
                KRATOS_WATCH(i_element->GetGeometry()[0]);
                KRATOS_WATCH(i_element->GetGeometry()[1]);
                KRATOS_WATCH(i_element->GetProperties().Id());
                KRATOS_WATCH(Var_Aux)
                //std::cout << "this element is" << this->i_element->Id() <<std::endl;
                KRATOS_ERROR(std::runtime_error, "Zero Nodal area found, Please check your model", "");
             }
        }
        KRATOS_CATCH("")
    }
};

}  // namespace Kratos.

#endif // KRATOS_ARTERY_TIME_INTEGRATOR_INCLUDED  defined


