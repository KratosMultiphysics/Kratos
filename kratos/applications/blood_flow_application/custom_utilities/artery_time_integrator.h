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
        const ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();

        //reset the acceleration on the nodes
        ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            in->FastGetSolutionStepValue(NODAL_MASS) = 0.0;
        }

        //compute the projection (first step
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                i!=ThisModelPart.ElementsEnd(); i++)
        {
       i->Initialize();
        }

        //compute the projection (first step
        for(ModelPart::ConditionIterator i = ThisModelPart.ConditionsBegin();
                i!=ThisModelPart.ConditionsEnd(); i++)
        {
       i->Initialize();
        }
    }


    void SolveStep(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

		double dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];

	InitializeWorkArray( ThisModelPart );

	//step1
	ComputeRHS( ThisModelPart );
        for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
        {
         double nodal_mass = in->GetSolutionStepValue(NODAL_MASS);

         array_1d<double,3> aux = in->FastGetSolutionStepValue(RHS);
         aux /= nodal_mass;
	    aux *= dt;

        int id = in->Id();
        double nodal_area = in->FastGetSolutionStepValue(NODAL_AREA);
        double nodal_flow = in->FastGetSolutionStepValue(FLOW);
        double lhs1=aux[0] / nodal_mass;
        double lhs2=aux[1] / nodal_mass;

        //NODAL_AREA is never prescribed
            if(in->IsFixed(NODAL_AREA) == false)
                in->FastGetSolutionStepValue(NODAL_AREA) += aux[0] ;

	    //FLOW
	    if(in->IsFixed(FLOW) == false)
          in->FastGetSolutionStepValue(FLOW) += aux[1];
        nodal_area = in->FastGetSolutionStepValue(NODAL_AREA);
         nodal_flow = in->FastGetSolutionStepValue(FLOW);

         in->FastGetSolutionStepValue(VELOCITY_X) = nodal_flow / nodal_area;

        }
//        std::cout << "area : ";
//        for(ModelPart::NodeIterator i = ThisModelPart.NodesBegin() ; i != ThisModelPart.NodesEnd() ; i++)
//        {
//            std::cout << i->FastGetSolutionStepValue(NODAL_AREA) << " , ";
//        }

//        std::cout << std::endl;
//        std::cout << "FLOW : ";
//        for(ModelPart::NodeIterator i = ThisModelPart.NodesBegin() ; i != ThisModelPart.NodesEnd() ; i++)
//        {
//            std::cout << i->FastGetSolutionStepValue(FLOW) << " , ";
//        }

//        std::cout << std::endl;

//        ActualizeWorkArray( ThisModelPart, dt/6.0 );
//	ActualizeSolution( ThisModelPart, dt*0.5 );
//
//		for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
//	  std::cout << in->Id() << " " << in->FastGetSolutionStepValue(RHS) << " " << in->FastGetSolutionStepValue(WORK) << " " << in->FastGetSolutionStepValue(FLOW) << " " << in->FastGetSolutionStepValue(NODAL_AREA) << std::endl;
//
//
//	//step2
//	ComputeRHS( ThisModelPart );
//	ActualizeWorkArray( ThisModelPart, dt/3.0 );
//	ActualizeSolution( ThisModelPart, dt*0.5 );
//
//	//step3
//	ComputeRHS( ThisModelPart );
//	ActualizeWorkArray( ThisModelPart, dt/3.0 );
//	ActualizeSolution( ThisModelPart, dt );
//
//	//step1
//	ComputeRHS( ThisModelPart );
//	ActualizeWorkArray( ThisModelPart, dt/6.0 );

	//assign work to solution
//       ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
//        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
//        {
//            array_1d<double,3>& aux = in->FastGetSolutionStepValue(WORK);
//
//	    //NODAL_AREA is never prescribed
//	    in->FastGetSolutionStepValue(NODAL_AREA) = aux[0];
//
//	    //FLOW
//	    if(in->IsFixed(FLOW) == false)
//	      in->FastGetSolutionStepValue(FLOW) = aux[1];
//        }

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

        //compute the projection (first step
	Vector rhs_el(4);

    for(ModelPart::ConditionsContainerType::iterator i = ThisModelPart.ConditionsBegin();
            i!=ThisModelPart.ConditionsEnd(); i++)
    {
        //KRATOS_WATCH(*i)
        i->CalculateRightHandSide(rhs_el,CurrentProcessInfo);

        Geometry< Node<3> >& geom = i->GetGeometry();

// 	    KRATOS_WATCH(geom);
        for(unsigned int j=0; j<geom.size(); j++)
        {
            int base = 2*j;
            array_1d<double,3>& rhs_node = geom[j].FastGetSolutionStepValue(RHS);

           rhs_node[0] += rhs_el[base];
           rhs_node[1] += rhs_el[base+1];
        }
    }

    for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
    {
    i->CalculateRightHandSide(rhs_el,CurrentProcessInfo);

// 	    KRATOS_WATCH(i->Id())
//
// 	    KRATOS_WATCH(rhs_el);


    Geometry< Node<3> >& geom = i->GetGeometry();

// 	    KRATOS_WATCH(geom);
    for(unsigned int j=0; j<geom.size(); j++)
    {
       int base = 2*j;
       array_1d<double,3>& rhs_node = geom[j].FastGetSolutionStepValue(RHS);

       rhs_node[0] += rhs_el[base];
       rhs_node[1] += rhs_el[base+1];
// 	       KRATOS_WATCH(rhs_node)
    }
}

//    std::cout << "RHS : ";
//    for(ModelPart::NodeIterator i = ThisModelPart.NodesBegin() ; i != ThisModelPart.NodesEnd() ; i++)
//    {
//        std::cout << i->FastGetSolutionStepValue(RHS)[0] / i->FastGetSolutionStepValue(NODAL_MASS) << " , " << i->FastGetSolutionStepValue(RHS)[1] / i->FastGetSolutionStepValue(NODAL_MASS) << ", ";
//    }

//    std::cout << std::endl;

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


    void ActualizeWorkArray(ModelPart& ThisModelPart, double dt)
    {
      //
      ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            noalias(in->FastGetSolutionStepValue(WORK)) += dt * in->FastGetSolutionStepValue(RHS);
        }

    }

    //takes in account BC
    void ActualizeSolution(ModelPart& ThisModelPart, double dt)
    {
      array_1d<double,3> aux;
	  ModelPart::NodesContainerType& rNodes = ThisModelPart.Nodes();

      //set to the value at the beginning of the step
      for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
	    if(in->IsFixed(NODAL_AREA) == false)
	      in->FastGetSolutionStepValue(NODAL_AREA) = in->FastGetSolutionStepValue(NODAL_AREA,1);

	    //FLOW
	    if(in->IsFixed(FLOW) == false)
	      in->FastGetSolutionStepValue(FLOW) = in->FastGetSolutionStepValue(FLOW,1);
        }

      //update the solution

        for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            double nodal_mass = in->GetSolutionStepValue(NODAL_MASS);
            aux = in->FastGetSolutionStepValue(RHS);
	    aux *= dt;
        double nodal_area = in->FastGetSolutionStepValue(NODAL_AREA);
        double nodal_flow = in->FastGetSolutionStepValue(FLOW);

	    if(in->IsFixed(NODAL_AREA) == false)
		in->FastGetSolutionStepValue(NODAL_AREA) += aux[0] / nodal_mass;

	    //FLOW
	    if(in->IsFixed(FLOW) == false)
          in->FastGetSolutionStepValue(FLOW) += aux[1] / nodal_mass;

        nodal_area = in->FastGetSolutionStepValue(NODAL_AREA);
         nodal_flow = in->FastGetSolutionStepValue(FLOW);
        }

    }

};

}  // namespace Kratos.

#endif // KRATOS_ARTERY_TIME_INTEGRATOR_INCLUDED  defined


