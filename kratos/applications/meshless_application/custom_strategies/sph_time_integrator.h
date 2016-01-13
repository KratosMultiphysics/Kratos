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


#if !defined(KRATOS_SPH_TIME_INTEGRATOR_INCLUDED )
#define  KRATOS_SPH_TIME_INTEGRATOR_INCLUDED



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
#include "meshless_application.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"



#include "geometries/point_3d.h"

#include "utilities/openmp_utils.h"
#include "custom_utilities/neighbours_calculator_SPH.h"
#include "custom_processes/node_and_element_erase_process.h"


#include "custom_elements/SPHparticle.h"
#include "meshless_application_variables.h"


namespace Kratos
{

//typedef  SPHparticle<KernelPoly6,KernelPoly6,KernelPoly6>                     ParticleType;
//typedef SPHparticle<KernelC2,KernelC2,KernelC2>                     ParticleType;

typedef  Neighbours_Calculator_SPH<Element> NeighboursCalculatorType;

typedef WeakPointerVector<Element> ParticleWeakVectorType;
typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;

typedef Node < 3 > PointType;
typedef PointerVector<PointType > PointVector;
typedef PointVector::iterator PointIterator;

typedef ModelPart::NodesContainerType::iterator NodeIterator;


class SPHTimeIntegrator
{
public:

    void Initialize(ModelPart& ThisModelPart)
    {
      //const ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();

        //need to compute neighbours
        NeighboursCalculatorType().Search_Neighbours(ThisModelPart);

        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            if (i->Id() == 1763){

                KRATOS_WATCH(i->GetValue(NEIGHBOUR_ELEMENTS).size());}
        }

        //compute the projection (first step
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            i->Initialize();

        }


        //DENSITY NORMALIZATION
//                        Normalize_density(ThisModelPart);

//                        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
//                            i!=ThisModelPart.ElementsEnd(); i++)
//                        {
//                            i->GetGeometry()(0)->GetValue(DENSITY) = i->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY);

//                        }



        ComputeRHS( ThisModelPart ); //here i compute a value RHS on every node


    }

    double EstimateDeltaTime(double time_step, ModelPart& ThisModelPart)
    {
        const ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();


        double dt = time_step;
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            double local_dt;
            i->Calculate(DELTA_TIME, local_dt,CurrentProcessInfo);
            if(local_dt < dt) { dt = local_dt;/* KRATOS_WATCH(i->Id());*/ }

        }
        return dt;
    }


    void Solve(ModelPart& ThisModelPart,const array_1d<double, 3 > & corner1, const array_1d<double, 3 > & corner2,
               ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY;
        ProcessInfo& process_info = ThisModelPart.GetProcessInfo();



        //Step 1 Update velocities
        double dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];

        for(ModelPart::NodesContainerType::iterator i = ThisModelPart.NodesBegin(); i!=ThisModelPart.NodesEnd(); i++)
        {
            //compute velocities
            if (i->FastGetSolutionStepValue(IS_STRUCTURE)==0.0){

                //const array_1d<double, 3 >& old_vel = i->FastGetSolutionStepValue(VELOCITY,1) ;
                array_1d<double, 3 >& vel = (i->FastGetSolutionStepValue(VELOCITY)) ;
                //array_1d<double, 3 >& disp = (i->FastGetSolutionStepValue(DISPLACEMENT)) ;
                const array_1d<double, 3 >& rhs = (i->FastGetSolutionStepValue(RHS)) ;

                noalias(vel) += dt*rhs;




            }
        }

        ////////////

        //Step 2 Apply velocity correction ( Obtain moving velocities )
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();i!=ThisModelPart.ElementsEnd(); i++){
            double tmp;
            i->Calculate(DUMMY_APPLY_XSPH,tmp,process_info);
        }


        //Step 3 Update positions
        for(ModelPart::NodesContainerType::iterator i = ThisModelPart.NodesBegin(); i!=ThisModelPart.NodesEnd(); i++)
        {
            //move_particles
            if (i->FastGetSolutionStepValue(IS_STRUCTURE)==0.0){
//                                                array_1d<double, 3 >& vel = (i->FastGetSolutionStepValue(VELOCITY)) ;

                array_1d<double, 3 >& vel = (i->FastGetSolutionStepValue(XSPH_VELOCITY)) ;
                array_1d<double, 3 >& disp = (i->FastGetSolutionStepValue(DISPLACEMENT)) ;

                //compute new positions

                //FORWARD EULER

                noalias(disp) += dt * vel  ;


                i->Coordinates() = i->GetInitialPosition();
                i->Coordinates() += disp;

            }
        }



        //Delete the particles out of the bounding box

        MarkOuterNodes(corner1,corner2,rNodes);

        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            if (i->GetGeometry()(0)->Is(TO_ERASE) == true ){
                ParticleWeakVectorType& neighbours = i->GetValue(NEIGHBOUR_ELEMENTS);

                for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
                {
                    neighbour_iterator->GetGeometry()(0)->Set(TO_ERASE, true) ;

                }
            }
        }

        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            if (i->GetGeometry()(0)->Is(TO_ERASE) == true ){
                ParticleWeakVectorType& neighbours = i->GetValue(NEIGHBOUR_ELEMENTS);

                for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
                {
                    neighbour_iterator->GetGeometry()(0)->Set(TO_ERASE, true) ;

                }
            }
        }

        NodeAndElementEraseProcess Eraser(ThisModelPart);
        Eraser.Execute();




        //Step 4 Get New neighbours
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();i!=ThisModelPart.ElementsEnd(); i++){
            i->GetValue(NEIGHBOUR_ELEMENTS).clear();}

        double t0 = OpenMPUtils::GetCurrentTime();
        mNeighbourCalculator.Search_Neighbours(ThisModelPart);

        double t1 = OpenMPUtils::GetCurrentTime();
        KRATOS_WATCH("**************");
        KRATOS_WATCH("NEIGHBOUR SEARCH TIME:");
        KRATOS_WATCH(t1-t0);
        KRATOS_WATCH("**************");


        //Step 5 Update Density and Pressure
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            double tmp;
            i->Calculate(DENSITY,tmp,process_info);
        }

//                        Normalize_density(ThisModelPart);



        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            double tmp;
            i->Calculate(PRESSURE,tmp,process_info);
        }

        // Assign pressure to the dummy particles (if they exist)
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            double tmp;
            i->Calculate(DUMMY_BOUNDARY_PRESSURES,tmp,process_info);
        }



        //Step 6 , Compute RHS for next step

        ComputeRHS( ThisModelPart ); //here i compute a value RHS on every node, to be used in next time step

        KRATOS_CATCH("");
    }


    void Normalize_density(ModelPart& ThisModelPart)
    {

        // DEnsity normalization, onyl to be used with collocative density approximation
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            ParticleWeakVectorType& neighbours = i->GetValue(NEIGHBOUR_ELEMENTS);

            double h=i->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS) ;
            double normalizer=0;


            for(ParticleWeakIteratorType neighbour_iterator =neighbours.begin(); neighbour_iterator != neighbours.end(); neighbour_iterator++)
            {
                double mass=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS) ;

                double density2=neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) ;
                array_1d<double,3> rvec = (neighbour_iterator->GetGeometry()[0].Coordinates());
                noalias(rvec) -= i->GetGeometry()(0)->Coordinates();
                double r = norm_2(rvec);

                double q;
                q=r/h;

                //Calculate the Value
                double Kernel=0.0;



//                if (q<0.0) {std::cout<<"WARNING: q is negative";}
//                else if (q>=0.0 && q<1.0) {Kernel=(7.0/(478.0*PI*h*h)) * (pow(3.0-q,5) - 6.0*pow(2.0-q,5) + 15.0*pow(1.0-q,5)); }
//                else if (q>=1.0 && q<2.0) {Kernel=(7.0/(478.0*PI*h*h)) * (pow(3.0-q,5) - 6.0*pow(2.0-q,5)); }
//                else if (q>=2.0 && q<3.0) {Kernel=(7.0/(478.0*PI*h*h)) * (pow(3.0-q,5)); }
//                else if (q>=3.0) {Kernel=0.0;}




                if (q<0.0) {std::cout<<"WARNING: q is negative";}
                else if (q>=0.0 && q<=2.0) {Kernel=(2.0/(PI*h*h)) * ((3.0/16.0)*q*q - 0.75*q + 0.75);}
                else if (q>2.0) {Kernel=0.0;}


                normalizer += (mass/density2) * Kernel;
            }

            //            if (normalizer > 1.0){
            //                KRATOS_WATCH(normalizer);
            //                KRATOS_WATCH(i->Id());
            //            }


            i->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY_NORM_PARAM)= normalizer;
        }



        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {


            i->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY) /= i->GetGeometry()(0)->FastGetSolutionStepValue(DENSITY_NORM_PARAM);


        }


    }


    void MarkOuterNodes(const array_1d<double, 3 > & corner1, const array_1d<double, 3 > & corner2,
                        ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY;

        //add a big number to the id of the nodes to be erased
        int n_erased = 0;
        double xmax, xmin, ymax, ymin, zmax, zmin;



        if (corner1[0] > corner2[0])
        {
            xmax = corner1[0];
            xmin = corner2[0];
        }
        else
        {
            xmax = corner2[0];
            xmin = corner1[0];
        }

        if (corner1[1] > corner2[1])
        {
            ymax = corner1[1];
            ymin = corner2[1];
        }
        else
        {
            ymax = corner2[1];
            ymin = corner1[1];
        }

        if (corner1[2] > corner2[2])
        {
            zmax = corner1[2];
            zmin = corner2[2];
        }
        else
        {
            zmax = corner2[2];
            zmin = corner1[2];
        }



        for (ModelPart::NodesContainerType::iterator in = rNodes.begin(); in != rNodes.end(); in++)
        {
            bool erase = false;
            double& x = in->X();
            double& y = in->Y();
            double& z = in->Z();

            if (x < xmin || x > xmax) erase = true;
            else if (y < ymin || y > ymax) erase = true;
            else if (z < zmin || z > zmax) erase = true;

            if (erase == true)
            {
                n_erased += 1;
                in->Set(TO_ERASE, true);
            }
        }

        KRATOS_WATCH(n_erased);

        KRATOS_CATCH("");
    }

private:

    NeighboursCalculatorType mNeighbourCalculator;

    void ComputeRHS(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;

        ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();

        // Assing 0 to the boundary forces of all particles
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            i->GetGeometry()(0)->FastGetSolutionStepValue(BOUNDARY_ACC) = ZeroVector(3);
        }


        Vector rhs_el(3);

        //Compute RHS of each element
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {
            i->CalculateRightHandSide(rhs_el,CurrentProcessInfo);
        }


        //Normalize RHS of each element
        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
            i!=ThisModelPart.ElementsEnd(); i++)
        {

            double tmp;
            //Normalize RHS
            i->Calculate(DUMMY_NORMALIZE_RHS,tmp,CurrentProcessInfo);


        }

        // Put the boundary acceleration if you want

//        for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
//            i!=ThisModelPart.ElementsEnd(); i++)
//        {

//            const array_1d<double,3>& bound_acc = i->GetGeometry()[0].FastGetSolutionStepValue(BOUNDARY_ACC);
//            array_1d<double,3>& rhs = i->GetGeometry()[0].FastGetSolutionStepValue(RHS);
//            noalias(rhs) += 2.5*bound_acc;

//        }


        KRATOS_CATCH("");
    }





};

}  // namespace Kratos.

#endif // KRATOS_SPH_TIME_INTEGRATOR_INCLUDED  defined


