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
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_ESTIMATE_TIME_STEP )
#define  KRATOS_ESTIMATE_TIME_STEP



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
#include "geometries/point.h"
#include "thermo_mechanical_application.h"
// #include "custom_conditions/environment_contact.h"
//#include "includes/variables.h"
#include "utilities/openmp_utils.h"




namespace Kratos
{
 	
	template<unsigned int TDim>
	class EstimateTimeStep
	{
	public:

		//**********************************************************************************************
		//**********************************************************************************************
		//

		double ComputeDt(ModelPart& ThisModelPart, const double dist_max, const double CFL, const double dt_min ,const double dt_max  )
		{			
		  KRATOS_TRY

        const unsigned int NumNodes = TDim +1;

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(ThisModelPart.NumberOfElements(),NumThreads,ElementPartition);

        std::vector<double> MaxProj(NumThreads,0.0);

        #pragma omp parallel shared(MaxProj)
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementIterator ElemBegin = ThisModelPart.ElementsBegin() + ElementPartition[k];
            ModelPart::ElementIterator ElemEnd = ThisModelPart.ElementsBegin() + ElementPartition[k+1];

            double& rMaxProj = MaxProj[k];

            double Area;
            array_1d<double, NumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, NumNodes, TDim> DN_DX;

            for( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
            {
                // Get the element's geometric parameters
                Geometry< Node<3> >& rGeom = itElem->GetGeometry();
			    double ele_dist = 0.0;
			    for (unsigned int kk = 0; kk < rGeom.size(); kk++){
				     double dist = rGeom[kk].FastGetSolutionStepValue(DISTANCE);
				  ele_dist +=  fabs(dist);
			    }

			    if(ele_dist <= dist_max){
                GeometryUtils::CalculateGeometryData(rGeom, DN_DX, N, Area);

                // Elemental Velocity
                array_1d<double,3> ElementVel = N[0]*itElem->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
                for (unsigned int i = 1; i < NumNodes; ++i)
                    ElementVel += N[i]*rGeom[i].FastGetSolutionStepValue(VELOCITY);

                // Velocity norm
                double VelNorm = ElementVel[0]*ElementVel[0];
                for (unsigned int d = 1; d < TDim; ++d)
                    VelNorm += ElementVel[d]*ElementVel[d];
                VelNorm = sqrt(VelNorm);

                // Maximum element size along the direction of velocity
                for (unsigned int i = 0; i < NumNodes; ++i)
                {
                    double Proj = 0.0;
                    for (unsigned int d = 0; d < TDim; ++d)
                        Proj += ElementVel[d]*DN_DX(i,d);
                    Proj = fabs(Proj);
                    if (Proj > rMaxProj) rMaxProj = Proj;
                }
			  }
            }
        }

        // Obtain the maximum projected element size (compare thread results)
        double Max = 0.0;
        for (int k = 0; k < NumThreads; ++k)
            if (Max < MaxProj[k]) Max = MaxProj[k];

        // Dt to obtain desired CFL
        double dt = CFL / Max;
        if(dt > dt_max)
            dt = dt_max;
		else if(dt < dt_min)
			dt = dt_min;

        //perform mpi sync if needed
        double global_dt = dt;
        ThisModelPart.GetCommunicator().MinAll(global_dt);
        dt = global_dt;

        return dt;

		  //array_1d<double, 3 > dx, dv;
		  //double deltatime = dt_max;
		  //double dvside, lside;

		  //for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
			 // i != ThisModelPart.ElementsEnd(); i++)
		  //{
		  //    //calculating shape functions values
		  //    Geometry< Node < 3 > >& geom = i->GetGeometry();
			 // double ele_dist = 0.0;
			 // for (unsigned int kk = 0; kk < geom.size(); kk++){
				//     double dist = geom[kk].FastGetSolutionStepValue(DISTANCE);
				//  ele_dist +=  fabs(dist);
			 // }

			 // if(ele_dist <= dist_max){
		  //    for (unsigned int i1 = 0; i1 < geom.size() - 1; i1++)
		  //    {
			 // for (unsigned int i2 = i1 + 1; i2 < geom.size(); i2++)
			 // {
			 //     dx[0] = geom[i2].X() - geom[i1].X();
			 //     dx[1] = geom[i2].Y() - geom[i1].Y();
			 //     dx[2] = geom[i2].Z() - geom[i1].Z();

			 //     lside = inner_prod(dx, dx);

			 //     noalias(dv) = geom[i2].FastGetSolutionStepValue(VELOCITY);
			 //     noalias(dv) -= geom[i1].FastGetSolutionStepValue(VELOCITY);

			 //     dvside = inner_prod(dx, dv);

			 //     double dt;
			 //     if (dvside < 0.0) //otherwise the side is "getting bigger" so the time step has not to be diminished
			 //     {
				//  dt = fabs(lside / dvside);
				//  if (dt < deltatime) deltatime = dt;
			 //     }

			 // }
		  //    }
			 //}
		  //}

		  //if (deltatime < dt_min)
		  //{
		  //    std::cout << "ATTENTION dt_min is being used" << std::endl;
		  //    deltatime = dt_min;
		  //}

		  //return deltatime;

		  KRATOS_CATCH("")	
		}



	private:


	};

}  // namespace Kratos.

#endif // ASSIGN_NO_SLIP_CONDITION  defined 


