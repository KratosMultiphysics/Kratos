/*
==============================================================================
KratosPFEMApplication 
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
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-11-19 15:38:01 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_EXPLICIT_DT_INCLUDED)
#define KRATOS_EXPLICIT_DT_INCLUDED


#include <string>
#include <iostream>
#include <algorithm>

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/triangle_2d_3.h"
#include "utilities/openmp_utils.h"

//#include "kratos/applications/MeshingApplication/meshing_application.h"

namespace Kratos
{ 
	class ExplicitDtProcess 
		: public Process 
	 {
	   public:

	      ExplicitDtProcess(double CFL,double min_dt,double max_dt, ModelPart& ThisModelPart )
			:Process(), cfl(CFL),Min_dt(min_dt),Max_dt(max_dt), mr_model_part(ThisModelPart)
		{
		}

	      /// Destructor.
	      virtual ~ExplicitDtProcess()
		{
		}
 

	      ///@}
	      ///@name Operators 
	      ///@{

	      void operator()()
		{
		  Execute();
		}
		
	    virtual void Execute()
		 {
		  KRATOS_TRY
		      
		    int NumThreads = OpenMPUtils::GetNumThreads();
		      std::vector< double > Threads_dt(NumThreads,10.0);

		      ModelPart::ElementsContainerType::iterator elem_bg = mr_model_part.ElementsBegin();
		      int n_elems = mr_model_part.Elements().size();	

		      #pragma omp parallel for firstprivate(n_elems, elem_bg)
		      for(int ii=0; ii<n_elems; ++ii)
			  {
			    //calculate min_dt
			    ModelPart::ElementsContainerType::iterator elem = elem_bg + ii;

			    double calc_dt = 1.0;
			    elem->Calculate(DELTA_TIME, calc_dt, mr_model_part.GetProcessInfo());

			    int k = OpenMPUtils::ThisThread();
			    if(calc_dt < Threads_dt[k])
				Threads_dt[k] = calc_dt;

			  }
		      #pragma omp barrier

	//KRATOS_WATCH(omp_get_thread_num());
	KRATOS_WATCH(NumThreads);

			  double DT = Max_dt;
			  for(int kk=0; kk<NumThreads; ++kk)
				if( Threads_dt[kk] < DT)
				      DT = Threads_dt[kk];


			    if(DT < Min_dt) DT = Min_dt;	

	// double DT = 0.00000001;
			    DT*=cfl;
			    mr_model_part.GetProcessInfo()[DELTA_TIME] = DT;
	KRATOS_WATCH("ExplicitDeltaT");
	KRATOS_WATCH(DT);
// 			return DT;

			KRATOS_WATCH("++++++++++++++++++++END OF ExplicitDtProcess PROCESS ^^^^^^^^^^^^^^^^^^^^^^");

		KRATOS_CATCH("")

		 }

		private:
                         double cfl,Min_dt,Max_dt;
                         ModelPart& mr_model_part;

	};

}//namespace kratos

#endif
