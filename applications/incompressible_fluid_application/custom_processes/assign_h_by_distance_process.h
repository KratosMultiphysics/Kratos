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
#if !defined(KRATOS_ASSIGN_H_BY_DISTANCE_PROCESS_INCLUDED)
#define KRATOS_ASSIGN_H_BY_DISTANCE_PROCESS_INCLUDED


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
	class AssignHByDistanceProcess 
		: public Process 
	 {
	   public:

	      AssignHByDistanceProcess(ModelPart& model_part, double min_H, double sec_min_H, double ref_dist, double air_H  )
			:Process(), mr_model_part(model_part), mr_min_H(min_H), mr_sec_min_H(sec_min_H), mr_ref_dist(ref_dist), mr_air_H(air_H)
		{
		}

	      /// Destructor.
	      virtual ~AssignHByDistanceProcess()
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

            double slope = (mr_sec_min_H - mr_min_H)/mr_ref_dist;
	    
            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector NodePartition;
            OpenMPUtils::DivideInPartitions(mr_model_part.Nodes().size(),NumThreads,NodePartition);	    
	    
            #pragma omp parallel   
            {	      
	        int k = OpenMPUtils::ThisThread();
                ModelPart::NodeIterator NodesBegin = mr_model_part.NodesBegin() + NodePartition[k];
                ModelPart::NodeIterator NodesEnd = mr_model_part.NodesBegin() + NodePartition[k+1];

                for (ModelPart::NodeIterator nd = NodesBegin; nd != NodesEnd; nd++)
                {	      	      
			    const double current_dist = nd->FastGetSolutionStepValue(DISTANCE);
			    if(current_dist <= mr_ref_dist && current_dist>0.0)
			      nd->FastGetSolutionStepValue(NODAL_H) = mr_min_H + slope*current_dist;
			    
			    if(current_dist <= 0.0)
			      nd->FastGetSolutionStepValue(NODAL_H) = mr_air_H;			  
		      
		 }            


			    KRATOS_WATCH("++++++++++++++++++++END OF SaveElementByFlagProcess PROCESS ^^^^^^^^^^^^^^^^^^^^^^");
		    }
            }
		private:
			ModelPart& mr_model_part;
		        double mr_min_H;
		        double mr_sec_min_H;
		        double mr_ref_dist;	
		        double mr_air_H;				
	};

}//namespace kratos

#endif
