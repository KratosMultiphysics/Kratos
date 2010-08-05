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


#if !defined(KRATOS_MARK_FOR_REFINEMENT )
#define  KRATOS_MARK_FOR_REFINEMENT



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
#include "incompressible_fluid_application.h"

namespace Kratos
{
	class RefinementUtilities
	{
	public:
	  
	  
	  RefinementUtilities(){}
	  ~RefinementUtilities(){}
	  
		//**********************************************************************************************
		//**********************************************************************************************
  		void MarkForRefinement(Variable<double>& rVariable, ModelPart& ThisModelPart, double admissible_ratio, unsigned int max_levels)
		{			
			KRATOS_TRY;
			
			//compute reference subscale energy
			double tot_mass = 0.0;
			double tot_kin_energy = 0.0;
			for(ModelPart::NodesContainerType::iterator it=ThisModelPart.NodesBegin(); it!=ThisModelPart.NodesEnd(); it++)
			{
			    const array_1d<double,3>& v = it->FastGetSolutionStepValue(VELOCITY);
			    const double& m = it->FastGetSolutionStepValue(NODAL_MASS);
			    
			    tot_kin_energy += 0.5*m*inner_prod(v,v);
			    tot_mass += m;		  
			}
			double avg_energy = tot_kin_energy / tot_mass;
			
			
			//reset the nodal values
			for(ModelPart::NodesContainerType::iterator it=ThisModelPart.NodesBegin(); it!=ThisModelPart.NodesEnd(); it++)
			{
			    it->GetValue(REFINEMENT_LEVEL) = 0;
			}
			
			//mark elements for splitting depending on the desired error ratio
			unsigned int number_of_splitted_elements = 0;			
			for(ModelPart::ElementsContainerType::iterator it=ThisModelPart.ElementsBegin(); it!=ThisModelPart.ElementsEnd(); it++)
			{
			    //compute the error ratio
			    double element_energy = 0.0;
			    it->Calculate(rVariable, element_energy, ThisModelPart.GetProcessInfo());
			    double ratio = element_energy/avg_energy;
			    
			    
			    if(ratio > admissible_ratio && it->GetValue(REFINEMENT_LEVEL) <= max_levels )
			    {
			      //mark for splitting
			      it->GetValue(SPLIT_ELEMENT) = true;
			      int& current_level = it->GetValue(REFINEMENT_LEVEL);
			      current_level += 1;
			      
			      number_of_splitted_elements++;
			      
			      //mark all of the nodes with the refinement level
			      Geometry< Node<3> >& geom = it->GetGeometry();
			      
			      for(unsigned int i=0; i<geom.size(); i++)
			      {
				 if(geom[i].GetValue(REFINEMENT_LEVEL) < current_level)
				 {
				   geom[i].GetValue(REFINEMENT_LEVEL) = current_level;
				 }
			      }
			    }
			}
			
			bool is_ok = false;
			unsigned int nit = 0;
			while(is_ok==false && nit<2*max_levels)
			{
			    is_ok = true;
			    
			    ModelPart::ElementsContainerType aux_elem_list;
			    
			    //fill a list with all of the elements that have to be refined to obtain a correct gradient
			    for(ModelPart::ElementsContainerType::iterator it=ThisModelPart.ElementsBegin(); it!=ThisModelPart.ElementsEnd(); it++)
			    {						  
				//mark all of the nodes with the refinement level
				Geometry< Node<3> >& geom = it->GetGeometry();
				
				//determine for each element the maximum and minimum level of refinement (basically if the neighbours have been refined)
				unsigned int level = geom[0].GetValue(REFINEMENT_LEVEL);
// 				unsigned int min_level = level;
				unsigned int max_level = level;
				for(unsigned int i=1; i<geom.size(); i++)
				{
				   level = geom[i].GetValue(REFINEMENT_LEVEL);
// 				   if(level < min_level) min_level = level;
				   if(level > max_level) max_level = level;
				}
				
				const int& current_level = it->GetValue(REFINEMENT_LEVEL);
				//if there is a difference of level of refinement greater than 1, then refine to have a smoother grading of elements
				if(max_level - current_level > 1 && current_level < max_levels)
				{
				  //more iterations of the overall algorithm are needed
				  is_ok = false;
				  
				  aux_elem_list.push_back( *(it.base()) );
				}
			    }
			    
			    //now signal such elements for spltting and color their nodes correctly    
			    for(ModelPart::ElementsContainerType::iterator it=aux_elem_list.begin(); it!=aux_elem_list.end(); it++)
			    {
				  int& current_level = it->GetValue(REFINEMENT_LEVEL);
				  //mark for splitting
				  it->GetValue(SPLIT_ELEMENT) = true;
				  current_level += 1; //here we increase the level
				  
				  number_of_splitted_elements++;
				  
				  //mark all of the nodes with the refinement level
				  Geometry< Node<3> >& geom = it->GetGeometry();
				  
				  for(unsigned int i=0; i<geom.size(); i++)
				  {
				    if(geom[i].GetValue(REFINEMENT_LEVEL) < current_level)
				    {
				      geom[i].GetValue(REFINEMENT_LEVEL) = current_level;
				    }
				  }
			    }
			    
			    //increase iterations of the grading algorithm
			    nit += 1;
			}
			
			KRATOS_WATCH("***********************************************************");
			std::cout << "total number of refined elements = " << number_of_splitted_elements << std::endl;

			

			KRATOS_CATCH("")
		}

	private:


	};

}  // namespace Kratos.

#endif // KRATOS_MARK_FOR_REFINEMENT  defined 


