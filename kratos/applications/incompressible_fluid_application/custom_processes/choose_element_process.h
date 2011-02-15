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
#if !defined(KRATOS_CHOOSE_ELEMENT_PROCESS_INCLUDED)
#define KRATOS_CHOOSE_ELEMENT_PROCESS_INCLUDED


#include <string>
#include <iostream>
#include <algorithm>

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/triangle_2d_3.h"
//#include "kratos/applications/MeshingApplication/meshing_application.h"

namespace Kratos
{ 
	class ChooseElementProcess 
		: public Process 
	 {
	   public:

	      ChooseElementProcess(ModelPart& model_part, unsigned int domain_size, char* water_element, char*  air_element)
			:Process(), mr_model_part(model_part), mdomain_size(domain_size), rElWater(KratosComponents<Element>::Get(water_element)), rElAir(KratosComponents<Element>::Get(air_element))
		{
		}

	      /// Destructor.
	      virtual ~ChooseElementProcess()
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
			KRATOS_WATCH("++++++++++++++++++++BEGIN OF CHOOSE PROCESS ^^^^^^^^^^^^^^^^^^^^^^");
			ModelPart::ElementsContainerType ElemPart;
			//KRATOS_WATCH(mr_model_part.Elements().size())

                        ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();
	                unsigned int n_elems = mr_model_part.Elements().size();	

			ElemPart.reserve(n_elems);		
			int water_num = 0;
			int air_num = 0;

			for(unsigned int ii = 0; ii< n_elems; ++ii)
			  {
                                ModelPart::ElementsContainerType::iterator Belem = el_begin + ii;
				Geometry< Node<3> >& geom = Belem->GetGeometry();
				
				Element::Pointer p_elem; 

				if(Belem->GetValue(IS_WATER_ELEMENT) == 1.0)
				  {
					p_elem = rElWater.Create(Belem->Id(),geom, Belem->pGetProperties());
					p_elem->GetValue(IS_WATER_ELEMENT) = 1.0;
					water_num++;						

				  }
			        else
				  {
					p_elem = rElAir.Create(Belem->Id(), geom, Belem->pGetProperties() );
					p_elem->GetValue(IS_WATER_ELEMENT) = 0.0;
					air_num++;

//KRATOS_WATCH("***********ASGSCompressible2D IS CHOSEN ******************");
				  }
                        ElemPart.push_back(p_elem);
				
			}

			mr_model_part.Elements() = ElemPart;

		 }

		private:
			ModelPart& mr_model_part;
			unsigned int mdomain_size;
			Element const& rElWater;
			Element const& rElAir;

	};

}//namespace kratos

#endif
