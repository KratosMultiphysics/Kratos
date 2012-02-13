
/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
 
/* *********************************************************   
*          
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2011-03-29 11:41:31 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_JOINTS_ELEMENTS_INCLUDED)
#define  KRATOS_JOINTS_ELEMENTS_INCLUDED
//System includes
//External includes
#include "boost/smart_ptr.hpp"
#include <cmath>

//Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/sd_math_utils.h"

#include "includes/model_part.h"
#include "includes/mesh.h"
#include "geometries/geometry.h"
#include "includes/element.h"



namespace Kratos
{
        template<const std::size_t Tdim>
        class Joint
        {     
	  public:
	    
	   /// Pointer definition 
           KRATOS_CLASS_POINTER_DEFINITION(Joint);
	    
	   Joint()
	   {
	     mfail = false;
	     mNodes.resize(Tdim);
	   }
	   
	   Joint(const std::vector<Node<3>::Pointer>& rNodes)
	   {
	     mfail   = false; 
	     mNodes.resize(Tdim);
	     mNodes  = rNodes;
	     //mlength = Length();
	   }
	   
	   ~Joint(){}
	  
	   void InsertNode(const int& pos, const Node<3>::Pointer& rNode)
	   {
	     mNodes[pos] = rNode;  
	     //mlength     = Length();
	   
	   }

//           double GetLength()
//           {
// 	    return mlength;
// 	  }
	  
          double Length()
          {
	     array_1d<double, 2> Node_One;
	     array_1d<double, 2> Node_Two;
	     array_1d<double, 2> L;
	     Node_One[0]  = (*this)[0]->X0();  
	     Node_One[1]  = (*this)[0]->Y0();  
	     Node_Two[0]  = (*this)[1]->X0();  
	     Node_Two[1]  = (*this)[1]->Y0();
	     noalias(L)   = Node_One-Node_Two;
	     return std::sqrt(inner_prod(L,L ));
	  }
          
          bool IsFail()
          {
	    return mfail;
	  }
          
          void SetFail()
          {
	    mfail = true;
	  }
          
	  Node<3>::Pointer GetNode(const int& pos)
	   {
	     return mNodes[pos];
	   }
	   
	  Node<3>::Pointer& operator[] (const int& pos)
	  {
	    return  mNodes[pos];
	  } 
	   
	  /// Assignment operator.
	  Joint& operator=(const Joint& rOther)
	  { 
	    mfail  = rOther.mfail;
	    mNodes = rOther.mNodes;     
	    return *this;
	  }

	  /// Copy constructor.
	  Joint(const Joint& rOther)
	  {
	    *this =  rOther;
	  }
	  
	  private:
	    bool mfail; 
	    //double mlength;
	    std::vector<Node<3>::Pointer> mNodes; 
	};
	
}//namespace Kratos.

#endif /*KRATOS_DISCONNECT_TRIANGLES*/

