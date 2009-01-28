/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last modified by:    $Author: pooyan $
//   Date:                $Date: 2009-01-14 14:43:15 $
//   Revision:            $Revision: 1.3 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "python/add_vector_to_python.h"
#include "python/vector_python_interface.h"
#include "python/vector_scalar_operator_python.h"
#include "python/vector_vector_operator_python.h"

namespace Kratos
{
	
namespace Python
{
		
  using namespace boost::python;
	
  template<class TContainerType> 
  struct UblasVectorModifier
  {
    typedef typename TContainerType::size_type index_type;
    static void Resize(TContainerType& ThisContainer, typename TContainerType::size_type NewSize)
    {
      ThisContainer.resize(NewSize, true);
    }
    static void MoveSlice(TContainerType& ThisContainer, index_type Index, index_type From, index_type To)
    {
      if(Index > From)
	{
	  ThisContainer.resize(ThisContainer.size() + Index - From, true);
	  std::copy_backward(ThisContainer.begin() + From, ThisContainer.begin() + To, ThisContainer.begin() + Index + To - From);
	}
      else
	{
	  std::copy(ThisContainer.begin() + From, ThisContainer.begin() + To, ThisContainer.begin() + Index);
	  ThisContainer.resize(ThisContainer.size() + Index - From, true);
	}
    }
  };
	
 	
  void  AddVectorToPython()
  {

    ReadonlyVectorPythonInterface<zero_vector<double> >::CreateInterface("ZeroVector")
      .def(init<zero_vector<double>::size_type>())
//       .def(VectorScalarOperatorPython<zero_vector<double>, double, vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, unit_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, scalar_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, mapped_vector<double>, mapped_vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, compressed_vector<double>, compressed_vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, coordinate_vector<double>, coordinate_vector<double> >())
       ;

    ReadonlyVectorPythonInterface<unit_vector<double> >::CreateInterface("UnitVector")
      .def(init<unit_vector<double>::size_type, vector<double>::size_type>())
//       .def(VectorScalarOperatorPython<unit_vector<double>, double, vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, zero_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, scalar_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, mapped_vector<double>, mapped_vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, compressed_vector<double>, compressed_vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, coordinate_vector<double>, coordinate_vector<double> >())
       ;

    ReadonlyVectorPythonInterface<scalar_vector<double> >::CreateInterface("ScalarVector")
      .def(init<scalar_vector<double>::size_type, scalar_vector<double>::value_type>())
//       .def(VectorScalarOperatorPython<scalar_vector<double>, double, vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, zero_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, unit_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, mapped_vector<double>, mapped_vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, compressed_vector<double>, compressed_vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, coordinate_vector<double>, coordinate_vector<double> >())
       ;

    VectorPythonInterface<vector<double>, UblasVectorModifier<vector<double> > >::CreateInterface("Vector")
      .def(init<vector<double>::size_type>())
      .def(init<vector_expression<vector<double> > >())
      .def(VectorScalarOperatorPython<vector<double>, double, vector<double> >())
//       .def(VectorVectorOperatorPython<vector<double>, zero_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<vector<double>, unit_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<vector<double>, scalar_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<vector<double>, mapped_vector<double>, vector<double> >())
       ;
 	}
	
}  // namespace Python.

} // Namespace Kratos

