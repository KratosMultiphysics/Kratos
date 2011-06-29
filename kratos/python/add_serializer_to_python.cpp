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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>
  

// Project includes
#include "includes/define.h"  
#include "includes/serializer.h"
#include "python/add_serializer_to_python.h"
#include "includes/model_part.h"

namespace Kratos
{
namespace Python 
{
 using namespace boost::python;

 template< class TObjectType >
 void BufferPushBack(Buffer& rBuffer, TObjectType& rObject)
        {
            return rBuffer.push_back(rObject);
        }

 template< class TObjectType >
 void BufferPopFront(Buffer& rBuffer, TObjectType& rObject)
        {
            return rBuffer.pop_front(rObject);
        }

 template< class TObjectType >
 void SerializerSave(Serializer& rSerializer, std::string const & rName, TObjectType& rObject)
        {
            return rSerializer.save(rName, rObject);
        }

 template< class TObjectType >
 void SerializerLoad(Serializer& rSerializer, std::string const & rName, TObjectType& rObject)
        {
            return rSerializer.load(rName, rObject);
        }

 void SerializerPrint(Serializer& rSerializer)
 {
     std::cout << "Serializer buffer:";
     std::cout << ((std::stringstream*)(rSerializer.pGetBuffer()))->str();
 }

  void  AddSerializerToPython()
  {  
	  class_<Buffer, Buffer::Pointer >("Buffer")
		  .def(init<>())
		  .def(init<Buffer::SizeType>())
		  .def("Size",&Buffer::size)
		  .def("Swap",&Buffer::swap)
		  .def("Clear",&Buffer::clear)
//                  .def("Resize"&Buffer::resize)
	      .def(self_ns::str(self))
	 ;

	  class_<Serializer, Serializer::Pointer, boost::noncopyable >("Serializer")
		  .def(init<>())
		  .def(init<std::string const&>())
		  .def(init<Serializer::TraceType>())
		  .def(init<std::string const&, Serializer::TraceType>())
		  .def("Load",SerializerLoad<ModelPart>)
		  .def("Save",SerializerSave<ModelPart>)
                  .def("Print", SerializerPrint)
		  //.def("",&Kernel::Initialize)
//	      .def(self_ns::str(self))
	 ;

              enum_<Serializer::TraceType>("SerializerTraceType")
                .value("SERIALIZER_NO_TRACE", Serializer::SERIALIZER_NO_TRACE)
                .value("SERIALIZER_TRACE_ERROR", Serializer::SERIALIZER_TRACE_ERROR)
                .value("SERIALIZER_TRACE_ALL", Serializer::SERIALIZER_TRACE_ALL)
                .export_values()
               ;

 	}
	
}  // namespace Python.

} // Namespace Kratos

