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
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2008-04-28 16:19:49 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "custom_utilities/trilinos_deactivation_utility.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/parallel_fill_communicator.h"
#include "custom_utilities/trilinos_refine_mesh.h"

namespace Kratos
{
    namespace Python
    {
        using namespace boost::python;
        
        void  AddCustomUtilitiesToPython()
        {
            class_<TrilinosDeactivationUtility, boost::noncopyable >
                    ("TrilinosDeactivationUtility",
                     init<>() )
                    .def("Deactivate", &TrilinosDeactivationUtility::Deactivate )
                    .def("Reactivate", &TrilinosDeactivationUtility::Reactivate )
                    .def("ReactivateStressFree", &TrilinosDeactivationUtility::ReactivateStressFree )
                    .def("ReactivateAll", &TrilinosDeactivationUtility::ReactivateAll )
                    .def("Initialize", &TrilinosDeactivationUtility::Initialize )
                    ;
		    
            class_<ParallelFillCommunicator, boost::noncopyable >
                    ("ParallelFillCommunicator",
                     init<ModelPart& >() )
                    .def("Execute", &ParallelFillCommunicator::Execute )
                    ;
		    
            class_<MPIRefineSimplexMesh, boost::noncopyable >
                    ("MPIRefineSimplexMesh",
                     init<ModelPart& >() )
                    .def("Local_Refine_Mesh", &MPIRefineSimplexMesh::Local_Refine_Mesh )
                    ;		    
		    
		    
		    
        }	
    }  // namespace Python.

} // Namespace Kratos

