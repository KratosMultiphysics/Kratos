/*
==============================================================================
KratosConvectionDiffusionApplication 
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
//   Last modified by:    $Author: julio.marti $
//   Date:                $Date:  $
//   Revision:            $Revision:  $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/rad_face_utilities.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
	
namespace Python
{
/*void GenerateModelPart(RadFaceUtilities& RadFaceUtilities,ModelPart& origin_model_part,ModelPart& destination_model_part,unsigned int domain_size )
{
    if(domain_size == 2)
    {
        RadFaceUtilities.GenerateModelPart(origin_model_part, destination_model_part, KratosComponents<Element>::Get("ConvDiffr2D"),KratosComponents<Condition>::Get("RadFace2D")	);
    }
    else if(domain_size == 3)
    {
        RadFaceUtilities.GenerateModelPart(origin_model_part, destination_model_part,KratosComponents<Element>::Get("ConvDiffr3D"),KratosComponents<Condition>::Get("RadFace3D")	);
    }
}*/
void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;


    class_<RadFaceUtilities>("RadFaceUtilities", init<>())
    .def("ConditionModelPart",&RadFaceUtilities::ConditionModelPart)
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


}
}  // namespace Python.

} // Namespace Kratos

