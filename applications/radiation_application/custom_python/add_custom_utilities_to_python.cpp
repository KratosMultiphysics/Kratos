//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//

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

