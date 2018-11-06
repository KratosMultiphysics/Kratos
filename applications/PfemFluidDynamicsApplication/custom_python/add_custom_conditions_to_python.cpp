//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes
#include <pybind11/pybind11.h>
//#include <pybind11/st1.h>
#include <cstring>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "python/add_mesh_to_python.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"

//Application includes
#include "custom_python/add_custom_conditions_to_python.h"

namespace Kratos
{

    namespace Python
    {
	namespace py = pybind11;

	typedef Condition                            ConditionBaseType;
	typedef Geometry<Node<3> >                        GeometryType;
	typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
	typedef GeometryType::PointsArrayType           NodesArrayType;

	void  AddCustomConditionsToPython(pybind11::module& m)
	{
	    // py::class_< FaceForce3D, FaceForce3D::Pointer, bases< ConditionBaseType > >
	    // ("FaceForce3D",
	    //  init<int, GeometryType::Pointer, Properties::Pointer>() )
	    // ;
	}

    }  // namespace Python.

}  // namespace Kratos.
