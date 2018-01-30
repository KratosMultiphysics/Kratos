//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <cstring>

// External includes
#include "boost/smart_ptr.hpp"

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
	using namespace boost::python;

	typedef Condition                            ConditionBaseType;
	typedef Geometry<Node<3> >                        GeometryType;
	typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
	typedef GeometryType::PointsArrayType           NodesArrayType;

	void  AddCustomConditionsToPython()
	{
	    // class_< FaceForce3D, FaceForce3D::Pointer, bases< ConditionBaseType > >
	    // ("FaceForce3D",
	    //  init<int, GeometryType::Pointer, Properties::Pointer>() )
	    // ;
	}

    }  // namespace Python.

}  // namespace Kratos.
