//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


// System includes
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

// Application includes
#include "custom_python/add_custom_conditions_to_python.h"



namespace Kratos{
namespace Python{


    typedef Condition                            ConditionBaseType;
    typedef Geometry<Node<3> >                        GeometryType;
    typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
    typedef GeometryType::PointsArrayType           NodesArrayType;

    void AddCustomConditionsToPython(pybind11::module& m)
    {

    }

}  // namespace Python.
}  // namespace Kratos.
