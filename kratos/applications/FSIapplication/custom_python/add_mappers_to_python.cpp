//
//   Project Name:        Kratos
//   Last modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:42 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_mappers_to_python.h"
#include "custom_utilities/AdvancedNMPointsMapper.hpp"
#include "custom_utilities/InterfacePreprocess.h"

#include "custom_utilities/shared_points_mapper.h"
#include "includes/node.h"

namespace Kratos
{

namespace Python
{
void  AddMappersToPython()
{

    using namespace boost::python;

    class_<SharedPointsMapper >("SharedPointsMapper",
                                init< const ModelPart::NodesContainerType&, const  ModelPart::NodesContainerType&, double>())
    .def("ScalarMap",&SharedPointsMapper::ScalarMap)
    .def("InverseScalarMap",&SharedPointsMapper::InverseScalarMap)
    .def("VectorMap",&SharedPointsMapper::VectorMap)
    .def("InverseVectorMap",&SharedPointsMapper::InverseVectorMap)
    ;

    class_<AdvancedNMPointsMapper>("AdvancedNMPointsMapper", init<const ModelPart&, ModelPart&>())
    .def("FindNeighbours",&AdvancedNMPointsMapper::FindNeighbours)
    .def("ScalarToNormalVectorMap",&AdvancedNMPointsMapper::ScalarToNormalVectorMap)
    .def("NormalVectorToScalarMap",&AdvancedNMPointsMapper::NormalVectorToScalarMap)
    .def("ScalarMap",&AdvancedNMPointsMapper::ScalarMap)
    .def("VectorMap",&AdvancedNMPointsMapper::VectorMap)
    ;

    class_<InterfacePreprocess>("InterfacePreprocess", init<>())
    .def("GenerateTriangleInterfacePart",&InterfacePreprocess::GenerateTriangleInterfacePart)
    .def("GenerateLineInterfacePart",&InterfacePreprocess::GenerateLineInterfacePart)
    ;
}

}  // namespace Python.

} // Namespace Kratos

