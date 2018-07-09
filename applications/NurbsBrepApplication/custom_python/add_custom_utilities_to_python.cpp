//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/NurbsBrepApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Riccardo Rossi
//                   Thomas Oberbichler
//

// System includes

// External includes
#include <ANurbs/Core>

// Project includes
#include "includes/define_python.h"
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_utilities/NurbsBrepModeler.h"
#include "custom_utilities/BrepModelGeometryReader.h"

#include "custom_utilities/geometries/NodeCurveGeometry3D.h"
#include "custom_utilities/geometries/NodeSurfaceGeometry3D.h"

namespace Kratos
{

namespace Python
{

    void  AddCustomUtilitiesToPython(pybind11::module& m)
    {
        using namespace pybind11;
        using namespace pybind11::literals;

        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


        class_<NurbsBrepModeler, typename NurbsBrepModeler::Pointer>(m, "NurbsBrepModeler")
            .def(init<ModelPart::Pointer>())
            .def("LoadGeometry", &NurbsBrepModeler::LoadGeometry)
            .def("CreateIntegrationDomain", &NurbsBrepModeler::CreateIntegrationDomain)
            .def("ApplyGeometryRefinement", &NurbsBrepModeler::ApplyGeometryRefinement)
            .def("ComputeArea", &NurbsBrepModeler::ComputeArea)
            .def("MapNode", &NurbsBrepModeler::MapNode)
            .def("GetInterfaceConditions", &NurbsBrepModeler::GetInterfaceConditions)
            ;

        class_<BrepModelGeometryReader, typename BrepModelGeometryReader::Pointer>(m, "BrepModelGeometryReader")
            .def(init<Parameters&>())
            .def("ReadGeometry", &BrepModelGeometryReader::ReadGeometry)
            .def("WriteGaussPoints", &BrepModelGeometryReader::WriteGaussPoints)
            .def("WriteGaussPointsJson", &BrepModelGeometryReader::WriteGaussPointsJson)
            ;


        // --- register geometries

        // register Point3d
        {
            using Type = ANurbs::Point3D;

            pybind11::class_<Type>(m, "Point3D")
                .def_property_readonly("X", &Type::X)
                .def_property_readonly("Y", &Type::Y)
                .def_property_readonly("Z", &Type::Z)
            ;
        }

        // register Interval
        {
            using Type = ANurbs::Interval<double>;

            pybind11::class_<Type>(m, "Interval")
                .def(pybind11::init<double, double>())
                .def_property_readonly("T0", &Type::T0)
                .def_property_readonly("T1", &Type::T1)
            ;
        }

        // register NodalCurveGeometry
        {
            using Type = NodeCurveGeometry3D;
            using Pointer = std::shared_ptr<Type>;

            pybind11::class_<Type, Pointer>(m, "NodeCurveGeometry3D")
                .def(pybind11::init<int, int>(),
                    "degree"_a,
                    "nbNodes"_a)
                .def("Knot", &Type::Knot,
                    "index"_a)
                .def("SetKnot", &Type::SetKnot,
                    "index"_a,
                    "value"_a)
                .def("Knots", &Type::Knots)
                .def("Node", &Type::Node,
                    "index"_a)
                .def("SetNode", &Type::SetNode,
                    "index"_a,
                    "value"_a)
                .def("Pole", &Type::Pole,
                    "index"_a)
                .def("SetPole", &Type::SetPole,
                    "index"_a,
                    "value"_a)
                .def("Weight", &Type::Weight,
                    "index"_a)
                .def("SetWeight", &Type::SetWeight,
                    "index"_a,
                    "value"_a)
                .def("PointAt", &Type::PointAt,
                    "t"_a)
                .def("Spans", &Type::Spans)
                .def("ValueAt", &Type::ValueAt<double>,
                    "variable"_a,
                    "t"_a)
                .def("ValueAt2", &Type::ValueAt2<double>,
                    "variable"_a,
                    "t"_a,
                    "order"_a)
            ;
        }

        // register NodalSurfaceGeometry
        {
            using Type = NodeSurfaceGeometry3D;
            using Pointer = std::shared_ptr<Type>;

            pybind11::class_<Type, Pointer>(m, "NodeSurfaceGeometry3D")
                .def(pybind11::init<int, int, int, int>(),
                    "degreeU"_a,
                    "degreeV"_a,
                    "nbNodesU"_a,
                    "nbNodesV"_a)
                .def("SetKnotU", &Type::SetKnotU,
                    "index"_a,
                    "value"_a)
                .def("SetKnotV", &Type::SetKnotV,
                    "index"_a,
                    "value"_a)
                .def("Node", &Type::Node,
                    "indexU"_a,
                    "indexV"_a)
                .def("SetNode", &Type::SetNode,
                    "indexU"_a,
                    "indexV"_a,
                    "value"_a)
                .def("Pole", &Type::Pole,
                    "indexU"_a,
                    "indexV"_a)
                .def("SetPole", &Type::SetPole,
                    "indexU"_a,
                    "indexV"_a,
                    "value"_a)
                .def("Weight", &Type::Weight,
                    "indexU"_a,
                    "indexV"_a)
                .def("SetWeight", &Type::SetWeight,
                    "indexU"_a,
                    "indexV"_a,
                    "value"_a)
                .def("PointAt", &Type::PointAt,
                    "u"_a,
                    "v"_a)
            ;
        }
    }

}  // namespace Python.

} // Namespace Kratos
