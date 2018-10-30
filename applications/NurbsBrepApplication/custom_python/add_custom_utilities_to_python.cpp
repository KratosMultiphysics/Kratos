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
            .def(init<ModelPart&>())
            .def("LoadGeometry", &NurbsBrepModeler::LoadGeometry)
            .def("CreateIntegrationDomain", &NurbsBrepModeler::CreateIntegrationDomain)
            .def("ApplyGeometryRefinement", &NurbsBrepModeler::ApplyGeometryRefinement)
            .def("ComputeArea", &NurbsBrepModeler::ComputeArea)
            .def("MapNode", &NurbsBrepModeler::MapNode)
            .def("GetInterfaceConditions", &NurbsBrepModeler::GetInterfaceConditions)
            .def("GetInterfaceConditionsAdvanced", &NurbsBrepModeler::GetInterfaceConditionsAdvanced)
            .def("GetUpdatedLocation", &NurbsBrepModeler::GetUpdatedLocation)
            .def("GetUpdatedLocationNewModelPart", &NurbsBrepModeler::GetUpdatedLocationNewModelPart)
            ;

        class_<BrepModelGeometryReader, typename BrepModelGeometryReader::Pointer>(m, "BrepModelGeometryReader")
            .def(init<Parameters&>())
            .def("ReadGeometry", &BrepModelGeometryReader::ReadGeometry)
            .def("WriteGaussPoints", &BrepModelGeometryReader::WriteGaussPoints)
            .def("WriteGaussPointsIteration", &BrepModelGeometryReader::WriteGaussPointsIteration)
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
                    "Degree"_a,
                    "NbNodes"_a)
                .def("Knot", &Type::Knot,
                    "Index"_a)
                .def("SetKnot", &Type::SetKnot,
                    "Index"_a,
                    "Value"_a)
                .def("Knots", &Type::Knots)
                .def("Node", &Type::Node,
                    "Index"_a)
                .def("SetNode", &Type::SetNode,
                    "Index"_a,
                    "Value"_a)
                .def("Pole", &Type::Pole,
                    "Index"_a)
                .def("SetPole", &Type::SetPole,
                    "Index"_a,
                    "Value"_a)
                .def("Weight", &Type::Weight,
                    "Index"_a)
                .def("SetWeight", &Type::SetWeight,
                    "Index"_a,
                    "Value"_a)
                .def("Spans", &Type::Spans)
                .def("PointAt", &Type::PointAt,
                    "T"_a)
                .def("DerivativesAt", &Type::DerivativesAt,
                    "T"_a,
                    "Order"_a)
                .def("ValueAt", &Type::ValueAt<double>,
                    "Variable"_a,
                    "T"_a)
                .def("ValueAt", &Type::ValueAt2<double>,
                    "Variable"_a,
                    "T"_a,
                    "Order"_a)
            ;
        }

        // register NodalSurfaceGeometry
        {
            using Type = NodeSurfaceGeometry3D;
            using Pointer = std::shared_ptr<Type>;

            pybind11::class_<Type, Pointer>(m, "NodeSurfaceGeometry3D")
                .def(pybind11::init<int, int, int, int>(),
                    "DegreeU"_a,
                    "DegreeV"_a,
                    "NbNodesU"_a,
                    "NbNodesV"_a)
                .def("SetKnotU", &Type::SetKnotU,
                    "Index"_a,
                    "Value"_a)
                .def("SetKnotV", &Type::SetKnotV,
                    "Index"_a,
                    "Value"_a)
                .def("Node", &Type::Node,
                    "IndexU"_a,
                    "IndexV"_a)
                .def("SetNode", &Type::SetNode,
                    "IndexU"_a,
                    "IndexV"_a,
                    "Value"_a)
                .def("Pole", &Type::Pole,
                    "IndexU"_a,
                    "IndexV"_a)
                .def("SetPole", &Type::SetPole,
                    "IndexU"_a,
                    "IndexV"_a,
                    "Value"_a)
                .def("Weight", &Type::Weight,
                    "IndexU"_a,
                    "IndexV"_a)
                .def("SetWeight", &Type::SetWeight,
                    "IndexU"_a,
                    "IndexV"_a,
                    "Value"_a)
                .def("PointAt", &Type::PointAt,
                    "U"_a,
                    "V"_a)
                .def("DerivativesAt", &Type::DerivativesAt,
                    "U"_a,
                    "V"_a,
                    "Order"_a)
                .def("ValueAt", &Type::ValueAt<double>,
                    "Variable"_a,
                    "U"_a,
                    "V"_a)
                .def("ValueAt", &Type::ValueAt2<double>,
                    "Variable"_a,
                    "U"_a,
                    "V"_a,
                    "Order"_a)
            ;
        }
    }

}  // namespace Python.

} // Namespace Kratos
