/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "iga_application_variables.h"

#include "spaces/ublas_space.h"
#include "custom_utilities/director_utilities.h"
#include "custom_utilities/iga_flags.h"
#include "geometries/nurbs_curve_geometry.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(
    pybind11::module& m)
{
    pybind11::class_< DirectorUtilities >(m, "DirectorUtilities")
        .def(pybind11::init<ModelPart&, Parameters>())
        .def("ComputeDirectors",
            &DirectorUtilities::ComputeDirectors)
        ;

    pybind11::class_< IgaFlags > iga_flags = pybind11::class_< IgaFlags >(m, "IgaFlags")
        .def(pybind11::init<>())
        ;

    iga_flags.attr("FIX_DISPLACEMENT_X") = IgaFlags::FIX_DISPLACEMENT_X;
    iga_flags.attr("FIX_DISPLACEMENT_Y") = IgaFlags::FIX_DISPLACEMENT_Y;
    iga_flags.attr("FIX_DISPLACEMENT_Z") = IgaFlags::FIX_DISPLACEMENT_Z;
    iga_flags.attr("FIX_ROTATION_X") = IgaFlags::FIX_ROTATION_X;
    iga_flags.attr("FIX_ROTATION_Y") = IgaFlags::FIX_ROTATION_Y;
    iga_flags.attr("FIX_ROTATION_Z") = IgaFlags::FIX_ROTATION_Z;

    m.def(
        "EvaluateNurbsCurveShapeFunctions",
        [](Geometry<Node>& rGeometry, const double Parameter) {

            using NurbsCurveType =
                NurbsCurveGeometry<3, PointerVector<Node>>;

            NurbsCurveType* p_nurbs_curve =
                dynamic_cast<NurbsCurveType*>(&rGeometry);

            NurbsCurveType::Pointer p_background_curve;

            if (p_nurbs_curve == nullptr) {
                p_background_curve =
                    std::dynamic_pointer_cast<NurbsCurveType>(
                        rGeometry.pGetGeometryPart(
                            Geometry<Node>::BACKGROUND_GEOMETRY_INDEX));

                KRATOS_ERROR_IF_NOT(p_background_curve)
                    << "Geometry does not contain a NurbsCurveGeometry."
                    << std::endl;

                p_nurbs_curve = p_background_curve.get();
            }

            array_1d<double, 3> local_coordinates =
                ZeroVector(3);

            local_coordinates[0] = Parameter;

            std::vector<std::size_t> control_point_ids;
            Vector shape_function_values;

            p_nurbs_curve->ShapeFunctionsValuesAndCPIndices(
                local_coordinates,
                control_point_ids,
                shape_function_values);

            return pybind11::make_tuple(
                control_point_ids,
                shape_function_values);
        });

    
    m.def(
        "GetNurbsCurveKnots",
        [](Geometry<Node>& rGeometry) {
            using NurbsCurveType =
                NurbsCurveGeometry<3, PointerVector<Node>>;

            NurbsCurveType* p_nurbs_curve =
                dynamic_cast<NurbsCurveType*>(&rGeometry);

            NurbsCurveType::Pointer p_background_curve;

            if (p_nurbs_curve == nullptr) {
                p_background_curve =
                    std::dynamic_pointer_cast<NurbsCurveType>(
                        rGeometry.pGetGeometryPart(
                            Geometry<Node>::BACKGROUND_GEOMETRY_INDEX));

                KRATOS_ERROR_IF_NOT(p_background_curve)
                    << "Geometry does not contain a NurbsCurveGeometry."
                    << std::endl;

                p_nurbs_curve = p_background_curve.get();
            }

            return p_nurbs_curve->Knots();
        });
}

} // namespace Python
} // Namespace Kratos
