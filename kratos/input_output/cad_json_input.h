#if !defined(KRATOS_CAD_JSON_INPUT_H_INCLUDED )
#define  KRATOS_CAD_JSON_INPUT_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/io.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/node.h"

// Geometries
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/brep_surface.h"
#include "geometries/geometry.h"

namespace Kratos
{

  ///@name Kratos Classes
  ///@{
  /// Short class definition.
  /** Gives IO capabilities for Nurbs based Brep models in the JSON format defined in 
  https://amses-journal.springeropen.com/articles/10.1186/s40323-018-0109-4.
  */
    template<class TNodeType, class TEmbeddedNodeType>
    class KRATOS_API(KRATOS_CORE) CadJsonInput : public IO
    {
    public:
        /// Pointer definition of CadJsonInput
        KRATOS_CLASS_POINTER_DEFINITION(CadJsonInput);

        /// Constructor.
        CadJsonInput(
            ModelPart & rModelPart,
            Parameters & rNurbsBrepGeometryJson,
            const int EchoLevel = 0)
            : mEchoLevel(EchoLevel)
        {};

        /// Destructor.
        virtual ~CadJsonInput() {};

    private:

        Vector ReadControlPointWeightVector(
            const Parameters& rParameters)
        {
            Vector control_point_weights = ZeroVector(rParameters.size());
            KRATOS_ERROR_IF(rParameters.size() == 0)
                << "Length of control point list is zero!" << std::endl;
            KRATOS_ERROR_IF(rParameters[0].size() != 4)
                << "Control points need to be provided in following structure: [[x, y, z, weight]] or [id, [x, y, z, weight]]"
                << std::endl;
                << "Size of inner vector incorrect!"
                << std::endl;

            SizeType number_of_entries = rParameters[0].size();
            KRATOS_ERROR_IF(number_of_entries != 1 || number_of_entries != 2)
                << "Control points need to be provided in following structure: [[x, y, z, weight]] or [id, [x, y, z, weight]]"
                << std::endl;

            for (IndexType cp_idx = 0; cp_idx < rParameters.size(); cp_idx++)
            {
                control_point_weights[cp_idx] = rParameters[cp_idx][number_of_entries - 1][3].GetDouble();
            }

            return control_point_weights;
        }

        template<class TThisNodeType>
        PointerVector<TThisNodeType> ReadControlPointVector(
            const Parameters& rParameters)
        {
            PointerVector<TThisNodeType> control_points(rParameters.size());

            for (IndexType cp_idx = 0; cp_idx < rParameters.size(); cp_idx++)
            {
                control_points[cp_idx] = rParameters[cp_idx][number_of_entries - 1][3].GetDouble();
            }

            return control_points;
        }

        template<class TThisNodeType>
        TThisNodeType ReadNode(
            const Parameters& rParameters)
        {
            SizeType number_of_entries = rParameters[0].size();
            KRATOS_ERROR_IF(number_of_entries != 1 || number_of_entries != 2)
                << "Control points need to be provided in following structure: [[x, y, z, weight]] or [id, [x, y, z, weight]]"
                << std::endl;
            if (number_of_entries == 1)
            {
                Vector cp = rParameters[0].GetVector();

                return TThisNodeType(0, cp[0], cp[1], cp[2]);
            }
            else
            {
                SizeType id = rParameters[0].GetInt();
                Vector cp = rParameters[1].GetVector();

                return TThisNodeType(id, cp[0], cp[1], cp[2]);
            }
        }

        template<int TWorkingSpaceDimension, class TThisNodeType>
        void ReadNurbsCurve(
            const Parameters& rParameters)
        {
            KRATOS_ERROR_IF_NOT(rParameters.Has("curve_direction"))
                << "Missing 'curve_direction' in nurbs curve" << std::endl;
            bool curve_direction = rTrimmingCurve["curve_direction"].GetBool();

            KRATOS_ERROR_IF_NOT(rParameters.Has("parameter_curve"))
                << "Missing 'parameter_curve' in nurbs curve" << std::endl;

            KRATOS_ERROR_IF_NOT(rParameters["parameter_curve"].Has("is_rational"))
                << "Missing 'is_rational' in nurbs parameter curve" << std::endl;
            bool is_rational = rTrimmingCurve["parameter_curve"]["is_rational"].GetBool();

            KRATOS_ERROR_IF_NOT(rParameters["parameter_curve"].Has("knot_vector"))
                << "Missing 'knot_vector' in nurbs parameter curve" << std::endl;
            Vector knot_vector = rTrimmingCurve["parameter_curve"]["knot_vector"].GetVector();

            KRATOS_ERROR_IF_NOT(rParameters["parameter_curve"].Has("degree"))
                << "Missing 'degree' in nurbs parameter curve" << std::endl;
            int polynomial_degree = rTrimmingCurve["parameter_curve"]["degree"].GetInt();

            PointerVector<TThisNodeType>> control_points = ReadControlPointVector<TThisNodeType>(
                rParameters["parameter_curve"]["control_points"]);

            Vector control_point_weights = ReadControlPointVector(
                rParameters["parameter_curve"]["control_points"]);

            NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>> nurbs_curve(
                control_points,
                polynomial_degree,
                knot_vector,
                control_point_weights);
        }

        template<int TWorkingSpaceDimension, class TThisNodeType>
        void ReadNurbsSurface(
            const Parameters& rParameters)
        {
            bool is_trimmed = brep_json["is_trimmed"].GetBool();
            bool is_rational = brep_json["is_rational"].GetBool();

            Vector knot_vector_u = brep_json["knot_vectors"][0].GetVector();
            Vector knot_vector_v = brep_json["knot_vectors"][1].GetVector();

            int p = brep_json["faces"][i]["surface"]["degrees"][0].GetInt();
            int q = brep_json["faces"][i]["surface"]["degrees"][1].GetInt();

            PointerVector<TThisNodeType >> control_points = ReadControlPointVector<TThisNodeType>(
                rParameters["control_points"]);

            Vector control_point_weights = ReadControlPointVector(
                rParameters["control_points"]);

            NurbsSurfaceGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>> nurbs_surface(
                control_points,
                p,
                q,
                knot_vector_u,
                knot_vector_v,
                control_point_weights);
        }

        template<class TThisNodeType>
        void ReadBoundaryLoops(
            const Parameters& rParameters)
        {

        }

        void ReadBrepSurfaces(
            const Parameters& rParameters)
        {
            for (int i = 0; i < rParameters.size(); i++)
            {
                KRATOS_ERROR_IF_NOT(rParameters.Has("brep_id") || rParameters.Has("brep_name"))
                    << "Missing 'brep_id' or 'brep_name' in brep face" << std::endl;

                KRATOS_ERROR_IF_NOT(rParameters.Has("surface"))
                    << "Missing 'surface' in brep face" << std::endl;

                ReadNurbsSurface<3, TNodeType>(rParameters["surface"]);

                KRATOS_ERROR_IF_NOT(rParameters.Has("boundary_loops"))
                    << "Missing 'boundary_loops' in brep face" << std::endl;
                ReadBoundaryLoops<TEmbeddedNodeType>(rParameters["boundary_loops"]);

                BrepSurface<PointerVector<TNodeType>, PointerVector<TEmbeddedNodeType>> brep_surface(
                    control_points,
                    polynomial_degree,
                    knot_vector,
                    control_point_weights);

                if (rParameters.Has("brep_id"))
                    brep_surface.SetId(rParameters["brep_id"].GetInt());
                if (rParameters.Has("brep_name"))
                    brep_surface.SetId(rParameters["brep_name"].GetString());
            }
        }

        void ReadBrepCurveOnSurfavce(
            const Parameters& rParameters)
        {

        }

        template<int TWorkingSpaceDimension, class TThisNodeType>
        void ReadBrep(
            const Parameters& rParameters)
        {
            if (rParameters.Has("faces"))
                ReadBrepSurfaces(rParameters["faces"]);

            if (rParameters.Has("edges"))
                ReadBrepCurveOnSurfavce(rParameters["edges"]);
        }

        int mEchoLevel;
    }; // Class CadJsonInput
}  // namespace Kratos.

#endif // KRATOS_CAD_JSON_INPUT_H_INCLUDED  defined