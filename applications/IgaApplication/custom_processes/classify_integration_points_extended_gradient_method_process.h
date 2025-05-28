//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Juan I. Camarotti

#if !defined(KRATOS_CLASSIFY_INTEGRATION_POINTS_EXTENDED_GRADIENT_METHOD_PROCESS_H_INCLUDED )
#define  KRATOS_CLASSIFY_INTEGRATION_POINTS_EXTENDED_GRADIENT_METHOD_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "utilities/function_parser_utility.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/geometry.h"
#include "geometries/point_3d.h"
#include "containers/array_1d.h" 

// Application includes
#include "custom_elements/laplacian_IGA_element.h"
#include "iga_application_variables.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/* @class ClassifyIntegrationPointsExtendedGradientMethod
 * @ingroup IgaApplication
 **/
class KRATOS_API(IGA_APPLICATION) ClassifyIntegrationPointsExtendedGradientMethodProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;
    typedef array_1d<double, 3> CoordinatesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;
    typedef typename GeometryType::PointsArrayType PointsArrayType;

    /// Pointer definition of ClassifyIntegrationPointsExtendedGradientMethodProcess
    KRATOS_CLASS_POINTER_DEFINITION(ClassifyIntegrationPointsExtendedGradientMethodProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ClassifyIntegrationPointsExtendedGradientMethodProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~ClassifyIntegrationPointsExtendedGradientMethodProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ExecuteBeforeSolutionLoop() override {
        Execute();
    };

    bool IsKnotSpanIntersected(
        const double u_min, const double u_max,
        const double v_min, const double v_max,
        Vector& unit_normal_vector);

    bool IsKnotSpanActive(
        const double u_min, const double u_max,
        const double v_min, const double v_max);

    bool LineIntersectsRectangle(
        const array_1d<double,3>& line_p1, 
        const array_1d<double,3>& line_p2,
        double u_min, double u_max, double v_min, double v_max);

    bool IsPointInsideRectangle(
        double x, double y, double u_min, 
        double u_max, double v_min, double v_max);

    bool DoSegmentsIntersect(
        const array_1d<double,3>& p1, const array_1d<double,3>& p2,
        const array_1d<double,3>& q1, const array_1d<double,3>& q2,
        array_1d<double,3>& line_knot_span_edge_intersection_point);

    Vector ComputeUnitNormal(
        const array_1d<double,3>& line_p1, 
        const array_1d<double,3>& line_p2);

    std::pair<bool, int> IsRectangleInsidePolygon(const array_1d<double, 3>& p1, const array_1d<double, 3>& p2,
                               const array_1d<double, 3>& p3, const array_1d<double, 3>& p4,
                               std::vector<ModelPart::ConditionType> polygon_segments);

    bool IsPointInsidePolygon(const array_1d<double, 3>& point, const std::vector<ModelPart::ConditionType>& polygon_segments, bool flag = false);

    void CreateBrepCurveOnSurfaceIntegrationPoints(std::shared_ptr<NurbsSurfaceGeometry<3, PointerVector<Node>>> p_nurbs_surface, 
                                               ModelPart& p_intersected_elements_sub_model_part, 
                                               ModelPart* p_loop);

    bool IsPointOnSegment(const array_1d<double, 3>& point, const array_1d<double, 3>& p1, const array_1d<double, 3>& p2, bool flag = false);


    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "background_model_part_name"                 : "please_specify_model_part_name",
            "skin_model_part_inner_name" : ["please_specify_skin_model_part_inner_name"],
            "skin_model_part_outer_name" : "please_specify_skin_model_part_outer_name",
            "number_of_inner_loops": 0
        })" );

        return default_parameters;
    }

    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ClassifyIntegrationPointsExtendedGradientMethodProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ClassifyIntegrationPointsExtendedGradientMethodProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Iga functionalities
    ///@{

    Model* mpModel = nullptr;
    Parameters mParameters;
    SizeType mEchoLevel;
    std::vector<array_1d<double, 3>> mOuterLoopPolygon;
    std::vector<ModelPart::ConditionType> mOuterLoopPolygonConditions;
    std::vector<std::vector<ModelPart::ConditionType>> mVectorInnerLoopsPolygonsConditions;

    ModelPart* mpBackgroundMeshModelPart = nullptr;
    ModelPart* mpSkinOuterLoopModelPart = nullptr;
    std::vector<ModelPart*> mpVectorSkinInnerLoops;
    IndexType mNumberOfInnerLoops = 0;

    ///@}
    ///@name Iga functionalities
    ///@{

    ///@}
    ///@}

    ///@}
    ///@name Utility
    ///@{


    ///@}
    ///@name Input and output
    ///@{


   

    ///@}

}; // Class ClassifyIntegrationPointsExtendedGradientMethodProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ClassifyIntegrationPointsExtendedGradientMethodProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ClassifyIntegrationPointsExtendedGradientMethodProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_CLASSIFY_INTEGRATION_POINTS_EXTENDED_GRADIENT_METHOD_PROCESS_H_INCLUDED 
