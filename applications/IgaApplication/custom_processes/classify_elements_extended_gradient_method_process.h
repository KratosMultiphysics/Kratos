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

#if !defined(KRATOS_CLASSIFY_ELEMENTS_EXTENDED_GRADIENT_METHOD_PROCESS_H_INCLUDED )
#define  KRATOS_CLASSIFY_ELEMENTS_EXTENDED_GRADIENT_METHOD_PROCESS_H_INCLUDED

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
#include "spatial_containers/spatial_containers.h" 

// Application includes
#include "custom_utilities/entity_point.h"
#include "custom_elements/laplacian_IGA_element.h"
#include "iga_application_variables.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/* @class ClassifyElementsExtendedGradientMethod
 * @ingroup IgaApplication
 **/
class KRATOS_API(IGA_APPLICATION) ClassifyElementsExtendedGradientMethodProcess
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
    
    using EntityType = ModelPart::ElementsContainerType::value_type;
    using EntityPointType = EntityPoint<EntityType>;
    using EntityPointVector = std::vector<typename EntityPointType::Pointer>;

    // Type definitions for tree-search
    using BucketType = Bucket<3, EntityPoint<EntityType>, EntityPointVector>;
    using KDTree = Tree<KDTreePartition<BucketType>>;


    /// Pointer definition of ClassifyElementsExtendedGradientMethodProcess
    KRATOS_CLASS_POINTER_DEFINITION(ClassifyElementsExtendedGradientMethodProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ClassifyElementsExtendedGradientMethodProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~ClassifyElementsExtendedGradientMethodProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ExecuteBeforeSolutionLoop() override {
        Execute();
    };

    bool LineIntersectsRectangle(
        const array_1d<double,3>& line_p1, 
        const array_1d<double,3>& line_p2,
        double u_min, double u_max, double v_min, double v_max);

    bool DoSegmentsIntersect(
        const array_1d<double,3>& p1, const array_1d<double,3>& p2,
        const array_1d<double,3>& q1, const array_1d<double,3>& q2);

    bool IsPointInsideRectangle(double x, double y, double u_min, double u_max, double v_min, double v_max);

    std::pair<bool, int> IsRectangleInsidePolygon(const array_1d<double, 3>& p1, const array_1d<double, 3>& p2,
                               const array_1d<double, 3>& p3, const array_1d<double, 3>& p4,
                               std::vector<array_1d<double, 3>> polygon_vertices);

    bool IsPointInsidePolygon(const array_1d<double, 3>& point, const std::vector<array_1d<double, 3>>& polygon_vertices);

    Vector ComputeUnitNormal(
        const array_1d<double,3>& line_p1, 
        const array_1d<double,3>& line_p2);

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "background_model_part_name"                 : "please_specify_model_part_name",
            "embedded_body_model_part_name"              : "please_specify_model_part_name",
            "keep_external_domain"                       : false
        })" );

        return default_parameters;
    }

    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ClassifyElementsExtendedGradientMethodProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ClassifyElementsExtendedGradientMethodProcess";
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
    EntityPointVector mpEntityPointVectorKDTree;

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

}; // Class ClassifyElementsExtendedGradientMethodProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ClassifyElementsExtendedGradientMethodProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ClassifyElementsExtendedGradientMethodProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_CLASSIFY_ELEMENTS_EXTENDED_GRADIENT_METHOD_PROCESS_H_INCLUDED 
