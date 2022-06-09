//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SPECIFICATIONS_UTILITIES)
#define KRATOS_SPECIFICATIONS_UTILITIES

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/constitutive_law.h"
#include "geometries/geometry_data.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

namespace
{
    template< class TContainerType>
    std::vector<std::string> GetDofsListFromGenericEntitiesSpecifications(const TContainerType& rContainer);

    // Definition of the map between the geometries in enum and string
    static std::unordered_map<std::string, GeometryData::KratosGeometryType> string_geometry_map = {
        {"generic_type",                          GeometryData::KratosGeometryType::Kratos_generic_type},
        {"Hexahedra3D20",                         GeometryData::KratosGeometryType::Kratos_Hexahedra3D20},
        {"Hexahedra3D27",                         GeometryData::KratosGeometryType::Kratos_Hexahedra3D27},
        {"Hexahedra3D8",                          GeometryData::KratosGeometryType::Kratos_Hexahedra3D8},
        {"Prism3D15",                             GeometryData::KratosGeometryType::Kratos_Prism3D15},
        {"Prism3D6",                              GeometryData::KratosGeometryType::Kratos_Prism3D6},
        {"Pyramid3D13",                           GeometryData::KratosGeometryType::Kratos_Pyramid3D13},
        {"Pyramid3D5",                            GeometryData::KratosGeometryType::Kratos_Pyramid3D5},
        {"Quadrilateral2D4",                      GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4},
        {"Quadrilateral2D8",                      GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8},
        {"Quadrilateral2D9",                      GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9},
        {"Quadrilateral3D4",                      GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4},
        {"Quadrilateral3D8",                      GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8},
        {"Quadrilateral3D9",                      GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9},
        {"Tetrahedra3D10",                        GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10},
        {"Tetrahedra3D4",                         GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4},
        {"Triangle2D3",                           GeometryData::KratosGeometryType::Kratos_Triangle2D3},
        {"Triangle2D6",                           GeometryData::KratosGeometryType::Kratos_Triangle2D6},
        {"Triangle3D3",                           GeometryData::KratosGeometryType::Kratos_Triangle3D3},
        {"Triangle3D6",                           GeometryData::KratosGeometryType::Kratos_Triangle3D6},
        {"Line2D2",                               GeometryData::KratosGeometryType::Kratos_Line2D2},
        {"Line2D3",                               GeometryData::KratosGeometryType::Kratos_Line2D3},
        {"Line3D2",                               GeometryData::KratosGeometryType::Kratos_Line3D2},
        {"Line3D3",                               GeometryData::KratosGeometryType::Kratos_Line3D3},
        {"Point2D",                               GeometryData::KratosGeometryType::Kratos_Point2D},
        {"Point3D",                               GeometryData::KratosGeometryType::Kratos_Point3D},
        {"Sphere3D1",                             GeometryData::KratosGeometryType::Kratos_Sphere3D1},
        {"NurbsCurve",                            GeometryData::KratosGeometryType::Kratos_Nurbs_Curve},
        {"NurbsSurface",                          GeometryData::KratosGeometryType::Kratos_Nurbs_Surface},
        {"NurbsVolume",                           GeometryData::KratosGeometryType::Kratos_Nurbs_Volume},
        {"NurbsCurveOnSurface",                   GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface},
        {"SurfaceInNurbsVolume",                  GeometryData::KratosGeometryType::Kratos_Surface_In_Nurbs_Volume},
        {"BrepCurve",                             GeometryData::KratosGeometryType::Kratos_Brep_Curve},
        {"BrepSurface",                           GeometryData::KratosGeometryType::Kratos_Brep_Surface},
        {"BrepCurveOnSurface",                    GeometryData::KratosGeometryType::Kratos_Brep_Curve_On_Surface},
        {"QuadraturePointGeometry",               GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Geometry},
        {"QuadraturePointCurveOnSurfaceGeometry", GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry},
        {"QuadraturePointSurfaceInVolumeGeometry",GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Surface_In_Volume_Geometry}
    };

    // Definition of the map between the dimension and integers
    static std::unordered_map<std::string, std::size_t> string_dimension_map = {
        {"2D",2},
        {"3D",3}
    };
}

/**
 * @class SpecificationsUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for evaluate specifications
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) SpecificationsUtilities
{
public:    
    ///@name Type Definitions
    ///@{

    /// Counted pointer of SpecificationsUtilities
    KRATOS_CLASS_POINTER_DEFINITION( SpecificationsUtilities );

    /**
     * @brief This enum defines a "hash" used to identify if implicit/explicit or static time integration is considered
     */
    enum class TimeIntegration
    {
        Static   = 0,
        Implicit = 1,
        Explicit = 2
    };

    /**
     * @brief This enum defines a "hash" used to identify if Lagrangian/Eulerian or ALE framework is considered
     */
    enum class Framework
    {
        Lagrangian = 0,
        Eulerian   = 1,
        ALE        = 2
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The default constructor
     */
    SpecificationsUtilities() = delete;
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method adds to the model part the missing variables
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    static void AddMissingVariables(ModelPart& rModelPart);

    /**
     * @brief This method adds to the model part the missing variables
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param EntitiesList List of entities check specifications
     */
    static void AddMissingVariablesFromEntitiesList(
        ModelPart& rModelPart,
        const Parameters EntitiesList
        );

    /**
     * @brief This method adds to the model part the missing variables from a given set of specifications
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param SpecificationsParameters The specification parameters
     * @param EntityName The name of the entity considered
     */
    static void AddMissingVariablesFromSpecifications(
        ModelPart& rModelPart,
        const Parameters SpecificationsParameters,
        const std::string EntityName = "NOT_DEFINED"
        );

    /**
     * @brief This method adds to the model part the missing dofs
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    static void AddMissingDofs(ModelPart& rModelPart);

    /**
     * @brief This method adds to the model part the missing dofs
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param EntitiesList List of entities check specifications
     */
    static void AddMissingDofsFromEntitiesList(
        ModelPart& rModelPart,
        const Parameters EntitiesList
        );

    /**
     * @brief This method adds to the model part the missing dofs from a given set of specifications
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param SpecificationsParameters The specification parameters
     * @param EntityName The name of the entity considered
     */
    static void AddMissingDofsFromSpecifications(
        ModelPart& rModelPart,
        const Parameters SpecificationsParameters,
        const std::string EntityName = "NOT_DEFINED"
        );

    /**
     * @brief This method gets dofs lists from specifications
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    static std::vector<std::string> GetDofsListFromSpecifications(const ModelPart& rModelPart);

    /**
     * @brief This method gets dofs lists from specifications (elements)
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    static std::vector<std::string> GetDofsListFromElementsSpecifications(const ModelPart& rModelPart);

    /**
     * @brief This method gets dofs lists from specifications (conditions)
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    static std::vector<std::string> GetDofsListFromConditionsSpecifications(const ModelPart& rModelPart);

    /**
     * @brief This method determine the flags used on the simulation
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    static void DetermineFlagsUsed(const ModelPart& rModelPart);

    /**
     * @brief This method detects the time integrations which are compatible. It throws a warning if incompatible time integration
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return The list of time integrations compatible
     */
    static std::vector<std::string> DetermineTimeIntegration(const ModelPart& rModelPart);

    /**
     * @brief This method detects the framework considered. It throws a warning if incompatible framework
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return The framework of the problem
     */
    static std::string DetermineFramework(const ModelPart& rModelPart);

    /**
     * @brief This method detects if the LHS is symmetric. It throws a warning if incompatible
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all the LHS are symmetric, false otherwise
     */
    static bool DetermineSymmetricLHS(const ModelPart& rModelPart);

    /**
     * @brief This method detects if the LHS is positive definite. It throws a warning if incompatible
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all the LHS are positive definite, false otherwise
     */
    static bool DeterminePositiveDefiniteLHS(const ModelPart& rModelPart);

    /**
     * @brief This method detects if the elements/conditions are compatible with its geometry
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all geometries are compatible, false otherwise
     */
    static bool DetermineIfCompatibleGeometries(const ModelPart& rModelPart);

    /**
     * @brief This method detects if all elements/conditions require time integration. It throws a warning if incompatible
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all the elements/conditions require time integration, false otherwise
     */
    static bool DetermineIfRequiresTimeIntegration(const ModelPart& rModelPart);

    /**
     * @brief This method detects if all elements/conditions are considering the proper CL
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all the elements/conditions are considering the adequate CL, false otherwise
     */
    static bool CheckCompatibleConstitutiveLaws(const ModelPart& rModelPart);

    /**
     * @brief This method detects if all elements/conditions are considering the proper geometrical polynomial degree
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return The highest geometrical polynomial degree
     */
    static int CheckGeometricalPolynomialDegree(const ModelPart& rModelPart);

    /**
     * @brief This method returns the documentation provided by the element/condition
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return Parameters containing a resume of all the documentation
     */
    static Parameters GetDocumention(const ModelPart& rModelPart);

    ///@}
}; // class SpecificationsUtilities
///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}
}  // namespace Kratos
#endif /* KRATOS_SPECIFICATIONS_UTILITIES defined */
