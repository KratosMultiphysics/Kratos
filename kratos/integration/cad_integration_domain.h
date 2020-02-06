//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Pooyan Dadvand
//
//


#if !defined(KRATOS_CAD_INTEGRATION_DOMAIN_H_INCLUDED )
#define  KRATOS_CAD_INTEGRATION_DOMAIN_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// CadIntegrationDomain for the the CAD entities.
/* Provides functionalities and processes to
* create a numerical model from a CAD model.
*/
class CadIntegrationDomain
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename GeometryType::PointsArrayType PointsArrayType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Creates integration points and applies elements and conditions
    static void CreateIntegrationDomain(
        ModelPart& rModelPart,
        ModelPart& rCadModelPart,
        const Parameters& rPhysicsParameters,
        int EchoLevel = 0);

private:
    /// Goes through the list of elements and conditions
    static void CreateIntegrationDomainElementConditionList(
        ModelPart& rModelPart,
        ModelPart& rCadModelPart,
        const Parameters& rElementConditionListParameters,
        int EchoLevel = 0);

    /// Creates list of elements/ condition per geometry
    static void CreateIntegrationDomainElementCondition(
        ModelPart& rModelPart,
        ModelPart& rCadModelPart,
        const Parameters& rParameters,
        int EchoLevel = 0);

    /// Creates list of rQuadraturePointGeometryList
    static void CreateQuadraturePointGeometries(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rCadSubModelPart,
        const Parameters& rParameters,
        int EchoLevel = 0);

    /// Creates elements from the quadrature_point_geometries
    static void CreateElements(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rCadSubModelPart,
        std::string& rElementName,
        int& rIdCounter,
        int EchoLevel = 0);

    /// Creates conditions from the rQuadraturePointGeometryList
    static void CreateConditions(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rCadSubModelPart,
        std::string& rConditionName,
        int& rIdCounter,
        int EchoLevel = 0);

    /// Searches for nodes with geometrical criteria
    static void GetGeometryPointsAt(
        GeometriesArrayType& rGeometryList,
        ModelPart& rCadSubModelPart,
        const Parameters& rParameters,
        IndexType SpecificationType,
        int EchoLevel = 0);

    /// Gets list of geometries from rModelPart
    static void GetGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters& rParameters,
        int EchoLevel = 0);

    ///@}

}; // namespace CadIntegrationDomain

}  // namespace Kratos.

#endif // KRATOS_CAD_INTEGRATION_DOMAIN_H_INCLUDED  defined