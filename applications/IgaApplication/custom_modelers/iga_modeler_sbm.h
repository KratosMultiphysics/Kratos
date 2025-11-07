//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//                   Andrea Gorgi
//


# pragma once

// System includes

// External includes

// Project includes
#include "modeler/modeler.h"
#include "includes/properties.h"

#include "integration/integration_info.h"
#include "spatial_containers/bins_dynamic.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(IGA_APPLICATION) IgaModelerSbm
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(IgaModelerSbm);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;
    typedef typename Properties::Pointer PropertiesPointerType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;


    using PointType = Node;
    using PointTypePointer = Node::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;
    using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    using PointerType = DynamicBins::PointerType;
    using NodePointerVector = GlobalPointersVector<NodeType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IgaModelerSbm()
        : Modeler()
    {
    }

    /// Constructor.
    IgaModelerSbm(
        Model & rModel,
        const Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
    }

    /// Destructor.
    virtual ~IgaModelerSbm() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<IgaModelerSbm>(rModel, ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupModelPart() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "IgaModelerSbm";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream & rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream & rOStream) const override
    {
    }

    ///@}

private:
    ///@name Iga functionalities
    ///@{

    Model* mpModel = nullptr;

    ///@}
    ///@name Iga functionalities
    ///@{

    /// Creates integration point geometries and applies elements and conditions
    void CreateIntegrationDomain(
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    void CreateIntegrationDomainPerUnit(
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    /// Creates list of rQuadraturePointGeometryList
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters,
        std::string geometry_type) const;
    
    /// Creates list of rQuadraturePointGeometryList for Sbm
    void CreateQuadraturePointGeometriesSbm(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters,
        std::string geometry_type) const;

    
    /// Creates list of rQuadraturePointGeometryList for Sbm 
    /// using the condition name of the associated projection on the skin boundary layer
    void CreateQuadraturePointGeometriesSbmByProjectionLayer(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters,
        std::string geometry_type) const;

    /// Creates list of rQuadraturePointGeometryList for Sbm 
    /// using a fixed condition name for the whole surrogate boundary loop
    void CreateQuadraturePointGeometriesSbmByFixedConditionName(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters,
        std::string geometry_type,
        std::string ConditionName) const;

    ///@}
    ///@name CAD functionalities
    ///@{

    /// Gets list of geometries from rModelPart
    void GetGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    ///@}
    ///@name Generate Elements and Conditions
    ///@{

    /// Creates elements from geometries
    void CreateElements(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rElementName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        const Vector KnotSpanSizes) const;

    /// Creates conditions from geometries
    void CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        const Vector KnotSpanSizes) const;

    /// Creates conditions from geometries
    void CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        ModelPart& rSkinModelPart,
        std::vector<int>& listIdClosestCondition,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        bool IsInner,
        const Vector KnotSpanSizes) const;

    void PrepareIntegrationOnTrueBoundary(ModelPart& analysis_model_part) const;
    
    /// Creates conditions from geometries (from skin projection)
    void CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rModelPart,
        ModelPart& rSkinModelPart,
        std::vector<int>& rListIdClosestCondition,
        std::vector<int>& rListIdSecondClosestCondition,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        const bool IsInner,
        const Vector KnotSpanSizes) const;

    ///@}
    ///@name Utility
    ///@{

    void CreateConditionsWithAdditionalData(
        GeometriesArrayType& rGeometries,
        ModelPart& rModelPart,
        const Parameters rAdditionalData,
        std::string& rDefaultConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        const Vector KnotSpanSizes) const;


    ///@}
    ///@name Serializer
    ///@{

    friend class Serializer;

    ///@}
}; // Class CadModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    IgaModelerSbm& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const IgaModelerSbm& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.
