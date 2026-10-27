//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "modeler/modeler.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) IgaModelerGsbm
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(IgaModelerGsbm);

    using SizeType = std::size_t;
    using IndexType = std::size_t;
    using GeometryType = Geometry<Node>;
    using GeometriesArrayType = GeometryType::GeometriesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IgaModelerGsbm()
        : Modeler()
    {
    }

    /// Constructor.
    IgaModelerGsbm(
        Model& rModel,
        const Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
        mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    }

    /// Destructor.
    ~IgaModelerGsbm() override = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<IgaModelerGsbm>(rModel, ModelParameters);
    }

    // Get the default parameters
    const Parameters GetDefaultParameters() const override;

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
        return "IgaModelerGsbm";
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr;

    ///@}
    ///@name Generate Integration Domain
    ///@{

    /// Creates integration point geometries and applies elements and conditions
    void CreateIntegrationDomain(
        ModelPart& rModelPart,
        const Parameters rPhysicsParameters) const;

    void CreateIntegrationDomainPerUnit(
        ModelPart& rModelPart,
        const Parameters rPhysicsParameters) const;

    /// Creates list of rQuadraturePointGeometryList
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    /// Creates elements with one patch per field (from "geometry_names"):
    /// one geometry per Gauss point, with a quadrature point of every patch
    void CreateQuadraturePointGeometriesMixed(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    /// Creates list of rQuadraturePointGeometryList for Gsbm
    void CreateQuadraturePointGeometriesGsbm(
        GeometryType& rBackground,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    /// Creates GSBM conditions on the skeleton of an earlier entry
    void CreateConditionsOnGsbmSkeleton(
        GeometriesArrayType& rFacets,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    ///@}
    ///@name CAD functionalities
    ///@{

    /// Gets list of geometries from rModelPart
    void GetGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rPhysicsParameters) const;


    ///@}
    ///@name Generate Elements and Conditions
    ///@{

    /// Creates elements from geometries
    void CreateElements(
        GeometriesArrayType::ptr_iterator GeometriesBegin,
        GeometriesArrayType::ptr_iterator GeometriesEnd,
        ModelPart& rModelPart,
        const std::string& rElementName) const;

    /// Creates conditions from geometries
    void CreateConditions(
        GeometriesArrayType::ptr_iterator GeometriesBegin,
        GeometriesArrayType::ptr_iterator GeometriesEnd,
        ModelPart& rModelPart,
        const std::string& rConditionName) const;

    ///@}
};

///@}

} // namespace Kratos
