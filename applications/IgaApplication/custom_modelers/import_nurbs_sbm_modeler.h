//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "modeler/modeler.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"
#include "geometries/nurbs_curve_geometry.h"

namespace Kratos {

class KRATOS_API(IGA_APPLICATION) ImportNurbsSbmModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ImportNurbsSbmModeler );

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceGeometryType;
    typedef typename NurbsSurfaceGeometryType::Pointer NurbsSurfaceGeometryPointerType;

    using NurbsCurveGeometryPointerType = NurbsCurveGeometry<2, PointerVector<Node>>::Pointer;

    using CoordinatesArrayType = Geometry<NodeType>::CoordinatesArrayType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ImportNurbsSbmModeler()
        : Modeler() {}

    /// Constructor.
    ImportNurbsSbmModeler(
        Model & rModel,
        const Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
    }

    /// Destructor.
    virtual ~ImportNurbsSbmModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<ImportNurbsSbmModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupGeometryModel() override;

    ///@}


    Parameters ReadParamatersFile(
        const std::string& rDataFileName) const;
    

private:
    ///@name Private Member Variables
    ///@{

    Model* mpModel;


};

} // End namesapce Kratos