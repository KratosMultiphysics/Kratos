//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti

#pragma once

// System includes

// External includes

// Project includes
#include "modeler/modeler.h"
#include "custom_utilities/mapping_intersection_utilities.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(MAPPING_APPLICATION) IgaFEMMappingGeometriesModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(IgaFEMMappingGeometriesModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    using NurbsCurveOnSurfacePointer = NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<NodeType>>::Pointer;
    using NurbsSurfacePointer = NurbsSurfaceGeometry<3, PointerVector<NodeType>>::Pointer;
    using BrepCurveOnSurfacePointer = BrepCurveOnSurface<PointerVector<NodeType>, false, PointerVector<Point>>::Pointer;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IgaFEMMappingGeometriesModeler()
        : Modeler()
    {
    }

    /// Constructor.
    IgaFEMMappingGeometriesModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
    {
        mpModels.resize(1);
        mpModels[0] = &rModel;
    }

    /// Destructor.
    virtual ~IgaFEMMappingGeometriesModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<IgaFEMMappingGeometriesModeler>(rModel, ModelParameters);
    }

    /// Adds the second model part to the modeler.
    void GenerateNodes(ModelPart& ThisModelPart) override
    {
        mpModels.push_back(&ThisModelPart.GetModel());
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupGeometryModel() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "IgaFEMMappingGeometriesModeler";
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
    std::vector<Model*> mpModels;

    void CopySubModelPart(ModelPart& rDestinationMP, ModelPart& rReferenceMP)
    {
        rDestinationMP.SetNodes(rReferenceMP.pNodes());
        rDestinationMP.SetNodalSolutionStepVariablesList(rReferenceMP.pGetNodalSolutionStepVariablesList());
        ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        rDestinationMP.SetConditions(coupling_conditions.pConditions());
    }

    void CreateIgaInterfaceBrepCurveOnSurface(ModelPart& rInterfaceModelPart);

    void CheckParameters();

    const Parameters GetDefaultParameters() const override
    {
        return Parameters( R"({
            "is_origin_iga"                : true,
            "is_surface_mapping"          : false,
            "use_initial_configuration"    : true,
            "echo_level"                   : 0,
        })");
    }

}; // Class IgaFEMMappingGeometriesModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    IgaFEMMappingGeometriesModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const IgaFEMMappingGeometriesModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.
