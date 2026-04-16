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
#include "geometries/brep_surface.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_geometry.h"
#include "custom_utilities/cad_geometry_reconstruction_utility.h"
#include <limits>

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(MAPPING_APPLICATION) IgaIgaMappingGeometriesModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(IgaIgaMappingGeometriesModeler);

    using SizeType = std::size_t;
    using IndexType = std::size_t;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using GeometryPointerType = typename GeometryType::Pointer; 
    using PointType = Point;
    using CoordinatesArrayType = typename PointType::CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IgaIgaMappingGeometriesModeler()
        : Modeler()
    {
    }

    /// Constructor.
    IgaIgaMappingGeometriesModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
    {
        mpModels.resize(1);
        mpModels[0] = &rModel;
    }

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<IgaIgaMappingGeometriesModeler>(rModel, ModelParameters);
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
        return "IgaIgaMappingGeometriesModeler";
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
    double mSearchRadius = std::numeric_limits<double>::max();

    void CopySubModelPartBrepSurfaceInterface(ModelPart& rDestinationMP, ModelPart& rReferenceMP)
    {
        rDestinationMP.SetNodes(rReferenceMP.pNodes());
        rDestinationMP.SetNodalSolutionStepVariablesList(rReferenceMP.pGetNodalSolutionStepVariablesList());
        rDestinationMP.SetConditions(rReferenceMP.pConditions());
        rDestinationMP.SetElements(rReferenceMP.pElements());
        rDestinationMP.SetProperties(rReferenceMP.pProperties());
        
        for(auto geometry_it = rReferenceMP.GetRootModelPart().GeometriesBegin(); geometry_it != rReferenceMP.GetRootModelPart().GeometriesEnd(); geometry_it++){
            IndexType geometry_id = geometry_it->Id();
            rDestinationMP.AddGeometry(rReferenceMP.GetRootModelPart().pGetGeometry(geometry_id));
        }
    }

    void ReconstructBrepGeometryFromFile(
        ModelPart& rInterfaceMP,
        const Parameters& rReconstructionSettings)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(rReconstructionSettings.Has("filename"))
            << "Missing \"filename\" in reconstruction settings." << std::endl;

        const std::string filename = rReconstructionSettings["filename"].GetString();

        const bool use_model_part_existing_nodes =
            rReconstructionSettings.Has("use_model_part_existing_nodes")
                ? rReconstructionSettings["use_model_part_existing_nodes"].GetBool()
                : true;

        if (use_model_part_existing_nodes) {
            CadGeometryReconstructionUtility::ReconstructModelPartBrepGeometryFromCadJson(
                filename,
                rInterfaceMP,
                mEchoLevel);
        } else {
            CadGeometryReconstructionUtility::ReconstructModelPartBrepGeometryFromCadJsonCreatingNodes(
                filename,
                rInterfaceMP,
                mEchoLevel);
        }

        KRATOS_CATCH("")
    }

    // Creates a coupling geometry connecting the origin and destination IGA surfaces
    void CreateIgaIgaSurfaceCouplingGeometry(ModelPart& rModelPartDomainA, ModelPart& rModelPartDomainB, ModelPart& rModelPartResult);

    void CheckParameters();

    const Parameters GetDefaultParameters() const override
    {
        return Parameters( R"({
            "is_surface_mapping"          : false,
            "search_radius"              : 1.0e+10,
            "use_initial_configuration"    : true,
            "echo_level"                   : 0,
        })");
    }

}; // Class IgaIgaMappingGeometriesModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    IgaIgaMappingGeometriesModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const IgaIgaMappingGeometriesModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.
