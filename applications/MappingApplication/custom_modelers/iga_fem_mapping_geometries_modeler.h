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

class KRATOS_API(MAPPING_APPLICATION) IgaFEMMappingGeometriesModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(IgaFEMMappingGeometriesModeler);

    using SizeType = std::size_t;
    using IndexType = std::size_t;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using GeometryPointerType = typename GeometryType::Pointer; 

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

    void CopySubModelPartIgaInterface(ModelPart& rDestinationMP, ModelPart& rReferenceMP)
    {
        // We just set the control points which have support on the interface (we are not interested in the CPs whose shape function is 0 on the interface)
        // For that, we loop over the quadrature points and check its shape functions
        
        // Set to keep track of already added node IDs (avoid duplicates)
        std::unordered_set<IndexType> added_node_ids;

        for (auto& r_cond : rReferenceMP.Conditions()) {
            auto& r_geom = r_cond.GetGeometry();
            auto& r_N = r_geom.ShapeFunctionsValues();

            // Loop over nodes in the geometry
            for (std::size_t i = 0; i < r_geom.size(); ++i) {
                const double N_i = r_N(0, i);

                if (N_i > 1e-5) { // Node has influence at quadrature point
                    const auto& node = r_geom[i];

                    if (added_node_ids.insert(node.Id()).second) {
                        // If node ID was not already added, add to destination
                        rDestinationMP.AddNode(r_geom.pGetPoint(i));
                    }
                }
            }
        }

        rDestinationMP.SetNodalSolutionStepVariablesList(rReferenceMP.pGetNodalSolutionStepVariablesList());
        ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        rDestinationMP.SetConditions(coupling_conditions.pConditions());
    }

    void CopySubModelPartFEMInterface(ModelPart& rDestinationMP, ModelPart& rReferenceMP)
    {
        rDestinationMP.SetNodes(rReferenceMP.pNodes());
        rDestinationMP.SetNodalSolutionStepVariablesList(rReferenceMP.pGetNodalSolutionStepVariablesList());
        ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        rDestinationMP.SetConditions(coupling_conditions.pConditions());
    }

    void CreateIgaInterfaceBrepCurveOnSurfaceConditions(ModelPart& rInterfaceModelPart);

    void CreateFEMInterfaceNurbsCurveConditions(ModelPart& rInterfaceModelPart);

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
