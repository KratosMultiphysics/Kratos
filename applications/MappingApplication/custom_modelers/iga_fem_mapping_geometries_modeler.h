//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan I. Camarotti
//                   Andrea Gorgi
//

#pragma once

// System includes

// External includes

// Project includes
#include "modeler/modeler.h"


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
        ModelPart& intermediate_model_part = rReferenceMP; 

        /*std::vector<IndexType> co_sim_conditions_id ;

        for (auto& r_cond : rReferenceMP.GetSubModelPart("PLEASE_SPECIFY").Conditions()) {
             co_sim_conditions_id.push_back(r_cond.Id());
        }

        for (IndexType i=0; i < co_sim_conditions_id.size(); i++){
            intermediate_model_part.RemoveCondition(co_sim_conditions_id[i]);
        }*/

        for (auto& r_cond : intermediate_model_part.Conditions()) {
            auto& r_geom = r_cond.GetGeometry();
            auto& r_N = r_geom.ShapeFunctionsValues();

            for (IndexType i = 0; i<r_N.size2();++i)
            {
                if(r_N(0,i) > 1e-5)
                {
                    rDestinationMP.AddNode(r_geom.pGetPoint(i));
                }
            }
        }
        rDestinationMP.SetNodalSolutionStepVariablesList(rReferenceMP.pGetNodalSolutionStepVariablesList());
        ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        rDestinationMP.SetConditions(coupling_conditions.pConditions());
        /*rDestinationMP.SetNodes(rReferenceMP.pNodes());
        rDestinationMP.SetNodalSolutionStepVariablesList(rReferenceMP.pGetNodalSolutionStepVariablesList());
        ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        rDestinationMP.SetConditions(coupling_conditions.pConditions());

        // Remove the control points which have all zero shape functions in the corresponding condition (quadrature point condition)
        for (auto& r_cond : rReferenceMP.Conditions()) {
            auto& r_geom = r_cond.GetGeometry();
            auto& r_N = r_geom.ShapeFunctionsValues();

            for (IndexType i = 0; i<r_N.size2();++i)
            {
                if(r_N(0,i) < 1e-8)
                {
                    rDestinationMP.RemoveNode(r_geom.pGetPoint(i));
                }
            }
        }*/
    }

    void CopySubModelPartFEM(ModelPart& rDestinationMP, ModelPart& rReferenceMP)
    {
        rDestinationMP.SetNodes(rReferenceMP.pNodes());
        rDestinationMP.SetNodalSolutionStepVariablesList(rReferenceMP.pGetNodalSolutionStepVariablesList());
        ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        rDestinationMP.SetConditions(coupling_conditions.pConditions());
    }

     void CopySubModelPartDestinationStrongSupport(ModelPart& rDestinationMP, ModelPart& rReferenceMP)
    {
        rDestinationMP.SetNodes(rReferenceMP.pNodes());
        rDestinationMP.SetNodalSolutionStepVariablesList(rReferenceMP.pGetNodalSolutionStepVariablesList());
        for(auto geometry_it = rReferenceMP.GeometriesBegin(); geometry_it != rReferenceMP.GeometriesEnd(); geometry_it++){
            IndexType geometry_id = geometry_it->Id();
            rDestinationMP.AddGeometry(rReferenceMP.pGetGeometry(geometry_id));
        }
    }

    void CreateInterfaceLineBrepCurveOnSurface(ModelPart& rInterfaceModelPart);

    void CreateInterfaceLineBrepCurveOnSurfaceStrongSupport(ModelPart& rInterfaceModelPart);

    void CreateInterfaceLineBrepCurve(ModelPart& rInterfaceModelPart);

    void CheckParameters();

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
