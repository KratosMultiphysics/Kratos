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
class KRATOS_API(MAPPING_APPLICATION) IgaIgaSurfaceMappingGeometriesModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(IgaIgaSurfaceMappingGeometriesModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IgaIgaSurfaceMappingGeometriesModeler()
        : Modeler()
    {
    }

    /// Constructor.
    IgaIgaSurfaceMappingGeometriesModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
    {
        mpModels.resize(1);
        mpModels[0] = &rModel;
    }

    /// Destructor.
    virtual ~IgaIgaSurfaceMappingGeometriesModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<IgaIgaSurfaceMappingGeometriesModeler>(rModel, ModelParameters);
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
        return "IgaIgaSurfaceMappingGeometriesModeler";
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
        rDestinationMP.SetElements(rReferenceMP.pElements());
        // ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        rDestinationMP.SetConditions(rReferenceMP.pConditions());

        for(auto geometry_it = rReferenceMP.GetRootModelPart().GeometriesBegin(); geometry_it != rReferenceMP.GetRootModelPart().GeometriesEnd(); geometry_it++){
            IndexType geometry_id = geometry_it->Id();
            rDestinationMP.AddGeometry(rReferenceMP.GetRootModelPart().pGetGeometry(geometry_id));
        }
    }

    void CheckParameters();

    void CreateIgaIgaSurfaceCouplingGeometry(ModelPart& rModelPartDomainA,
        ModelPart& rModelPartDomainB,
        ModelPart& rModelPartResult);

}; // Class IgaIgaSurfaceMappingGeometriesModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    IgaIgaSurfaceMappingGeometriesModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const IgaIgaSurfaceMappingGeometriesModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.
