//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//


#if !defined(KRATOS_IGA_COSIM_MODELER_H_INCLUDED )
#define  KRATOS_IGA_COSIM_MODELER_H_INCLUDED


// System includes

// External includes

// Project includes
#include "modeler/modeler.h"
#include "custom_utilities/mapping_intersection_iga_utilities.h"
#include "geometries/brep_surface.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/nurbs_curve_geometry.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(IGA_APPLICATION) IgaCosimModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(IgaCosimModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;

    typedef typename Properties::Pointer PropertiesPointerType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;


    typedef Point EmbeddedNodeType;
    typedef PointerVector<NodeType> ContainerNodeType;
    typedef PointerVector<EmbeddedNodeType> ContainerEmbeddedNodeType;
    typedef BrepSurface<ContainerNodeType, false, ContainerEmbeddedNodeType> BrepSurfaceType;
    typedef BrepCurveOnSurface<ContainerNodeType, false, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;

    typedef NurbsCurveGeometry<3, ContainerNodeType> NurbsCurveType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IgaCosimModeler()
        : Modeler()
    {
    }

    /// Constructor.
    IgaCosimModeler(
        Model & rModel,
        Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
    {
        mpModels.resize(1);
        mpModels[0] = &rModel;
    }

    /// Destructor.
    virtual ~IgaCosimModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<IgaCosimModeler>(rModel, ModelParameters);
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
        return "IgaCosimModeler";
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

    std::vector<Model*> mpModels;

    void CopySubModelPart(ModelPart& rDestinationMP, ModelPart& rReferenceMP)
    {
        KRATOS_WATCH(rReferenceMP)
        // rDestinationMP.SetNodes(rReferenceMP.pNodes());
        // rDestinationMP.SetConditions(rReferenceMP.pConditions()); //this doesn't work
        rDestinationMP.SetNodalSolutionStepVariablesList(rReferenceMP.pGetNodalSolutionStepVariablesList());
        KRATOS_WATCH(rReferenceMP.pGetNodalSolutionStepVariablesList())
        KRATOS_WATCH(rDestinationMP.pGetNodalSolutionStepVariablesList())
        ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        KRATOS_WATCH(coupling_conditions.pGetNodalSolutionStepVariablesList())
        rDestinationMP.SetNodes(coupling_conditions.pNodes());
        rDestinationMP.SetConditions(coupling_conditions.pConditions());
        rDestinationMP.AddGeometry(coupling_conditions.pGetGeometry(10));
        KRATOS_WATCH(rDestinationMP)
    }

    void CreateInterfaceLineCouplingConditions(ModelPart& rInterfaceModelPart, bool MasterToSlave);

    void CheckParameters();

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
    IgaCosimModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const IgaCosimModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_COSIM_MODELER_H_INCLUDED  defined