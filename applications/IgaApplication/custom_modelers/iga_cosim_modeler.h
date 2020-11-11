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

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;

    typedef typename Properties::Pointer PropertiesPointerType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;

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
        rDestinationMP.SetNodes(rReferenceMP.pNodes());
        ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        rDestinationMP.SetConditions(coupling_conditions.pConditions());
    }

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