//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//


#if !defined(KRATOS_MODELER_H_INCLUDED )
#define  KRATOS_MODELER_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/properties.h"
#include "spatial_containers/spatial_containers.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(Modeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef typename Properties::Pointer PropertiesPointerType;

    typedef Geometry<Node<3>> GeometryType;
    typedef PointerVector<GeometryType> GeometriesArrayType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Modeler() {}

    /// Destructor.
    virtual ~Modeler() {}

    ///@}
    ///@name Operations
    ///@{

    virtual void GenerateModelPart(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        Element const& rReferenceElement,
        Condition const& rReferenceBoundaryCondition);

    virtual void GenerateMesh(
        ModelPart& ThisModelPart,
        Element const& rReferenceElement,
        Condition const& rReferenceBoundaryCondition);

    virtual void GenerateNodes(
        ModelPart& ThisModelPart);

    ///@}
    ///@name Generate Elements and Conditions
    ///@{

    /// Creates elements from model part geometries
    static void CreateElements(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        std::string& rElementName,
        PropertiesPointerType pProperties,
        SizeType EchoLevel = 0);

    /// Creates elements from geometries
    template<class TContainerType>
    static void CreateElements(
        typename TContainerType::iterator& rGeometriesBegin,
        typename TContainerType::iterator& rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rElementName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        SizeType EchoLevel = 0);

    /// Creates conditions from model part geometries
    static void CreateConditions(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        std::string& rConditionName,
        PropertiesPointerType pProperties,
        SizeType EchoLevel = 0);

    /// Creates conditions from geometries
    template<class TContainerType>
    static void CreateConditions(
        typename TContainerType::iterator& rGeometriesBegin,
        typename TContainerType::iterator& rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        SizeType EchoLevel = 0);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Modeler";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Modeler& operator=(Modeler const& rOther);

    /// Copy constructor.
    Modeler(Modeler const& rOther);

    ///@}

}; // Class Modeler

///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Modeler& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Modeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_MODELER_H_INCLUDED  defined


