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


#if !defined(KRATOS_IGA_MODELER_H_INCLUDED )
#define  KRATOS_IGA_MODELER_H_INCLUDED


// System includes

// External includes

// Project includes
#include "modeler.h"
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) IgaModeler
    : Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(IgaModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IgaModeler()
        : Modeler()
    {}

    /// Destructor.
    virtual ~IgaModeler() = default;

    ///@}
    ///@name Stages
    ///@{

    void GenerateModelPart(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Parameters PrepareGeometryParameters,
        IndexType EchoLevel = 0) const override;

    ///@}
    ///@name Iga functionalities
    ///@{

    /// Creates integration point geometries and applies elements and conditions
    void CreateIntegrationDomain(
        ModelPart& rModelPart,
        ModelPart& rCadModelPart,
        const Parameters& rParameters,
        IndexType EchoLevel = 0) const;

    /// Creates list of rQuadraturePointGeometryList
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rCadSubModelPart,
        const Parameters& rParameters,
        IndexType EchoLevel = 0) const;

    ///@}
    ///@name CAD functionalities
    ///@{

    /// Gets list of geometries from rModelPart
    void GetCadGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters& rParameters,
        IndexType EchoLevel = 0) const;

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
    IgaModeler& operator=(IgaModeler const& rOther);

    /// Copy constructor.
    IgaModeler(IgaModeler const& rOther);

    ///@}

}; // Class CadModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    IgaModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const IgaModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_MODELER_H_INCLUDED  defined


