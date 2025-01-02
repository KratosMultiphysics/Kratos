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


#if !defined(KRATOS_IGA_IO_THB_MODELER_H_INCLUDED )
#define  KRATOS_IGA_IO_THB_MODELER_H_INCLUDED


// System includes

// External includes

// Project includes
#include "modeler.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) CadIoThbModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(CadIoThbModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CadIoThbModeler()
        : Modeler()
    {
    }

    /// Constructor.
    CadIoThbModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
    }

    /// Destructor.
    virtual ~CadIoThbModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<CadIoThbModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupGeometryModel() override;

    void SetupModelPart() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CadIoThbModeler";
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

    Model* mpModel = nullptr;

    ///@}
    ///@name Serializer
    ///@{

    friend class Serializer;

    ///@}
}; // Class CadIoThbModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    CadIoThbModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const CadIoThbModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_IO_THB_MODELER_H_INCLUDED  defined
