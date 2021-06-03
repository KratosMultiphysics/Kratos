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
//  Main authors:    Miguel Maso Sotomayor
//


#ifndef KRATOS_MDPA_IO_MODELER_H_INCLUDED
#define KRATOS_MDPA_IO_MODELER_H_INCLUDED


// System includes

// External includes

// Project includes
#include "modeler.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class MdpaIoModeler
 * @ingroup KratosCore
 * @brief Mdpa file reading
 * @details Import a model part from an mdpa file
 * @author Miguel Maso Sotomayor
 */
class KRATOS_API(KRATOS_CORE) MdpaIoModeler : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(MdpaIoModeler);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    MdpaIoModeler() : Modeler()
    {
    }

    /**
     * @brief Constructor with Model and Parameters
     */
    MdpaIoModeler(Model& rModel, Parameters ModelerParameters = Parameters());

    /**
     * @brief Destructor
     */
    virtual ~MdpaIoModeler() = default;

    /**
     * @brief Creates the Modeler Pointer
     */
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<MdpaIoModeler>(rModel, ModelParameters);
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
        return "MdpaIoModeler";
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
    ///@name Member variables
    ///@{

    Model* mpModel = nullptr;

    ///@}
    ///@name Private Operations
    ///@{

    const Parameters GetDefaultParameters() const;

    ///@}
    ///@name Serializer
    ///@{

    friend class Serializer;

    ///@}

}; // Class MdpaIoModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    MdpaIoModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const MdpaIoModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_MDPA_IO_MODELER_H_INCLUDED  defined
