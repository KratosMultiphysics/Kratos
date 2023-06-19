//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SERIAL_MODEL_PART_COMBINATOR_MODELER_INCLUDED_H)
#define KRATOS_SERIAL_MODEL_PART_COMBINATOR_MODELER_INCLUDED_H

// System includes

// External includes

// Project includes
#include "modeler/modeler.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"

namespace Kratos 
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
/**
 * @class SerialModelPartCombinatorModeler
 * @ingroup KratosCore
 * @brief This modeler combines several model parts in serial, reading its respective *.mdpa files (maybe others in the future)
 * @details Uses ModelPartCombinationUtilities to combine different ModelParts into one single ModelPart, with the corresponding sub ModelParts
 * @author Vicente Maataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) SerialModelPartCombinatorModeler 
    : public Modeler {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SerialModelPartCombinatorModeler
    KRATOS_CLASS_POINTER_DEFINITION(SerialModelPartCombinatorModeler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SerialModelPartCombinatorModeler() : Modeler()
    {
    }

    /**
     * @brief Constructor using a Model and Parameters
     * @param rModel Reference of the Model
     * @param ModelerParameters Parameters of the discretization
     */
    SerialModelPartCombinatorModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters()
        ) : Modeler(rModel, ModelerParameters),
            mpModel(&rModel),
            mParameters(ModelerParameters)
    {
    }

    /// Destructor.
    virtual ~SerialModelPartCombinatorModeler() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates the Modeler Pointer and returns a pointer to a new
     * SerialModelPartCombinatorModeler, created using the given input
     * @param rModel Reference of the Model
     * @param ModelerParameters Parameters of the discretization
     * @return a Pointer to the new Modeler
     */
    Modeler::Pointer Create(
        Model& rModel,
        const Parameters ModelParameters
        ) const override;

    /**
     * @brief Convert the geometry into an analysis suitable ModePart
     */
    void SetupModelPart() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SerialModelPartCombinatorModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr; /// The model of the problem

    Parameters mParameters;   /// The configuration parameters

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class SerialModelPartCombinatorModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream,
    SerialModelPartCombinatorModeler& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const SerialModelPartCombinatorModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_SERIAL_MODEL_PART_COMBINATOR_MODELER_INCLUDED_H  defined
