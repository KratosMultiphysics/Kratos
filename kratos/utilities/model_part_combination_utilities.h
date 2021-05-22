//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_MODEL_PART_COMBINATION_UTILITIES)
#define KRATOS_MODEL_PART_COMBINATION_UTILITIES

// System includes

// External includes

// Project includes
#include "containers/model.h"

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

/**
 * @class ModelPartCombinationUtilities
 * @ingroup KratosCore
 * @brief 
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) ModelPartCombinationUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ModelPartCombinationUtilities
    KRATOS_CLASS_POINTER_DEFINITION( ModelPartCombinationUtilities );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * The default constructor
     */
    ModelPartCombinationUtilities(Model& rModel):
        mrModel(rModel)
    {
    }

    virtual ~ModelPartCombinationUtilities()= default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method combines the model parts existing in the given Model, according to the provided parameters
     * @param ThisParameters The configuration parameters
     */
    void CombineModelParts(Parameters ThisParameters = Parameters(R"({})"));

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ModelPartCombinationUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

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

    Model& mrModel;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const;

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
};// class ModelPartCombinationUtilities

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
#endif /* KRATOS_MODEL_PART_COMBINATION_UTILITIES defined */
