//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ENTITY_ERASE_PROCESS_INCLUDED )
#define  KRATOS_ENTITY_ERASE_PROCESS_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"

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
 * @class EntitiesEraseProcess
 * @ingroup KratosCore
 * @brief This process removes the entities from a model part with the flag TO_ERASE
 * @details It does some things additionally:
 *  - RemoveFromAllLevels If the entities will be removed from all levels
 *  - AssignFlag If the flag will be assigned (this means removing all entities of the model part)
 *  - RemoveBelongings If removing the belogings of the removed entity
 * @tparam TEntity The type entity considered
 * @author Vicente Mataix Ferrandiz
*/
template<class TEntity>
class KRATOS_API(KRATOS_CORE) EntitiesEraseProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    // DEFINITION OF FLAGS TO CONTROL THE BEHAVIOUR
    KRATOS_DEFINE_LOCAL_FLAG(REMOVE_FROM_ALL_LEVELS); /// If the entities will be removed from all levels
    KRATOS_DEFINE_LOCAL_FLAG(ERASE_ALL_ENTITIES);     /// If the flag will be assigned (this means removing all entities of the model part)

    /// Pointer definition of EntitiesEraseProcess
    KRATOS_CLASS_POINTER_DEFINITION(EntitiesEraseProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The default parameters
     * @param rModelPart The model part where remove the entities
     * @param Options The flags definitions
     */
    explicit EntitiesEraseProcess(
        ModelPart& rModelPart,
        Flags Options = REMOVE_FROM_ALL_LEVELS.AsFalse() | ERASE_ALL_ENTITIES.AsFalse()
        );

    /**
     * @brief The default parameters (with parameters)
     * @param rModelPart The model part where remove the entities
     * @param ThisParameters The configuration parameters
     */
    explicit EntitiesEraseProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters
        );

    /// Destructor.
    ~EntitiesEraseProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "EntitiesEraseProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EntitiesEraseProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{


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

    ModelPart& mrModelPart; /// The model part to be computed

    Flags mrOptions;        /// Local flags

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters();

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    EntitiesEraseProcess& operator=(EntitiesEraseProcess const& rOther);

    /// Copy constructor.
    //EntitiesEraseProcess(EntitiesEraseProcess const& rOther);


    ///@}

}; // Class EntitiesEraseProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TEntity>
inline std::istream& operator >> (std::istream& rIStream,
                                  EntitiesEraseProcess<TEntity>& rThis);

/// output stream function
template<class TEntity>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const EntitiesEraseProcess<TEntity>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ENTITY_ERASE_PROCESS_INCLUDED  defined
