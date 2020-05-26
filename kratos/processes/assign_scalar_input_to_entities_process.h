//
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

#if !defined(KRATOS_assign_scalar_input_to_entities_process_H_INCLUDED )
#define  KRATOS_assign_scalar_input_to_entities_process_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class AssignScalarInputToEntitiesProcess
 * @ingroup KratosCore
 * @brief This function assigns a value from an input to a variable belonging to all of the entities in a given mesh
 * @details Can be used to any entities
 * @tparam TEntity The entity type
 * @author Vicente Mataix Ferrandiz
*/
template<class TEntity>
class KRATOS_API(KRATOS_CORE) AssignScalarInputToEntitiesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Node type definition
    typedef Node<3> NodeType;

    /// The container of the entities
    typedef PointerVectorSet<TEntity, IndexedObject> EntityContainerType;

    /// Pointer definition of AssignScalarInputToEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarInputToEntitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rModelPart The model part to be set
     * @param rParameters The configuration parameters
     */
    AssignScalarInputToEntitiesProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        );

    /// Destructor.
    ~AssignScalarInputToEntitiesProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the AssignScalarInputToEntitiesProcess algorithms.
     */
    void Execute() override;

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override
    {
        Execute();
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

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
        return "AssignScalarInputToEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignScalarInputToEntitiesProcess";
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

    /// Copy constructor.
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

    ModelPart& mrModelPart;                  /// The model part where to assign the values
    const Variable<double>* mpVariable = nullptr;  /// The pointer of the variable

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method assigns the value (with OMP)
     * @param rVar The variable to be assigned
     * @param Value The value to assign
     */
    template< class TVarType, class TDataType >
    void InternalAssignValue(TVarType& rVar, const TDataType Value)
    {
        auto& r_entities_array = GetEntitiesContainer();
        const int number_of_entities = static_cast<int>(r_entities_array.size());

        if(number_of_entities != 0) {
            const auto it_begin = r_entities_array.begin();

            #pragma omp parallel for
            for(int i = 0; i<number_of_entities; i++) {
                auto it_entity = it_begin + i;

                it_entity->SetValue(rVar, Value);
            }
        }
    }

    /**
     * @brief This method assigns the value (without OMP)
     * @param rVar The variable to be assigned
     * @param Value The value to assign
     */
    template< class TVarType, class TDataType >
    void InternalAssignValueSerial(TVarType& rVar, const TDataType Value)
    {
        auto& r_entities_array = GetEntitiesContainer();
        const int number_of_entities = static_cast<int>(r_entities_array.size());

        if(number_of_entities != 0) {
            const auto it_begin = r_entities_array.begin();

            for(int i = 0; i<number_of_entities; i++) {
                auto it_entity = it_begin + i;

                it_entity->SetValue(rVar, Value);
            }
        }
    }

    /**
     * @brief This method returns the current entity container
     * @return The current entity container
     */
    EntityContainerType& GetEntitiesContainer();

    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignScalarInputToEntitiesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TEntity>
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarInputToEntitiesProcess<TEntity>& rThis);

/// output stream function
template<class TEntity>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarInputToEntitiesProcess<TEntity>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_assign_scalar_input_to_entities_process_H_INCLUDED  defined
