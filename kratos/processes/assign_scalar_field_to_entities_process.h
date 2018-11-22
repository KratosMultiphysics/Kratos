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
//  Main authors:    Riccardo Rossi
//                   Josep Maria Carbonell
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ASSIGN_SCALAR_FIELD_TO_ENTITIES_PROCESS_H_INCLUDED )
#define  KRATOS_ASSIGN_SCALAR_FIELD_TO_ENTITIES_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "utilities/python_function_callback_utility.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class AssignScalarFieldToEntitiesProcess
 * @ingroup KratosCore
 * @brief This process is used in order to assign a function to a entity
 * @details This function assigns a value depending on a function to a variable in all of the conditions in a given mesh.
 * The behaviour is the following:
 * - Option 1 - Variable<double>: It is evaluated in the center of the entities
 * - Option 2 - array_1d_component_type: The same as Variable evaluated in the center of the element
 * - Option 3 - Variable< Vector > : The vector has to have a size equal to the number of nodes, and its values are computed per each entry of the vector using the coordinates of the nodes which occupies the same position in the geometry
 * @author Riccardo Rossi
 * @author Josep Maria Carbonell
 * @author Vicente Mataix Ferrandiz
*/
template<class TEntity>
class AssignScalarFieldToEntitiesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// The definitionof the node
    typedef Node<3> NodeType;

    /// The definitionof the geometry
    typedef Geometry<NodeType> GeometryType;

    /// The IndexType definition
    typedef std::size_t IndexType;

    /// The sizeType definition
    typedef std::size_t SizeType;

    /// The container of the entities
    typedef PointerVectorSet<TEntity, IndexedObject> EntityContainerType;

    /// The defition of a component of an array
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

    /// Pointer definition of AssignScalarFieldToEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarFieldToEntitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The default constructor
     * @param rModelPart The model part where the scalar field will be applied
     * @param rParameters The configuration parameters
     */
    AssignScalarFieldToEntitiesProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        );

    /// Destructor.
    ~AssignScalarFieldToEntitiesProcess() override {}

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
    ///@{.

    /**
     * @brief Execute method is used to execute the AssignScalarFieldToEntitiesProcess algorithms.
     */
    void Execute() override;

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override
    {
        Execute();
    }

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
        return "AssignScalarFieldToEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignScalarFieldToEntitiesProcess";
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
    AssignScalarFieldToEntitiesProcess(AssignScalarFieldToEntitiesProcess const& rOther);

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

    ModelPart& mrModelPart;                           /// The modelpart where compute

    PythonGenericFunctionUtility::Pointer mpFunction; /// The python function used, depends on X, Y, Z, and t

    std::string mVariableName;                        /// The name of the vaiable to assign

    std::size_t mMeshId = 0;                          /// The id of the mesh (0 by deafault)

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief It calls the function for a vector variable
     * @param pEntity The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     */
    void CallFunction(
        const typename TEntity::Pointer& pEntity,
        const double Time,
        Vector& rValue
        );

    /**
     * @brief It calls the function for components
     * @param pEntity The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     */
    void CallFunctionComponents(
        const typename TEntity::Pointer& pEntity,
        const double Time,
        double& rValue
        );

    /**
     * @brief It calls the function (local system)
     * @param pEntity The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     */
    void CallFunctionLocalSystem(
        const typename TEntity::Pointer& pEntity,
        const double Time,
        Vector& rValue
        );

    /**
     * @brief It calls the function (local system) for components
     * @param pEntity The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     */
    void CallFunctionLocalSystemComponents(
        const typename TEntity::Pointer& pEntity,
        const double Time,
        double& rValue
        );

    /**
     * @brief It assigns time dependency
     * @param pEntity The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     * @param Value The dependency value
     */
    void AssignTimeDependentValue(
        const typename TEntity::Pointer& pEntity,
        const double Time,
        Vector& rValue,
        const double Value
        );

    /**
     * @brief This is the methods that set the values globally (tries all the possible options)
     * @param rVar The variable to set
     * @param Time The current time
     */
    template< class TVarType >
    void InternalAssignValueVector(
        TVarType& rVar,
        const double Time
        )
    {
        auto& r_entities_array = GetEntitiesContainer();
        const int number_of_entities = static_cast<int>(r_entities_array.size());

        Vector Value;

        if(number_of_entities != 0) {
            auto it_begin = r_entities_array.begin();

            if(mpFunction->DependsOnSpace()) {
                if(mpFunction->UseLocalSystem()) {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(IndexType i = 0; i<number_of_entities; ++i) {
                        auto it_entity = it_begin + i;
                        this->CallFunctionLocalSystem(*(it_entity.base()), Time, Value);
                        it_entity->SetValue(rVar, Value);
                    }
                } else {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(IndexType i = 0; i<number_of_entities; ++i) {
                        auto it_entity = it_begin + i;
                        this->CallFunction(*(it_entity.base()), Time, Value);
                        it_entity->SetValue(rVar, Value);
                    }
                }
            } else { // only varies in time
                const double TimeValue = mpFunction->CallFunction(0.0, 0.0, 0.0,  Time);
                // WARNING: do not parallelize with openmp. python GIL prevents it
                for(IndexType i = 0; i<number_of_entities; ++i) {
                    auto it_entity = it_begin + i;
                    this->AssignTimeDependentValue(*(it_entity.base()), Time, Value,  TimeValue);
                    it_entity->SetValue(rVar, Value);
                }
            }
        }
    }

    /**
     * @brief This is the methods that set the values globally (tries all the possible options) for components
     * @param rVar The variable to set
     * @param Time The current time
     */
    template< class TVarType >
    void InternalAssignValueScalar(
        TVarType& rVar,
        const double Time
        )
    {
        auto& r_entities_array = GetEntitiesContainer();
        const int number_of_entities = static_cast<int>(r_entities_array.size());

        if(number_of_entities != 0) {
            auto it_begin = r_entities_array.begin();

            if(mpFunction->DependsOnSpace()) {
                double Value;

                if(mpFunction->UseLocalSystem()) {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(IndexType i = 0; i<number_of_entities; ++i) {
                        auto it_entity = it_begin + i;
                        this->CallFunctionLocalSystemComponents(*(it_entity.base()), Time, Value);
                        it_entity->SetValue(rVar, Value);
                    }
                } else {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(IndexType i = 0; i<number_of_entities; ++i) {
                        auto it_entity = it_begin + i;
                        this->CallFunctionComponents(*(it_entity.base()), Time, Value);
                        it_entity->SetValue(rVar, Value);
                    }
                }
            } else { // only varies in time
                const double TimeValue = mpFunction->CallFunction(0.0, 0.0, 0.0,  Time);
                // WARNING: do not parallelize with openmp. python GIL prevents it
                for(IndexType i = 0; i<number_of_entities; ++i) {
                    auto it_entity = it_begin + i;
                    it_entity->SetValue(rVar, TimeValue);
                }

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
    AssignScalarFieldToEntitiesProcess& operator=(AssignScalarFieldToEntitiesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignScalarFieldToEntitiesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TEntity>
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarFieldToEntitiesProcess<TEntity>& rThis);

/// output stream function
template<class TEntity>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarFieldToEntitiesProcess<TEntity>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_SCALAR_FIELD_TO_ENTITIES_PROCESS_H_INCLUDED  defined
