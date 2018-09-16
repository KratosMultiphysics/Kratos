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

#if !defined(KRATOS_ASSIGN_SCALAR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED )
#define  KRATOS_ASSIGN_SCALAR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED

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
 * @class AssignScalarFieldToConditionsProcess
 * @ingroup KratosCore
 * @brief This process is used in order to assign a function to a condition
 * @details This function assigns a value depending on a function to a variable in all of the conditions in a given mesh.
 * The behaviour is the following:
 * - Option 1 - Variable<double>: It is evaluated in the center of the condition
 * - Option 2 - array_1d_component_type: The same as Variable evaluated in the center of the element
 * - Option 3 - Variable< Vector > : The vector has to have a size equal to the number of nodes, and its values are computed per each entry of the vector using the coordinates of the nodes which occupies the same position in the geometry
 * @author Riccardo Rossi
 * @author Josep Maria Carbonell
 * @author Vicente Mataix Ferrandiz
*/
class AssignScalarFieldToConditionsProcess
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

    /// The defition of a component of an array
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

    /// Pointer definition of AssignScalarFieldToConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarFieldToConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The default constructor
     * @param rModelPart The model part where the scalar field will be applied
     * @param rParameters The configuration parameters
     */
    AssignScalarFieldToConditionsProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        );

    /// Destructor.
    ~AssignScalarFieldToConditionsProcess() override {}

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
     * @brief Execute method is used to execute the AssignScalarFieldToConditionsProcess algorithms.
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
        return "AssignScalarFieldToConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignScalarFieldToConditionsProcess";
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
    AssignScalarFieldToConditionsProcess(AssignScalarFieldToConditionsProcess const& rOther);

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
     * @param pCondition The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     */
    void CallFunction(
        const Condition::Pointer& pCondition,
        const double Time,
        Vector& rValue
        );

    /**
     * @brief It calls the function for components
     * @param pCondition The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     */
    void CallFunctionComponents(
        const Condition::Pointer& pCondition,
        const double Time,
        double& rValue
        );

    /**
     * @brief It calls the function (local system)
     * @param pCondition The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     */
    void CallFunctionLocalSystem(
        const Condition::Pointer& pCondition,
        const double Time,
        Vector& rValue
        );

    /**
     * @brief It calls the function (local system) for components
     * @param pCondition The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     */
    void CallFunctionLocalSystemComponents(
        const Condition::Pointer& pCondition,
        const double Time,
        double& rValue
        );

    /**
     * @brief It assigns time dependency
     * @param pCondition The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     * @param Value The dependency value
     */
    void AssignTimeDependentValue(
        const Condition::Pointer& pCondition,
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
        const SizeType nconditions = mrModelPart.GetMesh(mMeshId).Conditions().size();

        Vector Value;

        if(nconditions != 0) {
            auto it_begin = mrModelPart.GetMesh(mMeshId).ConditionsBegin();

            if(mpFunction->DependsOnSpace()) {
                if(mpFunction->UseLocalSystem()) {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(IndexType i = 0; i<nconditions; ++i) {
                        auto it_cond = it_begin + i;
                        this->CallFunctionLocalSystem(*(it_cond.base()), Time, Value);
                        it_cond->SetValue(rVar, Value);
                    }
                } else {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(IndexType i = 0; i<nconditions; ++i) {
                        auto it_cond = it_begin + i;
                        this->CallFunction(*(it_cond.base()), Time, Value);
                        it_cond->SetValue(rVar, Value);
                    }
                }
            } else { // only varies in time
                const double TimeValue = mpFunction->CallFunction(0.0, 0.0, 0.0,  Time);
                // WARNING: do not parallelize with openmp. python GIL prevents it
                for(IndexType i = 0; i<nconditions; ++i) {
                    auto it_cond = it_begin + i;
                    this->AssignTimeDependentValue(*(it_cond.base()), Time, Value,  TimeValue);
                    it_cond->SetValue(rVar, Value);
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
        const SizeType nconditions = mrModelPart.GetMesh(mMeshId).Conditions().size();

        if(nconditions != 0) {
            auto it_begin = mrModelPart.GetMesh(mMeshId).ConditionsBegin();

            if(mpFunction->DependsOnSpace()) {
                double Value;

                if(mpFunction->UseLocalSystem()) {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(IndexType i = 0; i<nconditions; ++i) {
                        auto it_cond = it_begin + i;
                        this->CallFunctionLocalSystemComponents(*(it_cond.base()), Time, Value);
                        it_cond->SetValue(rVar, Value);
                    }
                } else {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(IndexType i = 0; i<nconditions; ++i) {
                        auto it_cond = it_begin + i;
                        this->CallFunctionComponents(*(it_cond.base()), Time, Value);
                        it_cond->SetValue(rVar, Value);
                    }
                }
            } else { // only varies in time
                const double TimeValue = mpFunction->CallFunction(0.0, 0.0, 0.0,  Time);
                // WARNING: do not parallelize with openmp. python GIL prevents it
                for(IndexType i = 0; i<nconditions; ++i) {
                    auto it_cond = it_begin + i;
                    it_cond->SetValue(rVar, TimeValue);
                }

            }
        }
    }

    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    AssignScalarFieldToConditionsProcess& operator=(AssignScalarFieldToConditionsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignScalarFieldToConditionsProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarFieldToConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarFieldToConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_SCALAR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED  defined
