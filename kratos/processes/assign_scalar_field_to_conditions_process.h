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
//  Main authors:    Josep Maria Carbonell
//                   Riccardo Rossi
//

#if !defined(KRATOS_ASSIGN_SCALAR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED )
#define  KRATOS_ASSIGN_SCALAR_FIELD_TO_CONDITIONS_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/python_function_callback_utility.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class AssignScalarFieldToConditionsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;   
    
    /// Pointer definition of AssignScalarFieldToConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarFieldToConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AssignScalarFieldToConditionsProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : Process() , 
            mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "mesh_id": 0,
                "variable_name": "VARIABLE_NAME", 
                "interval"        : [0.0, 1e30],
                "value"           : "please give an expression in terms of the variable x, y, z, t",
                "local_axes" : {}
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mMeshId       = rParameters["mesh_id"].GetInt();
        mVariableName = rParameters["variable_name"].GetString();

        mpFunction = PythonGenericFunctionUtility::Pointer( new PythonGenericFunctionUtility(rParameters["value"].GetString(),  rParameters["local_axes"]));


        KRATOS_CATCH("");
    }



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


    /// Execute method is used to execute the AssignScalarFieldToConditionsProcess algorithms.
    void Execute() override
    {

        KRATOS_TRY;

        ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

        const double& rCurrentTime = rCurrentProcessInfo[TIME];


        if( KratosComponents< Variable<Vector> >::Has( mVariableName ) ) //case of double variable
        {
            InternalAssignValue<>(KratosComponents< Variable<Vector> >::Get(mVariableName), rCurrentTime);
        }
        else if( KratosComponents< array_1d_component_type >::Has( mVariableName ) ) //case of component variable
        {
            InternalAssignValueComponents<>(KratosComponents< array_1d_component_type >::Get(mVariableName), rCurrentTime);
        }
        else
        {
            KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mVariableName << std::endl;
        }

        KRATOS_CATCH("");

    }

    /// this function will be executed at every time step BEFORE performing the solve phase
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

    ModelPart& mrModelPart;
    PythonGenericFunctionUtility::Pointer mpFunction;
    std::string mVariableName;

    std::size_t mMeshId;

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * It calls the function 
     */
    void CallFunction(
        const Condition::Pointer& pCondition, 
        const double& time, 
        Vector& rValue
        )
    {

        Condition::GeometryType& rConditionGeometry = pCondition->GetGeometry();
        const unsigned int size = rConditionGeometry.size();

        if(rValue.size() !=  size)
        {
            rValue.resize(size,false);
        }

        for(unsigned int i=0; i<size; i++)
        {
            rValue[i] = mpFunction->CallFunction(rConditionGeometry[i].X(),rConditionGeometry[i].Y(),rConditionGeometry[i].Z(),time  );
        }

    }
    
    /**
     * It calls the function for components
     */
    void CallFunctionComponents(
        const Condition::Pointer& pCondition, 
        const double& time, 
        double& rValue
        )
    {
        Condition::GeometryType& rConditionGeometry = pCondition->GetGeometry();
        array_1d<double,3> center = rConditionGeometry.Center();

        rValue = mpFunction->CallFunction(center[0],center[1],center[2],time  );

    }

    /**
     * It calls the function (local system)
     */
    void CallFunctionLocalSystem(
        const Condition::Pointer& pCondition,
        const double& time, 
        Vector& rValue
        )
    {

        Condition::GeometryType& rConditionGeometry = pCondition->GetGeometry();
        const unsigned int size = rConditionGeometry.size();

        if(rValue.size() !=  size)
        {
            rValue.resize(size,false);
        }

        for(unsigned int i=0; i<size; i++)
        {
            rValue[i] = mpFunction->RotateAndCallFunction(rConditionGeometry[i].X(),rConditionGeometry[i].Y(),rConditionGeometry[i].Z(),time  );
        }
    }

    /**
     * It calls the function (local system) for components
     */
    void CallFunctionLocalSystemComponents(
        const Condition::Pointer& pCondition, 
        const double& time, 
        double& rValue
        )
    {
        Condition::GeometryType& rConditionGeometry = pCondition->GetGeometry();

        array_1d<double,3> center = rConditionGeometry.Center();

        rValue = mpFunction->RotateAndCallFunction(center[0],center[1],center[2],time  );
    }

    /**
     * It assigns time dependency
     */
    void AssignTimeDependentValue(
        const Condition::Pointer& pCondition, 
        const double& time, 
        Vector& rValue, 
        const double value
        )
    {

        Condition::GeometryType& rConditionGeometry = pCondition->GetGeometry();
        const unsigned int size = rConditionGeometry.size();

        if(rValue.size() !=  size)
        {
            rValue.resize(size,false);
        }

        for(unsigned int i=0; i<size; i++)
        {
            rValue[i] = value;
        }
    }

    template< class TVarType >
    void InternalAssignValue(
        TVarType& rVar, 
        const double& rTime
        )
    {
        const int nconditions = mrModelPart.GetMesh(mMeshId).Conditions().size();

        Vector Value;

        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator itBegin = mrModelPart.GetMesh(mMeshId).ConditionsBegin();

            if(mpFunction->DependsOnSpace())
            {
                if(mpFunction->UseLocalSystem())
                {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(int i = 0; i<nconditions; i++)
                    {
                        ModelPart::ConditionsContainerType::iterator it = itBegin + i;
                        this->CallFunctionLocalSystem(*(it.base()), rTime, Value);
                        it->SetValue(rVar, Value);
                    }
                }
                else
                {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(int i = 0; i<nconditions; i++)
                    {
                        ModelPart::ConditionsContainerType::iterator it = itBegin + i;
                        this->CallFunction(*(it.base()), rTime, Value);
                        it->SetValue(rVar, Value);
                    }
                }
            }
            else                                            // only varies in time
            {
                const double TimeValue = mpFunction->CallFunction(0.0, 0.0, 0.0,  rTime);
                // WARNING: do not parallelize with openmp. python GIL prevents it
                for(int i = 0; i<nconditions; i++)
                {
                    ModelPart::ConditionsContainerType::iterator it = itBegin + i;
                    this->AssignTimeDependentValue(*(it.base()), rTime, Value,  TimeValue);
                    it->SetValue(rVar, Value);
                }

            }
        }
    }

    template< class TVarType >
    void InternalAssignValueComponents(
        TVarType& rVar, 
        const double& rTime
        )
    {
        const int nconditions = mrModelPart.GetMesh(mMeshId).Conditions().size();
        
        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator itBegin = mrModelPart.GetMesh(mMeshId).ConditionsBegin();

            if(mpFunction->DependsOnSpace())
            {
                double Value;
                        
                if(mpFunction->UseLocalSystem())
                {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(int i = 0; i<nconditions; i++)
                    {
                        ModelPart::ConditionsContainerType::iterator it = itBegin + i;
                        this->CallFunctionLocalSystemComponents(*(it.base()), rTime, Value);
                        it->SetValue(rVar, Value);
                    }
                }
                else
                {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(int i = 0; i<nconditions; i++)
                    {
                        ModelPart::ConditionsContainerType::iterator it = itBegin + i;
                        this->CallFunctionComponents(*(it.base()), rTime, Value);
                        it->SetValue(rVar, Value);
                    }
                }
            }
            else // only varies in time
            {
                const double TimeValue = mpFunction->CallFunction(0.0, 0.0, 0.0,  rTime);
                // WARNING: do not parallelize with openmp. python GIL prevents it
                for(int i = 0; i<nconditions; i++)
                {
                    ModelPart::ConditionsContainerType::iterator it = itBegin + i;
                    it->SetValue(rVar, TimeValue);
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
