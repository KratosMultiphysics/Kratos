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

    /// Pointer definition of AssignScalarFieldToConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarFieldToConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AssignScalarFieldToConditionsProcess(ModelPart& model_part,
                                        Parameters rParameters
                                        ) : Process() , mr_model_part(model_part)
    {
        KRATOS_TRY
                        
        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "mesh_id": 0,
                "variable_name": "VARIABLE_NAME", 
                "interval"        : [0.0, 1e30],
                "spatially_varying": true, 
                "value"           : "please give an expression in terms of the variable x, y, z, t",
                "local_axes" : {}
            }  )" );

        
        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mmesh_id       = rParameters["mesh_id"].GetInt();
        mvariable_name = rParameters["variable_name"].GetString();

        mIsSpatialField = rParameters["spatially_varying"].GetBool();
        
        mpfunction = PythonGenericFunctionUtility::Pointer( new PythonGenericFunctionUtility(rParameters["value"].GetString(),  rParameters["local_axes"]));

//         if( KratosComponents< Variable<Vector> >::Has( mvariable_name ) == false ) //case of Vector variable
//         {
//         KRATOS_THROW_ERROR(std::runtime_error,"trying to set a variable that is not in the model_part - variable name is ",mvariable_name); 
//         }
        
        KRATOS_CATCH("");
    }



    /// Destructor.
    virtual ~AssignScalarFieldToConditionsProcess() {}


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
    virtual void Execute() override
    {

        KRATOS_TRY;

        ProcessInfo& rCurrentProcessInfo = mr_model_part.GetProcessInfo();

        const double& rCurrentTime = rCurrentProcessInfo[TIME];


        if( KratosComponents< Variable<Vector> >::Has( mvariable_name ) ) //case of double variable
        {
            InternalAssignValue<>(KratosComponents< Variable<Vector> >::Get(mvariable_name), rCurrentTime);
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not able to set the variable. Attempting to set variable:",mvariable_name);
        }

        KRATOS_CATCH("");

    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep() override
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
    virtual std::string Info() const
    {
        return "AssignScalarFieldToConditionsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AssignScalarFieldToConditionsProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    ModelPart& mr_model_part;
    PythonGenericFunctionUtility::Pointer mpfunction;
    std::string mvariable_name;

    bool mIsSpatialField = true;

    std::size_t mmesh_id;

    ///@}
    ///@name Private Operators
    ///@{

    void CallFunction(const Condition::Pointer& pCondition, const double& time, Vector& rValue)
    {
    
    Condition::GeometryType& rConditionGeometry = pCondition->GetGeometry();
    unsigned int size = rConditionGeometry.size();
    
    if(rValue.size() !=  size)
        rValue.resize(size,false);
        
    for(unsigned int i=0; i<size; i++)
        {
            rValue[i] = mpfunction->CallFunction(rConditionGeometry[i].X(),rConditionGeometry[i].Y(),rConditionGeometry[i].Z(),time  );
        }
        
    }

    void CallFunctionLocalSystem(const Condition::Pointer& pCondition, const double& time, Vector& rValue)
    {
    
    Condition::GeometryType& rConditionGeometry = pCondition->GetGeometry();
    unsigned int size = rConditionGeometry.size();
    
    if(rValue.size() !=  size)
        rValue.resize(size,false);
        
    for(unsigned int i=0; i<size; i++)
        {
            rValue[i] = mpfunction->RotateAndCallFunction(rConditionGeometry[i].X(),rConditionGeometry[i].Y(),rConditionGeometry[i].Z(),time  );
        }
    }
    
    void AssignTimeDependentValue(const Condition::Pointer& pCondition, const double& time, Vector& rValue, const double value)
    {
    
    Condition::GeometryType& rConditionGeometry = pCondition->GetGeometry();
    unsigned int size = rConditionGeometry.size();
    
    if(rValue.size() !=  size)
        rValue.resize(size,false);
        
    for(unsigned int i=0; i<size; i++)
        {
            rValue[i] = value;
        }
    }
    
    template< class TVarType >
    void InternalAssignValue(TVarType& rVar, const double& rTime)
    {
        const int nconditions = mr_model_part.GetMesh(mmesh_id).Conditions().size();

        Vector Value;
        
        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).ConditionsBegin();

            if(mIsSpatialField)
            {
                if(mpfunction->UseLocalSystem())
                {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(int i = 0; i<nconditions; i++)
                    {
                        ModelPart::ConditionsContainerType::iterator it = it_begin + i;
                        this->CallFunctionLocalSystem(*(it.base()), rTime, Value);
                        it->SetValue(rVar, Value);
                    }
                } 
                else
                {
                    // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(int i = 0; i<nconditions; i++)
                    {
                        ModelPart::ConditionsContainerType::iterator it = it_begin + i;
                        this->CallFunction(*(it.base()), rTime, Value);
                        it->SetValue(rVar, Value);
                    }
                }
            }
            else                                            // only varies in time
            {
                const double time_value = mpfunction->CallFunction(0.0, 0.0, 0.0,  rTime);
                // WARNING: do not parallelize with openmp. python GIL prevents it
                    for(int i = 0; i<nconditions; i++)
                    {
                        ModelPart::ConditionsContainerType::iterator it = it_begin + i;
                        this->AssignTimeDependentValue(*(it.base()), rTime, Value,  time_value);
                        it->SetValue(rVar, Value);
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
