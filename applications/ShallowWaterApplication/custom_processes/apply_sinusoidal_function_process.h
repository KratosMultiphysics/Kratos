//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_APPLY_SINUSOIDAL_FUNCTION_PROCESS_H_INCLUDED
#define KRATOS_APPLY_SINUSOIDAL_FUNCTION_PROCESS_H_INCLUDED


// System includes
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

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
 * @ingroup ShallowWaterApplication
 * @class ApplySinusoidalFunctionProcess
 * @brief The aim of this process is to generate sinusoidal waves
 */
template< class TVarType >
class KRATOS_API(SHALLOW_WATER_APPLICATION) ApplySinusoidalFunctionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t                     IndexType;
    typedef Node<3>                         NodeType;

    /// Pointer definition of ApplySinusoidalFunctionProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplySinusoidalFunctionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ApplySinusoidalFunctionProcess(
        ModelPart& rThisModelPart,
        TVarType& rThisVariable,
        Parameters& rThisParameters
    );

    /// Destructor.
    ~ApplySinusoidalFunctionProcess() override {};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief Perform a check with the parameters.
     */
    int Check() override;

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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ApplySinusoidalFunctionProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "ApplySinusoidalFunctionProcess";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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

    ModelPart& mrModelPart;
    TVarType& mrVariable;
    double mAmplitude;
    double mPeriod;
    double mAngularFrequency;
    double mPhase;
    double mVerticalShift;
    double mSmoothTime;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ValidateParameters(Parameters& rParameters);

    double Function(double& rTime);

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
    ApplySinusoidalFunctionProcess& operator=(ApplySinusoidalFunctionProcess const& rOther);

    /// Copy constructor.
    ApplySinusoidalFunctionProcess(ApplySinusoidalFunctionProcess const& rOther);


    ///@}

}; // Class ApplySinusoidalFunctionProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 ApplySinusoidalFunctionProcess& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const ApplySinusoidalFunctionProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_APPLY_SINUSOIDAL_FUNCTION_PROCESS_H_INCLUDED  defined
