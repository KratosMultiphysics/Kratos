//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

#if !defined(KRATOS_SET_INITIAL_STATE_H_INCLUDED )
#define  KRATOS_SET_INITIAL_STATE_H_INCLUDED

// System includes

// External includes

// Project includes

#include "processes/process.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/initial_state.h"
#include "includes/mat_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The SetInitialStateProcess.
/** This Operation is a derived class from the process.h
 *  
*/
template<SizeType TDim>
class SetInitialStateProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(SetInitialStateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SetInitialStateProcess(ModelPart& rModelPart) :
        mrModelPart(rModelPart);

    /// Full constructor.
    SetInitialStateProcess(
        ModelPart& rModelPart,
        const Vector& rInitialStrain,
        const Vector& rInitialStress,
        const Matrix& rInitialF) :
        mrModelPart(rModelPart), mInitialStrain(rInitialStrain), 
        mInitialStress(rInitialStress), mInitialF(rInitialF)
    {
    }

    /// Constructor with imposed vector.
    SetInitialStateProcess(
        ModelPart& rModelPart,
        const Vector& rInitialStateVector, 
        const int InitialStateType) :
        mrModelPart(rModelPart);

    /// Constructor with imposed F.
    SetInitialStateProcess(ModelPart& rModelPart,
        const Matrix& rInitialStateF) :
        mrModelPart(rModelPart);

    /// Destructor.
    ~SetInitialStateProcess() override {}


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


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override;

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

private:
    ///@name Static Member Variables
    ///@{


    ModelPart& mrModelPart;

    Vector mInitialStrain;
    Vector mInitialStress;
    Matrix mInitialF;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SetInitialStateProcess& operator=(Process const& rOther);

    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


}  // namespace Kratos.

#endif //  defined 


